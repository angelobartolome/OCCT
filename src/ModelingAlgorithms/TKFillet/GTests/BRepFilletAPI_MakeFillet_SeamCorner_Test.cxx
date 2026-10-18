// Copyright (c) 2026 OPEN CASCADE SAS
//
// This file is part of Open CASCADE Technology software library.
//
// This library is free software; you can redistribute it and/or modify it under
// the terms of the GNU Lesser General Public License version 2.1 as published
// by the Free Software Foundation, with special exception defined in the file
// OCCT_LGPL_EXCEPTION.txt. Consult the file LICENSE_LGPL_21.txt included in OCCT
// distribution for complete text of the license and disclaimer of any warranty.
//
// Alternatively, this file may be used under the terms of Open CASCADE
// commercial license or contractual agreement.

// Fillets and chamfers on an edge that ends at a vertex a periodic face's seam runs to.
//
// The model: a block with a hexagonal hole (apothem 16, vertices at 0, 60, ... degrees in the XZ
// plane) below y = -1, and a round hole (R16, along Y) above it. The circle touches the hexagon's
// sides at 30, 90, ... degrees, so the step at y = -1 is six faces that meet only at those points.
// The edge blended is the arc of the step from -30 to 30 degrees, between a step face and the
// round hole's cylinder. Only where the cylinder's seam lies changes between the tests.

#include <BRepAdaptor_Surface.hxx>
#include <BRepAlgoAPI_Cut.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakePolygon.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <BRepFilletAPI_MakeChamfer.hxx>
#include <BRepFilletAPI_MakeFillet.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepPrimAPI_MakePrism.hxx>
#include <BRep_Tool.hxx>
#include <NCollection_IndexedDataMap.hxx>
#include <NCollection_IndexedMap.hxx>
#include <NCollection_List.hxx>
#include <TopExp.hxx>
#include <TopTools_ShapeMapHasher.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Vertex.hxx>
#include <gp_Ax2.hxx>
#include <gp_Trsf.hxx>

#include <gtest/gtest.h>

#include <cmath>

namespace
{
constexpr double THE_SIZE = 0.2;

//! The stepped hole, its round hole's cylinder with its seam at theSeamDeg (measured from +X
//! towards +Z). theIndirect gives the cylinder a left-handed frame, as a face unified from
//! pieces can have.
TopoDS_Shape makeSteppedHole(const double theSeamDeg, const bool theIndirect)
{
  const TopoDS_Shape aBlock = BRepPrimAPI_MakeBox(gp_Pnt(-30, -11, -30), gp_Pnt(30, 0, 30)).Shape();

  const double               aCircum = 16.0 / std::cos(M_PI / 6.0);
  BRepBuilderAPI_MakePolygon aHexagon;
  for (int i = 0; i < 6; ++i)
  {
    aHexagon.Add(
      gp_Pnt(aCircum * std::cos(i * M_PI / 3.0), -12, aCircum * std::sin(i * M_PI / 3.0)));
  }
  aHexagon.Close();
  const TopoDS_Shape aHexHole =
    BRepPrimAPI_MakePrism(BRepBuilderAPI_MakeFace(aHexagon.Wire()).Face(), gp_Vec(0, 11, 0))
      .Shape();

  // Mirrored in z = 0, the cylinder's frame is left-handed; its x direction is mirrored first.
  const double aSeam = theSeamDeg * M_PI / 180.0;
  const gp_Ax2 anAxis(gp_Pnt(0, -1, 0),
                      gp_Dir(0, 1, 0),
                      gp_Dir(std::cos(aSeam), 0, theIndirect ? -std::sin(aSeam) : std::sin(aSeam)));
  TopoDS_Shape aRoundHole = BRepPrimAPI_MakeCylinder(anAxis, 16, 2).Shape();
  if (theIndirect)
  {
    gp_Trsf aMirror;
    aMirror.SetMirror(gp_Ax2(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1)));
    aRoundHole = BRepBuilderAPI_Transform(aRoundHole, aMirror, true).Shape();
  }

  return BRepAlgoAPI_Cut(BRepAlgoAPI_Cut(aBlock, aHexHole).Shape(), aRoundHole).Shape();
}

//! The arc of the step at y = -1 from -30 to 30 degrees, or a null edge.
TopoDS_Edge findStepArc(const TopoDS_Shape& theShape)
{
  const gp_Pnt aFrom(16.0 * std::cos(M_PI / 6.0), -1, -8);
  const gp_Pnt aTo(16.0 * std::cos(M_PI / 6.0), -1, 8);
  NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> anEdges;
  TopExp::MapShapes(theShape, TopAbs_EDGE, anEdges);
  for (int i = 1; i <= anEdges.Extent(); ++i)
  {
    TopoDS_Vertex aV1, aV2;
    TopExp::Vertices(TopoDS::Edge(anEdges(i)), aV1, aV2);
    const gp_Pnt aP1 = BRep_Tool::Pnt(aV1), aP2 = BRep_Tool::Pnt(aV2);
    if ((aP1.Distance(aFrom) < 1.e-4 && aP2.Distance(aTo) < 1.e-4)
        || (aP1.Distance(aTo) < 1.e-4 && aP2.Distance(aFrom) < 1.e-4))
    {
      return TopoDS::Edge(anEdges(i));
    }
  }
  return TopoDS_Edge();
}

//! The face of theShape on theEdge whose surface is of theType: the step face (a plane) or the
//! round hole's cylinder.
TopoDS_Face findFaceOf(const TopoDS_Shape&       theShape,
                       const TopoDS_Edge&        theEdge,
                       const GeomAbs_SurfaceType theType)
{
  NCollection_IndexedDataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>, TopTools_ShapeMapHasher>
    anEdgeFaces;
  TopExp::MapShapesAndAncestors(theShape, TopAbs_EDGE, TopAbs_FACE, anEdgeFaces);
  for (const TopoDS_Shape& aFace : anEdgeFaces.FindFromKey(theEdge))
  {
    if (BRepAdaptor_Surface(TopoDS::Face(aFace)).GetType() == theType)
    {
      return TopoDS::Face(aFace);
    }
  }
  return TopoDS_Face();
}

} // namespace

// The seam at 30 degrees runs to the arc's end vertex, here with the cylinder's frame left-handed:
// the corner code met the piece of extension left with no face on the cylinder's side, and the
// fillet was not done.
TEST(BRepFilletAPI_MakeFilletTest, SeamCorner_SeamAtArcEnd_IndirectCylinder)
{
  const TopoDS_Shape aShape = makeSteppedHole(30.0, true);
  const TopoDS_Edge  anArc  = findStepArc(aShape);
  ASSERT_FALSE(anArc.IsNull()) << "the step's arc from -30 to 30 degrees";

  BRepFilletAPI_MakeFillet aFillet(aShape);
  aFillet.Add(THE_SIZE, anArc);
  ASSERT_NO_THROW(aFillet.Build());
  EXPECT_TRUE(aFillet.IsDone());
}

// A distance-angle chamfer on the same corner, measured on the step face: it used to crash, on the
// piece of extension left with no face on the cylinder's side.
TEST(BRepFilletAPI_MakeChamferTest, SeamCorner_SeamAtArcEnd_DistAngle)
{
  for (const bool isIndirect : {false, true})
  {
    const TopoDS_Shape aShape = makeSteppedHole(30.0, isIndirect);
    const TopoDS_Edge  anArc  = findStepArc(aShape);
    ASSERT_FALSE(anArc.IsNull()) << "the step's arc from -30 to 30 degrees";

    BRepFilletAPI_MakeChamfer aChamfer(aShape);
    aChamfer.AddDA(THE_SIZE, M_PI / 4.0, anArc, findFaceOf(aShape, anArc, GeomAbs_Plane));
    ASSERT_NO_THROW(aChamfer.Build());
    EXPECT_TRUE(aChamfer.IsDone());
  }
}
