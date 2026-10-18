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
#include <BRepCheck_Analyzer.hxx>
#include <BRepFilletAPI_MakeChamfer.hxx>
#include <BRepFilletAPI_MakeFillet.hxx>
#include <BRepGProp.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepPrimAPI_MakePrism.hxx>
#include <BRep_Tool.hxx>
#include <GProp_GProps.hxx>
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

// The corner's faces are built to their tolerance only, and a corner whose connecting curve is cut
// at the seam integrates to a volume about 0.7% larger than the same corner uncut: the removed
// volumes are compared relatively.
constexpr double THE_VOLUME_TOLERANCE = 0.01;

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

double volumeOf(const TopoDS_Shape& theShape)
{
  GProp_GProps aProps;
  BRepGProp::VolumeProperties(theShape, aProps);
  return aProps.Mass();
}

//! What a THE_SIZE fillet on the step's arc removes, or -1 if it fails.
double filletRemoval(const TopoDS_Shape& theShape)
{
  BRepFilletAPI_MakeFillet aFillet(theShape);
  aFillet.Add(THE_SIZE, findStepArc(theShape));
  aFillet.Build();
  if (!aFillet.IsDone() || !BRepCheck_Analyzer(aFillet.Shape()).IsValid())
  {
    return -1.0;
  }
  return volumeOf(theShape) - volumeOf(aFillet.Shape());
}

//! What a THE_SIZE, 45 degree chamfer on the step's arc, measured on the step face, removes, or -1
//! if it fails.
double chamferRemoval(const TopoDS_Shape& theShape)
{
  const TopoDS_Edge         anArc = findStepArc(theShape);
  BRepFilletAPI_MakeChamfer aChamfer(theShape);
  aChamfer.AddDA(THE_SIZE, M_PI / 4.0, anArc, findFaceOf(theShape, anArc, GeomAbs_Plane));
  aChamfer.Build();
  if (!aChamfer.IsDone() || !BRepCheck_Analyzer(aChamfer.Shape()).IsValid())
  {
    return -1.0;
  }
  return volumeOf(theShape) - volumeOf(aChamfer.Shape());
}

void checkStepArcFillet(const double theSeamDeg, const bool theIndirect)
{
  const TopoDS_Shape aShape = makeSteppedHole(theSeamDeg, theIndirect);
  ASSERT_TRUE(BRepCheck_Analyzer(aShape).IsValid());
  const TopoDS_Edge anArc = findStepArc(aShape);
  ASSERT_FALSE(anArc.IsNull()) << "the step's arc from -30 to 30 degrees";

  BRepFilletAPI_MakeFillet aFillet(aShape);
  aFillet.Add(THE_SIZE, anArc);
  ASSERT_NO_THROW(aFillet.Build());
  ASSERT_TRUE(aFillet.IsDone());
  EXPECT_TRUE(BRepCheck_Analyzer(aFillet.Shape()).IsValid());
  // The same arc on the same faces: the seam's place must not change what the fillet removes.
  // The control has the seam at 90 degrees, where it can't matter.
  const double anExpected = filletRemoval(makeSteppedHole(90.0, false));
  ASSERT_GT(anExpected, 0.0);
  EXPECT_NEAR(volumeOf(aShape) - volumeOf(aFillet.Shape()),
              anExpected,
              THE_VOLUME_TOLERANCE * anExpected);
}

void checkStepArcChamfer(const double theSeamDeg, const bool theIndirect)
{
  const TopoDS_Shape aShape = makeSteppedHole(theSeamDeg, theIndirect);
  const TopoDS_Edge  anArc  = findStepArc(aShape);
  ASSERT_FALSE(anArc.IsNull()) << "the step's arc from -30 to 30 degrees";

  BRepFilletAPI_MakeChamfer aChamfer(aShape);
  aChamfer.AddDA(THE_SIZE, M_PI / 4.0, anArc, findFaceOf(aShape, anArc, GeomAbs_Plane));
  ASSERT_NO_THROW(aChamfer.Build());
  ASSERT_TRUE(aChamfer.IsDone());
  EXPECT_TRUE(BRepCheck_Analyzer(aChamfer.Shape()).IsValid());
  const double anExpected = chamferRemoval(makeSteppedHole(90.0, false));
  ASSERT_GT(anExpected, 0.0);
  EXPECT_NEAR(volumeOf(aShape) - volumeOf(aChamfer.Shape()),
              anExpected,
              THE_VOLUME_TOLERANCE * anExpected);
}
} // namespace

// Controls: the seam at 90 degrees, away from the arc.
TEST(BRepFilletAPI_MakeFilletTest, SeamCorner_SeamAwayFromArc)
{
  checkStepArcFillet(90.0, false);
  checkStepArcFillet(90.0, true);
  checkStepArcChamfer(90.0, false);
}

// The seam at 30 degrees runs to the arc's end vertex. It used to be made an edge of the corner and
// given a point beyond its end, so the cylinder face was rebuilt unorientable and the shell open.
TEST(BRepFilletAPI_MakeFilletTest, SeamCorner_SeamAtArcEnd_DirectCylinder)
{
  checkStepArcFillet(30.0, false);
}

// The same with the cylinder's frame left-handed: the corner code met the piece of extension left
// with no face on the cylinder's side, and the fillet was not done.
TEST(BRepFilletAPI_MakeFilletTest, SeamCorner_SeamAtArcEnd_IndirectCylinder)
{
  checkStepArcFillet(30.0, true);
}

// A distance-angle chamfer on the same corner, measured on the step face: it used to crash, on the
// piece of extension left with no face on the cylinder's side.
TEST(BRepFilletAPI_MakeChamferTest, SeamCorner_SeamAtArcEnd_DistAngle)
{
  checkStepArcChamfer(30.0, false);
  checkStepArcChamfer(30.0, true);
}
