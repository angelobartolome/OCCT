// Created on: 2025-08-26
// See header for license details

#include <BRepFilletAPI_MakeFilletNaive.hxx>

#include <BRepAlgoAPI_Common.hxx>
#include <Message_ProgressScope.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Shape.hxx>
#include <TopoDS_Vertex.hxx>

//-------------------------------------------------------------------------------------------------

BRepFilletAPI_MakeFilletNaive::BRepFilletAPI_MakeFilletNaive(const TopoDS_Shape&      theS,
                                                             const ChFi3d_FilletShape theFShape,
                                                             const Standard_Boolean   theNaiveMode)
    : BRepFilletAPI_MakeFillet(theS, theFShape),
      myBaseShape(theS),
      myNaiveMode(theNaiveMode)
{
}

//-------------------------------------------------------------------------------------------------

void BRepFilletAPI_MakeFilletNaive::SetNaiveMode(const Standard_Boolean theOn)
{
  myNaiveMode = theOn;
}

//-------------------------------------------------------------------------------------------------

void BRepFilletAPI_MakeFilletNaive::SetParams(const Standard_Real Tang,
                                              const Standard_Real Tesp,
                                              const Standard_Real T2d,
                                              const Standard_Real TApp3d,
                                              const Standard_Real TolApp2d,
                                              const Standard_Real Fleche)
{
  myHasParams = Standard_True;
  myTang      = Tang;
  myTesp      = Tesp;
  myT2d       = T2d;
  myTApp3d    = TApp3d;
  myTolApp2d  = TolApp2d;
  myFleche    = Fleche;
  // forward to base for current builder as well
  BRepFilletAPI_MakeFillet::SetParams(Tang, Tesp, T2d, TApp3d, TolApp2d, Fleche);
}

//-------------------------------------------------------------------------------------------------

void BRepFilletAPI_MakeFilletNaive::SetContinuity(const GeomAbs_Shape InternalContinuity,
                                                  const Standard_Real AngularTolerance)
{
  myHasContinuity      = Standard_True;
  myInternalContinuity = InternalContinuity;
  myAngularTolerance   = AngularTolerance;
  // forward to base for current builder as well
  BRepFilletAPI_MakeFillet::SetContinuity(InternalContinuity, AngularTolerance);
}

//-------------------------------------------------------------------------------------------------

void BRepFilletAPI_MakeFilletNaive::applyCapturedSettings(BRepFilletAPI_MakeFillet& theDst) const
{
  if (myHasParams)
  {
    theDst.SetParams(myTang, myTesp, myT2d, myTApp3d, myTolApp2d, myFleche);
  }
  if (myHasContinuity)
  {
    theDst.SetContinuity(myInternalContinuity, myAngularTolerance);
  }
  theDst.SetFilletShape(GetFilletShape());
}

//-------------------------------------------------------------------------------------------------

void BRepFilletAPI_MakeFilletNaive::Build(const Message_ProgressRange& theRange)
{
  if (!myNaiveMode)
  {
    // use standard behavior
    BRepFilletAPI_MakeFillet::Build(theRange);
    return;
  }

  // Naive mode: build each edge independently and intersect
  const Standard_Integer aNbContours = NbContours();
  if (aNbContours <= 0)
  {
    // nothing to do
    BRepFilletAPI_MakeFillet::Build(theRange);
    return;
  }

  // Count edges to setup progress
  Standard_Integer aTotalEdges = 0;
  for (Standard_Integer iC = 1; iC <= aNbContours; ++iC)
    aTotalEdges += NbEdges(iC);

  Message_ProgressScope aPS(theRange, "Naive fillet build", aTotalEdges > 0 ? aTotalEdges : 1);

  TopoDS_Shape aAccum = myBaseShape; // start from original

  // Process each edge separately
  for (Standard_Integer iC = 1; iC <= aNbContours; ++iC)
  {
    const Standard_Integer aNbEd = NbEdges(iC);
    for (Standard_Integer iE = 1; iE <= aNbEd; ++iE)
    {
      Message_ProgressRange aStep = aPS.Next();

      const TopoDS_Edge& aEdge = Edge(iC, iE);

      // Create a temporary fillet builder on the same base shape
      BRepFilletAPI_MakeFillet aTmp(myBaseShape, GetFilletShape());
      applyCapturedSettings(aTmp);

      // Add only this edge
      aTmp.Add(aEdge);

      // Propagate radius/law information if available for this edge
      if (IsConstant(iC, aEdge))
      {
        Standard_Real aRad = Radius(iC, aEdge);
        aTmp.SetRadius(aRad, aTmp.Contour(aEdge), aEdge);
      }
      else
      {
        Handle(Law_Function) aLaw = GetLaw(iC, aEdge);
        if (!aLaw.IsNull())
        {
          aTmp.SetLaw(aTmp.Contour(aEdge), aEdge, aLaw);
        }
      }

      // Try build this edge; on failure, skip but continue
      try
      {
        aTmp.Build(aStep);
      }
      catch (...)
      {
        // skip this edge on hard failure
        continue;
      }

      if (!aTmp.IsDone())
      {
        // If there is a partial result, use it; else skip
        if (!aTmp.HasResult())
          continue;
      }

      TopoDS_Shape aPart = aTmp.IsDone() ? aTmp.Shape() : aTmp.BadShape();
      if (aPart.IsNull())
        continue;

      // Intersect with accumulated result to apply combined material removals
      try
      {
        BRepAlgoAPI_Common aCommon(aAccum, aPart);
        aCommon.Build();
        if (aCommon.IsDone())
        {
          aAccum = aCommon.Shape();
        }
        // if intersect fails, keep previous accumulator and proceed
      }
      catch (...)
      {
        // ignore intersect errors and continue
      }
    }
  }

  // Store final result
  myShape = aAccum;
  Done();
}

//-------------------------------------------------------------------------------------------------
