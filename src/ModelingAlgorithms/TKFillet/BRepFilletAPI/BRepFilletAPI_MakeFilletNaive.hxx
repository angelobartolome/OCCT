// Created on: 2025-08-26
// Created by: Angelo Bartolome <angelo.m.bartolome@gmail.com
// Copyright (c) 2025 OPEN CASCADE SAS
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

#ifndef _BRepFilletAPI_MakeFilletNaive_HeaderFile
#define _BRepFilletAPI_MakeFilletNaive_HeaderFile

#include <BRepFilletAPI_MakeFillet.hxx>
#include <ChFi3d_FilletShape.hxx>
#include <GeomAbs_Shape.hxx>
#include <Message_ProgressRange.hxx>
#include <TopoDS_Shape.hxx>

//! BRepFilletAPI_MakeFillet extension that can build fillets per edge and
//! intersect them to tolerate overlapping fillets (naive mode).
class BRepFilletAPI_MakeFilletNaive : public BRepFilletAPI_MakeFillet
{
public:
  DEFINE_STANDARD_ALLOC

  //! Create fillet maker with optional naive mode.
  //! If naiveMode is true, Build() will process each edge independently
  //! and intersect results (common) to avoid aborting on overlapping fillets.
  Standard_EXPORT BRepFilletAPI_MakeFilletNaive(
    const TopoDS_Shape&      theS,
    const ChFi3d_FilletShape theFShape    = ChFi3d_Rational,
    const Standard_Boolean   theNaiveMode = Standard_False);

  //! Enable/disable naive mode at runtime.
  Standard_EXPORT void SetNaiveMode(const Standard_Boolean theOn);

  //! Query naive mode.
  Standard_EXPORT Standard_Boolean NaiveMode() const { return myNaiveMode; }

  //! Override Build to implement naive behavior when enabled.
  Standard_EXPORT virtual void Build(
    const Message_ProgressRange& theRange = Message_ProgressRange()) Standard_OVERRIDE;

  //! Store builder parameters so they can be mirrored to per-edge builds.
  Standard_EXPORT void SetParams(const Standard_Real Tang,
                                 const Standard_Real Tesp,
                                 const Standard_Real T2d,
                                 const Standard_Real TApp3d,
                                 const Standard_Real TolApp2d,
                                 const Standard_Real Fleche);

  //! Store internal continuity so it can be mirrored to per-edge builds.
  Standard_EXPORT void SetContinuity(const GeomAbs_Shape InternalContinuity,
                                     const Standard_Real AngularTolerance);

private:
  //! Apply parameters/continuity captured on this object to another fillet maker.
  void applyCapturedSettings(BRepFilletAPI_MakeFillet& theDst) const;

private:
  TopoDS_Shape     myBaseShape;
  Standard_Boolean myNaiveMode;

  // Captured settings (optional)
  Standard_Boolean myHasParams = Standard_False;
  Standard_Real    myTang = 0.0, myTesp = 0.0, myT2d = 0.0, myTApp3d = 0.0, myTolApp2d = 0.0,
                myFleche                = 0.0;
  Standard_Boolean myHasContinuity      = Standard_False;
  GeomAbs_Shape    myInternalContinuity = GeomAbs_C1;
  Standard_Real    myAngularTolerance   = 0.0;
};

#endif // _BRepFilletAPI_MakeFilletNaive_HeaderFile
