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

#include <Geom_BSplineCurve.hxx>
#include <Geom_BSplineSurface.hxx>
#include <GeomFill_AppSurf.hxx>
#include <GeomFill_Line.hxx>
#include <GeomFill_SectionGenerator.hxx>
#include <Precision.hxx>
#include <Standard_ConstructionError.hxx>
#include <Standard_NullObject.hxx>
#include <gtest/gtest.h>
#include <limits>

namespace
{
occ::handle<Geom_BSplineCurve> profile(const double z)
{
  NCollection_Array1<gp_Pnt> p(1, 3);
  p(1) = gp_Pnt(0, 0, z);
  p(2) = gp_Pnt(1, 1, z);
  p(3) = gp_Pnt(2, 0, z);
  NCollection_Array1<double> k(1, 2), w(1, 3);
  k(1) = 0;
  k(2) = 1;
  w(1) = 1;
  w(2) = 0.8;
  w(3) = 1.1;
  NCollection_Array1<int> m(1, 2);
  m.Init(3);
  return new Geom_BSplineCurve(p, w, k, m, 2);
}
}

TEST(GeomFill_SectionGenerator_Test, AppSurfPreservesEndpointTangentDirections)
{
  for (int n : {2, 3, 4, 5, 10})
    for (auto param : {Approx_IsoParametric, Approx_ChordLength, Approx_Centripetal})
      for (int ends : {1, 2, 3})
      {
        SCOPED_TRACE(testing::Message()
                     << "sections=" << n << " ends=" << ends << " param=" << param);
        GeomFill_SectionGenerator g;
        for (int i = 0; i < n; ++i)
          g.AddCurve(profile(2.0 * i + 0.2 * i * i));
        const gp_Vec d1(0.2, 0.1, 2), d2(-0.3, 0.2, 3);
        auto         first = occ::down_cast<Geom_Curve>(g.Curve(1)->Translated(d1));
        auto         last  = occ::down_cast<Geom_Curve>(g.Curve(n)->Translated(d2));
        g.SetTangents(ends & 1 ? first : nullptr, ends & 2 ? last : nullptr);
        g.Perform(Precision::PConfusion());
        GeomFill_AppSurf app(2, 5, 1.e-7, 1.e-7, 0);
        app.SetParType(param);
        app.Perform(new GeomFill_Line(n), g, true);
        ASSERT_TRUE(app.IsDone());
        Geom_BSplineSurface surface(app.SurfPoles(),
                                    app.SurfWeights(),
                                    app.SurfUKnots(),
                                    app.SurfVKnots(),
                                    app.SurfUMults(),
                                    app.SurfVMults(),
                                    app.UDegree(),
                                    app.VDegree());
        double              u0, u1, v0, v1;
        surface.Bounds(u0, u1, v0, v1);
        for (int end : {0, 1})
          for (double u : {0.0, 0.25, 0.5, 0.75, 1.0})
          {
            gp_Pnt p;
            gp_Vec du, dv;
            surface.D1(u, end ? v1 : v0, p, du, dv);
            EXPECT_LT(p.Distance(g.Curve(end ? n : 1)->Value(u)), 1.e-7);
            if (!(ends & (1 << end)))
              continue;
            const auto expected = end ? d2 : d1;
            EXPECT_GT(dv.Dot(expected), 0.0);
            EXPECT_LT(dv.Crossed(expected).Magnitude() / (dv.Magnitude() * expected.Magnitude()),
                      1.e-7);
          }
      }
}

TEST(GeomFill_SectionGenerator_Test, TangentCurvesShareProfileWithoutModifyingInputs)
{
  GeomFill_SectionGenerator g;
  auto                      first = profile(0), last = profile(5), tangent = profile(2);
  tangent->InsertKnot(0.4, 1, Precision::PConfusion());
  g.AddCurve(first);
  g.AddCurve(last);
  g.SetTangents(tangent, nullptr);
  GeomFill_Profiler& profiler = g;
  profiler.Perform(Precision::PConfusion());
  EXPECT_EQ(g.NbKnots(), 3);
  EXPECT_EQ(first->NbKnots(), 2);
  EXPECT_EQ(last->NbKnots(), 2);
  EXPECT_EQ(tangent->NbKnots(), 3);
  NCollection_Array1<gp_Pnt>   p(1, g.NbPoles());
  NCollection_Array1<gp_Vec>   d(1, g.NbPoles());
  NCollection_Array1<double>   w(1, g.NbPoles()), dw(1, g.NbPoles());
  NCollection_Array1<gp_Pnt2d> p2;
  NCollection_Array1<gp_Vec2d> d2;
  ASSERT_TRUE(g.Section(1, p, d, p2, d2, w, dw));
  for (const auto& v : d)
    EXPECT_LT(v.Crossed(gp_Vec(0, 0, 1)).Magnitude(), 1.e-10);
  g.SetTangents(nullptr, nullptr);
  g.Perform(Precision::PConfusion());
  EXPECT_FALSE(g.Section(1, p, d, p2, d2, w, dw));
}

TEST(GeomFill_SectionGenerator_Test, UnconstrainedTwoSectionsRemainLinear)
{
  GeomFill_SectionGenerator g;
  g.AddCurve(profile(0));
  g.AddCurve(profile(5));
  g.Perform(Precision::PConfusion());
  GeomFill_AppSurf app(2, 5, 1.e-7, 1.e-7, 0);
  app.Perform(new GeomFill_Line(2), g, true);
  ASSERT_TRUE(app.IsDone());
  EXPECT_EQ(app.VDegree(), 1);
}

TEST(GeomFill_SectionGenerator_Test, ZeroTangentFailsWithoutProducingInvalidSurface)
{
  GeomFill_SectionGenerator g;
  auto                      first = profile(0), last = profile(5);
  g.AddCurve(first);
  g.AddCurve(last);
  g.SetTangents(first, last);
  g.Perform(Precision::PConfusion());
  GeomFill_AppSurf app(2, 5, 1.e-7, 1.e-7, 0);
  app.Perform(new GeomFill_Line(2), g, true);
  EXPECT_FALSE(app.IsDone());
}

TEST(GeomFill_SectionGenerator, ParameterArrayOwnershipAndBounds)
{
  for (const int aLower : {-3, 0, 1, 4})
  {
    SCOPED_TRACE(aLower);
    occ::handle<NCollection_HArray1<double>> aParams =
      new NCollection_HArray1<double>(aLower, aLower + 2);
    for (int anIndex = 0; anIndex < 3; ++anIndex)
    {
      aParams->SetValue(aLower + anIndex, 2.0 + anIndex);
    }
    GeomFill_SectionGenerator aGenerator;
    aGenerator.SetParam(aParams);
    for (int anIndex = 0; anIndex < 3; ++anIndex)
    {
      EXPECT_DOUBLE_EQ(aParams->Value(aLower + anIndex), 2.0 + anIndex);
      EXPECT_DOUBLE_EQ(aGenerator.Parameter(anIndex + 1), 2.0 + anIndex);
    }
    aParams->SetValue(aLower, -99.0);
    EXPECT_DOUBLE_EQ(aGenerator.Parameter(1), 2.0);
    aGenerator.SetParam(aParams);
    EXPECT_DOUBLE_EQ(aGenerator.Parameter(1), -99.0);
  }
}

TEST(GeomFill_SectionGenerator, InvalidParametersDoNotReplacePreviousValues)
{
  GeomFill_SectionGenerator                aGenerator;
  occ::handle<NCollection_HArray1<double>> aParams = new NCollection_HArray1<double>(1, 2);
  aParams->SetValue(1, 0.0);
  aParams->SetValue(2, 1.0);
  aGenerator.SetParam(aParams);
  EXPECT_THROW(aGenerator.SetParam(nullptr), Standard_NullObject);
  aParams->SetValue(2, 0.0);
  EXPECT_THROW(aGenerator.SetParam(aParams), Standard_ConstructionError);
  aParams->SetValue(2, std::numeric_limits<double>::infinity());
  EXPECT_THROW(aGenerator.SetParam(aParams), Standard_ConstructionError);
  EXPECT_DOUBLE_EQ(aGenerator.Parameter(1), 0.0);
  EXPECT_DOUBLE_EQ(aGenerator.Parameter(2), 1.0);
}
