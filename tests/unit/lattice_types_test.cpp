//
// Created by genshen on 2018/11/13.
//

#include <gtest/gtest.h>
#include <lattice/lattice_types.h>

TEST(lattice_type_test_case1, lattice_type_test) {
  LatticeTypes type(LatticeTypes::Mo);
  EXPECT_EQ(type.isAtom(), true);
  EXPECT_EQ(type.isDumbbell(), false);

  LatticeTypes type2(LatticeTypes::MoMo);
  EXPECT_EQ(type2.isAtom(), false);
  EXPECT_EQ(type2.isDumbbell(), true);

  LatticeTypes type3(LatticeTypes::V);
  EXPECT_EQ(type3.isAtom(), false);
  EXPECT_EQ(type3.isDumbbell(), false);
}

TEST(lattice_type_combineToInter_test, lattice_type_test) {
  auto tp_new = LatticeTypes::combineToInter(LatticeTypes::Mo, LatticeTypes::Re);
  EXPECT_EQ(tp_new, LatticeTypes::MoRe);

  auto tp_new2 = LatticeTypes::combineToInter(LatticeTypes::Re, LatticeTypes::Mo);
  EXPECT_EQ(tp_new2, LatticeTypes::MoRe);

  auto tp_new3 = LatticeTypes::combineToInter(LatticeTypes::Mn, LatticeTypes::Mn);
  EXPECT_EQ(tp_new3, LatticeTypes::MnMn);

  auto tp_new4 = LatticeTypes::combineToInter(LatticeTypes::Re, LatticeTypes::Mn);
  EXPECT_EQ(tp_new4, LatticeTypes::MoMn);
}

// tests for special cases of combineToInter: combine v
TEST(lattice_type_combineToInter_V_test, lattice_type_test) {
  auto tp_new = LatticeTypes::combineToInter(LatticeTypes::V, LatticeTypes::Re);
  EXPECT_EQ(tp_new, LatticeTypes::Re);

  auto tp_new2 = LatticeTypes::combineToInter(LatticeTypes::V, LatticeTypes::V);
  EXPECT_EQ(tp_new2, LatticeTypes::V);

  auto tp_new3 = LatticeTypes::combineToInter(LatticeTypes::Mo, LatticeTypes::V);
  EXPECT_EQ(tp_new3, LatticeTypes::Mo);
}

// tests for special cases of diff: diff v
TEST(lattice_type_diff_V_test, lattice_type_test) {
  EXPECT_EQ(LatticeTypes{LatticeTypes::V}.diff(LatticeTypes{LatticeTypes::V}), LatticeTypes::V);
  EXPECT_EQ(LatticeTypes{LatticeTypes::MoRe}.diff(LatticeTypes{LatticeTypes::V}), LatticeTypes::MoRe);
  EXPECT_EQ(LatticeTypes{LatticeTypes::Mo}.diff(LatticeTypes{LatticeTypes::Re}),
            LatticeTypes::Mo); // not change
  EXPECT_EQ(LatticeTypes{LatticeTypes::V}.diff(LatticeTypes{LatticeTypes::Re}), LatticeTypes::V);
}

TEST(lattice_type_getHighLowEnd_test, lattice_type_test) {
  LatticeTypes type1(LatticeTypes::MoRe);
  EXPECT_EQ(type1.getHighEnd(), LatticeTypes::Mo);
  EXPECT_EQ(type1.getLowEnd(), LatticeTypes::Re);

  LatticeTypes type2(LatticeTypes::Re);
  EXPECT_EQ(type2.getHighEnd(), LatticeTypes::V);
  EXPECT_EQ(type2.getLowEnd(), LatticeTypes::Re);

  LatticeTypes type3(LatticeTypes::V);
  EXPECT_EQ(type3.getHighEnd(), LatticeTypes::V);
  EXPECT_EQ(type3.getLowEnd(), LatticeTypes::V);
}

TEST(lattice_type_isHighLowEnd_test, lattice_type_test) {
  LatticeTypes type1(LatticeTypes::MoRe);
  EXPECT_EQ(type1.isHighEnd(LatticeTypes::Mo), true);
  EXPECT_EQ(type1.isHighEnd(LatticeTypes::V), false);
  EXPECT_EQ(type1.isHighEnd(LatticeTypes::Re), false);
  EXPECT_EQ(type1.isLowEnd(LatticeTypes::Mo), false);
  EXPECT_EQ(type1.isLowEnd(LatticeTypes::V), false);
  EXPECT_EQ(type1.isLowEnd(LatticeTypes::Re), true);

  LatticeTypes type2(LatticeTypes::ReRe);
  EXPECT_EQ(type2.isHighEnd(LatticeTypes::Re), true);
  EXPECT_EQ(type2.isLowEnd(LatticeTypes::Re), true);
}

TEST(lattice_type_FirstSecondAtom_test, lattice_type_test) {
  LatticeTypes type1(LatticeTypes::MoRe);
  EXPECT_EQ(type1.getFirst(true), LatticeTypes::Mo);
  EXPECT_EQ(type1.getSecond(true), LatticeTypes::Re);

  EXPECT_EQ(type1.getFirst(false), LatticeTypes::Re);
  EXPECT_EQ(type1.getSecond(false), LatticeTypes::Mo);

  LatticeTypes type3(LatticeTypes::ReRe);
  EXPECT_EQ(type3.getFirst(true), LatticeTypes::Re);
  EXPECT_EQ(type3.getSecond(true), LatticeTypes::Re);
  EXPECT_EQ(type3.getFirst(false), LatticeTypes::Re);
  EXPECT_EQ(type3.getSecond(false), LatticeTypes::Re);
}

TEST(lattice_type_diff_test, lattice_type_test) {
  LatticeTypes type1(LatticeTypes::MoRe);
  LatticeTypes type2(LatticeTypes::Re);
  EXPECT_EQ(type1.diff(type2), LatticeTypes::Mo);

  LatticeTypes type3(LatticeTypes::V);
  EXPECT_EQ(type1.diff(type3), LatticeTypes::MoRe);

  LatticeTypes type4(LatticeTypes::Ni);
  EXPECT_EQ(type1.diff(type4), LatticeTypes::MoRe);
}
