//
// Created by genshen on 2019-03-08.
//

#include <counter.h>
#include <gtest/gtest.h>
#include <lattice/lattice_types.h>

TEST(counter_add_test, counter_test) {
  counter c;
  c.add(LatticeTypes::Mo);
  EXPECT_EQ(c.get(LatticeTypes::Mo), 1);
  c.add(LatticeTypes::Mo);
  EXPECT_EQ(c.get(LatticeTypes::Mo), 2);

  c.set(LatticeTypes::Mo, 1024);
  c.add(LatticeTypes::Mo);
  EXPECT_EQ(c.get(LatticeTypes::Mo), 1025);

  c.set(LatticeTypes::Re, 1024);
  c.add(LatticeTypes::Re);
  EXPECT_EQ(c.get(LatticeTypes::Re), 1025);
}

TEST(counter_sub_test, counter_test) {
  counter c;
  c.add(LatticeTypes::Mo);
  c.subtract(LatticeTypes::Mo);
  EXPECT_EQ(c.get(LatticeTypes::Mo), 0);

  c.set(LatticeTypes::Mo, 1024);
  c.subtract(LatticeTypes::Mo);
  EXPECT_EQ(c.get(LatticeTypes::Mo), 1023);

  c.set(LatticeTypes::Re, 1025);
  c.subtract(LatticeTypes::Re);
  EXPECT_EQ(c.get(LatticeTypes::Re), 1024);
}

TEST(counter_getAtomCount_test, counter_test) {
  counter c;
  c.set(LatticeTypes::V, 1);
  c.set(LatticeTypes::Mo, 2);
  c.set(LatticeTypes::Re, 4);
  c.set(LatticeTypes::ReRe, 8);
  c.set(LatticeTypes::ReNi, 16);
  c.set(LatticeTypes::MoMo, 32);

  EXPECT_EQ(c.getAtomCount(LatticeTypes::V), 1);
  EXPECT_EQ(c.getAtomCount(LatticeTypes::Mo), 2 + 2 * 32);
  EXPECT_EQ(c.getAtomCount(LatticeTypes::Re), 4 + 2 * 8 + 16);
  EXPECT_EQ(c.getAtomCount(LatticeTypes::Ni), 16);
}
