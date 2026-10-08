//
// Created by genshen on 2019-03-08.
//

#include "counter.h"
#include "type_define.h"
#include <cassert>
#include <iostream>

void counter::setLatTypeToStrFunc(const fn_lat_type_to_str fn_lat_type_to_str) {
  this->func_lat_type_to_str = fn_lat_type_to_str;
}

int counter::getAtomCount(const LatticeTypes::lat_type tp) {
#ifdef KMC_DEBUG_MODE
  assert(LatticeTypes{tp}.isAtom() || LatticeTypes{tp}.isVacancy());
#endif
  if (LatticeTypes{tp}.isVacancy()) {
    return data[tp];
  }
  int count = 0;
  LatticeTypes::lat_type high, low;
  for (const auto pair : data) {
    high = LatticeTypes{pair.first}.getHighEnd();
    low = LatticeTypes{pair.first}.getLowEnd();
    if (high == tp) {
      count += pair.second;
    }
    if (low == tp) {
      count += pair.second;
    }
  }
  return count;
}

counter counter::newCounter(LatticesList *p_list) {
  counter c;
  // for (_type_lattice_coord z = 0; z < p_list->meta.box_z; z++) {
  //   for (_type_lattice_coord y = 0; y < p_list->meta.box_y; y++) {
  //     // note: x is already doubled.
  //     for (_type_lattice_coord x = 0; x < p_list->meta.box_x; x++) {
  //       _type_lattice_coord latti_id = p_list->getId(p_list->meta.ghost_x + x, p_list->meta.ghost_y + y, p_list->meta.ghost_z + z);
  //       if (p_list->vac_hash.count(latti_id)) {
  //         c.add(LatticeTypes::V);
  //       } else if (p_list->re_hash.count(latti_id)){
  //         c.add(LatticeTypes::Re);
  //       } else if (p_list->mn_hash.count(latti_id)){
  //         c.add(LatticeTypes::Mn);
  //       } else if (p_list->ni_hash.count(latti_id)){
  //         c.add(LatticeTypes::Ni);
  //       } else if (p_list->si_hash.count(latti_id)){
  //         c.add(LatticeTypes::Si);
  //       } else if (p_list->momo_hash.count(latti_id)){
  //         c.add(LatticeTypes::MoMo);
  //       } else if (p_list->more_hash.count(latti_id)){
  //         c.add(LatticeTypes::MoRe);
  //       } else if (p_list->rere_hash.count(latti_id)){
  //         c.add(LatticeTypes::ReRe);
  //       } else {
  //         c.add(LatticeTypes::Mo);
  //       }
  //     }
  //   }
  // }
  return c;
}

std::ostream &operator<<(std::ostream &os, const counter &counter) {
  if (counter.func_lat_type_to_str) {
    for (const auto &c : counter.data) {
      os << "[" << counter.func_lat_type_to_str(c.first) << "]: " << c.second << "\n";
    }
  } else {
    for (const auto &c : counter.data) {
      os << "[" << c.first << "]: " << c.second << "\n";
    }
  }
  return os;
}