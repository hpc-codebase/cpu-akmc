//
// Created by genshen on 2019-01-24.
//

#include <cassert>

#include "../env.h"
#include "bonds/bonds_counter.h"
#include "vacancy_rates_solver.h"
#include <logs/logs.h>
#include <iostream>
VacRatesSolver::VacRatesSolver(LatticesList &lat_list, const double v, const double T) : RatesSolver(lat_list, v, T) {}

double VacRatesSolver::e0(const LatticeTypes::lat_type ghost_atom) const {
  switch (ghost_atom) {
    case LatticeTypes::Mo:
      return 0.65;
    case LatticeTypes::Re:
      return 0.54;
    case LatticeTypes::Mn:
      return 1.03;
    case LatticeTypes::Ni:
      return 0.68;
    case LatticeTypes::Si:
      return 0.42;
    default:
      return 0;
  }
}

double VacRatesSolver::deltaE(_type_lattice_id source_id, _type_lattice_id target_id, Lattice &source_lattice,
                              Lattice &target_lattice, const LatticeTypes::lat_type ghost_atom) {
#ifdef KMC_DEBUG_MODE
  {
    const bool debug_bool_1 = target_lattice.type.isDumbbell(); // always be false.
    assert(!debug_bool_1);
    const bool debug_bool_2 = source_lattice.type.isVacancy(); // always be true.
    assert(debug_bool_2);
    // in vacancy transition, ghost atom equals to the type of target lattice(target lattice will be a single atom)
    assert(target_lattice.type._type == ghost_atom);
  };
#endif

#ifdef EAM_POT_ENABLED
  return 0; // todo eam
#else
  // calculate system energy before transition.
  assert(source_id == source_lattice.id);
  assert(target_id == target_lattice.id);
  // assert(source_lattice.type.type != ghost_atom);
  double e_before = 0;
  {
    // bonds energy of src lattice contributed by its 1nn/2nn neighbour lattice.
    bonds::_type_pair_ia e_src = bonds::BondsCounter::count(&lattice_list, source_id, source_lattice.type);
    // bonds energy of des lattice contributed by its 1nn/2nn neighbour lattice(it is an atom).
    bonds::_type_pair_ia e_des = bonds::BondsCounter::count(&lattice_list, target_id, target_lattice.type);
    e_before = e_src + e_des;
  }
  // kiwi::logs::v(" ", " e_before is : {} .\n", e_before);

  // exchange atoms of vacancy and neighbour lattice atom/vacancy.
  // note: ghost_atom equals to origin target_lattice.type
  if (lattice_list.vac_hash.count(source_lattice.id)) {
    // lattice_list.vac_hash.erase(it_vac);
    lattice_list.vac_hash.emplace(std::make_pair(target_lattice.id, VacancyHash{}));
  } else assert(false);
  switch (ghost_atom)
  {
    case LatticeTypes::V:
      assert(false);
      break;
    case LatticeTypes::Re:
      if (lattice_list.re_hash.count(target_lattice.id)) {
        lattice_list.re_hash.erase(target_lattice.id);
        lattice_list.re_hash.emplace(source_lattice.id);
      } else assert(false);
      break;
    case LatticeTypes::Mn:
      if (lattice_list.mn_hash.count(target_lattice.id)) {
        lattice_list.mn_hash.erase(target_lattice.id);
        lattice_list.mn_hash.emplace(source_lattice.id);
      } else assert(false);
      break;
    case LatticeTypes::Ni:
      if (lattice_list.ni_hash.count(target_lattice.id)) {
        lattice_list.ni_hash.erase(target_lattice.id);
        lattice_list.ni_hash.emplace(source_lattice.id);
      } else assert(false);
      break;
    case LatticeTypes::Si:
      if (lattice_list.si_hash.count(target_lattice.id)) {
        lattice_list.si_hash.erase(target_lattice.id);
        lattice_list.si_hash.emplace(source_lattice.id);
      } else assert(false);
      break;
    case LatticeTypes::Mo:
      break;
    default:
      assert(false);
      break;
  }

  source_lattice.type = LatticeTypes{ghost_atom};
  target_lattice.type = LatticeTypes{LatticeTypes::V};

  // calculate system energy after transition.
  double e_after = 0;
  {
    bonds::_type_pair_ia e_src = bonds::BondsCounter::count(&lattice_list, source_id, source_lattice.type);
    bonds::_type_pair_ia e_des = bonds::BondsCounter::count_exchange(&lattice_list, target_id, source_id, source_lattice.type, target_lattice.type);
    e_after = e_src + e_des;
  }

  // exchange atoms back.
  if (lattice_list.vac_hash.count(target_lattice.id)) {
    lattice_list.vac_hash.erase(target_lattice.id);
    // lattice_list.vac_hash.emplace(std::make_pair(source_lattice.id, VacancyHash{}));
  } else assert(false);
  switch (ghost_atom)
  {
    case LatticeTypes::V:
      assert(false);
      break;
    case LatticeTypes::Re:
      if (lattice_list.re_hash.count(source_lattice.id)) {
        lattice_list.re_hash.erase(source_lattice.id);
        lattice_list.re_hash.emplace(target_lattice.id);
      } else assert(false);
      break;
    case LatticeTypes::Mn:
      if (lattice_list.mn_hash.count(source_lattice.id)) {
        lattice_list.mn_hash.erase(source_lattice.id);
        lattice_list.mn_hash.emplace(target_lattice.id);
      } else assert(false);
      break;
    case LatticeTypes::Ni:
      if (lattice_list.ni_hash.count(source_lattice.id)) {
        lattice_list.ni_hash.erase(source_lattice.id);
        lattice_list.ni_hash.emplace(target_lattice.id);
      } else assert(false);
      break;
    case LatticeTypes::Si:
      if (lattice_list.si_hash.count(source_lattice.id)) {
        lattice_list.si_hash.erase(source_lattice.id);
        lattice_list.si_hash.emplace(target_lattice.id);
      } else assert(false);
      break;
    case LatticeTypes::Mo:
      break;
    default:
      assert(false);
      break;
  }
  //kiwi::logs::v(" ", " e_after is : {} .\n", e_after);

  source_lattice.type = LatticeTypes{LatticeTypes::V};
  target_lattice.type = LatticeTypes{ghost_atom};

  return e_after - e_before;
#endif
  return 0;
}
