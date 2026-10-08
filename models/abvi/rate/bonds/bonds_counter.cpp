//
// Created by genshen on 2019-01-20.
//

#include "bonds_counter.h"
#include <logs/logs.h>
#include <iostream>
/**
 * table 7 in paper {Vincent, E., C. S. Becquart, and C. Domain. \
 * "Solute interaction with point defects in α Fe during thermal ageing: A combined ab initio and atomic kinetic Monte
 * Carlo approach."  \ Journal of nuclear materials 351.1-3 (2006): 88-99}.
 */
const std::map<bonds::PairBond::bond_type, bonds::_type_pair_ia> bonds::BondsCounter::_1nn_bonds = {
    {PairBond::VV, 0.199}, {PairBond::VMo, -0.755318}, {PairBond::VRe, -0.920318},
    {PairBond::VMn, -0.038}, {PairBond::VNi, -0.708818}, {PairBond::VSi, -0.344},
    
    {PairBond::MoMo, -1.991}, {PairBond::MoRe, -2.096}, {PairBond::MoMn, -0.446},
    {PairBond::MoNi, -1.8045}, {PairBond::MoSi, -0.683},
    
    {PairBond::ReRe, -2.201}, {PairBond::ReMn, -0.365},
    {PairBond::ReNi, -1.9095}, {PairBond::ReSi, -0.611},

    {PairBond::MnMn, -0.271}, {PairBond::MnNi, -0.366},
    {PairBond::MnSi, -0.488},

    {PairBond::NiNi, -1.618}, {PairBond::NiSi, -0.723},

    {PairBond::SiSi, -0.716},
};

const std::map<bonds::PairBond::bond_type, bonds::_type_pair_ia> bonds::BondsCounter::_2nn_bonds = {
    {PairBond::VV, -0.245}, {PairBond::VMo, -0.377659}, {PairBond::VRe, -0.460159},
    {PairBond::VMn, -0.203}, {PairBond::VNi, -0.354409}, {PairBond::VSi, -0.248},
    
    {PairBond::MoMo, -0.9955}, {PairBond::MoRe, -1.048}, {PairBond::MoMn, -0.631},
    {PairBond::MoNi, -0.90225}, {PairBond::MoSi, -0.469},

    {PairBond::ReRe, -1.1005}, {PairBond::ReMn, -0.621},
    {PairBond::ReNi, -0.95475}, {PairBond::ReSi, -0.566},

    {PairBond::MnMn, -0.611}, {PairBond::MnNi, -0.496},
    {PairBond::MnSi, -0.316},

    {PairBond::NiNi, -0.809}, {PairBond::NiSi, -0.521},

    {PairBond::SiSi, -0.389},
};

bonds::_type_pair_ia bonds::BondsCounter::count(LatticesList *lat_list, _type_lattice_id source_id,
                                                LatticeTypes src_atom_type) {
  Lattice _1nn_neighbour[LatticesList::MAX_1NN]; // todo new array many times.
  _type_neighbour_status _1nn_status = lat_list->get1nnStatus(source_id);
  _type_lattice_size x = source_id % lat_list->meta.size_x;
  _type_lattice_size y = (source_id / lat_list->meta.size_x) % lat_list->meta.size_y;
  _type_lattice_size z = source_id / (lat_list->meta.size_x * lat_list->meta.size_y);
  //lat_list->get1nn2(source_id, _1nn_neighbour);
  lat_list->get1nn2(x, y, z, _1nn_neighbour);
  // Traver all 1nn neighbour lattices, and calculate bond energy contribution.
  _type_pair_ia energy = 0.0;
  for (int b = 0; b < LatticesList::MAX_NEI_BITS; b++) {
    if ((_1nn_status >> b) & 0x1) {
      // we assume that src_atom_type is single atom or vacancy.
      if (_1nn_neighbour[b].type.isAtom() || _1nn_neighbour[b].type.isVacancy()) {
        // it is vacancy or single atom.
        const PairBond::bond_type bond = PairBond::makeBond(_1nn_neighbour[b].type, src_atom_type);
        energy += _1nn_bonds.at(bond);
      }
    }
  }
  //kiwi::logs::v(" ", " energy 1 is : {} .\n", energy);
  // Traver all 2nn neighbour lattices, and calculate bond energy contribution.
  Lattice _2nn_neighbour[LatticesList::MAX_2NN]; // todo new array many times.
  _type_neighbour_status _2nn_status = lat_list->get2nnStatus(source_id);
  // x = source_id % lat_list->meta.size_x;
  // y = (source_id / lat_list->meta.size_x) % box->lattice_list->meta.size_y;
  // z = source_id / (lat_list->meta.size_x * box->lattice_list->meta.size_y);
  lat_list->get2nn2(x, y, z, _2nn_neighbour);
  // kiwi::logs::v(" ", " x is : {} y is : {} z is : {} .\n", x, y, z);
  // for(int a = 0; a < 6; a++){
  //   kiwi::logs::v(" ", " {} is : {} .\n", a, _2nn_neighbour[a].type._type);
  // }
  // lat_list->get2nn(source_id, _2nn_neighbour);
  for (int b = 0; b < LatticesList::MAX_2NN; b++) {
    if ((_2nn_status >> b) & 0x1) {
      // we assume that src_atom_type is single atom or vacancy.
      if (_2nn_neighbour[b].type.isAtom() || _2nn_neighbour[b].type.isVacancy()) {
        // it is vacancy or single atom.
        const PairBond::bond_type bond = PairBond::makeBond(_2nn_neighbour[b].type, src_atom_type);
        energy += _2nn_bonds.at(bond);
      }
    }
  }
  //kiwi::logs::v(" ", " energy 2 is : {} .\n", energy);
  return energy;
}

bonds::_type_pair_ia bonds::BondsCounter::count_exchange(LatticesList *lat_list, _type_lattice_id now_id, _type_lattice_id source_id,
                                                         LatticeTypes atom, LatticeTypes src_atom_type) {
  Lattice _1nn_neighbour[LatticesList::MAX_1NN]; // todo new array many times.
  _type_neighbour_status _1nn_status = lat_list->get1nnStatus(now_id);
  _type_lattice_size x = now_id % lat_list->meta.size_x;
  _type_lattice_size y = (now_id / lat_list->meta.size_x) % lat_list->meta.size_y;
  _type_lattice_size z = now_id / (lat_list->meta.size_x * lat_list->meta.size_y);
  //lat_list->get1nn2(source_id, _1nn_neighbour);
  lat_list->get1nn3(x, y, z, _1nn_neighbour, source_id, atom);
  // Traver all 1nn neighbour lattices, and calculate bond energy contribution.
  _type_pair_ia energy = 0.0;
  for (int b = 0; b < LatticesList::MAX_NEI_BITS; b++) {
    if ((_1nn_status >> b) & 0x1) {
      // we assume that src_atom_type is single atom or vacancy.
      if (_1nn_neighbour[b].type.isAtom() || _1nn_neighbour[b].type.isVacancy()) {
        // it is vacancy or single atom.
        const PairBond::bond_type bond = PairBond::makeBond(_1nn_neighbour[b].type, src_atom_type);
        energy += _1nn_bonds.at(bond);
      }
    }
  }
  //kiwi::logs::v(" ", " energy 1 is : {} .\n", energy);
  // Traver all 2nn neighbour lattices, and calculate bond energy contribution.
  Lattice _2nn_neighbour[LatticesList::MAX_2NN]; // todo new array many times.
  _type_neighbour_status _2nn_status = lat_list->get2nnStatus(now_id);
  // x = source_id % lat_list->meta.size_x;
  // y = (source_id / lat_list->meta.size_x) % box->lattice_list->meta.size_y;
  // z = source_id / (lat_list->meta.size_x * box->lattice_list->meta.size_y);
  lat_list->get2nn3(x, y, z, _2nn_neighbour, source_id, atom);
  // kiwi::logs::v(" ", " x is : {} y is : {} z is : {} .\n", x, y, z);
  // for(int a = 0; a < 6; a++){
  //   kiwi::logs::v(" ", " {} is : {} .\n", a, _2nn_neighbour[a].type._type);
  // }
  // lat_list->get2nn(source_id, _2nn_neighbour);
  for (int b = 0; b < LatticesList::MAX_2NN; b++) {
    if ((_2nn_status >> b) & 0x1) {
      // we assume that src_atom_type is single atom or vacancy.
      if (_2nn_neighbour[b].type.isAtom() || _2nn_neighbour[b].type.isVacancy()) {
        // it is vacancy or single atom.
        const PairBond::bond_type bond = PairBond::makeBond(_2nn_neighbour[b].type, src_atom_type);
        energy += _2nn_bonds.at(bond);
      }
    }
  }
  //kiwi::logs::v(" ", " energy 2 is : {} .\n", energy);
  return energy;
}