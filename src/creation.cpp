//
// Created by genshen on 2019/10/6.
//

#include "creation.h"
#include "type_define.h"
#include "utils/macros.h"
#include "utils/random/random.h"
#include "utils/random/rng_type.h"
#include "utils/simulation_domain.h"
#include <cassert>
#include <fstream>
#include <utils/bundle.h>
#include <utils/mpi_utils.h>
#include <iostream>
#include <logs/logs.h>

void creation::createAtomsRandom(int64_t seed_create_types, LatticesList *lats,
                                 const std::vector<LatticeTypes::lat_type> types,
                                 const std::vector<int64_t> types_ratio, const comm::ColoredDomain *p_domain) {
  int64_t ratio_total = 0;
  for (int64_t i = 0; i < types.size(); i++) {
    ratio_total += types_ratio[i];
  }
  std::uniform_int_distribution<> types_dis(1, ratio_total);
  r::type_rng types_rng(seed_create_types);
  // todo skip ghost lattices.
  int64_t n = 0;

  // for (_type_lattice_size z = 0; z < lats->meta.size_z; z++) {
  //   for (_type_lattice_size y = 0; y < lats->meta.size_y; y++) {
  //     for (_type_lattice_size x = 0; x < lats->meta.size_x; x++) {
  //       assert(n == x + lats->meta.size_x * (y + z * lats->meta.size_y));
  //       const int64_t rand_hit = types_dis(types_rng);
  //       LatticeTypes::lat_type latti_type = LatticeTypes::randomAtomsType(types.data(), types_ratio.data(), types.size(), rand_hit);
  //       // lats->_lattices[z][y][x].type._type = latti_type;
  //       switch (latti_type)
  //       {
  //         case LatticeTypes::Cu:
  //           lats->cu_hash.emplace(n);
  //           break;
  //         case LatticeTypes::Ni:
  //           lats->ni_hash.emplace(n);
  //           break;
  //         case LatticeTypes::Mn:
  //           lats->mn_hash.emplace(n);
  //           break;
  //         case LatticeTypes::Fe:
  //           // 不存储Fe
  //           break;
  //         default:
  //           break;
  //       }
  //       n++;
  //     }
  //   }
  // }
  // unsigned long cu_count = 0;
  // unsigned long ni_count = 0;
  // unsigned long mn_count = 0;
  // kiwi::logs::v(" ", " Number of elements in cu_hash is : {} .\n", lats->cu_hash.size());
  // kiwi::logs::v(" ", " Number of elements in ni_hash is : {} .\n", lats->ni_hash.size());
  // kiwi::logs::v(" ", " Number of elements in mn_hash is : {} .\n", lats->mn_hash.size());
  // for (_type_lattice_size z = 0; z < lats->meta.size_z; z++) {
  //   for (_type_lattice_size y = 0; y < lats->meta.size_y; y++) {
  //     for (_type_lattice_size x = 0; x < lats->meta.size_x; x++) {
  //       switch (lats->_lattices[z][y][x].type._type)
  //       {
  //       case LatticeTypes::Cu:
  //         cu_count++;
  //         break;
  //       case LatticeTypes::Ni:
  //         ni_count++;
  //         break;
  //       case LatticeTypes::Mn:
  //         mn_count++;
  //         break;
  //       default:
  //         break;
  //       }
  //     }
  //   }
  // }
  // kiwi::logs::v(" ", " Number of cu elements in latti is : {} .\n", cu_count);
  // kiwi::logs::v(" ", " Number of ni elements in latti is : {} .\n", ni_count);
  // kiwi::logs::v(" ", " Number of mn elements in latti is : {} .\n", mn_count);
  // lats->forAllLattices(
  //     [&](const _type_lattice_coord x, const _type_lattice_coord y, const _type_lattice_coord z, Lattice &lattice) {
  //       // set lattice types randomly.
  //       const unsigned int rand_hit = types_dis(types_rng);
  //       lattice.type._type = LatticeTypes::randomAtomsType(types.data(), types_ratio.data(), types.size(), rand_hit);
  //       return true;
  //     });
}

void creation::createRandom(int64_t seed_create_types, int64_t seed_create_vacancy, LatticesList *lats,
                            VacancyList *va_list, const std::vector<LatticeTypes::lat_type> types,
                            const std::vector<int64_t> types_ratio, const int64_t va_count,
                            const comm::ColoredDomain *p_domain) {
  _type_lattice_count va_local =
      va_count / SimulationDomain::comm_sim_pro.all_ranks +
      (SimulationDomain::comm_sim_pro.own_rank < (va_count % SimulationDomain::comm_sim_pro.all_ranks) ? 1 : 0);
  _type_lattice_count re_local =
      types_ratio[1] / SimulationDomain::comm_sim_pro.all_ranks +
      (SimulationDomain::comm_sim_pro.own_rank < (types_ratio[1] % SimulationDomain::comm_sim_pro.all_ranks) ? 1 : 0);
  _type_lattice_count mn_local =
      types_ratio[2] / SimulationDomain::comm_sim_pro.all_ranks +
      (SimulationDomain::comm_sim_pro.own_rank < (types_ratio[2] % SimulationDomain::comm_sim_pro.all_ranks) ? 1 : 0);
  _type_lattice_count ni_local =
      types_ratio[3] / SimulationDomain::comm_sim_pro.all_ranks +
      (SimulationDomain::comm_sim_pro.own_rank < (types_ratio[3] % SimulationDomain::comm_sim_pro.all_ranks) ? 1 : 0);
  _type_lattice_count si_local =
      types_ratio[4] / SimulationDomain::comm_sim_pro.all_ranks +
      (SimulationDomain::comm_sim_pro.own_rank < (types_ratio[4] % SimulationDomain::comm_sim_pro.all_ranks) ? 1 : 0);
  _type_lattice_count momo_local =
      types_ratio[5] / SimulationDomain::comm_sim_pro.all_ranks +
      (SimulationDomain::comm_sim_pro.own_rank < (types_ratio[5] % SimulationDomain::comm_sim_pro.all_ranks) ? 1 : 0);
  _type_lattice_count more_local =
      types_ratio[6] / SimulationDomain::comm_sim_pro.all_ranks +
      (SimulationDomain::comm_sim_pro.own_rank < (types_ratio[6] % SimulationDomain::comm_sim_pro.all_ranks) ? 1 : 0);
  _type_lattice_count rere_local =
      types_ratio[10] / SimulationDomain::comm_sim_pro.all_ranks +
      (SimulationDomain::comm_sim_pro.own_rank < (types_ratio[10] % SimulationDomain::comm_sim_pro.all_ranks) ? 1 : 0);
  lats->vac_hash.reserve(static_cast<size_t>(va_local * 1.2));
  lats->re_hash.reserve(static_cast<size_t>(re_local * 1.2));
  lats->mn_hash.reserve(static_cast<size_t>(mn_local * 1.2));
  lats->ni_hash.reserve(static_cast<size_t>(ni_local * 1.2));
  lats->si_hash.reserve(static_cast<size_t>(si_local * 1.2));
  lats->momo_hash.reserve(static_cast<size_t>(momo_local * 1.2));
  lats->more_hash.reserve(static_cast<size_t>(more_local * 1.2));
  lats->rere_hash.reserve(static_cast<size_t>(rere_local * 1.2));
  // create lattice types
  // createAtomsRandom(seed_create_types, lats, types, types_ratio, p_domain);

  // create vacancy
  // to use SimulationDomain, make sure this function is called after
  // simulation::createDomain()
  //const int64_t max_lattice_size = lats->meta.box_x * lats->meta.box_y * lats->meta.box_z;
  const int64_t max_lattice_size = (lats->meta.box_x - 2 * lats->meta.ghost_x) * (lats->meta.box_y - 2 * lats->meta.ghost_y) * (lats->meta.box_z - 2 * lats->meta.ghost_z);
  std::uniform_int_distribution<long long> types_dis(0, max_lattice_size - 1);
  r::type_rng types_rng(seed_create_types);
  for(_type_lattice_count re_i = 0; re_i < re_local;) {
    _type_lattice_id local_id = types_dis(types_rng);
    _type_lattice_coord x = local_id % (lats->meta.box_x - 2 * lats->meta.ghost_x);
    local_id = local_id / (lats->meta.box_x - 2 * lats->meta.ghost_x);
    _type_lattice_coord y = local_id % (lats->meta.box_y - 2 * lats->meta.ghost_y);
    _type_lattice_coord z = local_id / (lats->meta.box_y - 2 * lats->meta.ghost_y);
    _type_lattice_id latti_id = lats->getId(x + 2 * lats->meta.ghost_x, y + 2 * lats->meta.ghost_y, z + 2 * lats->meta.ghost_z);

    if(!lats->re_hash.count(latti_id)){
      lats->re_hash.emplace(latti_id);
      re_i++;
    }
  }

  for(_type_lattice_count mn_i = 0; mn_i < mn_local;) {
    _type_lattice_id local_id = types_dis(types_rng);
    _type_lattice_coord x = local_id % (lats->meta.box_x - 2 * lats->meta.ghost_x);
    local_id = local_id / (lats->meta.box_x - 2 * lats->meta.ghost_x);
    _type_lattice_coord y = local_id % (lats->meta.box_y - 2 * lats->meta.ghost_y);
    _type_lattice_coord z = local_id / (lats->meta.box_y - 2 * lats->meta.ghost_y);
    _type_lattice_id latti_id = lats->getId(x + 2 * lats->meta.ghost_x, y + 2 * lats->meta.ghost_y, z + 2 * lats->meta.ghost_z);

    if(!lats->re_hash.count(latti_id) && !lats->mn_hash.count(latti_id)) {
      lats->mn_hash.emplace(latti_id);
      mn_i++;
    }
  }

  for(_type_lattice_count ni_i = 0; ni_i < ni_local;) {
    _type_lattice_id local_id = types_dis(types_rng);
    _type_lattice_coord x = local_id % (lats->meta.box_x - 2 * lats->meta.ghost_x);
    local_id = local_id / (lats->meta.box_x - 2 * lats->meta.ghost_x);
    _type_lattice_coord y = local_id % (lats->meta.box_y - 2 * lats->meta.ghost_y);
    _type_lattice_coord z = local_id / (lats->meta.box_y - 2 * lats->meta.ghost_y);
    _type_lattice_id latti_id = lats->getId(x + 2 * lats->meta.ghost_x, y + 2 * lats->meta.ghost_y, z + 2 * lats->meta.ghost_z);

    if(!lats->re_hash.count(latti_id) && !lats->mn_hash.count(latti_id) && !lats->ni_hash.count(latti_id)) {
      lats->ni_hash.emplace(latti_id);
      ni_i++;
    }
  }

  for(_type_lattice_count si_i = 0; si_i < si_local;) {
    _type_lattice_id local_id = types_dis(types_rng);
    _type_lattice_coord x = local_id % (lats->meta.box_x - 2 * lats->meta.ghost_x);
    local_id = local_id / (lats->meta.box_x - 2 * lats->meta.ghost_x);
    _type_lattice_coord y = local_id % (lats->meta.box_y - 2 * lats->meta.ghost_y);
    _type_lattice_coord z = local_id / (lats->meta.box_y - 2 * lats->meta.ghost_y);
    _type_lattice_id latti_id = lats->getId(x + 2 * lats->meta.ghost_x, y + 2 * lats->meta.ghost_y, z + 2 * lats->meta.ghost_z);

    if(!lats->re_hash.count(latti_id) && !lats->mn_hash.count(latti_id) && !lats->ni_hash.count(latti_id) && !lats->si_hash.count(latti_id)) {
      lats->si_hash.emplace(latti_id);
      si_i++;
    }
  }

  for(_type_lattice_count momo_i = 0; momo_i < momo_local;) {
    _type_lattice_id local_id = types_dis(types_rng);
    _type_lattice_coord x = local_id % (lats->meta.box_x - 2 * lats->meta.ghost_x);
    local_id = local_id / (lats->meta.box_x - 2 * lats->meta.ghost_x);
    _type_lattice_coord y = local_id % (lats->meta.box_y - 2 * lats->meta.ghost_y);
    _type_lattice_coord z = local_id / (lats->meta.box_y - 2 * lats->meta.ghost_y);
    _type_lattice_id latti_id = lats->getId(x + 2 * lats->meta.ghost_x, y + 2 * lats->meta.ghost_y, z + 2 * lats->meta.ghost_z);

    if(!lats->re_hash.count(latti_id) && !lats->mn_hash.count(latti_id) && !lats->ni_hash.count(latti_id) && !lats->si_hash.count(latti_id) && !lats->momo_hash.count(latti_id)) {
      lats->momo_hash.emplace(latti_id);
      momo_i++;
    }
  }

  for(_type_lattice_count more_i = 0; more_i < more_local;) {
    _type_lattice_id local_id = types_dis(types_rng);
    _type_lattice_coord x = local_id % (lats->meta.box_x - 2 * lats->meta.ghost_x);
    local_id = local_id / (lats->meta.box_x - 2 * lats->meta.ghost_x);
    _type_lattice_coord y = local_id % (lats->meta.box_y - 2 * lats->meta.ghost_y);
    _type_lattice_coord z = local_id / (lats->meta.box_y - 2 * lats->meta.ghost_y);
    _type_lattice_id latti_id = lats->getId(x + 2 * lats->meta.ghost_x, y + 2 * lats->meta.ghost_y, z + 2 * lats->meta.ghost_z);

    if(!lats->re_hash.count(latti_id) && !lats->mn_hash.count(latti_id) && !lats->ni_hash.count(latti_id) && !lats->si_hash.count(latti_id) && !lats->momo_hash.count(latti_id) && !lats->more_hash.count(latti_id)) {
      lats->more_hash.emplace(latti_id);
      more_i++;
    }
  }

  for(_type_lattice_count rere_i = 0; rere_i < rere_local;) {
    _type_lattice_id local_id = types_dis(types_rng);
    _type_lattice_coord x = local_id % (lats->meta.box_x - 2 * lats->meta.ghost_x);
    local_id = local_id / (lats->meta.box_x - 2 * lats->meta.ghost_x);
    _type_lattice_coord y = local_id % (lats->meta.box_y - 2 * lats->meta.ghost_y);
    _type_lattice_coord z = local_id / (lats->meta.box_y - 2 * lats->meta.ghost_y);
    _type_lattice_id latti_id = lats->getId(x + 2 * lats->meta.ghost_x, y + 2 * lats->meta.ghost_y, z + 2 * lats->meta.ghost_z);

    if(!lats->re_hash.count(latti_id) && !lats->mn_hash.count(latti_id) && !lats->ni_hash.count(latti_id) && !lats->si_hash.count(latti_id) && !lats->momo_hash.count(latti_id) && !lats->more_hash.count(latti_id) && !lats->rere_hash.count(latti_id)) {
      lats->rere_hash.emplace(latti_id);
      rere_i++;
    }
  }
  // const comm::_type_lattice_size max_lattice_size = lats->meta.box_x * lats->meta.box_y * lats->meta.box_z;
  std::uniform_int_distribution<long long> va_dis(0, max_lattice_size - 1);
  r::type_rng va_rng(seed_create_vacancy);
#ifdef KMC_DEBUG_MODE
  assert(max_lattice_size >= 1);
#endif
  for (_type_lattice_id va_i = 0; va_i < va_local;) {
    _type_lattice_id local_id = va_dis(va_rng);
    _type_lattice_coord x = local_id % (lats->meta.box_x - 2 * lats->meta.ghost_x);
    local_id = local_id / (lats->meta.box_x - 2 * lats->meta.ghost_x);
    _type_lattice_coord y = local_id % (lats->meta.box_y - 2 * lats->meta.ghost_y);
    _type_lattice_coord z = local_id / (lats->meta.box_y - 2 * lats->meta.ghost_y);
    _type_lattice_id latti_id = lats->getId(x + 2 * lats->meta.ghost_x, y + 2 * lats->meta.ghost_y, z + 2 * lats->meta.ghost_z);
    // _type_lattice_coord x = local_id % lats->meta.box_x;
    // local_id = local_id / lats->meta.box_x;
    // _type_lattice_coord y = local_id % lats->meta.box_y;
    // _type_lattice_coord z = local_id / lats->meta.box_y;
    // int64_t long latti_id = lats->getId(x + lats->meta.ghost_x, y + lats->meta.ghost_y, z + lats->meta.ghost_z);
    if (!lats->vac_hash.count(latti_id) && !lats->re_hash.count(latti_id) && !lats->mn_hash.count(latti_id) && !lats->ni_hash.count(latti_id) && !lats->si_hash.count(latti_id) && !lats->momo_hash.count(latti_id) && !lats->more_hash.count(latti_id) && !lats->rere_hash.count(latti_id)) {
      lats->vac_hash.emplace(std::make_pair(latti_id, VacancyHash{}));
      va_i++;
    }
  }
  if(SimulationDomain::comm_sim_pro.own_rank == 0) kiwi::logs::v(" ", " Vac count is : {} Re count is : {} Mn count is : {} Ni count is : {} Si count is : {} MoMo count is : {} MoRe count is : {} ReRe count is : {}.\n", va_local, re_local, mn_local, ni_local, si_local, momo_local, more_local, rere_local);
}

void creation::setGlobalId(LatticesList *lats_list, const comm::Region<comm::_type_lattice_coord> lbr,
                           const comm::Region<comm::_type_lattice_coord> gbr,
                           std::array<int64_t, comm::DIMENSION_SIZE> phase_space) {
  for (comm::_type_lattice_coord z = lbr.z_low; z < lbr.z_high; z++) {
    for (comm::_type_lattice_coord y = lbr.y_low; y < lbr.y_high; y++) {
      for (comm::_type_lattice_coord x = BCC_DBX * lbr.x_low; x < BCC_DBX * lbr.x_high; x++) {
        Lattice &lattice = lats_list->getLat(x, y, z);
        const comm::_type_lattice_coord gx = (x - lats_list->meta.ghost_x) + BCC_DBX * gbr.x_low;
        const comm::_type_lattice_coord gy = (y - lats_list->meta.ghost_y) + gbr.y_low;
        const comm::_type_lattice_coord gz = (z - lats_list->meta.ghost_z) + gbr.z_low;
        // convert xyz to global id.
        const comm::_type_lattice_coord gid = gx + BCC_DBX * phase_space[0] * (gy + gz * phase_space[1]);
        lattice.id = gid;
      }
    }
  }
}

void creation::createFromPile(const std::string pipe_file, int64_t seed_create_types, LatticesList *lats,
                              VacancyList *va_list, const std::vector<LatticeTypes::lat_type> types,
                              const std::vector<int64_t> types_ratio, const comm::ColoredDomain *p_domain) {
  createAtomsRandom(seed_create_types, lats, types, types_ratio, p_domain);

  // open vacancies file
  std::vector<std::array<_type_lattice_coord, 3>> vac_positions{};
  if (kiwi::mpiUtils::global_process.own_rank == MASTER_PROCESSOR) {
    std::ifstream vac_file;
    vac_file.open(pipe_file);
    _type_lattice_coord x, y, z;
    while (vac_file >> x >> y >> z) {
      vac_positions.push_back({x, y, z});
    }
  }

  kiwi::Bundle bundle = kiwi::Bundle();
  bundle.newPackBuffer(1024 * 1024); // make sure it is large enough

  if (kiwi::mpiUtils::global_process.own_rank == MASTER_PROCESSOR) {
    bundle.put(vac_positions);
  }
  MPI_Bcast(bundle.getPackedData(), bundle.getPackedDataCap(), MPI_BYTE, MASTER_PROCESSOR,
            MPI_COMM_WORLD); // synchronize data

  if (kiwi::mpiUtils::global_process.own_rank != MASTER_PROCESSOR) { // unpack
    int cursor = 0;
    bundle.get(cursor, vac_positions);
  }
  bundle.freePackBuffer();

  // place vacancies
  for (auto vac : vac_positions) {
    _type_lattice_coord x = vac[0]; // vac[0] is doubled.
    _type_lattice_coord y = vac[1];
    _type_lattice_coord z = vac[2];
    if (p_domain->local_sub_box_lattice_region.isIn(x / 2, y, z)) {
      Lattice &lat = lats->getLat(x + lats->meta.ghost_x, y + lats->meta.ghost_y, z + lats->meta.ghost_z);
      if (!lat.type.isVacancy()) {
        lat.type._type = LatticeTypes::V;
      }
    }
  }

  // update vacancies list
  va_list->reindex(lats, p_domain->local_sub_box_lattice_region);
}