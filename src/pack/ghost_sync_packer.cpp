//
// Created by genshen on 2019/10/14.
//

#include "ghost_sync_packer.h"
#include "utils/macros.h"
#include <iostream>
#include <omp.h>

GhostSyncPacker::GhostSyncPacker(LatticesList *lats_list) : lats(lats_list) {}

const unsigned long GhostSyncPacker::sendLength(const std::vector<comm::Region<pack_region_type>> send_regions,
                                                const int dimension, const int direction) {
  unsigned long size_ghost = 0;
  for (auto r : send_regions) {
    size_ghost += r.volume();
  }
  return BCC_DBX * size_ghost;
}

const int GhostSyncPacker::sendLength2(const int dimension,
                                       std::unordered_set<_type_lattice_id>& now_exchange_surface_x, 
                                       std::unordered_set<_type_lattice_id>& now_exchange_surface_y,
                                       std::unordered_set<_type_lattice_id>& now_exchange_surface_z) {
  int size_send = 0;
  // for (auto r : send_regions) {
  //   // 返回体积 既发送总量
  //   size_send += r.volume();
  // }
  // return BCC_DBX * size_send;
  switch (dimension)
  {
    case 0:
      size_send = now_exchange_surface_x.size();
      break;
    case 1:
      size_send = now_exchange_surface_y.size();
      break;
    case 2:
      size_send = now_exchange_surface_z.size();
      break;
    default:
      assert(false);
      break;
  }
  return size_send;
}

void GhostSyncPacker::onSend(pack_date_type *buffer, const std::vector<comm::Region<pack_region_type>> send_regions,
                             const unsigned long send_len, const int dimension, const int direction) {
  unsigned long len = 0;
  for (auto &r : send_regions) {
    for (int z = r.z_low; z < r.z_high; z++) {
      for (int y = r.y_low; y < r.y_high; y++) {
        for (int x = r.x_low; x < r.x_high; x++) {
          _type_lattice_id latti_id = lats->getId(BCC_DBX * x, y, z);
          buffer[len].id = latti_id;

          if (lats->vac_hash.count(latti_id)) {
            buffer[len].type = LatticeTypes{LatticeTypes::V};
          } else if (lats->re_hash.count(latti_id)) {
            buffer[len].type = LatticeTypes{LatticeTypes::Re};
          } else if (lats->mn_hash.count(latti_id)) {
            buffer[len].type = LatticeTypes{LatticeTypes::Mn};
          } else if (lats->ni_hash.count(latti_id)) {
            buffer[len].type = LatticeTypes{LatticeTypes::Ni};
          } else if (lats->si_hash.count(latti_id)) {
            buffer[len].type = LatticeTypes{LatticeTypes::Si};
          } else if (lats->momo_hash.count(latti_id)) {
            buffer[len].type = LatticeTypes{LatticeTypes::MoMo};
          } else if (lats->more_hash.count(latti_id)) {
            buffer[len].type = LatticeTypes{LatticeTypes::MoRe};
          } else if (lats->rere_hash.count(latti_id)) {
            buffer[len].type = LatticeTypes{LatticeTypes::ReRe};
          } else {
            buffer[len].type = LatticeTypes{LatticeTypes::Mo};
          }
          len++;

          _type_lattice_id latti_id2 = lats->getId(BCC_DBX * x + 1, y, z);
          buffer[len].id = latti_id2;
          if (lats->vac_hash.count(latti_id2)) {
            buffer[len].type = LatticeTypes{LatticeTypes::V};
          } else if (lats->re_hash.count(latti_id2)) {
            buffer[len].type = LatticeTypes{LatticeTypes::Re};
          } else if (lats->mn_hash.count(latti_id2)) {
            buffer[len].type = LatticeTypes{LatticeTypes::Mn};
          } else if (lats->ni_hash.count(latti_id2)) {
            buffer[len].type = LatticeTypes{LatticeTypes::Ni};
          } else if (lats->si_hash.count(latti_id2)) {
            buffer[len].type = LatticeTypes{LatticeTypes::Si};
          } else if (lats->momo_hash.count(latti_id2)) {
            buffer[len].type = LatticeTypes{LatticeTypes::MoMo};
          } else if (lats->more_hash.count(latti_id2)) {
            buffer[len].type = LatticeTypes{LatticeTypes::MoRe};
          } else if (lats->rere_hash.count(latti_id2)) {
            buffer[len].type = LatticeTypes{LatticeTypes::ReRe};
          } else {
            buffer[len].type = LatticeTypes{LatticeTypes::Mo};
          }
          // buffer[len++] = lats->_lattices[z][y][BCC_DBX * x];
          // buffer[len++] = lats->_lattices[z][y][BCC_DBX * x + 1];
          len++;
        }
      }
    }
  }
}

void GhostSyncPacker::onSend2(ChangeLattice *buffer, const int send_len, const int dimension,
                              std::unordered_set<_type_lattice_id>& now_exchange_surface_x, 
                              std::unordered_set<_type_lattice_id>& now_exchange_surface_y,
                              std::unordered_set<_type_lattice_id>& now_exchange_surface_z) {

  int len = 0;

  if(dimension == 0) {
    for (const auto& pair : now_exchange_surface_x) {
      buffer[len].x = pair % lats->meta.size_x;
      buffer[len].y = (pair / lats->meta.size_x) % lats->meta.size_y;
      buffer[len].z = pair / (lats->meta.size_x * lats->meta.size_y);
      buffer[len].type._type = lats->getType(pair);
      len++;
    }
  } else if(dimension == 1) {
    for (const auto& pair : now_exchange_surface_y) {
      buffer[len].x = pair % lats->meta.size_x;
      buffer[len].y = (pair / lats->meta.size_x) % lats->meta.size_y;
      buffer[len].z = pair / (lats->meta.size_x * lats->meta.size_y);
      buffer[len].type._type = lats->getType(pair);
      len++;
    }
  }else {
    for (const auto& pair : now_exchange_surface_z) {
      buffer[len].x = pair % lats->meta.size_x;
      buffer[len].y = (pair / lats->meta.size_x) % lats->meta.size_y;
      buffer[len].z = pair / (lats->meta.size_x * lats->meta.size_y);
      buffer[len].type._type = lats->getType(pair);
      len++;
    }
  }
  assert(len == send_len);
}

void GhostSyncPacker::onReceive(pack_date_type *buffer, const std::vector<comm::Region<pack_region_type>> recv_regions,
                                const unsigned long receive_len, const int dimension, const int direction) {
  unsigned long len = 0;
  for (auto &r : recv_regions) {
    for (int z = r.z_low; z < r.z_high; z++) {
      for (int y = r.y_low; y < r.y_high; y++) {
        for (int x = r.x_low; x < r.x_high; x++) {
          _type_lattice_id latti_id = lats->getId(BCC_DBX * x, y, z);
          if (lats->vac_hash.count(latti_id)) {
            lats->vac_hash.erase(latti_id);
          } else if (lats->re_hash.count(latti_id)) {
            lats->re_hash.erase(latti_id);
          } else if (lats->mn_hash.count(latti_id)) {
            lats->mn_hash.erase(latti_id);
          } else if (lats->ni_hash.count(latti_id)) {
            lats->ni_hash.erase(latti_id);
          } else if (lats->si_hash.count(latti_id)) {
            lats->si_hash.erase(latti_id);
          } else if (lats->momo_hash.count(latti_id)) {
            lats->momo_hash.erase(latti_id);
          } else if (lats->more_hash.count(latti_id)) {
            lats->more_hash.erase(latti_id);
          } else if (lats->rere_hash.count(latti_id)) {
            lats->rere_hash.erase(latti_id);
          }

          switch (buffer[len].type._type)
          {
            case LatticeTypes::V:
              lats->vac_hash.emplace(std::make_pair(latti_id, VacancyHash{}));
              break;
            case LatticeTypes::Re:
              lats->re_hash.emplace(latti_id);
              break;
            case LatticeTypes::Mn:
              lats->mn_hash.emplace(latti_id);
              break;
            case LatticeTypes::Ni:
              lats->ni_hash.emplace(latti_id);
              break;
            case LatticeTypes::Si:
              lats->si_hash.emplace(latti_id);
              break;
            case LatticeTypes::MoMo:
              lats->momo_hash.emplace(latti_id);
              break;
            case LatticeTypes::MoRe:
              lats->more_hash.emplace(latti_id);
              break;
            case LatticeTypes::ReRe:
              lats->rere_hash.emplace(latti_id);
              break;
            case LatticeTypes::Mo:
              // 不做操作
              break;
            default:
              break;
          }
          len++;

          _type_lattice_id latti_id2 = lats->getId(BCC_DBX * x + 1, y, z);
          if (lats->vac_hash.count(latti_id2)) {
            lats->vac_hash.erase(latti_id2);
          } else if (lats->re_hash.count(latti_id2)) {
            lats->re_hash.erase(latti_id2);
          } else if (lats->mn_hash.count(latti_id2)) {
            lats->mn_hash.erase(latti_id2);
          } else if (lats->ni_hash.count(latti_id2)) {
            lats->ni_hash.erase(latti_id2);
          } else if (lats->si_hash.count(latti_id2)) {
            lats->si_hash.erase(latti_id2);
          } else if (lats->momo_hash.count(latti_id2)) {
            lats->momo_hash.erase(latti_id2);
          } else if (lats->more_hash.count(latti_id2)) {
            lats->more_hash.erase(latti_id2);
          } else if (lats->rere_hash.count(latti_id2)) {
            lats->rere_hash.erase(latti_id2);
          }

          switch (buffer[len].type._type)
          {
            case LatticeTypes::V:
              lats->vac_hash.emplace(std::make_pair(latti_id2, VacancyHash{}));
              break;
            case LatticeTypes::Re:
              lats->re_hash.emplace(latti_id2);
              break;
            case LatticeTypes::Mn:
              lats->mn_hash.emplace(latti_id2);
              break;
            case LatticeTypes::Ni:
              lats->ni_hash.emplace(latti_id2);
              break;
            case LatticeTypes::Si:
              lats->si_hash.emplace(latti_id2);
              break;
            case LatticeTypes::MoMo:
              lats->momo_hash.emplace(latti_id2);
              break;
            case LatticeTypes::MoRe:
              lats->more_hash.emplace(latti_id2);
              break;
            case LatticeTypes::ReRe:
              lats->rere_hash.emplace(latti_id2);
              break;
            case LatticeTypes::Mo:
              // 不做操作
              break;
            default:
              break;
          }
          // lats->_lattices[z][y][BCC_DBX * x].type = buffer[len++].type;
          // lats->_lattices[z][y][BCC_DBX * x + 1].type = buffer[len++].type;
          len++;
        }
      }
    }
  }
}


void GhostSyncPacker::onReceive2(ChangeLattice *buffer, const int receive_len, const int dimension,
                                 std::array<std::unordered_set<_type_lattice_id>, 8>& exchange_surface_x,
                                 std::array<std::unordered_set<_type_lattice_id>, 8>& exchange_surface_y,
                                 std::array<std::unordered_set<_type_lattice_id>, 8>& exchange_surface_z,
                                 const _type_lattice_count *sub_box_lattice_size, 
                                 const _type_lattice_count *neighbour_local_sub_box, unsigned int next_id,
                                 const comm::ColoredDomain *p_domain) {
  
  int len = 0;

  _type_lattice_id x, y, z;
  _type_lattice_id latti_id;
  if(dimension == 2){
    for(len = 0; len < receive_len; len++){
      x = buffer[len].x;
      y = buffer[len].y;
      if(next_id > 3){
        buffer[len].z += sub_box_lattice_size[dimension];
      } else {
        buffer[len].z -= sub_box_lattice_size[dimension] + (neighbour_local_sub_box[dimension] - sub_box_lattice_size[dimension]);
      }
      z = buffer[len].z;
      latti_id = lats->getId(x, y, z);
      if (lats->vac_hash.count(latti_id)) {
        lats->vac_hash.erase(latti_id);
      } else if (lats->re_hash.count(latti_id)) {
        lats->re_hash.erase(latti_id);
      } else if (lats->mn_hash.count(latti_id)) {
        lats->mn_hash.erase(latti_id);
      } else if (lats->ni_hash.count(latti_id)) {
        lats->ni_hash.erase(latti_id);
      } else if (lats->si_hash.count(latti_id)) {
        lats->si_hash.erase(latti_id);
      } else if (lats->momo_hash.count(latti_id)) {
        lats->momo_hash.erase(latti_id);
      } else if (lats->more_hash.count(latti_id)) {
        lats->more_hash.erase(latti_id);
      } else if (lats->rere_hash.count(latti_id)) {
        lats->rere_hash.erase(latti_id);
      }

      switch (buffer[len].type._type)
      {
        case LatticeTypes::V:
          lats->vac_hash.emplace(std::make_pair(latti_id, VacancyHash{}));
          break;
        case LatticeTypes::Re:
          lats->re_hash.emplace(latti_id);
          break;
        case LatticeTypes::Mn:
          lats->mn_hash.emplace(latti_id);
          break;
        case LatticeTypes::Ni:
          lats->ni_hash.emplace(latti_id);
          break;
        case LatticeTypes::Si:
          lats->si_hash.emplace(latti_id);
          break;
        case LatticeTypes::MoMo:
          lats->momo_hash.emplace(latti_id);
          break;
        case LatticeTypes::MoRe:
          lats->more_hash.emplace(latti_id);
          break;
        case LatticeTypes::ReRe:
          lats->rere_hash.emplace(latti_id);
          break;
        case LatticeTypes::Mo:
          // 不做操作
          break;
        default:
          break;
      }
    }
  } else if(dimension == 1) {
    for(len = 0; len < receive_len; len++) {
      x = buffer[len].x;
      if(next_id == 2 || next_id == 3 || next_id == 6 || next_id == 7){
        buffer[len].y += sub_box_lattice_size[dimension];
      }else{
        buffer[len].y -= sub_box_lattice_size[dimension] + (neighbour_local_sub_box[dimension] - sub_box_lattice_size[dimension]);
      }
      y = buffer[len].y;
      z = buffer[len].z;
      latti_id = lats->getId(x, y, z);
      if (lats->vac_hash.count(latti_id)) {
        lats->vac_hash.erase(latti_id);
      } else if (lats->re_hash.count(latti_id)) {
        lats->re_hash.erase(latti_id);
      } else if (lats->mn_hash.count(latti_id)) {
        lats->mn_hash.erase(latti_id);
      } else if (lats->ni_hash.count(latti_id)) {
        lats->ni_hash.erase(latti_id);
      } else if (lats->si_hash.count(latti_id)) {
        lats->si_hash.erase(latti_id);
      } else if (lats->momo_hash.count(latti_id)) {
        lats->momo_hash.erase(latti_id);
      } else if (lats->more_hash.count(latti_id)) {
        lats->more_hash.erase(latti_id);
      } else if (lats->rere_hash.count(latti_id)) {
        lats->rere_hash.erase(latti_id);
      }

      switch (buffer[len].type._type)
      {
        case LatticeTypes::V:
          lats->vac_hash.emplace(std::make_pair(latti_id, VacancyHash{}));
          break;
        case LatticeTypes::Re:
          lats->re_hash.emplace(latti_id);
          break;
        case LatticeTypes::Mn:
          lats->mn_hash.emplace(latti_id);
          break;
        case LatticeTypes::Ni:
          lats->ni_hash.emplace(latti_id);
          break;
        case LatticeTypes::Si:
          lats->si_hash.emplace(latti_id);
          break;
        case LatticeTypes::MoMo:
          lats->momo_hash.emplace(latti_id);
          break;
        case LatticeTypes::MoRe:
          lats->more_hash.emplace(latti_id);
          break;
        case LatticeTypes::ReRe:
          lats->rere_hash.emplace(latti_id);
          break;
        case LatticeTypes::Mo:
          // 不做操作
          break;
        default:
          break;
      }
    }
  } else if(dimension == 0) {
    for(len = 0; len < receive_len; len++) {
      if(next_id % 2 == 1) {
        buffer[len].x += BCC_DBX * (sub_box_lattice_size[dimension]);
      } else {
        buffer[len].x -= BCC_DBX * (sub_box_lattice_size[dimension] + (neighbour_local_sub_box[dimension] - sub_box_lattice_size[dimension]));
      }
      x = buffer[len].x;
      y = buffer[len].y;
      z = buffer[len].z;
      latti_id = lats->getId(x, y, z);
      if (lats->vac_hash.count(latti_id)) {
        lats->vac_hash.erase(latti_id);
      } else if (lats->re_hash.count(latti_id)) {
        lats->re_hash.erase(latti_id);
      } else if (lats->mn_hash.count(latti_id)) {
        lats->mn_hash.erase(latti_id);
      } else if (lats->ni_hash.count(latti_id)) {
        lats->ni_hash.erase(latti_id);
      } else if (lats->si_hash.count(latti_id)) {
        lats->si_hash.erase(latti_id);
      } else if (lats->momo_hash.count(latti_id)) {
        lats->momo_hash.erase(latti_id);
      } else if (lats->more_hash.count(latti_id)) {
        lats->more_hash.erase(latti_id);
      } else if (lats->rere_hash.count(latti_id)) {
        lats->rere_hash.erase(latti_id);
      }

      switch (buffer[len].type._type)
      {
        case LatticeTypes::V:
          lats->vac_hash.emplace(std::make_pair(latti_id, VacancyHash{}));
          break;
        case LatticeTypes::Re:
          lats->re_hash.emplace(latti_id);
          break;
        case LatticeTypes::Mn:
          lats->mn_hash.emplace(latti_id);
          break;
        case LatticeTypes::Ni:
          lats->ni_hash.emplace(latti_id);
          break;
        case LatticeTypes::Si:
          lats->si_hash.emplace(latti_id);
          break;
        case LatticeTypes::MoMo:
          lats->momo_hash.emplace(latti_id);
          break;
        case LatticeTypes::MoRe:
          lats->more_hash.emplace(latti_id);
          break;
        case LatticeTypes::ReRe:
          lats->rere_hash.emplace(latti_id);
          break;
        case LatticeTypes::Mo:
          // 不做操作
          break;
        default:
          break;
      }
    }
  } else {
    assert(false);
  }

  // 存储转发
  const int dims[comm::DIMENSION_SIZE] = {comm::DIM_X, comm::DIM_Y, comm::DIM_Z};
  std::array<std::vector<comm::Region<comm::_type_lattice_coord>>, comm::DIMENSION_SIZE> send_regions; // send regions in each dimension
  for (int sect = 0; sect < 8; sect++) {
    for (int d = 1; d < comm::DIMENSION_SIZE; d++) {
      send_regions[d] = comm::fwCommSectorSendRegion(sect, dims[d], p_domain->lattice_size_ghost,
                                                     p_domain->local_split_coord, p_domain->local_sub_box_lattice_region);
      for (auto &r : send_regions[d]) {
        _type_lattice_size x_low = BCC_DBX * r.x_low;
        _type_lattice_size y_low = r.y_low;
        _type_lattice_size z_low = r.z_low;
        _type_lattice_size x_high = BCC_DBX * r.x_high;
        _type_lattice_size y_high = r.y_high;
        _type_lattice_size z_high = r.z_high;

        for(len = 0; len < receive_len; len++) {
          x = buffer[len].x;
          y = buffer[len].y;
          z = buffer[len].z;
          if(x_low <= x && x < x_high && y_low <= y && y < y_high && z_low <= z && z < z_high){
            _type_lattice_id latti_id = lats->getId(x, y, z);
            if(d == 0) {
              auto it_from = exchange_surface_x[sect].find(latti_id);
              if(it_from == exchange_surface_x[sect].end()) exchange_surface_x[sect].emplace(latti_id);
            } else if(d == 1) {
              auto it_from = exchange_surface_y[sect].find(latti_id);
              if(it_from == exchange_surface_y[sect].end()) exchange_surface_y[sect].emplace(latti_id);
            } else {
              auto it_from = exchange_surface_z[sect].find(latti_id);
              if(it_from == exchange_surface_z[sect].end()) exchange_surface_z[sect].emplace(latti_id);
            }
          }
        }
      }
    }
  }
 
}