//
// Created by runchu on 2019-10-16.
//

#include "ghost_init_packer.h"
#include "utils/macros.h"
#include <comm/preset/comm_forwarding_region.h>
#include <iostream>
#include <omp.h>

GhostInitPacker::GhostInitPacker(const comm::ColoredDomain *p_domain, LatticesList *lats_list)
    : p_domain(p_domain), lats(lats_list) {}

const unsigned long GhostInitPacker::sendLength(const int dimension, const int direction) {
  comm::Region<comm::_type_lattice_size> send_region = comm::fwCommLocalSendRegion(
      p_domain->lattice_size_ghost, p_domain->local_sub_box_lattice_region, dimension, direction);
  return BCC_DBX * send_region.volume();
}

int GhostInitPacker::sendLength2(const int dimension, const int direction) {
  comm::Region<comm::_type_lattice_size> send_region = comm::fwCommLocalSendRegion(
      p_domain->lattice_size_ghost, p_domain->local_sub_box_lattice_region, dimension, direction);

  int len = 0;

  for (_type_lattice_id z = send_region.z_low; z < send_region.z_high; z++) {
    for (_type_lattice_id y = send_region.y_low; y < send_region.y_high; y++) {
      for (_type_lattice_id x = send_region.x_low; x < send_region.x_high; x++) {
        _type_lattice_id latti_id = lats->getId(BCC_DBX * x, y, z);
        // lats->isUniqueLattiId(latti_id);
        if (lats->vac_hash.count(latti_id) || lats->re_hash.count(latti_id) || lats->mn_hash.count(latti_id) || 
            lats->ni_hash.count(latti_id) || lats->si_hash.count(latti_id) || lats->momo_hash.count(latti_id) || 
            lats->more_hash.count(latti_id) || lats->rere_hash.count(latti_id)) {
          len++;
        }
        
        _type_lattice_id latti_id2 = lats->getId(BCC_DBX * x + 1, y, z);
        // lats->isUniqueLattiId(latti_id2);
        // if (!(it_vac2 == lats->vac_hash.end() && it_cu2 == lats->cu_hash.end() && it_mn2 == lats->mn_hash.end() && it_ni2 == lats->ni_hash.end())) {
        //   len++;
        // }
        if (lats->vac_hash.count(latti_id2) || lats->re_hash.count(latti_id2) || lats->mn_hash.count(latti_id2) || 
            lats->ni_hash.count(latti_id2) || lats->si_hash.count(latti_id2) || lats->momo_hash.count(latti_id2) || 
            lats->more_hash.count(latti_id2) || lats->rere_hash.count(latti_id2)) {
          len++;
        }
      }
    }
  }

  return len;

}
void GhostInitPacker::onSend(Lattice *buffer, const unsigned long send_len, const int dimension,
              const int direction){
  const comm::Region<comm::_type_lattice_size> send_region = comm::fwCommLocalSendRegion(
      p_domain->lattice_size_ghost, p_domain->local_sub_box_lattice_region, dimension, direction);
  unsigned long len = 0;
  for (int z = send_region.z_low; z < send_region.z_high; z++) {
    for (int y = send_region.y_low; y < send_region.y_high; y++) {
      for (int x = send_region.x_low; x < send_region.x_high; x++) {
        unsigned long latti_id = lats->getId(BCC_DBX * x, y, z);
        buffer[len].id = latti_id;
        buffer[len].type._type = lats->getType(latti_id);
        // if (lats->vac_hash.count(latti_id)) {
        //   buffer[len].type = LatticeTypes{LatticeTypes::V};
        // } else if (lats->re_hash.count(latti_id)) {
        //   buffer[len].type = LatticeTypes{LatticeTypes::Re};
        // } else if (lats->mn_hash.count(latti_id)) {
        //   buffer[len].type = LatticeTypes{LatticeTypes::Mn};
        // } else if (lats->ni_hash.count(latti_id)) {
        //   buffer[len].type = LatticeTypes{LatticeTypes::Ni};
        // } else if (lats->si_hash.count(latti_id)) {
        //   buffer[len].type = LatticeTypes{LatticeTypes::Si};
        // } else if (lats->momo_hash.count(latti_id)) {
        //   buffer[len].type = LatticeTypes{LatticeTypes::MoMo};
        // } else if (lats->more_hash.count(latti_id)) {
        //   buffer[len].type = LatticeTypes{LatticeTypes::MoRe};
        // } else if (lats->rere_hash.count(latti_id)) {
        //   buffer[len].type = LatticeTypes{LatticeTypes::ReRe};
        // } else {
        //   buffer[len].type = LatticeTypes{LatticeTypes::Mo};
        // }
        len++;

        unsigned long latti_id2 = lats->getId(BCC_DBX * x + 1, y, z);
        buffer[len].id = latti_id2;
        buffer[len].type._type = lats->getType(latti_id2);
        // if (lats->vac_hash.count(latti_id2)) {
        //   buffer[len].type = LatticeTypes{LatticeTypes::V};
        // } else if (lats->re_hash.count(latti_id2)) {
        //   buffer[len].type = LatticeTypes{LatticeTypes::Re};
        // } else if (lats->mn_hash.count(latti_id2)) {
        //   buffer[len].type = LatticeTypes{LatticeTypes::Mn};
        // } else if (lats->ni_hash.count(latti_id2)) {
        //   buffer[len].type = LatticeTypes{LatticeTypes::Ni};
        // } else if (lats->si_hash.count(latti_id2)) {
        //   buffer[len].type = LatticeTypes{LatticeTypes::Si};
        // } else if (lats->momo_hash.count(latti_id2)) {
        //   buffer[len].type = LatticeTypes{LatticeTypes::MoMo};
        // } else if (lats->more_hash.count(latti_id2)) {
        //   buffer[len].type = LatticeTypes{LatticeTypes::MoRe};
        // } else if (lats->rere_hash.count(latti_id2)) {
        //   buffer[len].type = LatticeTypes{LatticeTypes::ReRe};
        // } else {
        //   buffer[len].type = LatticeTypes{LatticeTypes::Mo};
        // }
        // buffer[len++] = lats->_lattices[z][y][BCC_DBX * x];
        // buffer[len++] = lats->_lattices[z][y][BCC_DBX * x + 1];
        len++;
      }
    }
  }              
}

void GhostInitPacker::onSend2(Lattice *buffer, int send_len, const int dimension,
                             const int direction) {
  const comm::Region<comm::_type_lattice_size> send_region = comm::fwCommLocalSendRegion(
      p_domain->lattice_size_ghost, p_domain->local_sub_box_lattice_region, dimension, direction);

  int len = 0;
  _type_lattice_count index = 0;

  for (_type_lattice_id z = send_region.z_low; z < send_region.z_high; z++) {
    for (_type_lattice_id y = send_region.y_low; y < send_region.y_high; y++) {
      for (_type_lattice_id x = send_region.x_low; x < send_region.x_high; x++) {
        _type_lattice_id latti_id = lats->getId(BCC_DBX * x, y, z);
        
        // lats->isUniqueLattiId(latti_id);
        // 这里的 id 存的是 本次传输的数组下标，而不是原子的局部 id
        if (lats->vac_hash.count(latti_id)) {
          buffer[index].id = len;
          buffer[index].type = LatticeTypes{LatticeTypes::V};
          index++;
        } else if (lats->re_hash.count(latti_id)) {
          buffer[index].id = len;
          buffer[index].type = LatticeTypes{LatticeTypes::Re};
          index++;
        } else if (lats->mn_hash.count(latti_id)) {
          buffer[index].id = len;
          buffer[index].type = LatticeTypes{LatticeTypes::Mn};
          index++;
        } else if (lats->ni_hash.count(latti_id)) {
          buffer[index].id = len;
          buffer[index].type = LatticeTypes{LatticeTypes::Ni};
          index++;
        } else if (lats->si_hash.count(latti_id)) {
          buffer[index].id = len;
          buffer[index].type = LatticeTypes{LatticeTypes::Si};
          index++;
        } else if (lats->momo_hash.count(latti_id)) {
          buffer[index].id = len;
          buffer[index].type = LatticeTypes{LatticeTypes::MoMo};
          index++;
        } else if (lats->more_hash.count(latti_id)) {
          buffer[index].id = len;
          buffer[index].type = LatticeTypes{LatticeTypes::MoRe};
          index++;
        } else if (lats->rere_hash.count(latti_id)) {
          buffer[index].id = len;
          buffer[index].type = LatticeTypes{LatticeTypes::ReRe};
          index++;
        } else {
          // Mo
        }
        len++;
        
        _type_lattice_id latti_id2 = lats->getId(BCC_DBX * x + 1, y, z);
        if (lats->vac_hash.count(latti_id2)) {
          buffer[index].id = len;
          buffer[index].type = LatticeTypes{LatticeTypes::V};
          index++;
        } else if (lats->re_hash.count(latti_id2)) {
          buffer[index].id = len;
          buffer[index].type = LatticeTypes{LatticeTypes::Re};
          index++;
        } else if (lats->mn_hash.count(latti_id2)) {
          buffer[index].id = len;
          buffer[index].type = LatticeTypes{LatticeTypes::Mn};
          index++;
        } else if (lats->ni_hash.count(latti_id2)) {
          buffer[index].id = len;
          buffer[index].type = LatticeTypes{LatticeTypes::Ni};
          index++;
        } else if (lats->si_hash.count(latti_id2)) {
          buffer[index].id = len;
          buffer[index].type = LatticeTypes{LatticeTypes::Si};
          index++;
        } else if (lats->momo_hash.count(latti_id2)) {
          buffer[index].id = len;
          buffer[index].type = LatticeTypes{LatticeTypes::MoMo};
          index++;
        } else if (lats->more_hash.count(latti_id2)) {
          buffer[index].id = len;
          buffer[index].type = LatticeTypes{LatticeTypes::MoRe};
          index++;
        } else if (lats->rere_hash.count(latti_id2)) {
          buffer[index].id = len;
          buffer[index].type = LatticeTypes{LatticeTypes::ReRe};
          index++;
        } else {
          // Mo
        }
        len++;
        // buffer[len++] = lats->_lattices[z][y][BCC_DBX * x];
        // buffer[len++] = lats->_lattices[z][y][BCC_DBX * x + 1];
      }
    }
  }
  assert(index == send_len);
}

void GhostInitPacker::onReceive(Lattice buffer[], const unsigned long receive_len, const int dimension,
                 const int direction){
  const comm::Region<comm::_type_lattice_size> recv_region = comm::fwCommLocalRecvRegion(
      p_domain->lattice_size_ghost, p_domain->local_sub_box_lattice_region, dimension, direction);
  unsigned long len = 0;
  for (int z = recv_region.z_low; z < recv_region.z_high; z++) {
    for (int y = recv_region.y_low; y < recv_region.y_high; y++) {
      for (int x = recv_region.x_low; x < recv_region.x_high; x++) {
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

void GhostInitPacker::onReceive2(Lattice *buffer, const int receive_len, const int dimension, const int direction) {
  const comm::Region<comm::_type_lattice_size> recv_region = comm::fwCommLocalRecvRegion(
      p_domain->lattice_size_ghost, p_domain->local_sub_box_lattice_region, dimension, direction);

  int len = 0;
  _type_lattice_count index = 0;

  if(receive_len == 0){
    // 说明全是Fe原子
    for (_type_lattice_id z = recv_region.z_low; z < recv_region.z_high; z++) {
      for (_type_lattice_id y = recv_region.y_low; y < recv_region.y_high; y++) {
        for (_type_lattice_id x = recv_region.x_low; x < recv_region.x_high; x++) {
          _type_lattice_id latti_id = lats->getId(BCC_DBX * x, y, z);
          // lats->isUniqueLattiId(latti_id);
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
          _type_lattice_id latti_id2 = lats->getId(BCC_DBX * x + 1, y, z);
          // lats->isUniqueLattiId(latti_id2);
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
        }
      }
    }
  } else {
    
    // kiwi::logs::v(" ", " DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD .\n");
    for (_type_lattice_id z = recv_region.z_low; z < recv_region.z_high; z++) {
      for (_type_lattice_id y = recv_region.y_low; y < recv_region.y_high; y++) {
        for (_type_lattice_id x = recv_region.x_low; x < recv_region.x_high; x++) {
          
          _type_lattice_id latti_id = lats->getId(BCC_DBX * x, y, z);
          // lats->isUniqueLattiId(latti_id);
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

          if(buffer[index].id == len) {
            // 在 buffer 中找到了
            switch (buffer[index].type._type)
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
            index++;
          }
          len++;

          _type_lattice_id latti_id2 = lats->getId(BCC_DBX * x + 1, y, z);
          // lats->isUniqueLattiId(latti_id);
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

          if(buffer[index].id == len){
            switch (buffer[index].type._type)
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
            index++;
          }
          len++;
        }
      }
    }
  }
  assert(index == receive_len);
  // std::cout << "3333333333333333333333333333333333" << std::endl;
  // todo set more information if the lattice is vacancy or dumbbell.
}
