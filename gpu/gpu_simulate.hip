#include "gpu_simulate.h"

// Global Variables
namespace gpuSimulate {

    // device variable
    double *d_sum_rates;
    double *d_config;
    dev_event *d_event;

    dev_Vacancy *deviceVacancyInfo;
    dev_nnLattice *deviceNeighbourInfo;

    // host variable
    _type_lattice_count currentVacCount;
    double h_config[2];

    dev_Vacancy *hostVacancyInfo;
    dev_nnLattice *hostNeighbourInfo;

    int blockSize;
    int numBlocks;
    _type_lattice_count vacHostArrayCapacity;

    double *sum_rates;
    dev_event* h_event;

    _type_lattice_id to_x;
    _type_lattice_id to_y;
    _type_lattice_id to_z;

    LatticesList *lattice_list;
    comm::ColoredDomain *p_domain;

}

__device__ double _1nnMapToSwitch(bonds::PairBond::bond_type bondType){
    switch (bondType) {
      case bonds::PairBond::MoMo:
        return -1.991;
      case bonds::PairBond::VV:
        return 0.199;
      case bonds::PairBond::ReRe:
        return -2.201;
      case bonds::PairBond::NiNi:
        return -1.618;
      case bonds::PairBond::MnMn:
        return -0.271;
      case bonds::PairBond::SiSi:
        return -0.716;
      case bonds::PairBond::VMo:
        return -0.755318;
      case bonds::PairBond::MoRe:
        return -2.096;
      case bonds::PairBond::MoNi:
        return -1.8045;
      case bonds::PairBond::MoMn:
        return -0.446;
      case bonds::PairBond::MoSi:
        return -0.683;
      case bonds::PairBond::VRe:
        return -0.920318;
      case bonds::PairBond::VNi:
        return -0.708818;
      case bonds::PairBond::VMn:
        return -0.038;
      case bonds::PairBond::VSi:
        return -0.344;
      case bonds::PairBond::ReNi:
        return -1.9095;
      case bonds::PairBond::ReMn:
        return -0.365;
      case bonds::PairBond::ReSi:
        return -0.611;
      case bonds::PairBond::MnNi:
        return -0.366;
      case bonds::PairBond::MnSi:
        return -0.488;
      case bonds::PairBond::NiSi:
        return -0.723;
      default:
        // 处理默认情况
        return 0.0;
    }
}

__device__ double _2nnMapToSwitch(bonds::PairBond::bond_type bondType){
    switch (bondType) {
      case bonds::PairBond::MoMo:
        return -0.9955;
      case bonds::PairBond::VV:
        return -0.245;
      case bonds::PairBond::ReRe:
        return -1.1005;
      case bonds::PairBond::NiNi:
        return -0.809;
      case bonds::PairBond::MnMn:
        return -0.611;
      case bonds::PairBond::SiSi:
        return -0.389;
      case bonds::PairBond::VMo:
        return -0.377659;
      case bonds::PairBond::MoRe:
        return -1.048;
      case bonds::PairBond::MoNi:
        return -0.90225;
      case bonds::PairBond::MoMn:
        return -0.631;
      case bonds::PairBond::MoSi:
        return -0.469;
      case bonds::PairBond::VRe:
        return -0.460159;
      case bonds::PairBond::VNi:
        return -0.354409;
      case bonds::PairBond::VSi:
        return -0.248;
      case bonds::PairBond::VMn:
        return -0.203;
      case bonds::PairBond::ReNi:
        return -0.95475;
      case bonds::PairBond::ReMn:
        return -0.621;
      case bonds::PairBond::ReSi:
        return -0.566;
      case bonds::PairBond::MnNi:
        return -0.496;
      case bonds::PairBond::MnSi:
        return -0.316;
      case bonds::PairBond::NiSi:
        return -0.521;
      default:
        // 处理默认情况
        return 0.0;
    }
}

__global__ void calcRatesGPU(dev_Vacancy *deviceVacancyInfo, dev_nnLattice *deviceNeighbourInfo,
                             double* d_sum_rates, double attempt_freq, double temperature, _type_lattice_count vac_total) {
  // 获取线程的全局ID
  _type_lattice_count tid = blockIdx.x * blockDim.x + threadIdx.x;

  // 洗牌函数、寄存器和共享内存优化原子加
  if (tid == 0) *d_sum_rates = 0.0;
  __shared__ double tmp_sh;
  if (threadIdx.x == 0) tmp_sh = 0.0;
  __syncthreads();

  // 判断当前tid是否可用
  if(tid < vac_total) {
    double totalRatesOfCurVac = 0.0;
    
    deviceVacancyInfo[tid].rates[0] = 0;
    deviceVacancyInfo[tid].rates[1] = 0;
    deviceVacancyInfo[tid].rates[2] = 0;
    deviceVacancyInfo[tid].rates[3] = 0;
    deviceVacancyInfo[tid].rates[4] = 0;
    deviceVacancyInfo[tid].rates[5] = 0;
    deviceVacancyInfo[tid].rates[6] = 0;
    deviceVacancyInfo[tid].rates[7] = 0;

    double vac_energy = 0;
    int first_ir = 1;
    unsigned char atom_nei_status = 0;
    for (int b = 0; b < 8; b++) {
      if (deviceNeighbourInfo[tid * NN_TOTAL + b].type._type <= LatticeTypes::Si && deviceNeighbourInfo[tid * NN_TOTAL + b].type._type != LatticeTypes::V) { // can trans
        atom_nei_status |= 1 << b;
      }
    }

    for (unsigned char b = 0; b < 8; b++) {
      if ((atom_nei_status >> b) & 1) { // the neighbour lattice is available
        double active_energy;
        switch (deviceNeighbourInfo[tid * NN_TOTAL + b].type._type) {
          case LatticeTypes::Mo:
            active_energy = 0.65;
            break;
          case LatticeTypes::Re:
            active_energy = 0.54;
            break;
          case LatticeTypes::Mn:
            active_energy = 1.03;
            break;
          case LatticeTypes::Ni:
            active_energy = 0.68;
            break;
          case LatticeTypes::Si:
            active_energy = 0.42;
            break;
          default:
            active_energy = 0;
            break;
        }
        // calculate system energy before transition.
        // const LatticeTypes::lat_type ghost_atom = d_lat_list[tid * 8 + b].type._type;
        double e_before = 0;
        {
          // bonds energy of src lattice contributed by its 1nn/2nn neighbour lattice.
          double e_src, e_des;
          unsigned char _1nn_status = 0xFF;
          unsigned char _2nn_status = 0x2F;
          double energy = 0.0;
          if(first_ir == 1) {
            // Traver all 1nn neighbour lattices, and calculate bond energy contribution.
            for (int b2 = 0; b2 < 8; b2++) {
              if ((_1nn_status >> b2) & 0x1) {
                // we assume that src_atom_type is single atom or vacancy.
                if (deviceNeighbourInfo[tid * NN_TOTAL + b2].type._type <= LatticeTypes::Si) {
                  // it is vacancy or single atom.
                  LatticeTypes::lat_type temp_1nn = deviceNeighbourInfo[tid * NN_TOTAL + b2].type._type;
                  LatticeTypes::lat_type temp_src = LatticeTypes::V;
                  if (temp_1nn > temp_src) {
                    LatticeTypes::lat_type temp = temp_1nn;
                    temp_1nn = temp_src;
                    temp_src = temp;
                  }
                  // PairBond::bond_type bond = PairBond::makeBond(_1nn_neighbour[b]->type, src_atom_type);
                  const bonds::PairBond::bond_type bond = static_cast<bonds::PairBond::bond_type>((temp_1nn << 8) | temp_src);
                  vac_energy += _1nnMapToSwitch(bond);
                }
              }
            }
        
            // Traver all 2nn neighbour lattices, and calculate bond energy contribution.
            for (int b2 = 0; b2 < 6; b2++) {
              if ((_2nn_status >> b2) & 0x1) {
                // we assume that src_atom_type is single atom or vacancy.

                if (deviceNeighbourInfo[tid * NN_TOTAL + 8 + b2].type._type <= LatticeTypes::Si) {
                  // it is vacancy or single atom.
                  LatticeTypes::lat_type temp_2nn = deviceNeighbourInfo[tid * NN_TOTAL + 8 + b2].type._type;
                  LatticeTypes::lat_type temp_src = LatticeTypes::V;
                  if (temp_2nn > temp_src) {
                    LatticeTypes::lat_type temp = temp_2nn;
                    temp_2nn = temp_src;
                    temp_src = temp;
                  }
                  // PairBond::bond_type bond = PairBond::makeBond(_2nn_neighbour[b3]->type, src_atom_type);
                  const bonds::PairBond::bond_type bond = static_cast<bonds::PairBond::bond_type>((temp_2nn << 8) | temp_src);
                  vac_energy += _2nnMapToSwitch(bond);
                }
              }
            }
            first_ir = 0;
          }
          e_src = vac_energy;
//_______________________________________
          energy = 0.0;

          for (int b2 = 0; b2 < 8; b2++) {
            if ((_1nn_status >> b2) & 0x1) {
              // we assume that src_atom_type is single atom or vacancy.
              if (deviceNeighbourInfo[tid * NN_TOTAL + 14 + b * 8 + b2].type._type <= LatticeTypes::Si) {
                // it is vacancy or single atom.
                LatticeTypes::lat_type temp_1nn = deviceNeighbourInfo[tid * NN_TOTAL + 14 + b * 8 + b2].type._type;
                LatticeTypes::lat_type temp_src = deviceNeighbourInfo[tid * NN_TOTAL + b].type._type;
                if (temp_1nn > temp_src) {
                  LatticeTypes::lat_type temp = temp_1nn;
                  temp_1nn = temp_src;
                  temp_src = temp;
                }
                // PairBond::bond_type bond = PairBond::makeBond(_1nn_neighbour[b]->type, src_atom_type);
                const bonds::PairBond::bond_type bond = static_cast<bonds::PairBond::bond_type>((temp_1nn << 8) | temp_src);
                energy += _1nnMapToSwitch(bond);
              }
            }
          }

          for (int b2 = 0; b2 < 6; b2++) {
            if ((_2nn_status >> b2) & 0x1) {
              // we assume that src_atom_type is single atom or vacancy.
              if (deviceNeighbourInfo[tid * NN_TOTAL + 78 + b * 6 + b2].type._type <= LatticeTypes::Si) {
                // it is vacancy or single atom.
                LatticeTypes::lat_type temp_2nn = deviceNeighbourInfo[tid * NN_TOTAL + 78 + b * 6 + b2].type._type;
                LatticeTypes::lat_type temp_src = deviceNeighbourInfo[tid * NN_TOTAL + b].type._type;
                if (temp_2nn > temp_src) {
                  LatticeTypes::lat_type temp = temp_2nn;
                  temp_2nn = temp_src;
                  temp_src = temp;
                }
                // PairBond::bond_type bond = PairBond::makeBond(_2nn_neighbour[b3]->type, src_atom_type);
                const bonds::PairBond::bond_type bond = static_cast<bonds::PairBond::bond_type>((temp_2nn << 8) | temp_src);
                energy += _2nnMapToSwitch(bond);
              }
            }
          }
          // if(is_first) *energy_temp2 = energy;
          // is_first = false;
          e_des = energy;

          // double e_des = bonds::BondsCounter::count(&lattice_list, target_id, target_lattice.type);
          e_before = e_src + e_des;
          // energy_temp2[b] = energy;
        }

// ______________________________________________________________________________________________________________________________________



        // calculate system energy after transition.

        double e_after = 0;
        {
          // bonds energy of src lattice contributed by its 1nn/2nn neighbour lattice.
          double e_src, e_des;
          unsigned char _1nn_status = 0xFF;
              
          double energy = 0.0;

          for (int b2 = 0; b2 < 8; b2++) {
            if ((_1nn_status >> b2) & 0x1) {
              // we assume that src_atom_type is single atom or vacancy.
              LatticeTypes::lat_type temp_type = deviceNeighbourInfo[tid * NN_TOTAL + b2].type._type;
              if (b2 == b){
                temp_type = LatticeTypes::V;
              }
              if (temp_type <= LatticeTypes::Si) {
                // it is vacancy or single atom.
                LatticeTypes::lat_type temp_1nn = temp_type; 
                LatticeTypes::lat_type temp_src = deviceNeighbourInfo[tid * NN_TOTAL + b].type._type;
                if (temp_1nn > temp_src) {
                  LatticeTypes::lat_type temp = temp_1nn;
                  temp_1nn = temp_src;
                  temp_src = temp;
                }
                // PairBond::bond_type bond = PairBond::makeBond(_1nn_neighbour[b]->type, src_atom_type);
                bonds::PairBond::bond_type bond = static_cast<bonds::PairBond::bond_type>((temp_1nn << 8) | temp_src);
                energy += _1nnMapToSwitch(bond);
              }
            }
          }

          // Traver all 2nn neighbour lattices, and calculate bond energy contribution.
          unsigned char _2nn_status = 0x2F;
          
          // lat_list->get2nn(source_id, _2nn_neighbour);
          for (int b2 = 0; b2 < 6; b2++) {
            if ((_2nn_status >> b2) & 0x1) {
              // we assume that src_atom_type is single atom or vacancy.
              if (deviceNeighbourInfo[tid * NN_TOTAL + 8 + b2].type._type <= LatticeTypes::Si) {
                // it is vacancy or single atom.
                LatticeTypes::lat_type temp_2nn = deviceNeighbourInfo[tid * NN_TOTAL + 8 + b2].type._type;
                LatticeTypes::lat_type temp_src = deviceNeighbourInfo[tid * NN_TOTAL + b].type._type;
                if (temp_2nn > temp_src) {
                  LatticeTypes::lat_type temp = temp_2nn;
                  temp_2nn = temp_src;
                  temp_src = temp;
                }
                // PairBond::bond_type bond = PairBond::makeBond(_2nn_neighbour[b3]->type, src_atom_type);
                const bonds::PairBond::bond_type bond = static_cast<bonds::PairBond::bond_type>((temp_2nn << 8) | temp_src);
                energy += _2nnMapToSwitch(bond);
              }
            }
          }
          e_src = energy;
//_______________________________________
          energy = 0.0;
          for (int b2 = 0; b2 < 8; b2++) {
            if ((_1nn_status >> b2) & 0x1) {
              // we assume that src_atom_type is single atom or vacancy.
              LatticeTypes::lat_type temp_type = deviceNeighbourInfo[tid * NN_TOTAL + 14 + b * 8 + b2].type._type;
              if (b2 + b == 7) {
                temp_type = deviceNeighbourInfo[tid * NN_TOTAL + b].type._type;
              }
              // if (d_target_lattice_1neighbour[tid * 64 + b * 8 + b2].id == d_lattices[tid].id){
              //   temp_type = d_lat_list[tid * 8 + b].type._type;
              // }
              if (temp_type <= LatticeTypes::Si) {
                // it is vacancy or single atom.
                LatticeTypes::lat_type temp_1nn = temp_type;
                LatticeTypes::lat_type temp_src = LatticeTypes::V;
                if (temp_1nn > temp_src) {
                  LatticeTypes::lat_type temp = temp_1nn;
                  temp_1nn = temp_src;
                  temp_src = temp;
                }
                // PairBond::bond_type bond = PairBond::makeBond(_1nn_neighbour[b]->type, src_atom_type);
                const bonds::PairBond::bond_type bond = static_cast<bonds::PairBond::bond_type>((temp_1nn << 8) | temp_src);
                energy += _1nnMapToSwitch(bond);
              }
            }
          }
          for (int b2 = 0; b2 < 6; b2++) {
            if ((_2nn_status >> b2) & 0x1) {
              // we assume that src_atom_type is single atom or vacancy.
              if (deviceNeighbourInfo[tid * NN_TOTAL + 78 + b * 6 + b2].type._type <= LatticeTypes::Si) {
                // it is vacancy or single atom.
                LatticeTypes::lat_type temp_2nn = deviceNeighbourInfo[tid * NN_TOTAL + 78 + b * 6 + b2].type._type;;
                LatticeTypes::lat_type temp_src = LatticeTypes::V;
                if (temp_2nn > temp_src) {
                  LatticeTypes::lat_type temp = temp_2nn;
                  temp_2nn = temp_src;
                  temp_src = temp;
                }
                // PairBond::bond_type bond = PairBond::makeBond(_2nn_neighbour[b3]->type, src_atom_type);
                const bonds::PairBond::bond_type bond = static_cast<bonds::PairBond::bond_type>((temp_2nn << 8) | temp_src);
                energy += _2nnMapToSwitch(bond);
              }
            }
          }
          e_des = energy;
          e_after = e_src + e_des;
        }
        active_energy += (e_after - e_before) / 2;
        double rate = attempt_freq * __expf(-active_energy / (8.617330350E-5 * temperature));
        deviceVacancyInfo[tid].rates[b] = rate;
        totalRatesOfCurVac += rate;
      }
    }
    int delta = 1;
    for (int j = 0; j < 6; j++) {
      totalRatesOfCurVac += __shfl_down(totalRatesOfCurVac, delta, 64);
      delta += delta;
    }
    if (threadIdx.x % 64 == 0) atomicAdd(&tmp_sh, totalRatesOfCurVac);
    __syncthreads();
    if (threadIdx.x == 0) atomicAdd(d_sum_rates, tmp_sh);
  }
}

__global__ void selectEventGPU(dev_Vacancy *deviceVacancyInfo, dev_nnLattice *deviceNeighbourInfo, 
                               dev_event *d_event, _type_lattice_count vac_total, double excepted_rate) {

  // 获取线程的全局ID
  _type_lattice_count tid = blockIdx.x * blockDim.x + threadIdx.x;
  if(tid == 0) {
    double rate_accumulator = 0;
    int tag = 1;
    for(_type_lattice_count i = 0; i < vac_total; i++) {
      if(tag == 1) {
        for (int rate_index = 0; rate_index < 8; rate_index++) {
          rate_accumulator += deviceVacancyInfo[i].rates[rate_index];
          if(rate_accumulator > excepted_rate) {
            if(tag == 1) {
              d_event->from_id = deviceVacancyInfo[i].id;
              d_event->to_id = deviceNeighbourInfo[i * NN_TOTAL + rate_index].id;
              d_event->to_type._type = deviceNeighbourInfo[i * NN_TOTAL + rate_index].type._type;
              tag = 0;
              break;
            }
          }
        } 
      } else break;
    }
  }
}

void selectAndPerformEventGPU(double excepted_rand, int rank, _type_lattice_count step, int sect,
                              std::array<std::unordered_set<_type_lattice_id>, 7>& exchange_ghost, const unsigned int sector_id,
                              std::array<std::unordered_set<_type_lattice_id>, 8>& exchange_surface_x, 
                              std::array<std::unordered_set<_type_lattice_id>, 8>& exchange_surface_y,
                              std::array<std::unordered_set<_type_lattice_id>, 8>& exchange_surface_z,
                              std::array<comm::Region<comm::_type_lattice_size>, 56>& ghost_region){
  using namespace gpuSimulate;

  blockSize = 1;
  numBlocks = 1;
  hipLaunchKernelGGL(selectEventGPU, numBlocks, blockSize, 0, 0, 
                     deviceVacancyInfo,
                     deviceNeighbourInfo,
                     d_event,
                     currentVacCount,
                     *sum_rates * excepted_rand);

  CHECK_FOR_ERROR();

  HANDLE_HIP(hipMemcpy(h_event, d_event, sizeof(dev_event), hipMemcpyDeviceToHost));
  // kiwi::logs::v(" ", " excepted_rate is : {} step is : {} sect is : {} from_id is : {} to_id is : {}. \n", *sum_rates * excepted_rand, step, sect, h_event->from_id, h_event->to_id);
  to_x = h_event->to_id % lattice_list->meta.size_x;
  to_y = (h_event->to_id / lattice_list->meta.size_x) % lattice_list->meta.size_y;
  to_z = h_event->to_id / (lattice_list->meta.size_x * lattice_list->meta.size_y);

  exchange_lattices(h_event);
  if(lattice_list->meta.isGhostLat(h_event->to_id)) {
    for(int i = 0; i < 7; i++) {
      if(ghost_region[sector_id * 7 + i].x_low <= to_x && to_x < ghost_region[sector_id * 7 + i].x_high && 
         ghost_region[sector_id * 7 + i].y_low <= to_y && to_y < ghost_region[sector_id * 7 + i].y_high && 
         ghost_region[sector_id * 7 + i].z_low <= to_z && to_z < ghost_region[sector_id * 7 + i].z_high) {
          
        if(!exchange_ghost[i].count(h_event->to_id)) exchange_ghost[i].emplace(h_event->to_id);
        break;
      }
    }
  }
  if (lattice_list->meta.isSurfaceLat(h_event->from_id)) add_surface(h_event->from_id, exchange_surface_x, exchange_surface_y, exchange_surface_z);
  if (lattice_list->meta.isSurfaceLat(h_event->to_id)) add_surface(h_event->to_id, exchange_surface_x, exchange_surface_y, exchange_surface_z);
}

double calculate_GPU(_type_lattice_id *vac_idArray,  const _type_lattice_count& vac_count) {
  using namespace gpuSimulate;
  
  *sum_rates = 0;
  currentVacCount = vac_count;
  // 如果超出预分配数组大小则重新分配锁页内存和显存
  if(vac_count > vacHostArrayCapacity) {
    vacHostArrayCapacity = vac_count * RE_SCALE_SIZE;
    hipHostFree(hostVacancyInfo);
    hipHostFree(hostNeighbourInfo);
    HANDLE_HIP(hipFree(deviceVacancyInfo));
    HANDLE_HIP(hipFree(deviceNeighbourInfo));

    hipHostMalloc((void**)&hostVacancyInfo, vacHostArrayCapacity * sizeof(dev_Vacancy));
    hipHostMalloc((void**)&hostNeighbourInfo, vacHostArrayCapacity * NN_TOTAL * sizeof(dev_nnLattice));
    HANDLE_HIP(hipMalloc((void**)&deviceVacancyInfo, vacHostArrayCapacity * sizeof(dev_Vacancy)));
    HANDLE_HIP(hipMalloc((void**)&deviceNeighbourInfo, vacHostArrayCapacity * NN_TOTAL * sizeof(dev_nnLattice)));
  }

  lattice_list->initGpuInfo(vac_count, vac_idArray, hostVacancyInfo, hostNeighbourInfo);
  HANDLE_HIP(hipMemcpy(deviceVacancyInfo, hostVacancyInfo, vac_count * sizeof(dev_Vacancy), hipMemcpyHostToDevice));
  HANDLE_HIP(hipMemcpy(deviceNeighbourInfo, hostNeighbourInfo, vac_count * NN_TOTAL * sizeof(dev_nnLattice), hipMemcpyHostToDevice));
  // for (_type_lattice_count j = 0; j < vac_count; j++) {
  //   kiwi::logs::v(" ", " latti_id is : {}. \n", hostVacancyInfo[j].id);
  //   for (int i = 0; i < 8; i++) {
  //     kiwi::logs::v(" ", " nn_index is : {} nn_id is : {} nn_type is : {} . \n", i, deviceNeighbourInfo[j * NN_TOTAL + i].id, deviceNeighbourInfo[j * NN_TOTAL + i].type._type);
  //   }
  // }
  blockSize = 128;
  numBlocks = (vac_count + blockSize - 1) / blockSize;

  // 调用GPU内核函数
  hipLaunchKernelGGL(calcRatesGPU, numBlocks, blockSize, 0, 0, 
                     deviceVacancyInfo,
                     deviceNeighbourInfo,
                     d_sum_rates,
                     h_config[0],
                     h_config[1], 
                     vac_count);

  free(vac_idArray);

  CHECK_FOR_ERROR();
  HANDLE_HIP(hipMemcpy(sum_rates, d_sum_rates, sizeof(double), hipMemcpyDeviceToHost));
  return *sum_rates;
}

void initialize_gpu(int process_rank){
  int gpu_count;
  HANDLE_HIP(hipGetDeviceCount(&gpu_count));
  HANDLE_HIP(hipSetDevice(process_rank % gpu_count));
  HANDLE_HIP(hipDeviceReset());
}

void gpu_prepare(double& v, double& T, LatticesList *p_list, comm::ColoredDomain *_p_domain){
  using namespace gpuSimulate;

  lattice_list = p_list;
  p_domain = _p_domain;
  h_config[0] = v;   // v
  h_config[1] = T;   // T

  HANDLE_HIP(hipMalloc((void**)&d_sum_rates, sizeof(double)));
  HANDLE_HIP(hipMalloc((void**)&d_event, sizeof(dev_event)));
  
  vacHostArrayCapacity = lattice_list->vac_hash.size() / 6 + 1;

  HANDLE_HIP(hipMalloc((void**)&deviceVacancyInfo, vacHostArrayCapacity * sizeof(dev_Vacancy)));
  HANDLE_HIP(hipMalloc((void**)&deviceNeighbourInfo, vacHostArrayCapacity * NN_TOTAL * sizeof(dev_nnLattice)));

  hipHostMalloc((void**)&sum_rates, sizeof(double));
  hipHostMalloc((void**)&h_event, sizeof(dev_event));

  hipHostMalloc((void**)&hostVacancyInfo, vacHostArrayCapacity * sizeof(dev_Vacancy));
  hipHostMalloc((void**)&hostNeighbourInfo, vacHostArrayCapacity * NN_TOTAL * sizeof(dev_nnLattice));
}

void gpu_final(){
  using namespace gpuSimulate;

  // 释放GPU内存
  HANDLE_HIP(hipFree(d_sum_rates));
  HANDLE_HIP(hipFree(d_event));

  HANDLE_HIP(hipFree(deviceVacancyInfo));
  HANDLE_HIP(hipFree(deviceNeighbourInfo));

  hipHostFree(h_event);
  hipHostFree(sum_rates);

  hipHostFree(hostVacancyInfo);
  hipHostFree(hostNeighbourInfo);
}

void exchange_lattices(dev_event* h_event){
  using namespace gpuSimulate;

  if (lattice_list->vac_hash.count(h_event->from_id)) {
    lattice_list->vac_hash.erase(h_event->from_id);
    lattice_list->vac_hash.emplace(std::make_pair(h_event->to_id, VacancyHash{}));
  } else assert(false);

  switch (h_event->to_type._type)
  {
    case LatticeTypes::Re:
      lattice_list->re_hash.erase(h_event->to_id);
      lattice_list->re_hash.emplace(h_event->from_id);
      break;
    case LatticeTypes::Mn:
      lattice_list->mn_hash.erase(h_event->to_id);
      lattice_list->mn_hash.emplace(h_event->from_id);
      break;
    case LatticeTypes::Ni:
      lattice_list->ni_hash.erase(h_event->to_id);
      lattice_list->ni_hash.emplace(h_event->from_id);
      break;
    case LatticeTypes::Si:
      lattice_list->si_hash.erase(h_event->to_id);
      lattice_list->si_hash.emplace(h_event->from_id);
      break;
    case LatticeTypes::Mo:
      break;
    default:
      assert(false);
      break;
  }
}

void add_surface(_type_lattice_id& surface_id,
                 std::array<std::unordered_set<_type_lattice_id>, 8>& exchange_surface_x, 
                 std::array<std::unordered_set<_type_lattice_id>, 8>& exchange_surface_y,
                 std::array<std::unordered_set<_type_lattice_id>, 8>& exchange_surface_z) {
  using namespace gpuSimulate;

  _type_lattice_id surface_x = surface_id % lattice_list->meta.size_x;
  _type_lattice_id surface_y = (surface_id / lattice_list->meta.size_x) % lattice_list->meta.size_y;
  _type_lattice_id surface_z = surface_id / (lattice_list->meta.size_x * lattice_list->meta.size_y);
  const int dims[comm::DIMENSION_SIZE] = {comm::DIM_X, comm::DIM_Y, comm::DIM_Z};
  // todo the regions can be static.
  std::array<std::vector<comm::Region<comm::_type_lattice_coord>>, comm::DIMENSION_SIZE> send_regions; // send regions in each dimension
  // 应该可以不循环 8 次的 ？ 如果发送区没有重叠的话
  for (int sect = 0; sect < 8; sect++) {
    for (int d = 0; d < comm::DIMENSION_SIZE; d++) {
      send_regions[d] = comm::fwCommSectorSendRegion(sect, dims[d], p_domain->lattice_size_ghost,
                                                     p_domain->local_split_coord, p_domain->local_sub_box_lattice_region);
    }
    for(int d = 0; d < comm::DIMENSION_SIZE; d++) {
      for (auto &r : send_regions[d]) {
        _type_lattice_size x_low = 2 * r.x_low;
        _type_lattice_size y_low = r.y_low;
        _type_lattice_size z_low = r.z_low;
        _type_lattice_size x_high = 2 * r.x_high;
        _type_lattice_size y_high = r.y_high;
        _type_lattice_size z_high = r.z_high;
    
        if(x_low <= surface_x && surface_x < x_high && y_low <= surface_y && surface_y < y_high && z_low <= surface_z && surface_z < z_high) {
          if(d == 0){
            if(!exchange_surface_x[sect].count(surface_id)) exchange_surface_x[sect].emplace(surface_id);
          } else if(d == 1) {
            if(!exchange_surface_y[sect].count(surface_id)) exchange_surface_y[sect].emplace(surface_id);
          } else {
            if(!exchange_surface_z[sect].count(surface_id)) exchange_surface_z[sect].emplace(surface_id);
          }
        }
      }
    }
  }
}

void output_time() {
  using namespace gpuSimulate;
  // kiwi::logs::v(" ", " time_event is : {} s. \n", time_event);
  // kiwi::logs::v(" ", " time_event_trans is : {} s. \n", time_event_trans);
  // kiwi::logs::v(" ", " time_event_other : {} s. \n", time_event_other);
  // kiwi::logs::v(" ", " time_update_construct : {} s. \n", time_update_construct);
  // kiwi::logs::v(" ", " time_update_trans : {} s. \n", time_update_trans);
}