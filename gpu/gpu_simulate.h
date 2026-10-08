#ifndef GPU_SIMULATE_H
#define GPU_SIMULATE_H

#include "hip/hip_runtime.h"
#include <hip/hip_runtime.h>
#include "gpuError.h"
#include "DeviceVacRatesSolver.h"
#include <cassert>
#include <iostream>

#include <logs/logs.h>
#include "../models/abvi/kmc.h"
#include "../models/abvi/box.h"
#include "../models/abvi/rate/vacancy_rates_solver.h"
#include "../models/abvi/defect/vacancy_list.h"
#include "../models/abvi/rate/bonds/pair_bond.hpp"
#include "../models/abvi/rate/rates_types.h"
#include "../models/abvi/event.h"

#include "../src/type_define.h"
#include "../src/lattice/lattice.h"
#include "../src/lattice/lattices_list.h"
#include "../src/lattice/lattice_types.h"
#include "../src/lattice/lattice_list_meta.h"

double calculate_GPU(_type_lattice_id *vac_idArray,  const _type_lattice_count& vac_count);

void selectAndPerformEventGPU(double excepted_rand, int rank, _type_lattice_count step, int sect,
                              std::array<std::unordered_set<_type_lattice_id>, 7>& exchange_ghost, const unsigned int sector_id,
                              std::array<std::unordered_set<_type_lattice_id>, 8>& exchange_surface_x, 
                              std::array<std::unordered_set<_type_lattice_id>, 8>& exchange_surface_y,
                              std::array<std::unordered_set<_type_lattice_id>, 8>& exchange_surface_z,
                              std::array<comm::Region<comm::_type_lattice_size>, 56>& ghost_region);

void initialize_gpu(int process_rank);

void gpu_final();

void gpu_prepare(double& v, double& T, LatticesList *p_list, comm::ColoredDomain *_p_domain);

void exchange_lattices(dev_event* h_event);

void add_surface(_type_lattice_id& surface_id,
                 std::array<std::unordered_set<_type_lattice_id>, 8>& exchange_surface_x, 
                 std::array<std::unordered_set<_type_lattice_id>, 8>& exchange_surface_y,
                 std::array<std::unordered_set<_type_lattice_id>, 8>& exchange_surface_z);
                 
void output_time();
#endif /*GPU_SIMULATE_H*/