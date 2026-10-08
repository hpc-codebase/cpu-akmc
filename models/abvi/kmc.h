//
// Created by genshen on 2018-12-12.
//

#ifndef MISA_KMC_KMC_H
#define MISA_KMC_KMC_H

#include <vector>
#include <array>
#include <cstdint>
#include <random>
#include "utils/random/rng_type.h"
#include "lattice/lattice.h"
#include "box.h"
#include "event.h"
#include "plugin/event_listener.h"
#include "type_define.h"
#include <models/model_adapter.h>
#include <comm/domain/colored_domain.h>
#include <unordered_map>
#include "abvi/defect/vac_hash.h"
#include <comm/preset/sector_forwarding_region.h>
#include "utils/simulation_domain.h"
#include <comm/comm.hpp>
#include <cmath>

// #include "../../src/algorithms/sl/sublattice.h"

// #include "../gpu/gpu_simulate.h"
/*!
 * \brief the model routine of KMC simulation, including rate calculation, event selecting
 * and execution implementation.
 */

/*!
 * \brief CPU 版 recb 验证统计，字段与 akmc GPU 版 RecbVerifyStats 一一对应，
 * 保证两侧 [recb-verify]/[recb-detail]/[recb-ops]/[recb-directions] 日志
 * 统计口径一致、可逐项对照。每次 recb launch 由 recb_checki 重置并汇总打印。
 */
struct RecbVerifyStatsCpu {
  unsigned int sector_pairs = 0;               // 本次 launch 本扇区参与游走的对数（只计 MoMo/MoRe，与 GPU 同口径）
  unsigned long long attempts = 0;             // 所有 walker 随机选向的总尝试次数
  unsigned long long successful_hops = 0;      // 实际完成换位的跳跃次数
  unsigned int max_successful_hops = 0;        // 单个 walker 一次 launch 内最多成功跳跃数
  unsigned int multi_hop_pairs = 0;            // 成功跳跃 >=2 次的 walker 数
  unsigned int captured_walks = 0;             // 因近邻哑铃聚集(numSIA>=2 / MoRe 且 numRe>=5)被捕获终止的 walker 数
  unsigned int captured_no_move = 0;           // 一步未走即被捕获的 walker 数（captured 子集）
  unsigned int exited_walks = 0;               // 走出本扇区边界而终止的 walker 数
  unsigned int capped_walks = 0;               // 触发跳跃上限的 walker 数（CPU 无保险丝，恒 0）
  unsigned int stale_source_walkers = 0;       // 出发位置已不是哑铃的 walker 数
  unsigned int unclassified_walkers = 0;       // 未归入任何终止类别的 walker 数（非 0 即统计/游走 bug）
  unsigned long long transaction_failures = 0; // 锁事务失败次数（CPU 串行无锁事务，恒 0）
  unsigned int roster_live_count = 0;          // 游走结束后存活对数（花名册口径 = momo+more 哈希大小）
  unsigned int hash_pair_count = 0;            // 哈希中的对数（CPU 哈希即花名册，与上项恒等）
  // ---- detail 级（对应 GPU [recb-detail]）----
  unsigned long long target_mo = 0, target_re = 0, target_other = 0; // 随机目标位类型分布
  unsigned int roster_type_mismatches = 0;     // 出发位置类型与对身份不符次数（= stale_source）
  unsigned int direction_counts[14] = {};      // 14 个随机跳跃方向的选取次数（均匀性检验）
  // ---- 四类操作（对应 GPU [recb-ops]）：0=MoMo->Mo, 1=MoMo->Re, 2=MoRe->Mo, 3=MoRe->Re ----
  unsigned long long op_attempts[4] = {};
  unsigned long long op_successes[4] = {};
  unsigned long long op_failures[4] = {};      // CPU 串行提交恒成功，恒 0
};

class ABVIModel : public ModelAdapter<event::SelectedEvent> {

public:
  /**
   * \brief initialize kmc with simulation box.
   */
  explicit ABVIModel(Box *box, double v, double T);

  /**
   * \brief calculate the transition rates of each defect
   *  in a region at current process
   *
   * The calculating methods depends on the lattice type(single atom, vacancy and dumbbell).
   * Different types of lattice have different methods or formulas to calculate the rate.
   * see the implementation for more details.
   *
   * \param region a region, this function will calculate the transition rates in this region
   * \note the x lattice size is not doubled in \param region parameter
   * After this step, the rate of every transition direction of each lattice will be set.
   * \return return the sum of rates of all KMC events in this region,
   *  including dumbbell transition and vacancy transition and defect generation.
   */
  _type_rate calcRates(const lat_region region) override; // todo

  _type_rate calcRatesGPU(const lat_region region, int sect) override;
  // _type_rate calcRatesGPU(const lat_region region, bool *sector_first) override; // todo

  /*
   * \brief select an event randomly from rates list in a given region.  从给定区域的机率列表中随机选择一个事件。
   *
   * \param excepted_rate which equals to total rate* random number between 0-1. 等于总机率*0-1之间的随机数。
   * \param total_rates the sum rates 机率之和
   * \note the x lattice size is not doubled in \param region parameter \param region参数中的x晶格大小没有加倍
   * \return the selected event.  返回所选事件。
   */
  event::SelectedEvent select(const lat_region region, const _type_rate excepted_rate,
                              const _type_rate sum_rates) override;

  /**
   * \brief perform the selected KMC event.
   *
   */
  void perform(const event::SelectedEvent event, const lat_region region, int rank, _type_lattice_count step, int sect, const unsigned int sector_id, const unsigned int next_sector_id) override;

  void selectAndPerformOnGPU(const _type_rate rate, int rank, _type_lattice_count step, int sect, const unsigned int sector_id, const unsigned int next_sector_id) override;

  void recb_checki(const lat_region region, const unsigned int sector_id) override;

  void recb_solver(_type_lattice_id id, const lat_region& region, const unsigned int& sector_id) override;

  void reindex(const lat_region region) override;

  /**
   * \brief set kmc event listener.
   * \param p_listener pointer to the event listener.
   */
  void setEventListener(EventListener *p_listener);

  void setColoredDomain(comm::ColoredDomain *_p_domain);

  void addExchange_ghost(const _type_lattice_id& latticeId, const unsigned int sector_id);

  void addExchange_surface(_type_lattice_id surface_id);

  void clear_exchange_ghost() override;

  void clear_exchange_surface(const unsigned int  next_sect) override;

  unsigned long defectSize() override;

  void set_ghost_region() override;
  
protected:
  double time = 0;

public:
  Box *box = nullptr; // todo init box pointer
  /**
   * \brief [recb-verify] CPU 版 recb 验证统计（每次 recb launch 由 recb_checki 重置），
   * 字段口径与 akmc GPU 版 RecbVerifyStats 一致。
   */
  RecbVerifyStatsCpu recb_stats;
  /**
   * \brief recb 专用随机数生成器：与事件选择（全局 r::random()）、时间推进（time_inc_rng）完全独立。
   * recb_solver 的方向选择与 MoRe 接受/拒绝判断都从这里取随机数，保证 recb 步数变化不会挤占事件/时间流。
   */
  r::type_rng recb_rng;
  std::uniform_real_distribution<double> recb_dis{0.0, 1.0};
  void setRecbSeed(uint32_t seed) { recb_rng.seed(seed); }
  double recbRand() { return recb_dis(recb_rng); }
  /**
   * \brief attempt frequency.
   */
  const double v;
  /**
   * \brief temperature
   */
  const double T; 

  /**
   * \brief pointer to event listener.
   * event callback function will be called when executing a kmc event.
   */
  EventListener *p_event_listener = nullptr;

  /**
   * \brief it returns the rate of defect generation.
   * \return
   */
  _type_rate defectGenRate();

  // const comm::ColoredDomain *p_domain = nullptr;


};

#endif // MISA_KMC_KMC_H
