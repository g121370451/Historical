#pragma once

#include <map>
#include <utility>
#include "utils/ThreadPool.h"
#include "entity/two_hop_label.h"
#include "utils/ExecutionTimer.h"
#include "entity/graph.h"
#include "entity/two_hop_label.h"
#include "entity/nonhop_global_params.h"

namespace experiment::nonhop::algorithm2021::decrease {
    template<typename weight_type, typename hop_weight_type>
    class StrategyA2021NonHopDecrease {
    public:
        void operator()(graph<weight_type> &instance_graph, two_hop_case_info<hop_weight_type> &mm,
                        std::vector<std::pair<int, int> > &v, std::vector<weight_type> &w_new,
                        ThreadPool &pool_dynamic, std::vector<std::future<int> > &results_dynamic, int time);
#ifdef _DEBUG
        std::vector<affected_label<hop_weight_type>> CL_globals;
#endif
        std::vector<record_in_increase<hop_weight_type>> list;
    private:
        void ProDecreasep_batch(graph<weight_type> &instance_graph,
                                std::vector<std::vector<two_hop_label<hop_weight_type> > > *L, PPR_TYPE::PPR_type *PPR,
                                std::vector<affected_label<hop_weight_type>> &CL_curr, std::vector<affected_label<hop_weight_type>> *CL_next,
                                ThreadPool &pool_dynamic, std::vector<std::future<int> > &results_dynamic,
                                int time);
    };

    template<typename weight_type, typename hop_weight_type>
    inline void
    StrategyA2021NonHopDecrease<weight_type, hop_weight_type>::ProDecreasep_batch(graph<weight_type> &instance_graph,
                                                                                  std::vector<std::vector<two_hop_label<hop_weight_type> > > *L,
                                                                                  PPR_TYPE::PPR_type *PPR,
                                                                                  std::vector<affected_label<hop_weight_type>> &CL_curr,
                                                                                  std::vector<affected_label<hop_weight_type>> *CL_next,
                                                                                  ThreadPool &pool_dynamic,
                                                                                  std::vector<std::future<int> > &results_dynamic,
                                                                                  int time) {
        for (const affected_label<hop_weight_type> &it: CL_curr) {
            results_dynamic.emplace_back(pool_dynamic.enqueue([time, it, L, PPR, CL_next, &instance_graph,this] {
                mtx_595_1.lock();
                int current_tid = Qid_595.front();
                Qid_595.pop();
                mtx_595_1.unlock();
                auto &counter = experiment::result::global_csv_config.old_counter;
                auto &shard = counter.get_thread_maintain_shard(current_tid);

                const int v = it.first;
                int u = it.second;

                L_lock[u].lock();
                std::vector<two_hop_label<hop_weight_type> > Lu = (*L)[u]; // to avoid interlocking
                L_lock[u].unlock();

                for (auto nei: instance_graph[v]) {
                    int vnei = nei.first;
                    long long dnew = it.dis + nei.second;
#ifdef _DEBUG
                    if (dnew < 0) {
                        std::cout << "overflow happen in maintain decrease 2021 with nonhop" << " it_dis is "
                                  << it.dis
                                  << " nei.second is " << nei.second << std::endl;
                    }
#endif
                    if (u < vnei) {
                        L_lock[vnei].lock();
                        auto [query_dis, query_hub] = graph_weighted_two_hop_extract_distance_and_hub_in_current_with_csv(
                                (*L)[vnei], Lu, vnei, u, shard); // query_result is {distance, common hub}
                        L_lock[vnei].unlock();
                        if (query_dis > dnew) {
                            L_lock[vnei].lock();
                            insert_sorted_two_hop_label_with_csv((*L)[vnei], u, dnew, time, shard);
                            L_lock[vnei].unlock();
                            mtx_list_check.lock();
                            this->list.emplace_back(vnei, u, dnew, query_dis, time);
                            mtx_list_check.unlock();
                            mtx_595_1.lock();
                            CL_next->emplace_back(vnei, u, dnew);
                            if (status::currentTimeMode == status::MaintainTimeMode::SLOT1) {
                                shard.diffuse_count_slot1++;
                            } else {
                                shard.diffuse_count_slot2++;
                            }
                            mtx_595_1.unlock();
                        } else {
                            L_lock[vnei].lock();
                            hop_weight_type search_result = search_sorted_two_hop_label_in_current_with_csv(
                                    (*L)[vnei], u, shard);
                            L_lock[vnei].unlock();
                            if (search_result < std::numeric_limits<hop_weight_type>::max() && search_result > dnew) {
                                L_lock[vnei].lock();
                                insert_sorted_two_hop_label_with_csv((*L)[vnei], u,
                                                                     dnew, time, shard);
                                L_lock[vnei].unlock();
                                mtx_list_check.lock();
                                this->list.emplace_back(vnei, u, dnew, query_dis, time);
                                mtx_list_check.unlock();
                                mtx_595_1.lock();
                                CL_next->emplace_back(vnei, u, dnew);
                                if (status::currentTimeMode == status::MaintainTimeMode::SLOT1) {
                                    shard.diffuse_count_slot1++;
                                } else {
                                    shard.diffuse_count_slot2++;
                                }
                                mtx_595_1.unlock();
                            }
                            if (query_hub != -1 && query_hub != u) {
                                ppr_lock[vnei].lock();
                                PPR_TYPE::PPR_insert_with_csv(PPR, vnei, query_hub, u, shard);
                                ppr_lock[vnei].unlock();
                            }
                            if (query_hub != -1 && query_hub != vnei) {
                                ppr_lock[u].lock();
                                PPR_TYPE::PPR_insert_with_csv(PPR, u, query_hub, vnei, shard);
                                ppr_lock[u].unlock();
                            }
                        }
                    }
                }
                mtx_595_1.lock();
                Qid_595.push(current_tid);
                mtx_595_1.unlock();
                return 1;
            }));
        }

        for (auto &&result: results_dynamic) {
            result.get();
        }
        std::vector<std::future<int> >().swap(results_dynamic);
    }

    template<typename weight_type, typename hop_weight_type>
    void
    StrategyA2021NonHopDecrease<weight_type, hop_weight_type>::operator()(graph<weight_type> &instance_graph,
                                                                          two_hop_case_info<hop_weight_type> &mm,
                                                                          std::vector<std::pair<int, int> > &v,
                                                                          std::vector<weight_type> &w_new,
                                                                          ThreadPool &pool_dynamic,
                                                                          std::vector<std::future<int> > &
                                                                          results_dynamic,
                                                                          int time) {
        int current_tid = Qid_595.front();
        auto &counter = experiment::result::global_csv_config.old_counter;
        auto &shard = counter.get_thread_maintain_shard(current_tid);
        std::map<std::pair<int, int>, weight_type> w_new_map;
        size_t batch_size = v.size();
        for (size_t i = 0; i < batch_size; i++) {
            if (v[i].first > v[i].second) {
                std::swap(v[i].first, v[i].second);
            }
            if (w_new_map.count(v[i]) == 0 || w_new_map[v[i]] > w_new[i]) {
                w_new_map[v[i]] = w_new[i];
            }
        }

        std::vector<affected_label<hop_weight_type>> CL_curr, CL_next;

        std::vector<std::vector<two_hop_label<hop_weight_type> > > &L = mm.L;
        /*
        the following part does not suit parallel computation:
        the reason is that L is changed below, and as a result, in each following loop, L[v2] or L[v1] is locked at each step,
        which means that following loops cannot be actually parallelized
        */

        for (auto &w_new_item: w_new_map) {
            int v1 = w_new_item.first.first, v2 = w_new_item.first.second;
            weight_type w_item_dis = w_new_item.second;
            for (int sl = 0; sl < 2; sl++) {
                if (sl == 1) {
                    std::swap(v1, v2);
                }
                for (two_hop_label<hop_weight_type> &it: L[v1]) {
                    int _v_item = it.vertex;
                    long long dis = it.distance + w_item_dis;
                    if (it.distance == std::numeric_limits<weight_type>::max()) {
                        continue;
                    }
#ifdef _DEBUG
                    if (dis < 0) {
                        std::cout << "overflow happen in 2021 maintain decrease step1 it dis is " << it.distance
                                  << " _w_new is " << w_item_dis << std::endl;
                    }
#endif
                    if (_v_item <= v2 && it.t_e == std::numeric_limits<int>::max()) {
                        auto [query_dis,query_hub] = graph_weighted_two_hop_extract_distance_and_hub_in_current_with_csv(
                                    L[_v_item],
                                    L[v2],
                                    _v_item, v2,
                                    shard); // query_result is {distance, common hub}
                        if (query_dis > dis) {
                            insert_sorted_two_hop_label_with_csv(L[v2], _v_item, dis, time, shard);
                            this->list.emplace_back(v2, _v_item, dis, query_dis, time);
                            CL_curr.emplace_back(v2, _v_item, dis);
                            if (status::currentTimeMode == status::MaintainTimeMode::SLOT1) {
                                shard.diffuse_count_slot1++;
                            } else {
                                shard.diffuse_count_slot2++;
                            }
                        } else {
                            auto search_result = search_sorted_two_hop_label_in_current_with_csv(L[v2], _v_item, shard);
                            if (search_result < std::numeric_limits<hop_weight_type>::max() && search_result > dis) {
                                insert_sorted_two_hop_label_with_csv(L[v2], _v_item, dis, time,
                                                                     shard);
                                this->list.emplace_back(v2, _v_item, dis, search_result, time);
                                CL_curr.emplace_back(v2, _v_item, dis);
                                if (status::currentTimeMode == status::MaintainTimeMode::SLOT1) {
                                    shard.diffuse_count_slot1++;
                                } else {
                                    shard.diffuse_count_slot2++;
                                }
                            }
                            if (query_hub != -1 && query_hub != _v_item) {
                                PPR_TYPE::PPR_insert_with_csv(&(mm.PPR), v2, query_hub, _v_item, shard);
                            }
                            if (query_hub != -1 && query_hub != v2) {
                                PPR_TYPE::PPR_insert_with_csv(&(mm.PPR), _v_item, query_hub, v2, shard);
                            }
                        }
                    }
                }
            }
        }
        while (!CL_curr.empty()) {
#ifdef _DEBUG
            CL_globals.insert(CL_globals.end(), CL_curr.begin(), CL_curr.end());
#endif
            ProDecreasep_batch(instance_graph, &mm.L, &mm.PPR, CL_curr, &CL_next, pool_dynamic, results_dynamic, time);
            CL_curr = CL_next;
            std::vector<affected_label<hop_weight_type>>().swap(CL_next);
        }
    }
}
