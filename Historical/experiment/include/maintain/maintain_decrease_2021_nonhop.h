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
                        ThreadPool &pool_dynamic, std::vector<std::future<int> > &results_dynamic, int time) const;

    private:
        void ProDecreasep_batch(graph<weight_type> &instance_graph,
                                std::vector<std::vector<two_hop_label<hop_weight_type> > > *L, PPR_TYPE::PPR_type *PPR,
                                std::vector<affected_label> &CL_curr, std::vector<affected_label> *CL_next,
                                ThreadPool &pool_dynamic, std::vector<std::future<int> > &results_dynamic,
                                int time) const;
    };

    template<typename weight_type, typename hop_weight_type>
    void
    StrategyA2021NonHopDecrease<weight_type, hop_weight_type>::operator()(graph<weight_type> &instance_graph,
                                                                          two_hop_case_info<hop_weight_type> &mm,
                                                                          std::vector<std::pair<int, int> > &v,
                                                                          std::vector<weight_type> &w_new,
                                                                          ThreadPool &pool_dynamic,
                                                                          std::vector<std::future<int> > &
                                                                          results_dynamic,
                                                                          int time) const {
        int current_tid = Qid_595.front();
        auto &counter = experiment::result::global_csv_config.old_counter;
        auto &shard = counter.get_thread_maintain_shard(current_tid);
        std::map<std::pair<int, int>, weight_type> w_new_map;
        int batch_size = v.size();
        for (int i = 0; i < batch_size; i++) {
            if (v[i].first > v[i].second) {
                std::swap(v[i].first, v[i].second);
            }
            if (w_new_map.count(v[i]) == 0) {
                w_new_map[v[i]] = w_new[i];
            } else if (w_new_map[v[i]] > w_new[i]) {
                w_new_map[v[i]] = w_new[i];
            }
        }

        std::vector<affected_label> CL_curr, CL_next;

        auto &L = mm.L;
        /*
        the following part does not suit parallel computation:
        the reason is that L is changed below, and as a result, in each following loop, L[v2] or L[v1] is locked at each step,
        which means that following loops cannot be actually parallized
        */

        for (auto &w_new_item: w_new_map) {
            int v1 = w_new_item.first.first, v2 = w_new_item.first.second;
            int w_new = w_new_item.second;
            for (int sl = 0; sl < 2; sl++) {
                if (sl == 1) {
                    std::swap(v1, v2);
                }
                for (auto it: L[v1]) {
                    int v = it.vertex;
                    long long dis = it.distance + w_new;
                    if (v <= v2) {
                        auto query_result = graph_weighted_two_hop_extract_distance_and_hub_in_current_with_csv(L[v],
                            L[v2],
                            v, v2,
                            shard); // query_result is {distance, common hub}
                        if (query_result.first > dis) {
                            mtx_595[v2].lock();
                            insert_sorted_two_hop_label_with_csv(L[v2], v, dis, time, shard);
                            mtx_595[v2].unlock();
                            CL_curr.emplace_back(v2, v, dis);
                            if (status::currentTimeMode == status::MaintainTimeMode::SLOT1) {
                                shard.diffuse_count_slot1++;
                            } else {
                                shard.diffuse_count_slot2++;
                            }
                        } else {
                            auto search_result = search_sorted_two_hop_label_in_current_with_csv(L[v2], v, shard);
                            if (search_result.first < 1e7 && search_result.first > dis) {
                                mtx_595[v2].lock();
                                insert_sorted_two_hop_label_with_csv(L[v2], search_result.second, dis, time,
                                                                     shard);
                                mtx_595[v2].unlock();
                                // L[v2][search_result.second].distance = dis;
                                CL_curr.emplace_back(v2, v, dis);
                                if (status::currentTimeMode == status::MaintainTimeMode::SLOT1) {
                                    shard.diffuse_count_slot1++;
                                } else {
                                    shard.diffuse_count_slot2++;
                                }
                            }
                            if (query_result.second != v) {
                                PPR_TYPE::PPR_insert_with_csv(&(mm.PPR), v2, query_result.second, v, shard);
                            }
                            if (query_result.second != v2) {
                                PPR_TYPE::PPR_insert_with_csv(&(mm.PPR), v, query_result.second, v2, shard);
                            }
                        }
                    }
                }
            }
        }
        while (!CL_curr.empty()) {
            // std::cout << "2021 decrease cl_curr size is " << CL_curr.size() << std::endl;
            ProDecreasep_batch(instance_graph, &mm.L, &mm.PPR, CL_curr, &CL_next, pool_dynamic, results_dynamic, time);
            CL_curr = CL_next;
            std::vector<affected_label>().swap(CL_next);
        }
    }

    template<typename weight_type, typename hop_weight_type>
    inline void
    StrategyA2021NonHopDecrease<weight_type, hop_weight_type>::ProDecreasep_batch(graph<weight_type> &instance_graph,
        std::vector<std::vector<two_hop_label<hop_weight_type> > > *L,
        PPR_TYPE::PPR_type *PPR,
        std::vector<affected_label> &CL_curr,
        std::vector<affected_label> *CL_next,
        ThreadPool &pool_dynamic,
        std::vector<std::future<int> > &results_dynamic,
        int time) const {
        bool is_debug = false;
        if (CL_next->size() > 100000) {
            is_debug = true;
        }
        for (const affected_label &it: CL_curr) {
            results_dynamic.emplace_back(pool_dynamic.enqueue([time, it, L, PPR, CL_next, &instance_graph, is_debug] {
                int v = it.first, u = it.second;
                mtx_595_1.lock();
                int current_tid = Qid_595.front();
                Qid_595.pop();
                mtx_595_1.unlock();
                auto &counter = experiment::result::global_csv_config.old_counter;
                auto &shard = counter.get_thread_maintain_shard(current_tid);
                mtx_595[u].lock();
                std::vector<two_hop_label<hop_weight_type> > Lu = (*L)[u]; // to avoid interlocking
                //std::cout << "lu address is " << &Lu << std::endl;
                //std::cout << "L[u] address is " << &((*(L))[u]) << std::endl;
                mtx_595[u].unlock();

                for (auto nei: instance_graph[v]) {
                    int vnei = nei.first;
                    long long dnew = it.dis + nei.second;

                    if (u < vnei) {
                        mtx_595[vnei].lock();
                        auto query_result = graph_weighted_two_hop_extract_distance_and_hub_in_current_with_csv(
                            (*L)[vnei], Lu, vnei, u, shard); // query_result is {distance, common hub}
                        mtx_595[vnei].unlock();
                        if (query_result.first > dnew) {
                            mtx_595[vnei].lock();
                            //if (is_debug) {
                            //	std::cout << "2021 decrease function two hop label1 is from " << query_result.first.t_s << " cost " << query_result.first.distance
                            //		<< " hop label2 is from " << query_result.second.t_s << " cost " << query_result.second.distance
                            //		<< " and all time is " << query_result.first.distance + query_result.second.distance << " greater than " << dnew << std::endl;
                            //}

                            insert_sorted_two_hop_label_with_csv((*L)[vnei], u, dnew, time, shard);
                            mtx_595[vnei].unlock();
                            mtx_595_1.lock();
                            CL_next->emplace_back(vnei, u, dnew);
                            mtx_595_1.unlock();
                            if (status::currentTimeMode == status::MaintainTimeMode::SLOT1) {
                                shard.diffuse_count_slot1++;
                            } else {
                                shard.diffuse_count_slot2++;
                            }
                        } else {
                            mtx_595[vnei].lock();
                            std::pair<int, int> search_result = search_sorted_two_hop_label_in_current_with_csv(
                                (*L)[vnei], u, shard);
                            mtx_595[vnei].unlock();
                            if (search_result.first < 1e7 && search_result.first > dnew) {
                                mtx_595[vnei].lock();
                                // std::cout << "decrease label has better answer : old label is " << vnei << " to " << search_result.vertex << " old value is " << search_result.distance << " to " << dnew << " t_s is " << search_result.t_s << std::endl;
                                insert_sorted_two_hop_label_with_csv((*L)[vnei], search_result.second,
                                                                     dnew, time, shard);
                                // (*L)[vnei][search_result.second].distance = dnew;
                                mtx_595[vnei].unlock();
                                mtx_595_1.lock();
                                CL_next->emplace_back(vnei, u, dnew);
                                mtx_595_1.unlock();
                                if (status::currentTimeMode == status::MaintainTimeMode::SLOT1) {
                                    shard.diffuse_count_slot1++;
                                } else {
                                    shard.diffuse_count_slot2++;
                                }
                            }
                            if (query_result.second != u) {
                                mtx_5952[vnei].lock();
                                PPR_TYPE::PPR_insert_with_csv(PPR, vnei, query_result.second, u, shard);
                                mtx_5952[vnei].unlock();
                            }
                            if (query_result.second != vnei) {
                                mtx_5952[u].lock();
                                PPR_TYPE::PPR_insert_with_csv(PPR, u, query_result.second, vnei, shard);
                                mtx_5952[u].unlock();
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
}
