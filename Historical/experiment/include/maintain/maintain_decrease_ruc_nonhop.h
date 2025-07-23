#pragma once

#include <map>
#include <utility>
#include "utils/ThreadPool.h"
#include "entity/two_hop_label.h"
#include "utils/ExecutionTimer.h"
#include "entity/graph.h"
#include "entity/two_hop_label.h"
#include "entity/nonhop_global_params.h"


namespace experiment::nonhop::ruc::decrease {
    template<typename weight_type, typename hop_weight_type>
    class Strategy2024NonHopDecrease {
    public:
        void operator()(graph<weight_type> &instance_graph, two_hop_case_info<hop_weight_type> &mm,
                        std::vector<std::pair<int, int>> &v, std::vector<weight_type> &w_new,
                        ThreadPool &pool_dynamic, std::vector<std::future<int>> &results_dynamic,
                        int time);
#ifdef _DEBUG
        std::vector<affected_label<hop_weight_type>> CL_globals;
#endif
        std::vector<record_in_increase<hop_weight_type>> list;
    private:
        void decrease_maintain_step1_batch(std::map<std::pair<int, int>, weight_type> &v_map,
                                           std::vector<std::vector<two_hop_label<hop_weight_type>>> *L,
                                           PPR_TYPE::PPR_type *PPR, std::vector<affected_label<hop_weight_type>> *CL,
                                           ThreadPool &pool_dynamic,
                                           std::vector<std::future<int>> &results_dynamic);

        void DIFFUSE_batch(graph<weight_type> &instance_graph,
                           std::vector<std::vector<two_hop_label<hop_weight_type>>> *L,
                           PPR_TYPE::PPR_type *PPR, std::vector<affected_label<hop_weight_type>> &CL,
                           ThreadPool &pool_dynamic, std::vector<std::future<int>> &results_dynamic,
                           int t);
    };

    template<typename weight_type, typename hop_weight_type>
    inline void
    Strategy2024NonHopDecrease<weight_type, hop_weight_type>::operator()(graph<weight_type> &instance_graph,
                                                                         two_hop_case_info<hop_weight_type> &mm,
                                                                         std::vector<std::pair<int, int>> &v,
                                                                         std::vector<weight_type> &w_new,
                                                                         ThreadPool &pool_dynamic,
                                                                         std::vector<std::future<int>> &results_dynamic,
                                                                         int time) {
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
        std::vector<affected_label<hop_weight_type>> CL;
        decrease_maintain_step1_batch(w_new_map, &mm.L, &mm.PPR, &CL, pool_dynamic, results_dynamic);
#ifdef _DEBUG
        CL_globals.insert(CL_globals.end(), CL.begin(), CL.end());
#endif
        DIFFUSE_batch(instance_graph, &mm.L, &mm.PPR, CL, pool_dynamic, results_dynamic, time);
    }

    template<typename weight_type, typename hop_weight_type>
    inline void Strategy2024NonHopDecrease<weight_type, hop_weight_type>::decrease_maintain_step1_batch(
            std::map<std::pair<int, int>, weight_type> &v_map,
            std::vector<std::vector<two_hop_label<hop_weight_type>>> *L, PPR_TYPE::PPR_type *PPR,
            std::vector<affected_label<hop_weight_type>> *CL, ThreadPool &pool_dynamic,
            std::vector<std::future<int>> &results_dynamic) {
        for (std::pair<std::pair<int, int>, weight_type> v_item: v_map) {
            results_dynamic.emplace_back(pool_dynamic.enqueue([&v_item, L, PPR, CL] {
                mtx_595_1.lock();
                int current_tid = Qid_595.front();
                Qid_595.pop();
                mtx_595_1.unlock();

                auto &counter = experiment::result::global_csv_config.ruc_counter;
                auto &shard = counter.get_thread_maintain_shard(current_tid);

                int v1 = v_item.first.first, v2 = v_item.first.second;
                weight_type w_new = v_item.second;

                for (int sl = 0; sl < 2; sl++) {
                    if (sl == 1) {
                        std::swap(v1, v2);
                    }
                    for (auto &it: (*L)[v1]) {
                        if (it.distance != std::numeric_limits<hop_weight_type>::max()
                            && it.vertex <= v2
                            && it.t_e == std::numeric_limits<int>::max()) {
                            hop_weight_type dnew = it.distance + w_new;
#ifdef _DEBUG
                            if (dnew < 0) {
                                std::cout << "overflow happen in maintain decrease ruc with hop" << std::endl;
                            }
#endif
                            auto [query_dis,query_hub] = graph_weighted_two_hop_extract_distance_and_hub_in_current_with_csv(
                                (*L)[it.vertex], (*L)[v2], it.vertex, v2,
                                shard); // query_result is {distance, common hub}
                            if (query_dis > dnew) {
                                mtx_595_1.lock();
                                CL->push_back(affected_label{v2, it.vertex, dnew});
                                if (status::currentTimeMode == status::MaintainTimeMode::SLOT1) {
                                    shard.diffuse_count_slot1++;
                                } else {
                                    shard.diffuse_count_slot2++;
                                }
                                mtx_595_1.unlock();
                            } else {
                                hop_weight_type search_result = search_sorted_two_hop_label_in_current_with_csv(
                                        (*L)[v2], it.vertex, shard);
                                if (search_result > dnew &&
                                    search_result < std::numeric_limits<hop_weight_type>::max()) {
                                    mtx_595_1.lock();
                                    CL->push_back(affected_label{v2, it.vertex, dnew});
                                    if (status::currentTimeMode == status::MaintainTimeMode::SLOT1) {
                                        shard.diffuse_count_slot1++;
                                    } else {
                                        shard.diffuse_count_slot2++;
                                    }
                                    mtx_595_1.unlock();
                                }
                                if (query_hub != -1 && query_hub != it.vertex) {
                                    mtx_5952[v2].lock();
                                    PPR_TYPE::PPR_insert_with_csv(PPR, v2, query_hub, it.vertex,
                                                                  shard);
                                    mtx_5952[v2].unlock();
                                }
                                if (query_hub != -1 && query_hub != v2) {
                                    mtx_5952[it.vertex].lock();
                                    PPR_TYPE::PPR_insert_with_csv(PPR, it.vertex, query_hub, v2,
                                                                  shard);
                                    mtx_5952[it.vertex].unlock();
                                }
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
        std::vector<std::future<int>>().swap(results_dynamic);
    }

    template<typename weight_type, typename hop_weight_type>
    inline void Strategy2024NonHopDecrease<weight_type, hop_weight_type>::DIFFUSE_batch(
            graph<weight_type> &instance_graph, std::vector<std::vector<two_hop_label<hop_weight_type>>> *L,
            PPR_TYPE::PPR_type *PPR, std::vector<affected_label<hop_weight_type>> &CL, ThreadPool &pool_dynamic,
            std::vector<std::future<int>> &results_dynamic, int t) {
        // Deduplication
        std::map<std::pair<int, int>, hop_weight_type> CL_edge_map;
        for (auto &it: CL) {
            if (CL_edge_map.count({it.first, it.second}) == 0 || CL_edge_map[{it.first, it.second}] > it.dis) {
                CL_edge_map[{it.first, it.second}] = it.dis;
            }
        }

        // extract each unique hub v and its (u,dis) list
        std::map<int, std::vector<std::pair<int, int>>> CL_map; // CL_map[v]=(u1,dis1),(u2,dis2)...
        for (auto &CL_edgeItem: CL_edge_map) {
            int u = CL_edgeItem.first.first;
            int v = CL_edgeItem.first.second;
            hop_weight_type dis = CL_edgeItem.second;
            if (CL_map.count(v) == 0) {
                std::vector<std::pair<int, int>> vec_with_hub_v;
                vec_with_hub_v.emplace_back(u, dis);
                CL_map[v] = vec_with_hub_v;
            } else {
                std::vector<std::pair<int, int>> &vec_with_hub_v = CL_map[v];
                vec_with_hub_v.emplace_back(u, dis);
            }
        }

        std::vector<std::pair<int, std::vector<std::pair<int, int>>>> CL_map_vec(CL_map.begin(),
                                                                                 CL_map.end());
        sort(CL_map_vec.begin(), CL_map_vec.end(),
             [](const std::pair<int, std::vector<std::pair<int, int>>> &a,
                const std::pair<int, std::vector<std::pair<int, int>>> &b) { return a.first < b.first; });
        // each thread processes one unique hub
        for (auto &cl_item: CL_map_vec) {
            results_dynamic.emplace_back(pool_dynamic.enqueue([t, cl_item, L, &instance_graph, PPR, this] {
                mtx_595_1.lock();
                int current_tid = Qid_595.front();
                Qid_595.pop();
                mtx_595_1.unlock();
                auto &counter = experiment::result::global_csv_config.ruc_counter;
                auto &shard = counter.get_thread_maintain_shard(current_tid);

                int v = cl_item.first;
                std::vector<std::pair<int, int>> vec_with_hub_v = cl_item.second;

                mtx_595[v].lock();
                auto Lv = (*L)[v]; // to avoid interlocking
                mtx_595[v].unlock();

                std::vector<int> Dis_changed;
                std::vector<std::pair<hop_weight_type, int>> &DIS = Dis<hop_weight_type>[current_tid];
                std::vector<handle_t_for_DIFFUSE> &Q_HANDLES = Q_handles[current_tid];
                auto &Q_VALUE = Q_value<hop_weight_type>[current_tid];

                boost::heap::fibonacci_heap<node_for_DIFFUSE> Q;

                for (auto &vec_with_hub_v_item: vec_with_hub_v) {
                    int u = vec_with_hub_v_item.first;
                    int du = vec_with_hub_v_item.second;
                    DIS[u] = {du, v}; // <distance, hub responsible for this distance>
                    Dis_changed.push_back(u);
                    Q_HANDLES[u] = Q.push(node_for_DIFFUSE(u, du));
                    Q_VALUE[u] = du;
                }

                while (!Q.empty()) {

                    node_for_DIFFUSE temp2 = Q.top();
                    int x = temp2.index;
                    int dx = temp2.disx;
                    Q.pop();
                    Q_VALUE[x] = 1e7;

                    mtx_595[x].lock();
                    hop_weight_type d_old = search_sorted_two_hop_label_in_current_with_csv((*L)[x], v,
                                                                                            shard);
                    mtx_595[x].unlock();
                    if (d_old > dx) {
                        mtx_595[x].lock();
                        insert_sorted_two_hop_label_with_csv((*L)[x], v, dx, t, shard);
                        mtx_595[x].unlock();
                        mtx_list_check.lock();
                        this->list.emplace_back(x, v, dx, d_old, time);
                        mtx_list_check.unlock();
                    } else {
                        continue;
                    }

                    for (auto &nei: instance_graph[x]) {
                        int xnei = nei.first;
                        int d_new = dx + nei.second;
#ifdef _DEBUG
                        if (d_new < 0) {
                            std::cout << "overflow happen in maintain decrease ruc diffuse with nonhop"
                                      << " dx is "
                                      << dx << " nei.second is " << nei.second << std::endl;
                        }
#endif
                        if (v < xnei) {
                            if (DIS[xnei].first == -1) {
                                mtx_595[xnei].lock();
                                DIS[xnei] = graph_weighted_two_hop_extract_distance_and_hub_in_current_with_csv(
                                        (*L)[xnei], Lv, xnei, v, shard);
                                mtx_595[xnei].unlock();
                                Dis_changed.push_back(xnei);
                            }
                            if (DIS[xnei].first > d_new) {
                                DIS[xnei] = {d_new, v};
                                if (Q_VALUE[xnei] >= 1e7) {
                                    Q_HANDLES[xnei] = Q.push(node_for_DIFFUSE(xnei, d_new));
                                    if (status::currentTimeMode == status::MaintainTimeMode::SLOT1) {
                                        shard.diffuse_count_slot1++;
                                    } else {
                                        shard.diffuse_count_slot2++;
                                    }
                                } else {
                                    Q.update(Q_HANDLES[xnei], node_for_DIFFUSE(xnei, d_new));
                                }
                                Q_VALUE[xnei] = d_new;
                            } else {
                                mtx_595[xnei].lock();
                                hop_weight_type search_result = search_sorted_two_hop_label_in_current_with_csv(
                                        (*L)[xnei], v, shard);
                                mtx_595[xnei].unlock();
                                if (search_result < std::numeric_limits<hop_weight_type>::max() &&
                                    std::min(search_result, Q_VALUE[xnei]) > d_new) {
                                    if (Q_VALUE[xnei] >= 1e7) {
                                        Q_HANDLES[xnei] = Q.push(node_for_DIFFUSE(xnei, d_new));
                                        if (status::currentTimeMode == status::MaintainTimeMode::SLOT1) {
                                            shard.diffuse_count_slot1++;
                                        } else {
                                            shard.diffuse_count_slot2++;
                                        }
                                    } else {
                                        Q.update(Q_HANDLES[xnei], node_for_DIFFUSE(xnei, d_new));
                                    }
                                    Q_VALUE[xnei] = d_new;
                                }
                                if (DIS[xnei].second != v) {
                                    mtx_5952[xnei].lock();
                                    PPR_TYPE::PPR_insert_with_csv(PPR, xnei, DIS[xnei].second, v, shard);
                                    mtx_5952[xnei].unlock();
                                }
                                if (DIS[xnei].second != xnei) {
                                    mtx_5952[v].lock();
                                    PPR_TYPE::PPR_insert_with_csv(PPR, v, DIS[xnei].second, xnei, shard);
                                    mtx_5952[v].unlock();
                                }
                            }
                        }
                    }
                }
                for (int i: Dis_changed) {
                    DIS[i] = {-1, -1};
                }
                // Q_HANDLES.clear();
                // Q_HANDLES.resize(Q_HANDLES.size());
                Q_VALUE.resize(Q_VALUE.size(), 1e7);
                mtx_595_1.lock();
                Qid_595.push(current_tid);
                mtx_595_1.unlock();
                return 1;
            }));
        }

        for (auto &&result: results_dynamic) {
            result.get();
        }
        std::vector<std::future<int>>().swap(results_dynamic);
    }

}



