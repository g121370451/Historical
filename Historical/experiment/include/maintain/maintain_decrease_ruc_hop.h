#pragma once

#include "entity/graph.h"
#include "utils/ThreadPool.h"
#include "entity/two_hop_label.h"
#include "entity/hop_global_params.h"
#include <map>

namespace experiment::hop::ruc::decrease {
    template<typename weight_type, typename hop_weight_type>
    class Strategy2024HopDecrease {
    public:
        void operator()(graph<weight_type> &instance_graph, two_hop_case_info<hop_weight_type> &mm,
                        std::vector<std::pair<int, int> > &v, std::vector<weight_type> &w_new,
                        ThreadPool &pool_dynamic, std::vector<std::future<int> > &results_dynamic, int time);
#ifdef _DEBUG
        std::vector<hop_constrained_affected_label<hop_weight_type>> CL_globals;
#endif
        std::vector<record_in_increase_with_hop<hop_weight_type>> list;
    private:

        void decrease_maintain_step1_batch(std::map<std::pair<int, int>, hop_weight_type> &v_map,
                                           std::vector<std::vector<two_hop_label<hop_weight_type> > > *L,
                                           PPR_TYPE::PPR_type *PPR,
                                           std::vector<hop_constrained_affected_label<hop_weight_type> > *CL,
                                           ThreadPool &pool_dynamic, std::vector<std::future<int> > &results_dynamic,
                                           int upper_k);

        void DIFFUSE_batch(graph<weight_type> &instance_graph,
                           std::vector<std::vector<two_hop_label<hop_weight_type> > > *L, PPR_TYPE::PPR_type *PPR,
                           std::vector<hop_constrained_affected_label<hop_weight_type> > &CL,
                           ThreadPool &pool_dynamic, std::vector<std::future<int> > &results_dynamic, int upper_k,
                           int time);
    };
    template<typename weight_type, typename hop_weight_type>
    void Strategy2024HopDecrease<weight_type, hop_weight_type>::operator()(graph<weight_type> &instance_graph,
                                                                           two_hop_case_info<hop_weight_type> &mm,
                                                                           std::vector<std::pair<int, int> > &v,
                                                                           std::vector<weight_type> &w_new,
                                                                           ThreadPool &pool_dynamic,
                                                                           std::vector<std::future<int> > &
                                                                           results_dynamic, int time) {
        std::map<std::pair<int, int>, hop_weight_type> w_new_map;
        const size_t batch_size = v.size();
        for (size_t i = 0; i < batch_size; i++) {
            if (v[i].first > v[i].second) {
                std::swap(v[i].first, v[i].second);
            }
            if (!w_new_map.contains(v[i]) || w_new_map[v[i]] > w_new[i]) {
                w_new_map[v[i]] = w_new[i];
            }
        }
        std::vector<hop_constrained_affected_label<hop_weight_type> > CL;
        // 每次传入的要是和线程数大小相等
        decrease_maintain_step1_batch(w_new_map, &mm.L, &mm.PPR, &CL, pool_dynamic, results_dynamic, mm.upper_k);
#ifdef _DEBUG
        CL_globals.insert(CL_globals.end(), CL.begin(), CL.end());
#endif
        DIFFUSE_batch(instance_graph, &mm.L, &mm.PPR, CL, pool_dynamic, results_dynamic, mm.upper_k, time);
    }
    template<typename weight_type, typename hop_weight_type>
    void Strategy2024HopDecrease<weight_type, hop_weight_type>::decrease_maintain_step1_batch(
            std::map<std::pair<int, int>, hop_weight_type> &v_map,
            std::vector<std::vector<two_hop_label<hop_weight_type> > > *L, PPR_TYPE::PPR_type *PPR,
            std::vector<hop_constrained_affected_label<hop_weight_type> > *CL, ThreadPool &pool_dynamic,
            std::vector<std::future<int> > &results_dynamic, int upper_k) {
        for (const auto &v_map_item: v_map) {
            results_dynamic.emplace_back(pool_dynamic.enqueue([&v_map_item, L, PPR, CL, upper_k] {
                mtx_599_1.lock();
                int current_tid = Qid_599.front();
                Qid_599.pop();
                mtx_599_1.unlock();

                auto &counter = experiment::result::global_csv_config.ruc_counter;
                auto &shard = counter.get_thread_maintain_shard(current_tid);

                int v1 = v_map_item.first.first, v2 = v_map_item.first.second;
                hop_weight_type w_new = v_map_item.second;
                for (int sl = 0; sl < 2; sl++) {
                    if (sl == 1) {
                        std::swap(v1, v2);
                    }
                    for (auto &it: (*L)[v1]) {
                        if (it.distance != std::numeric_limits<hop_weight_type>::max()
                            && it.hub_vertex <= v2
                            && it.hop + 1 < upper_k
                            && it.t_e == std::numeric_limits<int>::max()) {
                            hop_weight_type dnew = it.distance + w_new;
#ifdef _DEBUG
                            if (dnew < 0) {
                                std::cout << "overflow happen in maintain decrease ruc with hop" << std::endl;
                            }
#endif
                            auto [query_dis, query_hop, query_hub] =
                                    graph_weighted_two_hop_extract_distance_and_hop_and_hub_in_current_with_csv(
                                            (*L)[it.hub_vertex], (*L)[v2], it.hub_vertex, v2, it.hop + 1,
                                            shard); // query_result is {distance, common hub}
                            if (query_dis > dnew) {
                                mtx_599_1.lock();
                                CL->push_back(hop_constrained_affected_label<hop_weight_type>{
                                        v2, it.hub_vertex, it.hop + 1, dnew
                                });
                                if (status::currentTimeMode == status::MaintainTimeMode::SLOT1) {
                                    shard.diffuse_count_slot1++;
                                } else {
                                    shard.diffuse_count_slot2++;
                                }
                                mtx_599_1.unlock();
                            } else {
                                auto [label_dis, label_hop] =
                                        search_sorted_two_hop_label_in_current_with_less_than_k_limit_with_csv(
                                                (*L)[v2], it.hub_vertex, it.hop + 1, shard);
                                if (label_dis < std::numeric_limits<hop_weight_type>::max() && label_dis > dnew) {
                                    mtx_599_1.lock();
                                    CL->push_back(hop_constrained_affected_label<hop_weight_type>{
                                            v2, it.hub_vertex, it.hop + 1, dnew
                                    });
                                    if (status::currentTimeMode == status::MaintainTimeMode::SLOT1) {
                                        shard.diffuse_count_slot1++;
                                    } else {
                                        shard.diffuse_count_slot2++;
                                    }
                                    mtx_599_1.unlock();
                                }
                                if (query_hub != -1 && query_hub != it.hub_vertex) {
                                    ppr_lock[v2].lock();
                                    PPR_TYPE::PPR_insert_with_csv(PPR, v2, query_hub, it.hub_vertex);
                                    ppr_lock[v2].unlock();
                                }
                                if (query_hub != -1 && query_hub != v2) {
                                    ppr_lock[it.hub_vertex].lock();
                                    PPR_TYPE::PPR_insert_with_csv(PPR, it.hub_vertex, query_hub, v2);
                                    ppr_lock[it.hub_vertex].unlock();
                                }
                            }
                        }
                    }
                }
                mtx_599_1.lock();
                Qid_599.push(current_tid);
                mtx_599_1.unlock();
                return 1;
            }));
        }

        for (auto &&result: results_dynamic) {
            result.get();
        }
        std::vector<std::future<int> >().swap(results_dynamic);
    }

    template<typename weight_type, typename hop_weight_type>
    void Strategy2024HopDecrease<weight_type, hop_weight_type>::DIFFUSE_batch(graph<weight_type> &instance_graph,
                                                                              std::vector<std::vector<two_hop_label<
                                                                                      hop_weight_type> > > *L,
                                                                              PPR_TYPE::PPR_type *PPR,
                                                                              std::vector<hop_constrained_affected_label
                                                                                      <hop_weight_type> > &CL,
                                                                              ThreadPool &pool_dynamic,
                                                                              std::vector<std::future<int> > &
                                                                              results_dynamic,
                                                                              int upper_k, int time) {
        std::map<hop_constrained_pair_label, hop_weight_type> CL_edge_map;
        for (auto &it: CL) {
            if (!CL_edge_map.contains({it.first, it.second, it.hop}) ||
                CL_edge_map[{it.first, it.second, it.hop}] > it.dis) {
                CL_edge_map[{it.first, it.second, it.hop}] = it.dis;
            }
        }

        // extract each unique hub v and its (u,hop,dis) list
        std::map<int, std::vector<hop_constrained_label_v2<hop_weight_type> > > CL_map;
        for (auto &it: CL_edge_map) {
            int u = it.first.first;
            int v = it.first.second;
            int hop = it.first.hop;
            hop_weight_type dis = it.second;
            if (!CL_map.contains(v)) {
                std::vector<hop_constrained_label_v2<hop_weight_type> > vec_with_hub_v;
                hop_constrained_label_v2 tmp(u, hop, dis);
                vec_with_hub_v.emplace_back(tmp);
                CL_map[v] = vec_with_hub_v;
            } else {
                std::vector<hop_constrained_label_v2<hop_weight_type> > &vec_with_hub_v = CL_map[v];
                hop_constrained_label_v2 tmp(u, hop, dis);
                vec_with_hub_v.emplace_back(tmp);
            }
        }
        long long size = 0;
        for (auto &cl_item: CL_map) {
            // results_dynamic.emplace_back(pool_dynamic.enqueue([time, &cl_item, L, &instance_graph, PPR, upper_k] {
            results_dynamic.emplace_back(
                    pool_dynamic.enqueue([time, &cl_item, L, &instance_graph, PPR, upper_k, &size, this] {
                        mtx_599_1.lock();
                        int current_tid = Qid_599.front();
                        Qid_599.pop();
                        mtx_599_1.unlock();
                        auto &counter = experiment::result::global_csv_config.ruc_counter;
                        auto &shard = counter.get_thread_maintain_shard(current_tid);

                        int v = cl_item.first;
                        std::vector<hop_constrained_label_v2<hop_weight_type>> &vec_with_hub_v = cl_item.second;

                        L_lock[v].lock();
                        auto Lv = (*L)[v]; // to avoid interlocking
                        L_lock[v].unlock();

                        std::vector<std::pair<int, int>> dist_hop_changes;
                        std::vector<std::vector<hop_weight_type>> &dist_hop = dist_hop_599_v2<hop_weight_type>[current_tid];
                        boost::heap::fibonacci_heap<hop_constrained_node_for_DIFFUSE<hop_weight_type> > pq;
                        std::map<std::pair<int, int>, std::pair<hop_constrained_handle_t_for_DIFFUSE<hop_weight_type>,
                                hop_weight_type> > Q_handle;
                        std::vector<int> hubs(instance_graph.size(), -1);
                        auto &Q_VALUE = Q_value<hop_weight_type>[current_tid];

                        for (auto &it: vec_with_hub_v) {
                            int u = it.hub_vertex;
                            int h_v = it.hop;
                            hop_weight_type du = it.distance;
                            L_lock[u].lock();
                            auto [query_dis, query_hop, query_hub] =
                                    graph_weighted_two_hop_extract_distance_and_hop_and_hub_in_current_with_csv(
                                            (*L)[u], Lv,
                                            u, v,
                                            h_v, shard);
                            L_lock[u].unlock();
                            if (query_dis <= du && query_dis != std::numeric_limits<hop_weight_type>::max()) {
                                if (query_hub != -1 && query_hub != v) {
                                    ppr_lock[u].lock();
                                    PPR_TYPE::PPR_insert(PPR, u, query_hub, v);
                                    ppr_lock[u].unlock();
                                }
                                if (query_hub != -1 && query_hub != u) {
                                    ppr_lock[v].lock();
                                    PPR_TYPE::PPR_insert(PPR, v, query_hub, u);
                                    ppr_lock[v].unlock();
                                }
                            } else {
                                hop_constrained_node_for_DIFFUSE<hop_weight_type> tmp;
                                tmp.index = u;
                                tmp.hop = h_v;
                                tmp.disx = du;
                                Q_handle[{u, h_v}] = {pq.push({tmp}), du};
                                Q_VALUE[u][h_v] = du;
                            }

                        }
                        while (!pq.empty()) {
                            int x = pq.top().index;
                            int xhv = pq.top().hop;
                            hop_weight_type dx = pq.top().disx;
                            ++size;
                            pq.pop();
                            if (Q_VALUE[x][xhv - 1] != -1 && Q_VALUE[x][xhv - 1] <= dx) {
                                continue;
                            }

                            L_lock[x].lock();
                            auto [label_dis, label_hop] = search_sorted_two_hop_label_in_current_with_less_than_k_limit_with_csv(
                                    (*L)[x], v, xhv, shard);
                            L_lock[x].unlock();
                            if (dx >= 0 && dx < label_dis) {
                                L_lock[x].lock();
                                insert_sorted_hop_constrained_two_hop_label_with_csv((*L)[x], v, xhv, dx, time, shard);
                                L_lock[x].unlock();
                                Q_VALUE[x][xhv] = dx;
                                mtx_list_check.lock();
                                this->list.emplace_back(x, v, xhv, dx, label_dis, time);
                                mtx_list_check.unlock();
                            } else {
                                continue;
                            }
                            if (xhv + 1 > upper_k)
                                continue;

                            for (const auto &nei: instance_graph[x]) {
                                int xnei = nei.first;
                                if (v < xnei) {
                                    int hop_nei = xhv + 1;
                                    hop_weight_type d_new = dx + static_cast<hop_weight_type>(nei.second);
#ifdef _DEBUG
                                    if (d_new < 0) {
                                        std::cout << "overflow happen in maintain decrease ruc diffuse with hop"
                                                  << " dx is "
                                                  << dx << " nei.second is " << nei.second << std::endl;
                                    }
#endif

                                    hop_constrained_node_for_DIFFUSE<hop_weight_type> node = {xnei, hop_nei, d_new};
                                    if (dist_hop[xnei][hop_nei] == -1) {
                                        // Q_handle[{xnei, hop_nei}] = {pq.push(node), d_new};
                                        // Q_VALUE[xnei][hop_nei] = d_new;
                                        L_lock[xnei].lock();
                                        auto allDistances =
                                                graph_weighted_two_hop_extract_all_distances_under_hop_limit(
                                                        (*L)[xnei], Lv, xnei, v, upper_k, shard);
                                        //std::pair<int, int> temp_dis = hop_constrained_extract_distance_and_hop(*L, xnei, v, xhv + 1);
                                        L_lock[xnei].unlock();
                                        for (int i = 1; i <= upper_k; ++i) {
                                            auto [query_dis, query_hop, query_hub] = allDistances[i];
                                            if (query_dis != std::numeric_limits<hop_weight_type>::max()) {
                                                hubs[xnei] = query_hub;
                                            }
                                            dist_hop[xnei][i] = query_dis;
                                            dist_hop_changes.emplace_back(xnei, i);
                                        }

                                    }

                                    if (d_new < dist_hop[xnei][upper_k]) {
                                        //if (Q_handle.find({xnei, hop_nei}) != Q_handle.end())
                                        if (Q_handle.contains({xnei, hop_nei})) {
                                            if (Q_handle[{xnei, hop_nei}].second > d_new) {
                                                pq.update(Q_handle[{xnei, hop_nei}].first, node);
                                                Q_handle[{xnei, hop_nei}].second = d_new;
                                            }
                                        } else {
                                            Q_handle[{xnei, hop_nei}] = {pq.push(node), d_new};
                                        }
                                        for (int index = hop_nei; index <= upper_k; index++) {
                                            dist_hop[xnei][hop_nei] = d_new;
                                        }
                                        hubs[xnei] = v;
                                    } else if (d_new < dist_hop[xnei][hop_nei]) {
                                        for (int index = hop_nei; index <= upper_k; index++) {
                                            if (dist_hop[xnei][index] > d_new) {
                                                dist_hop[xnei][index] = d_new;
                                            }
                                        }
                                        if (Q_VALUE[xnei][hop_nei] == -1) {
                                            L_lock[xnei].lock();
                                            // 查询当前的单侧值 如果更小则加入
                                            auto [best_dis, best_hop] =
                                                    search_sorted_two_hop_label_in_current_with_less_than_k_limit_with_csv(
                                                            (*L)[xnei], v, hop_nei, shard);
                                            L_lock[xnei].unlock();
                                            for (int index = best_hop; index <= hop_nei; ++index) {
                                                if (Q_VALUE[xnei][index] == -1) {
                                                    Q_VALUE[xnei][index] = best_dis;
                                                }
                                            }
                                        }
                                        if (Q_VALUE[xnei][hop_nei] > d_new) {
                                            if (Q_handle.contains({xnei, hop_nei})) {
                                                if (Q_handle[{xnei, hop_nei}].second > d_new) {
                                                    pq.update(Q_handle[{xnei, hop_nei}].first, node);
                                                    Q_handle[{xnei, hop_nei}].second = d_new;
                                                }
                                            } else {
                                                Q_handle[{xnei, hop_nei}] = {pq.push(node), d_new};
                                            }
                                            Q_VALUE[xnei][hop_nei] = d_new;
                                        }
                                    }


                                    if (dist_hop[xnei][hop_nei] < d_new) {
                                        if (hubs[xnei] != -1 && hubs[xnei] != v) {
                                            ppr_lock[xnei].lock();
                                            PPR_TYPE::PPR_insert(PPR, xnei, hubs[xnei], v);
                                            ppr_lock[xnei].unlock();
                                        }
                                        if (hubs[xnei] != -1 && hubs[xnei] != xnei) {
                                            ppr_lock[v].lock();
                                            PPR_TYPE::PPR_insert(PPR, v, hubs[xnei], xnei);
                                            ppr_lock[v].unlock();
                                        }
                                    }
                                }
                            }
                        }
                        for (const auto &[fst, snd]: dist_hop_changes) {
                            dist_hop[fst][snd] = {-1};
                        }
                        Q_VALUE.resize(Q_VALUE.size(),
                                       std::vector<hop_weight_type>(upper_k + 1,
                                                                    std::numeric_limits<hop_weight_type>::max()));
                        mtx_599_1.lock();
                        Qid_599.push(current_tid);
                        mtx_599_1.unlock();
                        return 1;
                    }));
        }
        for (auto &&result: results_dynamic) {
            result.get();
        }
        std::cout << "ruc size is " << size << std::endl;
        std::vector<std::future<int> >().swap(results_dynamic);
    }
}
