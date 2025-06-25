#pragma once

#include <map>
#include <utility>
#include "utils/ThreadPool.h"
#include "entity/two_hop_label.h"
#include "utils/ExecutionTimer.h"
#include "entity/graph.h"
#include "entity/two_hop_label.h"
#include "entity/hop_global_params.h"

namespace experiment::hop::algorithm2021::increase {
    template<typename weight_type, typename hop_weight_type>
    class StrategyA2021HopIncrease {
    public:
        void operator()(graph<weight_type> &instance_graph, two_hop_case_info<hop_weight_type> &mm,
                        std::vector<std::pair<int, int> > &v, std::vector<weight_type> &w_old_vec,
                        ThreadPool &pool_dynamic, std::vector<std::future<int> > &results_dynamic, int time);

    private:
        void PI11(graph<int> &instance_graph, std::vector<std::vector<two_hop_label<hop_weight_type> > > *L,
                  std::vector<hop_constrained_affected_label<hop_weight_type> > &al1_curr,
                  std::vector<hop_constrained_affected_label<hop_weight_type> > *al1_next,
                  std::map<std::pair<int, int>, weight_type> &w_old_map,
                  ThreadPool &pool_dynamic, std::vector<std::future<int> > &results_dynamic, int time) const;

        void PI12(graph<int> &instance_graph, std::vector<std::vector<two_hop_label<hop_weight_type> > > *L,
                  PPR_TYPE::PPR_type *PPR,
                  std::vector<hop_constrained_affected_label<hop_weight_type> > &al1_curr,
                  std::vector<hop_constrained_pair_label> *al2_next, ThreadPool &pool_dynamic,
                  std::vector<std::future<int> > &results_dynamic, int upper_k, int time) const;

        void PI22(graph<int> &instance_graph, std::vector<std::vector<two_hop_label<hop_weight_type> > > *L,
                  PPR_TYPE::PPR_type *PPR,
                  std::vector<hop_constrained_pair_label> &al2_curr, std::vector<hop_constrained_pair_label> *al2_next,
                  ThreadPool &pool_dynamic, std::vector<std::future<int> > &results_dynamic, int upper_k,
                  int time) const;
    };

    template<typename weight_type, typename hop_weight_type>
    void StrategyA2021HopIncrease<weight_type, hop_weight_type>::operator()(
            graph<weight_type> &instance_graph, two_hop_case_info<hop_weight_type> &mm,
            std::vector<std::pair<int, int> > &v, std::vector<weight_type> &w_old_vec,
            ThreadPool &pool_dynamic, std::vector<std::future<int> > &results_dynamic, int time) {
        std::map<std::pair<int, int>, weight_type> w_old_map;
        const size_t batch_size = v.size();
        for (size_t i = 0; i < batch_size; i++) {
            if (v[i].first > v[i].second) {
                std::swap(v[i].first, v[i].second);
            }
            if (!w_old_map.contains(v[i])) {
                w_old_map[v[i]] = w_old_vec[i];
            }
        }

        const int current_tid = Qid_599.front();

        auto &counter = result::global_csv_config.old_counter;
        auto &shard = counter.get_thread_maintain_shard(current_tid);

        std::vector<hop_constrained_affected_label<hop_weight_type> > al1_curr, al1_next;
        std::vector<hop_constrained_pair_label> al2_curr, al2_next;
        for (auto &iter: w_old_map) {
            int v1 = iter.first.first;
            int v2 = iter.first.second;
            weight_type w_old = iter.second;
            for (const auto &it: mm.L[v1]) {
                if (it.hub_vertex <= v2 && it.t_e == std::numeric_limits<int>::max()) {
                    hop_weight_type search_weight = search_sorted_two_hop_label_in_current_with_equal_k_limit_with_csv(
                            mm.L[v2], it.hub_vertex,
                            it.hop + 1, shard).first;
                    if (search_weight >= it.distance + w_old &&
                        search_weight < std::numeric_limits<hop_weight_type>::max()) {
                        if (search_weight > it.distance + w_old) {
                            std::cout << "judge the affected label :search_weight is " << search_weight
                                      << " old_length is " << it.distance + w_old << std::endl;
                        }
                        al1_curr.push_back(
                                hop_constrained_affected_label<hop_weight_type>(v2, it.hub_vertex, it.hop + 1,
                                                                                it.distance + w_old));
                    }
                }
            }
            for (auto it: mm.L[v2]) {
                if (it.hub_vertex <= v1 && it.t_e == std::numeric_limits<int>::max()) {
                    hop_weight_type search_weight = search_sorted_two_hop_label_in_current_with_equal_k_limit_with_csv(
                            mm.L[v1], it.hub_vertex,
                            it.hop + 1, shard).first;
                    if (search_weight >= it.distance + w_old &&
                        search_weight < std::numeric_limits<hop_weight_type>::max()) {
                        if (search_weight > it.distance + w_old) {
                            std::cout << "judge the affected label :search_weight is " << search_weight
                                      << " old_length is " << it.distance + w_old << std::endl;
                        }
                        al1_curr.push_back(
                                hop_constrained_affected_label<hop_weight_type>(v1, it.hub_vertex, it.hop + 1,
                                                                                it.distance + w_old));
                    }
                }
            }
        }
        while (al1_curr.size() || !al2_curr.empty()) {
            PI11(instance_graph, &mm.L, al1_curr, &al1_next, w_old_map, pool_dynamic, results_dynamic, time);
            PI12(instance_graph, &mm.L, &mm.PPR, al1_curr, &al2_next, pool_dynamic, results_dynamic, mm.upper_k, time);
            PI22(instance_graph, &mm.L, &mm.PPR, al2_curr, &al2_next, pool_dynamic, results_dynamic, mm.upper_k, time);

            al1_curr = al1_next;
            al2_curr = al2_next;
            std::vector<hop_constrained_affected_label<hop_weight_type> >().swap(al1_next);
            std::vector<hop_constrained_pair_label>().swap(al2_next);
        }
        std::cout << "2021 ppr size " << experiment::PPR_TYPE::getSize(mm.PPR) << std::endl;
    };

    template<typename weight_type, typename hop_weight_type>
    void StrategyA2021HopIncrease<weight_type, hop_weight_type>::PI11(graph<int> &instance_graph,
                                                                      std::vector<std::vector<two_hop_label<
                                                                              hop_weight_type> > > *L,
                                                                      std::vector<hop_constrained_affected_label<
                                                                              hop_weight_type> > &al1_curr,
                                                                      std::vector<hop_constrained_affected_label<
                                                                              hop_weight_type> > *al1_next,
                                                                      std::map<std::pair<int, int>, weight_type> &
                                                                      w_old_map,
                                                                      ThreadPool &pool_dynamic,
                                                                      std::vector<std::future<int> > &results_dynamic,
                                                                      int time) const {
        for (auto it: al1_curr) {
            results_dynamic.emplace_back(pool_dynamic.enqueue([time, it, L, al1_next, &instance_graph, &w_old_map] {
                mtx_599_1.lock();
                int current_tid = Qid_599.front();
                Qid_599.pop();
                mtx_599_1.unlock();
                auto &counter = experiment::result::global_csv_config.old_counter;
                auto &shard = counter.get_thread_maintain_shard(current_tid);

                for (auto nei: instance_graph[it.first]) {
                    L_lock[nei.first].lock();
                    hop_weight_type search_weight = search_sorted_two_hop_label_in_current_with_equal_k_limit_with_csv(
                            (*L)[nei.first],
                            it.second,
                            it.hop + 1, shard).first;
                    L_lock[nei.first].unlock();
                    weight_type w_old = nei.second;
                    if (w_old_map.count(std::pair<int, int>(it.first, nei.first)) > 0) {
                        w_old = w_old_map[std::pair<int, int>(it.first, nei.first)];
                    } else if (w_old_map.count(std::pair<int, int>(nei.first, it.first)) > 0) {
                        w_old = w_old_map[std::pair<int, int>(nei.first, it.first)];
                    } else {
                        w_old = nei.second;
                    }

                    if (it.dis + w_old <= search_weight && search_weight < std::numeric_limits<hop_weight_type>::max()) {
                        mtx_599_1.lock();
                        al1_next->push_back(
                                hop_constrained_affected_label<hop_weight_type>{nei.first, it.second, it.hop + 1,
                                                                                it.dis + w_old});
                        mtx_599_1.unlock();
                    }
                }
                L_lock[it.first].lock();
                insert_sorted_hop_constrained_two_hop_label_with_csv<hop_weight_type>((*L)[it.first], it.second, it.hop,
                                                                     std::numeric_limits<hop_weight_type>::max(),
                                                                     time, shard);
                // this does not change the size of L[it->first] here, so does not need to lock here
                L_lock[it.first].unlock();
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
    void StrategyA2021HopIncrease<weight_type, hop_weight_type>::PI12(graph<int> &instance_graph,
                                                                      std::vector<std::vector<two_hop_label<
                                                                              hop_weight_type> > > *L,
                                                                      PPR_TYPE::PPR_type *PPR,
                                                                      std::vector<hop_constrained_affected_label<
                                                                              hop_weight_type> > &al1_curr,
                                                                      std::vector<hop_constrained_pair_label> *al2_next,
                                                                      ThreadPool &pool_dynamic,
                                                                      std::vector<std::future<int> > &results_dynamic,
                                                                      int upper_k, int time) const {
        for (auto it: al1_curr) {
            results_dynamic.emplace_back(pool_dynamic.enqueue([time, it, L, PPR, al2_next, &instance_graph, upper_k] {
                mtx_599_1.lock();
                int current_tid = Qid_599.front();
                Qid_599.pop();
                mtx_599_1.unlock();

                auto &counter = experiment::result::global_csv_config.old_counter;
                auto &shard = counter.get_thread_maintain_shard(current_tid);

                // v is diffuse origin vertex. u is hub.
                int v = it.first, u = it.second;
                ppr_lock[v].lock();
                std::vector<int> temp = PPR_TYPE::PPR_retrieve(*PPR, v, u);
                ppr_lock[v].unlock();
                temp.push_back(u);

                L_lock[v].lock();
                auto Lv = (*L)[v]; // to avoid interlocking
                L_lock[v].unlock();

                for (auto t: temp) {
                    if (v < t) {
                        hop_weight_type d1 = std::numeric_limits<hop_weight_type>::max();
                        int hop_vn = 0;
                        for (auto nei: instance_graph[t]) {
                            L_lock[nei.first].lock();
                            std::pair<hop_weight_type, int> dis_hop = search_sorted_two_hop_label_in_current_with_csv(
                                    (*L)[nei.first], v, shard);
                            L_lock[nei.first].unlock();
                            if (d1 > dis_hop.first + nei.second) {
                                d1 = dis_hop.first + nei.second;
                                hop_vn = dis_hop.second;
                            }
                        }
                        for (int hop_i = 1; hop_i <= hop_vn + 1; hop_i++) {
                            if (hop_i > upper_k)
                                break;
                            hop_weight_type di = std::numeric_limits<hop_weight_type>::max();
                            for (auto nei: instance_graph[t]) {
                                L_lock[nei.first].lock();
                                di = std::min(
                                        di, search_sorted_two_hop_label_in_current_with_less_than_k_limit_with_csv(
                                                (*L)[nei.first], v,
                                                hop_i - 1, shard).first + nei.second);
                                L_lock[nei.first].unlock();
                            }

                            L_lock[t].lock();
                            auto [query_dis, query_hop, query_hub] =
                                    graph_weighted_two_hop_extract_distance_and_hop_and_hub_in_current_with_csv(
                                            (*L)[t], Lv,
                                            t, v,
                                            hop_i, shard);
                            L_lock[t].unlock();

                            if (query_dis > di) {
                                L_lock[t].lock();
                                insert_sorted_hop_constrained_two_hop_label_with_csv((*L)[t], v, hop_i, di, time,
                                                                                     shard);
                                L_lock[t].unlock();
                                mtx_599_1.lock();
                                al2_next->emplace_back(t, v, hop_i);
                                mtx_599_1.unlock();
                            } else {
                                if (query_hub != -1 && query_hub != v) {
                                    ppr_lock[t].lock();
                                    PPR_TYPE::PPR_insert(PPR, t, query_hub, v);
                                    ppr_lock[t].unlock();
                                }
                                if (query_hub != -1 && query_hub != t) {
                                    ppr_lock[v].lock();
                                    PPR_TYPE::PPR_insert(PPR, v, query_hub, t);
                                    ppr_lock[v].unlock();
                                }
                            }
                        }
                    }
                    if (t < v) {
                        hop_weight_type d1 = std::numeric_limits<hop_weight_type>::max();
                        int hop_vn = 0;
                        for (auto nei: instance_graph[v]) {
                            L_lock[nei.first].lock();
                            std::pair<hop_weight_type, int> dis_hop = search_sorted_two_hop_label_in_current_with_csv(
                                    (*L)[nei.first], t, shard);
                            L_lock[nei.first].unlock();
                            if (d1 > dis_hop.first + nei.second) {
                                d1 = dis_hop.first + nei.second;
                                hop_vn = dis_hop.second;
                            }
                        }

                        for (int hop_i = 1; hop_i <= hop_vn + 1; hop_i++) {
                            if (hop_i > upper_k)
                                break;
                            hop_weight_type di = std::numeric_limits<hop_weight_type>::max();
                            for (auto nei: instance_graph[v]) {
                                L_lock[nei.first].lock();
                                di = std::min(
                                        di, search_sorted_two_hop_label_in_current_with_less_than_k_limit_with_csv(
                                                (*L)[nei.first], t,
                                                hop_i - 1, shard).first + nei.second);
                                L_lock[nei.first].unlock();
                            }
                            L_lock[t].lock();
                            auto [query_dis, query_hop, query_hub] =
                                    graph_weighted_two_hop_extract_distance_and_hop_and_hub_in_current_with_csv(
                                            (*L)[t], Lv,
                                            t, v,
                                            hop_i, shard);
                            L_lock[t].unlock();

                            if (query_dis > di) {
                                L_lock[v].lock();
                                insert_sorted_hop_constrained_two_hop_label_with_csv((*L)[v], t, hop_i, di, time,
                                                                                     shard);
                                L_lock[v].unlock();
                                mtx_599_1.lock();
                                al2_next->emplace_back(v, t, hop_i);
                                mtx_599_1.unlock();
                            } else {
                                if (query_hub != -1 && query_hub != v) {
                                    ppr_lock[t].lock();
                                    PPR_TYPE::PPR_insert(PPR, t, query_hub, v);
                                    ppr_lock[t].unlock();
                                }
                                if (query_hub != -1 && query_hub != t) {
                                    ppr_lock[v].lock();
                                    PPR_TYPE::PPR_insert(PPR, v, query_hub, t);
                                    ppr_lock[v].unlock();
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
    void StrategyA2021HopIncrease<weight_type, hop_weight_type>::PI22(graph<int> &instance_graph,
                                                                      std::vector<std::vector<two_hop_label<
                                                                              hop_weight_type> > > *L,
                                                                      PPR_TYPE::PPR_type *PPR,
                                                                      std::vector<hop_constrained_pair_label> &al2_curr,
                                                                      std::vector<hop_constrained_pair_label> *al2_next,
                                                                      ThreadPool &pool_dynamic,
                                                                      std::vector<std::future<int> > &results_dynamic,
                                                                      int upper_k, int time) const {
        for (auto &it: al2_curr) {
            results_dynamic.emplace_back(pool_dynamic.enqueue([time, &it, L, PPR, al2_next, &instance_graph, upper_k] {
                mtx_599_1.lock();
                int current_tid = Qid_599.front();
                Qid_599.pop();
                mtx_599_1.unlock();

                auto &counter = experiment::result::global_csv_config.old_counter;
                auto &shard = counter.get_thread_maintain_shard(current_tid);
                L_lock[it.second].lock();
                auto Lxx = (*L)[it.second]; // to avoid interlocking
                L_lock[it.second].unlock();

                if (it.hop + 1 > upper_k) {
                    mtx_599_1.lock();
                    Qid_599.push(current_tid);
                    mtx_599_1.unlock();
                    return 0;
                }

                for (auto nei: instance_graph[it.first]) {
                    if (nei.first > it.second) {
                        L_lock[it.first].lock();
                        hop_weight_type search_result =
                                search_sorted_two_hop_label_in_current_with_less_than_k_limit_with_csv((*L)[it.first],
                                                                                                       it.second,
                                                                                                       it.hop,
                                                                                                       shard).first +
                                nei.second;
                        L_lock[it.first].unlock();
                        L_lock[nei.first].lock();
                        auto [query_dis, query_hop, query_hub] =
                                graph_weighted_two_hop_extract_distance_and_hop_and_hub_in_current_with_csv(
                                        (*L)[nei.first], Lxx, nei.first, it.second, it.hop + 1, shard);
                        L_lock[nei.first].unlock();
                        if (query_dis > search_result) {
                            L_lock[nei.first].lock();
                            insert_sorted_hop_constrained_two_hop_label_with_csv((*L)[nei.first], it.second, it.hop + 1,
                                                                                 search_result, time, shard);
                            L_lock[nei.first].unlock();
                            mtx_599_1.lock();
                            al2_next->emplace_back(nei.first, it.second, it.hop + 1);
                            mtx_599_1.unlock();
                        } else {
                            if (query_hub != -1 && query_hub != it.second) {
                                ppr_lock[nei.first].lock();
                                PPR_TYPE::PPR_insert(PPR, nei.first, query_hub, it.second);
                                ppr_lock[nei.first].unlock();
                            }
                            if (query_hub != -1 && query_hub != nei.first) {
                                ppr_lock[it.second].lock();
                                PPR_TYPE::PPR_insert(PPR, it.second, query_hub, nei.first);
                                ppr_lock[it.second].unlock();
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
}
