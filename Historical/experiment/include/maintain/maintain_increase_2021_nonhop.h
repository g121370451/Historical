#pragma once

#include <map>
#include <utility>
#include "utils/ThreadPool.h"
#include "entity/two_hop_label.h"
#include "utils/ExecutionTimer.h"
#include "entity/graph.h"
#include "entity/two_hop_label.h"
#include "entity/nonhop_global_params.h"

namespace experiment::nonhop::algorithm2021::increase {
    template<typename weight_type, typename hop_weight_type>
    class StrategyA2021NonHopIncrease {
    public:
        void operator()(graph<weight_type> &instance_graph, two_hop_case_info<hop_weight_type> &mm,
                        std::vector<std::pair<int, int>> &v, std::vector<weight_type> &w_old_vec,
                        ThreadPool &pool_dynamic, std::vector<std::future<int>> &results_dynamic, int time);

    private:
        void PI11(graph<weight_type> &instance_graph,
                  std::vector<std::vector<two_hop_label<hop_weight_type>>> *L,
                  std::vector<affected_label> &al1_curr, std::vector<affected_label> *al1_next,
                  std::map<std::pair<int, int>, weight_type> &w_old_map,
                  ThreadPool &pool_dynamic, std::vector<std::future<int>> &results_dynamic, int time);

        void PI12(graph<weight_type> &instance_graph,
                  std::vector<std::vector<two_hop_label<hop_weight_type>>> *L, PPR_TYPE::PPR_type *PPR,
                  std::vector<affected_label> &al1_curr, std::vector<pair_label> *al2_next,
                  ThreadPool &pool_dynamic, std::vector<std::future<int>> &results_dynamic, int time);

        void PI22(graph<weight_type> &instance_graph,
                  std::vector<std::vector<two_hop_label<hop_weight_type>>> *L, PPR_TYPE::PPR_type *PPR,
                  std::vector<pair_label> &al2_curr, std::vector<pair_label> *al2_next,
                  ThreadPool &pool_dynamic, std::vector<std::future<int>> &results_dynamic, int time);
    };

    template<typename weight_type, typename hop_weight_type>
    inline void StrategyA2021NonHopIncrease<weight_type, hop_weight_type>::operator()(
            graph<weight_type> &instance_graph, two_hop_case_info<hop_weight_type> &mm,
            std::vector<std::pair<int, int>> &v, std::vector<weight_type> &w_old_vec,
            ThreadPool &pool_dynamic, std::vector<std::future<int>> &results_dynamic, int time) {
        int current_tid = Qid_595.front();
        auto &counter = experiment::result::global_csv_config.old_counter;
        auto &shard = counter.get_thread_maintain_shard(current_tid);

        std::map<std::pair<int, int>, weight_type> w_old_map;
        size_t batch_size = v.size();
        for (size_t i = 0; i < batch_size; i++) {
            if (v[i].first > v[i].second) {
                std::swap(v[i].first, v[i].second);
            }
            if (w_old_map.count(v[i]) == 0) {
                w_old_map[v[i]] = w_old_vec[i];
            }
        }

        std::vector<affected_label> al1_curr, al1_next;
        std::vector<pair_label> al2_curr, al2_next;

        for (auto &iter: w_old_map) {
            int v1 = iter.first.first;
            int v2 = iter.first.second;
            weight_type w_old = iter.second;

            for (auto it: mm.L[v1]) {
                mtx_595[v2].lock();
                hop_weight_type search_weight = search_sorted_two_hop_label_in_current_with_csv(mm.L[v2], it.vertex,
                                                                                                shard);
                mtx_595[v2].unlock();
                if (it.vertex <= v2 && search_weight >= it.distance + w_old &&
                    search_weight < std::numeric_limits<int>::max()) {
                    al1_curr.push_back(affected_label(v2, it.vertex, it.distance + w_old));
                }
            }
            for (auto it: mm.L[v2]) {
                mtx_595[v1].lock();
                hop_weight_type search_weight = search_sorted_two_hop_label_in_current_with_csv(mm.L[v1], it.vertex,
                                                                                                shard);
                mtx_595[v1].unlock();
                if (it.vertex <= v1 && search_weight >= it.distance + w_old &&
                    search_weight < std::numeric_limits<int>::max()) {
                    al1_curr.push_back(affected_label(v1, it.vertex, it.distance + w_old));
                }
            }
        }

        while (!al1_curr.empty() || !al2_curr.empty()) {

            // if (std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin_time).count() > max_run_time_nanosec) {
            //	throw reach_limit_time_string;
            // }
            PI11(instance_graph, &mm.L, al1_curr, &al1_next, w_old_map, pool_dynamic, results_dynamic,
                 time);
            PI12(instance_graph, &mm.L, &mm.PPR, al1_curr, &al2_next, pool_dynamic, results_dynamic, time);
            PI22(instance_graph, &mm.L, &mm.PPR, al2_curr, &al2_next, pool_dynamic, results_dynamic, time);
            // std::cout << "increase 2021 al1_cuur size is " << al1_curr.size() << " al2_curr size is " << al2_curr.size() << " al2_next size is " << al2_next.size() << std::endl;
            al1_curr = al1_next;
            al2_curr = al2_next;
            std::vector<affected_label>().swap(al1_next);
            std::vector<pair_label>().swap(al2_next);
        }
    };

    template<typename weight_type, typename hop_weight_type>
    inline void
    StrategyA2021NonHopIncrease<weight_type, hop_weight_type>::PI11(graph<weight_type> &instance_graph,
                                                                    std::vector<std::vector<two_hop_label<hop_weight_type>>> *L,
                                                                    std::vector<affected_label> &al1_curr,
                                                                    std::vector<affected_label> *al1_next,
                                                                    std::map<std::pair<int, int>, weight_type> &w_old_map,
                                                                    ThreadPool &pool_dynamic,
                                                                    std::vector<std::future<int>> &results_dynamic,
                                                                    int time) {
        for (auto it: al1_curr) {
            results_dynamic.emplace_back(
                    pool_dynamic.enqueue([it, L, al1_next, &instance_graph, &w_old_map, time] {
                        mtx_595_1.lock();
                        int current_tid = Qid_595.front();
                        Qid_595.pop();
                        mtx_595_1.unlock();
                        auto &counter = experiment::result::global_csv_config.old_counter;
                        auto &shard = counter.get_thread_maintain_shard(current_tid);

                        for (auto nei: instance_graph[it.first]) {
                            mtx_595[nei.first].lock();
                            hop_weight_type search_weight = search_sorted_two_hop_label_in_current_with_csv(
                                    (*L)[nei.first], it.second, shard);
                            mtx_595[nei.first].unlock();
                            weight_type w_old;
                            if (w_old_map.count(std::pair<int, int>(it.first, nei.first)) > 0) {
                                w_old = w_old_map[std::pair<int, int>(it.first, nei.first)];
                            } else if (w_old_map.count(std::pair<int, int>(nei.first, it.first)) > 0) {
                                w_old = w_old_map[std::pair<int, int>(nei.first, it.first)];
                            } else {
                                w_old = nei.second;
                            }
                            if (it.dis + w_old == search_weight) {
                                mtx_595_1.lock();
                                al1_next->push_back(
                                        affected_label(nei.first, it.second, search_weight));
                                mtx_595_1.unlock();
                            }
                        }
                        mtx_595[it.first].lock();
                        insert_sorted_two_hop_label_with_csv(
                                (*L)[it.first], it.second,
                                MAX_VALUE, time,
                                shard); // this does not change the size of L[it->first] here, so does not need to lock here
                        mtx_595[it.first].unlock();
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
    };

    template<typename weight_type, typename hop_weight_type>
    inline void
    StrategyA2021NonHopIncrease<weight_type, hop_weight_type>::PI12(graph<weight_type> &instance_graph,
                                                                    std::vector<std::vector<two_hop_label<hop_weight_type>>> *L,
                                                                    PPR_TYPE::PPR_type *PPR,
                                                                    std::vector<affected_label> &al1_curr,
                                                                    std::vector<pair_label> *al2_next,
                                                                    ThreadPool &pool_dynamic,
                                                                    std::vector<std::future<int>> &results_dynamic,
                                                                    int time) {
        for (auto it: al1_curr) {
            results_dynamic.emplace_back(
                    pool_dynamic.enqueue([it, L, PPR, al2_next, &instance_graph, time] {
                        mtx_595_1.lock();
                        int current_tid = Qid_595.front();
                        Qid_595.pop();
                        mtx_595_1.unlock();
                        auto &counter = experiment::result::global_csv_config.old_counter;
                        auto &shard = counter.get_thread_maintain_shard(current_tid);

                        int v = it.first, u = it.second;
                        mtx_5952[v].lock();
                        std::vector<int> temp = PPR_TYPE::PPR_retrieve(*PPR, v, u);
                        mtx_5952[v].unlock();
                        temp.push_back(u);

                        mtx_595[v].lock();
                        auto Lv = (*L)[v]; // to avoid interlocking
                        mtx_595[v].unlock();

                        for (auto t: temp) {
                            if (v < t) {
                                hop_weight_type d1 = MAX_VALUE;
                                for (auto nei: instance_graph[t]) {
                                    mtx_595[nei.first].lock();
                                    d1 = std::min(d1, search_sorted_two_hop_label_in_current_with_csv(
                                            (*L)[nei.first], v, shard) + nei.second);
                                    mtx_595[nei.first].unlock();
                                }
                                mtx_595[t].lock();
                                auto query_result = graph_weighted_two_hop_extract_distance_and_hub_in_current_with_csv(
                                        (*L)[t], Lv, t, v, shard);
                                mtx_595[t].unlock();
                                if (query_result.first > d1) {
                                    mtx_595[t].lock();
                                    insert_sorted_two_hop_label_with_csv((*L)[t], v, d1, time, shard);
                                    mtx_595[t].unlock();
                                    mtx_595_1.lock();
                                    al2_next->emplace_back(t, v);
                                    mtx_595_1.unlock();
                                    if (status::currentTimeMode == status::MaintainTimeMode::SLOT1) {
                                        shard.diffuse_count_slot1++;
                                    } else {
                                        shard.diffuse_count_slot2++;
                                    }
                                } else {
                                    if (query_result.second != v) {
                                        mtx_5952[t].lock();
                                        PPR_TYPE::PPR_insert_with_csv(PPR, t, query_result.second, v,
                                                                      shard);
                                        mtx_5952[t].unlock();
                                    }
                                    if (query_result.second != t) {
                                        mtx_5952[v].lock();
                                        PPR_TYPE::PPR_insert_with_csv(PPR, v, query_result.second, t,
                                                                      shard);
                                        mtx_5952[v].unlock();
                                    }
                                }
                            }
                            if (t < v) {
                                hop_weight_type d1 = MAX_VALUE;
                                for (auto nei: instance_graph[v]) {
                                    mtx_595[nei.first].lock();
                                    d1 = std::min(d1, search_sorted_two_hop_label_in_current_with_csv(
                                            (*L)[nei.first], t, shard) + nei.second);
                                    mtx_595[nei.first].unlock();
                                }
                                mtx_595[t].lock();
                                auto query_result = graph_weighted_two_hop_extract_distance_and_hub_in_current_with_csv(
                                        (*L)[t], Lv, t, v, shard);
                                mtx_595[t].unlock();
                                if (query_result.first > d1) {
                                    mtx_595[v].lock();
                                    insert_sorted_two_hop_label_with_csv((*L)[v], t, d1, time, shard);
                                    mtx_595[v].unlock();
                                    mtx_595_1.lock();
                                    al2_next->emplace_back(v, t);
                                    mtx_595_1.unlock();
                                    if (status::currentTimeMode == status::MaintainTimeMode::SLOT1) {
                                        shard.diffuse_count_slot1++;
                                    } else {
                                        shard.diffuse_count_slot2++;
                                    }
                                } else {
                                    if (query_result.second != v) {
                                        mtx_5952[t].lock();
                                        PPR_TYPE::PPR_insert_with_csv(PPR, t, query_result.second, v,
                                                                      shard);
                                        mtx_5952[t].unlock();
                                    }
                                    if (query_result.second != t) {
                                        mtx_5952[v].lock();
                                        PPR_TYPE::PPR_insert_with_csv(PPR, v, query_result.second, t,
                                                                      shard);
                                        mtx_5952[v].unlock();
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
    };

    template<typename weight_type, typename hop_weight_type>
    inline void
    StrategyA2021NonHopIncrease<weight_type, hop_weight_type>::PI22(graph<weight_type> &instance_graph,
                                                                    std::vector<std::vector<two_hop_label<hop_weight_type>>> *L,
                                                                    PPR_TYPE::PPR_type *PPR,
                                                                    std::vector<pair_label> &al2_curr,
                                                                    std::vector<pair_label> *al2_next,
                                                                    ThreadPool &pool_dynamic,
                                                                    std::vector<std::future<int>> &results_dynamic,
                                                                    int time) {

        for (auto &it: al2_curr) {
            results_dynamic.emplace_back(
                    pool_dynamic.enqueue([&it, L, PPR, al2_next, &instance_graph, time] {
                        mtx_595_1.lock();
                        int current_tid = Qid_595.front();
                        Qid_595.pop();
                        mtx_595_1.unlock();
                        auto &counter = experiment::result::global_csv_config.old_counter;
                        auto &shard = counter.get_thread_maintain_shard(current_tid);

                        mtx_595[it.second].lock();
                        auto Lxx = (*L)[it.second]; // to avoid interlocking
                        mtx_595[it.second].unlock();

                        for (auto nei: instance_graph[it.first]) {
                            if (nei.first > it.second) {
                                mtx_595[it.first].lock();
                                hop_weight_type search_result =
                                        search_sorted_two_hop_label_in_current_with_csv((*L)[it.first],
                                                                                        it.second,
                                                                                        shard) +
                                        nei.second;
                                mtx_595[it.first].unlock();
                                mtx_595[nei.first].lock();
                                auto query_result = graph_weighted_two_hop_extract_distance_and_hub_in_current_with_csv(
                                        (*L)[nei.first], Lxx, nei.first, it.second, shard);
                                mtx_595[nei.first].unlock();
                                if (query_result.first > search_result) {
                                    mtx_595[nei.first].lock();
                                    insert_sorted_two_hop_label_with_csv((*L)[nei.first], it.second,
                                                                         search_result, time, shard);
                                    mtx_595[nei.first].unlock();
                                    mtx_595_1.lock();
                                    al2_next->push_back(pair_label(nei.first, it.second));
                                    mtx_595_1.unlock();
                                    if (status::currentTimeMode == status::MaintainTimeMode::SLOT1) {
                                        shard.diffuse_count_slot1++;
                                    } else {
                                        shard.diffuse_count_slot2++;
                                    }
                                } else {
                                    if (query_result.second != it.second) {
                                        mtx_5952[nei.first].lock();
                                        PPR_TYPE::PPR_insert_with_csv(PPR, nei.first, query_result.second,
                                                                      it.second, shard);
                                        mtx_5952[nei.first].unlock();
                                    }
                                    if (query_result.second != nei.first) {
                                        mtx_5952[it.second].lock();
                                        PPR_TYPE::PPR_insert_with_csv(PPR, it.second, query_result.second,
                                                                      nei.first, shard);
                                        mtx_5952[it.second].unlock();
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
}


