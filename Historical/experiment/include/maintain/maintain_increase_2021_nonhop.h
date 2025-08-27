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

#ifdef _DEBUG
        std::vector<record_in_increase<hop_weight_type> > list_infinite;
        std::vector<pair_label> global_al2;
#endif
        std::vector<record_in_increase<hop_weight_type> > list;
    private:
        void PI11(graph<weight_type> &instance_graph,
                  std::vector<std::vector<two_hop_label<hop_weight_type>>> *L,
                  std::vector<affected_label<hop_weight_type>> &al1_curr,
                  std::vector<affected_label<hop_weight_type>> *al1_next,
                  std::map<std::pair<int, int>, weight_type> &w_old_map,
                  ThreadPool &pool_dynamic, std::vector<std::future<int>> &results_dynamic, int time);

        void PI12(graph<weight_type> &instance_graph,
                  std::vector<std::vector<two_hop_label<hop_weight_type>>> *L, PPR_TYPE::PPR_type *PPR,
                  std::vector<affected_label<hop_weight_type>> &al1_curr, std::vector<pair_label> *al2_next,
                  ThreadPool &pool_dynamic, std::vector<std::future<int>> &results_dynamic, int time);

        void PI22(graph<weight_type> &instance_graph,
                  std::vector<std::vector<two_hop_label<hop_weight_type>>> *L, PPR_TYPE::PPR_type *PPR,
                  std::vector<pair_label> &al2_curr, std::vector<pair_label> *al2_next,
                  ThreadPool &pool_dynamic, std::vector<std::future<int>> &results_dynamic, int time);
    };

    template<typename weight_type, typename hop_weight_type>
    inline void StrategyA2021NonHopIncrease<weight_type, hop_weight_type>::operator()(
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
        int current_tid = Qid_595.front();

        auto &counter = experiment::result::global_csv_config.old_counter;
        auto &shard = counter.get_thread_maintain_shard(current_tid);

        std::vector<affected_label<hop_weight_type>> al1_curr, al1_next;
        std::vector<pair_label> al2_curr, al2_next;

        for (const std::pair<const std::pair<int, int>, weight_type> &iter: w_old_map) {
            int v1 = iter.first.first;
            int v2 = iter.first.second;
            weight_type w_old = iter.second;

            for (const two_hop_label<hop_weight_type> &it: mm.L[v1]) {
                if (it.distance == std::numeric_limits<hop_weight_type>::max()) {
                    continue;
                }
                hop_weight_type d_new = it.distance + w_old;
#ifdef _DEBUG
                if (d_new < 0) {
                    std::cout << "overflow happen in increase 2021 with hop" << std::endl;
                }
#endif
                if (it.vertex <= v2 && it.t_e == std::numeric_limits<int>::max()) {
                    hop_weight_type search_weight = search_sorted_two_hop_label_in_current_with_csv(mm.L[v2], it.vertex,
                        shard);
                    if (search_weight == d_new &&
                        search_weight < std::numeric_limits<hop_weight_type>::max()) {
                        al1_curr.push_back(affected_label(v2, it.vertex, d_new));
                    }
                }
            }
            for (auto it: mm.L[v2]) {
                if (it.distance == std::numeric_limits<hop_weight_type>::max()) {
                    continue;
                }
                hop_weight_type d_new = it.distance + w_old;
#ifdef _DEBUG
                if (d_new < 0) {
                    std::cout << "overflow happen in increase 2021 with hop" << std::endl;
                }
#endif
                if (it.vertex <= v1 && it.t_e == std::numeric_limits<int>::max()) {
                    hop_weight_type search_weight = search_sorted_two_hop_label_in_current_with_csv(mm.L[v1], it.vertex,
                        shard);
                    if (search_weight == d_new
                        && search_weight < std::numeric_limits<hop_weight_type>::max()) {
                        al1_curr.push_back(affected_label(v1, it.vertex, d_new));
                    }
                }
            }
        }
        unsigned long long int al1Size = 0;
        unsigned long long int al2Size = 0;
#ifdef _DEBUG
        double cost1 = 0;
        double cost2 = 0;
        double cost3 = 0;
#endif
        while (!al1_curr.empty() || !al2_curr.empty()) {
            std::cout << " al1 size is " << al1_curr.size() << " al2_size is " << al2_curr.size() << std::endl;
            al1Size += al1_curr.size();
            al2Size += al2_curr.size();
#ifdef _DEBUG
            auto time1 = std::chrono::steady_clock::now();
#endif
            PI11(instance_graph, &mm.L, al1_curr, &al1_next, w_old_map, pool_dynamic, results_dynamic,
                 time);
#ifdef _DEBUG
            auto time2 = std::chrono::steady_clock::now();
#endif
            PI12(instance_graph, &mm.L, &mm.PPR, al1_curr, &al2_next, pool_dynamic, results_dynamic, time);
#ifdef _DEBUG
            auto time3 = std::chrono::steady_clock::now();
#endif
            PI22(instance_graph, &mm.L, &mm.PPR, al2_curr, &al2_next, pool_dynamic, results_dynamic, time);
#ifdef _DEBUG
            auto time4 = std::chrono::steady_clock::now();
            cost1 += std::chrono::duration_cast<std::chrono::duration<double> >(time2 - time1).count();
            cost2 += std::chrono::duration_cast<std::chrono::duration<double> >(time3 - time2).count();
            cost3 += std::chrono::duration_cast<std::chrono::duration<double> >(time4 - time3).count();
#endif
            al1_curr = al1_next;
            al2_curr = al2_next;
#ifdef _DEBUG
            for (auto& item : al2_next) {
                this->global_al2.push_back(std::move(item));
            }
#endif
            std::vector<affected_label<hop_weight_type>>().swap(al1_next);
            std::vector<pair_label>().swap(al2_next);
        }
#ifdef _DEBUG
        std::cout << "2021 " << cost1 << " " << cost2 << " " << cost3 << std::endl;
#endif
        std::cout << "2021 algorithm al1Size is " << al1Size << " al2Size is " << al2Size << std::endl;
    }

    template<typename weight_type, typename hop_weight_type>
    inline void
    StrategyA2021NonHopIncrease<weight_type, hop_weight_type>::PI11(graph<weight_type> &instance_graph,
                                                                    std::vector<std::vector<two_hop_label<
                                                                        hop_weight_type>>> *L,
                                                                    std::vector<affected_label<hop_weight_type> > &
                                                                    al1_curr,
                                                                    std::vector<affected_label<hop_weight_type> > *
                                                                    al1_next,
                                                                    std::map<std::pair<int, int>, weight_type> &
                                                                    w_old_map,
                                                                    ThreadPool &pool_dynamic,
                                                                    std::vector<std::future<int> > &results_dynamic,
                                                                    int time) {
#ifdef _DEBUG
        for (auto &it: al1_curr) {
            this->list_infinite.emplace_back(it.first, it.second,
                                             std::numeric_limits<hop_weight_type>::max(), it.dis, time);
        }
#endif
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
                        if (it.second < nei.first) {
                            L_lock[nei.first].lock();
                            hop_weight_type search_weight = search_sorted_two_hop_label_in_current_with_csv(
                                (*L)[nei.first], it.second, shard);
                            L_lock[nei.first].unlock();
                            weight_type w_old = nei.second;
                            if (w_old_map.count(std::pair<int, int>(it.first, nei.first)) > 0) {
                                w_old = w_old_map[std::pair<int, int>(it.first, nei.first)];
                            } else if (w_old_map.count(std::pair<int, int>(nei.first, it.first)) > 0) {
                                w_old = w_old_map[std::pair<int, int>(nei.first, it.first)];
                            } else {
                                w_old = nei.second;
                            }
                            hop_weight_type d_old = it.dis + w_old;
#ifdef _DEBUG
                            if (d_old < 0) {
                                std::cout << "overflow happen in maintain increase 2021 with hop pi11"
                                        << " it dis is "
                                        << it.dis << " w_old is " << w_old << std::endl;
                            }
#endif
                            if (d_old == search_weight &&
                                search_weight < std::numeric_limits<hop_weight_type>::max()) {
                                mtx_595_1.lock();
                                al1_next->push_back(
                                    affected_label(nei.first, it.second, search_weight));
                                mtx_595_1.unlock();
                            }
                        }
                    }
                    L_lock[it.first].lock();
                    insert_sorted_two_hop_label_with_csv(
                        (*L)[it.first], it.second,
                        std::numeric_limits<hop_weight_type>::max(), time,
                        shard); // this does not change the size of L[it->first] here, so does not need to lock here
                    L_lock[it.first].unlock();
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
    };

    template<typename weight_type, typename hop_weight_type>
    inline void
    StrategyA2021NonHopIncrease<weight_type, hop_weight_type>::PI12(graph<weight_type> &instance_graph,
                                                                    std::vector<std::vector<two_hop_label<
                                                                        hop_weight_type> > > *L,
                                                                    PPR_TYPE::PPR_type *PPR,
                                                                    std::vector<affected_label<hop_weight_type> > &
                                                                    al1_curr,
                                                                    std::vector<pair_label> *al2_next,
                                                                    ThreadPool &pool_dynamic,
                                                                    std::vector<std::future<int> > &results_dynamic,
                                                                    int time) {
        for (auto it: al1_curr) {
            results_dynamic.emplace_back(
                pool_dynamic.enqueue([it, L, PPR, al2_next, &instance_graph, time, this] {
                    mtx_595_1.lock();
                    int current_tid = Qid_595.front();
                    Qid_595.pop();
                    mtx_595_1.unlock();
                    auto &counter = experiment::result::global_csv_config.old_counter;
                    auto &shard = counter.get_thread_maintain_shard(current_tid);

                    int v = it.first, u = it.second;
                    ppr_lock[v].lock();
                    std::vector<int> temp = PPR_TYPE::PPR_retrieve(*PPR, v, u);
                    ppr_lock[v].unlock();
                    temp.push_back(u);

                    for (auto t: temp) {
                        int diffuseVertex = std::max(v, t);
                        int targetVertex = std::min(v, t);
                        hop_weight_type d1 = std::numeric_limits<hop_weight_type>::max();
                        for (auto nei: instance_graph[diffuseVertex]) {
                            L_lock[nei.first].lock();
                            hop_weight_type dis = search_sorted_two_hop_label_in_current_with_csv(
                                (*L)[nei.first], targetVertex, shard);
                            L_lock[nei.first].unlock();
                            if (dis == std::numeric_limits<hop_weight_type>::max()) {
                                continue;
                            }
                            hop_weight_type d_new = dis + nei.second;
#ifdef _DEBUG
                            if (d_new < 0) {
                                std::cout << "overflow happen in pi12 maintain increase 2021 with nonhop"
                                        << std::endl;
                            }
#endif
                            d1 = std::min(d1, d_new);
                        }
                        L_lock[targetVertex].lock();
                        L_lock[diffuseVertex].lock();
                        auto [query_dis, query_hub] =
                                graph_weighted_two_hop_extract_distance_and_hub_in_current_with_csv(
                                    (*L)[diffuseVertex], (*L)[targetVertex],
                                    diffuseVertex, targetVertex, shard);
                        L_lock[diffuseVertex].unlock();
                        L_lock[targetVertex].unlock();
                        if (query_dis > d1) {
                            L_lock[diffuseVertex].lock();
                            insert_sorted_two_hop_label_with_csv((*L)[diffuseVertex], targetVertex, d1, time, shard);
                            L_lock[diffuseVertex].unlock();
                            mtx_list_check.lock();
                            this->list.emplace_back(diffuseVertex, targetVertex, d1, query_dis, time);
                            mtx_list_check.unlock();
                            mtx_595_1.lock();
                            al2_next->emplace_back(diffuseVertex, targetVertex);
                            mtx_595_1.unlock();
                            if (status::currentTimeMode == status::MaintainTimeMode::SLOT1) {
                                shard.diffuse_count_slot1++;
                            } else {
                                shard.diffuse_count_slot2++;
                            }
                        } else if (query_dis != std::numeric_limits<hop_weight_type>::max() && query_dis <= d1) {
                            if (query_hub != -1 && query_hub != v) {
                                ppr_lock[t].lock();
                                PPR_TYPE::PPR_insert_with_csv(PPR, t, query_hub, v,
                                                              shard);
                                ppr_lock[t].unlock();
                            }
                            if (query_hub != -1 && query_hub != t) {
                                ppr_lock[v].lock();
                                PPR_TYPE::PPR_insert_with_csv(PPR, v, query_hub, t,
                                                              shard);
                                ppr_lock[v].unlock();
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
    };

    template<typename weight_type, typename hop_weight_type>
    inline void
    StrategyA2021NonHopIncrease<weight_type, hop_weight_type>::PI22(graph<weight_type> &instance_graph,
                                                                    std::vector<std::vector<two_hop_label<
                                                                        hop_weight_type> > > *L,
                                                                    PPR_TYPE::PPR_type *PPR,
                                                                    std::vector<pair_label> &al2_curr,
                                                                    std::vector<pair_label> *al2_next,
                                                                    ThreadPool &pool_dynamic,
                                                                    std::vector<std::future<int> > &results_dynamic,
                                                                    int time) {
        for (auto &it: al2_curr) {
            results_dynamic.emplace_back(
                pool_dynamic.enqueue([&it, L, PPR, al2_next, &instance_graph, time, this] {
                    mtx_595_1.lock();
                    int current_tid = Qid_595.front();
                    Qid_595.pop();
                    mtx_595_1.unlock();

                    auto &counter = experiment::result::global_csv_config.old_counter;
                    auto &shard = counter.get_thread_maintain_shard(current_tid);
                    L_lock[it.second].lock();
                    auto Lxx = (*L)[it.second]; // to avoid interlocking
                    L_lock[it.second].unlock();

                    L_lock[it.first].lock();
                    hop_weight_type search_result =
                            search_sorted_two_hop_label_in_current_with_csv((*L)[it.first],
                                                                            it.second,
                                                                            shard);
                    L_lock[it.first].unlock();
                    if (search_result != std::numeric_limits<hop_weight_type>::max()) {
                        for (auto nei: instance_graph[it.first]) {
                            if (nei.first > it.second) {
                                hop_weight_type d_new = search_result + nei.second;
#ifdef _DEBUG
                                if (d_new < 0) {
                                    std::cout << "overflow happen in maintain increase 2021 with hop pi22"
                                            << std::endl;
                                }
#endif
                                L_lock[nei.first].lock();
                                auto [query_dis, query_hub] =
                                        graph_weighted_two_hop_extract_distance_and_hub_in_current_with_csv(
                                            (*L)[nei.first], Lxx, nei.first, it.second, shard);
                                L_lock[nei.first].unlock();
                                if (query_dis > d_new) {
                                    L_lock[nei.first].lock();
                                    insert_sorted_two_hop_label_with_csv((*L)[nei.first], it.second,
                                                                         d_new, time, shard);
                                    L_lock[nei.first].unlock();
                                    mtx_list_check.lock();
                                    this->list.emplace_back(nei.first, it.second, d_new, query_dis,
                                                            time);
                                    mtx_list_check.unlock();
                                    mtx_595_1.lock();
                                    al2_next->push_back(pair_label(nei.first, it.second));
                                    mtx_595_1.unlock();
                                    if (status::currentTimeMode == status::MaintainTimeMode::SLOT1) {
                                        shard.diffuse_count_slot1++;
                                    } else {
                                        shard.diffuse_count_slot2++;
                                    }
                                } else {
                                    if (query_hub != -1 && query_hub != it.second) {
                                        ppr_lock[nei.first].lock();
                                        PPR_TYPE::PPR_insert_with_csv(PPR, nei.first, query_hub,
                                                                      it.second, shard);
                                        ppr_lock[nei.first].unlock();
                                    }
                                    if (query_hub != -1 && query_hub != nei.first) {
                                        ppr_lock[it.second].lock();
                                        PPR_TYPE::PPR_insert_with_csv(PPR, it.second, query_hub,
                                                                      nei.first, shard);
                                        ppr_lock[it.second].unlock();
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
        std::vector<std::future<int> >().swap(results_dynamic);
    }
}
