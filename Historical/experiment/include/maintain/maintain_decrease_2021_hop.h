#pragma once

#include <map>
#include "entity/graph.h"
#include "entity/two_hop_label.h"
#include "utils/ThreadPool.h"
#include "entity/hop_global_params.h"

namespace experiment::hop::algorithm2021::decrease {
    template<typename weight_type, typename hop_weight_type>
    class StrategyA2021HopDecrease {
    public:
        void operator()(graph<weight_type> &instance_graph, two_hop_case_info<hop_weight_type> &mm,
                        std::vector<std::pair<int, int> > v, std::vector<weight_type> w_new,
                        ThreadPool &pool_dynamic, std::vector<std::future<int> > &results_dynamic, int time) const;

    private:
        void ProDecreasep_batch(graph<weight_type> &instance_graph,
                                std::vector<std::vector<two_hop_label<hop_weight_type> > > *L,
                                PPR_TYPE::PPR_type *PPR,
                                std::vector<hop_constrained_affected_label<hop_weight_type> > &CL_curr,
                                std::vector<hop_constrained_affected_label<hop_weight_type> > *CL_next,
                                ThreadPool &pool_dynamic, std::vector<std::future<int> > &results_dynamic, int upper_k,
                                int t) const;
    };

    template<typename weight_type, typename hop_weight_type>
    void StrategyA2021HopDecrease<weight_type, hop_weight_type>::ProDecreasep_batch(graph<weight_type> &instance_graph,
                                                                                    std::vector<std::vector<two_hop_label<hop_weight_type> > > *L,
                                                                                    PPR_TYPE::PPR_type *PPR,
                                                                                    std::vector<hop_constrained_affected_label<hop_weight_type> > &CL_curr,
                                                                                    std::vector<hop_constrained_affected_label<hop_weight_type> > *CL_next,
                                                                                    ThreadPool &pool_dynamic,
                                                                                    std::vector<std::future<int> > &results_dynamic,
                                                                                    int upper_k, int time) const {
        for (const auto &it: CL_curr) {
            results_dynamic.emplace_back(pool_dynamic.enqueue([time, it, L, PPR, CL_next, &instance_graph, upper_k] {
                mtx_599_1.lock();
                int current_tid = Qid_599.front();
                Qid_599.pop();
                mtx_599_1.unlock();

                auto &counter = experiment::result::global_csv_config.old_counter;
                auto &shard = counter.get_thread_maintain_shard(current_tid);

                const int v = it.first;
                int u = it.second;

                L_lock[u].lock();
                auto Lu = (*L)[u]; // to avoid interlocking
                L_lock[u].unlock();

                if (it.hop + 1 > upper_k) {
                    mtx_599_1.lock();
                    Qid_599.push(current_tid);
                    mtx_599_1.unlock();
                    return 0;
                }
                for (auto nei: instance_graph[v]) {
                    int vnei = nei.first;
                    int hop_u = it.hop;
                    hop_weight_type dnew = it.dis + nei.second;
                    if(dnew <0){
                        std::cout <<"overflow happen in maintain decrease 2021 with hop" << " it_dis is " << it.dis <<" nei.second is " << nei.second <<std::endl;
                        exit(0);
                    }
                    if (u < vnei) {
                        L_lock[vnei].lock();
                        auto [query_dis, query_hop, query_hub] = graph_weighted_two_hop_extract_distance_and_hop_and_hub_in_current_with_csv(
                                (*L)[vnei], Lu,
                                vnei, u,
                                hop_u + 1, shard); // query_result is {distance, common hub}
                        L_lock[vnei].unlock();
                        if (query_dis > dnew) {
                            L_lock[vnei].lock();
                            insert_sorted_hop_constrained_two_hop_label_with_csv((*L)[vnei], u, hop_u + 1, dnew, time,
                                                                                 shard);
                            L_lock[vnei].unlock();
                            mtx_599_1.lock();
                            CL_next->push_back(
                                    hop_constrained_affected_label<hop_weight_type>{vnei, u, hop_u + 1, dnew});
                            mtx_599_1.unlock();
                        } else {
                            L_lock[vnei].lock();
                            auto [label_dis, label_hop] = search_sorted_two_hop_label_in_current_with_less_than_k_limit_with_csv(
                                    (*L)[vnei], u,
                                    hop_u + 1, shard);
                            L_lock[vnei].unlock();
                            if (label_dis < std::numeric_limits<hop_weight_type>::max() && label_dis > dnew) {
                                L_lock[vnei].lock();
                                insert_sorted_hop_constrained_two_hop_label_with_csv((*L)[vnei], u, hop_u + 1, dnew, time,
                                                                            shard);
                                L_lock[vnei].unlock();
                                mtx_599_1.lock();
                                CL_next->push_back(
                                        hop_constrained_affected_label<hop_weight_type>{vnei, u, hop_u + 1, dnew});
                                mtx_599_1.unlock();
                            }
                            if (query_hub != u) {
                                ppr_lock[vnei].lock();
                                PPR_TYPE::PPR_insert(PPR, vnei, query_hub, u);
                                ppr_lock[vnei].unlock();
                            }
                            if (query_hub != vnei) {
                                ppr_lock[u].lock();
                                PPR_TYPE::PPR_insert(PPR, u, query_hub, vnei);
                                ppr_lock[u].unlock();
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
    void StrategyA2021HopDecrease<weight_type, hop_weight_type>::operator()(graph<weight_type> &instance_graph,
                                                                            two_hop_case_info<hop_weight_type> &mm,
                                                                            std::vector<std::pair<int, int> > v,
                                                                            std::vector<weight_type> w_new,
                                                                            ThreadPool &pool_dynamic,
                                                                            std::vector<std::future<int> > &
                                                                            results_dynamic,
                                                                            int time) const {
        int current_tid = Qid_599.front();
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

        std::vector<hop_constrained_affected_label<hop_weight_type> > CL_curr, CL_next;

        auto &L = mm.L;
        /*
        the following part does not suit parallel computation:
        the reason is that L is changed below, and as a result, in each following loop, L[v2] or L[v1] is locked at each step,
        which means that following loops cannot be actually parallized
        */
        for (auto &w_new_item: w_new_map) {
            int v1 = w_new_item.first.first, v2 = w_new_item.first.second;
            weight_type _w_new = w_new_item.second;
            for (int sl = 0; sl < 2; sl++) {
                if (sl == 1) {
                    std::swap(v1, v2);
                }
                for (const auto &it: L[v1]) {
                    int _v = it.hub_vertex;
                    int hop_v = it.hop;
                    hop_weight_type dis = it.distance + _w_new;
                    if (it.distance == std::numeric_limits<weight_type>::max()) {
                        continue;
                    }
                    if (dis < 0) {
                        mm.print_L_vk(v1);
                        std::cout << it << std::endl;
                        std::cout << "overflow happen in 2021 maintain decrease step1 it dis is " << it.distance <<" _w_new is "<< _w_new << std::endl;
                    }
                    if (_v <= v2 && hop_v + 1 < mm.upper_k && it.t_e == std::numeric_limits<int>::max()) {
                        auto [query_dis, query_hop, query_hub] =
                                graph_weighted_two_hop_extract_distance_and_hop_and_hub_in_current_with_csv(
                                        L[_v], L[v2], _v,
                                        v2, hop_v + 1,
                                        shard);
                        //                        auto query_result = hop_constrained_extract_distance_and_hub(L, _v, v2, hop_v + 1); // query_result is {distance, common hub}

                        if (query_dis > dis) {
                            insert_sorted_hop_constrained_two_hop_label_with_csv(L[v2], _v, hop_v + 1, dis, time,
                                                                                 shard);
                            CL_curr.push_back(hop_constrained_affected_label<hop_weight_type>(v2, _v, hop_v + 1, dis));
                        } else {
                            auto [label_dis, label_hop] =
                                    search_sorted_two_hop_label_in_current_with_less_than_k_limit_with_csv(L[v2], _v,
                                                                                                           hop_v + 1,
                                                                                                           shard);
                            if (label_dis < std::numeric_limits<hop_weight_type>::max() && label_dis > dis) {
                                insert_sorted_hop_constrained_two_hop_label_with_csv((L)[v2], _v, hop_v + 1, dis, time,
                                                                                     shard);
                                CL_curr.push_back(
                                        hop_constrained_affected_label<hop_weight_type>(v2, _v, hop_v + 1, dis));
                            }
                            if (query_hub != _v) {
                                PPR_TYPE::PPR_insert(&(mm.PPR), v2, query_hub, _v);
                            }
                            if (query_hub != v2) {
                                PPR_TYPE::PPR_insert(&(mm.PPR), _v, query_hub, v2);
                            }
                        }
                    }
                }
            }
        }
        while (!CL_curr.empty()) {
            ProDecreasep_batch(instance_graph, &mm.L, &mm.PPR, CL_curr, &CL_next, pool_dynamic, results_dynamic,
                               mm.upper_k, time);
            CL_curr = CL_next;
            std::vector<hop_constrained_affected_label<hop_weight_type> >().swap(CL_next);
        }
    }
}
