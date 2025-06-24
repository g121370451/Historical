#pragma once
#include <boost/heap/fibonacci_heap.hpp>
#include <algorithm>
#include "utils/ExecutionTimer.h"
#include "entity/graph.h"
#include "entity/two_hop_label.h"
#include "utils/ThreadPool.h"
#include "entity/nonhop_global_params.h"

namespace experiment::nonhop {
    inline ExecutionTimer generator_timer;

    template<typename weight_type, typename hop_weight_type>
    class GeneratorNonHop {
    public:
        void operator()(graph<weight_type> &input_graph, two_hop_case_info<hop_weight_type> &case_info) const;

    private:
        void PLL_dij_function(int v_k, graph<weight_type> &input_graph) const;

        std::vector<std::vector<two_hop_label<hop_weight_type> > > sortL(int num_of_threads) const;

        void clean_L(two_hop_case_info<hop_weight_type> &case_info, int thread_num) const;

        void PLL_clear_global_values() const;
    };

    template<typename weight_type, typename hop_weight_type>
    void GeneratorNonHop<weight_type, hop_weight_type>::operator()(
        graph<weight_type> &input_graph, two_hop_case_info<hop_weight_type> &case_info) const {
        //----------------------------------- step 1: initialization ------------------------------------------------------------------
        generator_timer.startSubtask("step 1: initialization");
        int num_of_threads = case_info.thread_num;
        int N = input_graph.ADJs.size();

        experiment::nonhop::L_temp_595<hop_weight_type>.resize(N);
        PPR_595.resize(N);
        generator_timer.endSubtask();
        //---------------------------------------------------------------------------------------------------------------------------------------

        //----------------------------------------------- step 2: generate labels ---------------------------------------------------------------
        generator_timer.startSubtask("step 2: generate labels"); {
            // to save RAM of ThreadPool
            /*seaching shortest paths*/
            ThreadPool pool(num_of_threads);
            std::vector<std::future<int> > results; // return typename: xxx
            P_dij_595<hop_weight_type>.resize(num_of_threads);
            T_dij_595<hop_weight_type>.resize(num_of_threads);
            Q_handles_595<hop_weight_type>.resize(num_of_threads);
            std::queue<int>().swap(Qid_595);
            for (int i = 0; i < num_of_threads; i++) {
                P_dij_595<hop_weight_type>[i].resize(N, std::numeric_limits<hop_weight_type>::max());
                T_dij_595<hop_weight_type>[i].resize(N, std::numeric_limits<hop_weight_type>::max());
                Q_handles_595<hop_weight_type>[i].resize(N);
                Qid_595.push(i);
            }
            int last_check_vID = N - 1;
            for (int v_k = 0; v_k <= last_check_vID; v_k++) {
                results.emplace_back(
                    pool.enqueue([this, v_k, &input_graph, last_check_vID] {
                        // pass const type value j to thread; [] can be empty
                        this->PLL_dij_function(v_k, input_graph);
                        return 1; // return to results; the return type must be the same with results
                    }));
            }
            for (auto &&result: results)
                result.get(); // all threads finish here
            results.clear();
        }
        generator_timer.endSubtask();
        //---------------------------------------------------------------------------------------------------------------------------------------

        //----------------------------------------------- step 3: sortL ---------------------------------------------------------------
        generator_timer.startSubtask("step 3: sortL");
        case_info.L = this->sortL(num_of_threads);
        case_info.PPR = PPR_595;
        generator_timer.endSubtask();
        //---------------------------------------------------------------------------------------------------------------------------------------

        //----------------------------------------------- step 3: canonical_repair ---------------------------------------------------------------
        generator_timer.startSubtask("step 4: canonical_repair");
        this->clean_L(case_info, num_of_threads);
        generator_timer.endSubtask();
        //---------------------------------------------------------------------------------------------------------------------------------------
        this->PLL_clear_global_values();
    };

    template<typename weight_type, typename hop_weight_type>
    void GeneratorNonHop<weight_type, hop_weight_type>::PLL_dij_function(
        int v_k, graph<weight_type> &input_graph) const {
        mtx_595[max_N_ID_for_mtx_595 - 1].lock();
        int used_id = Qid_595.front();
        Qid_595.pop();
        mtx_595[max_N_ID_for_mtx_595 - 1].unlock();

        std::vector<int> P_changed_vertices, T_changed_vertices;
        std::vector<hop_weight_type> &T_dij = T_dij_595<hop_weight_type>[used_id], P_dij = P_dij_595<hop_weight_type>[
            used_id];

        std::vector<PLL_handle_t_for_sp<hop_weight_type> > &Q_handles = Q_handles_595<hop_weight_type>[used_id];

        boost::heap::fibonacci_heap<two_hop_label<hop_weight_type> > Q;
        two_hop_label<hop_weight_type> node(0);
        node.vertex = v_k;
        node.distance = 0;
        Q_handles[v_k] = Q.push(node);
        P_dij[v_k] = 0;
        P_changed_vertices.push_back(v_k);
        mtx_595[v_k].lock();
        for (const auto &xx: L_temp_595<hop_weight_type>[v_k]) {
            int L_v_k_i_vertex = xx.vertex;
            T_dij[L_v_k_i_vertex] = xx.distance;
            T_changed_vertices.push_back(L_v_k_i_vertex);
        }
        mtx_595[v_k].unlock();
        int new_label_num = 0;
        int count = 0;
        while (!Q.empty()) {
            count++;
            node = Q.top();
            Q.pop();
            int u = node.vertex;

            int P_u = node.distance;

            if (v_k > u) {
                continue;
            }

            int query_v_k_u = std::numeric_limits<hop_weight_type>::max();
            int common_hub_for_query_v_k_u = 0;
            mtx_595[u].lock(); // put lock in for loop is very slow
            for (const auto &xx: L_temp_595<hop_weight_type>[u]) {
                if (T_dij[xx.vertex] == std::numeric_limits<hop_weight_type>::max()) {
                    continue;
                }
                hop_weight_type dis = xx.distance + T_dij[xx.vertex];
                if (dis < 0) {
                    std::cout << "overflow happen in generation with nonhop" << std::endl;
                }
                if (query_v_k_u > dis) {
                    query_v_k_u = dis;
                    common_hub_for_query_v_k_u = xx.vertex;
                }
            }
            mtx_595[u].unlock();

            if (P_u < query_v_k_u) {
                node.vertex = v_k;
                node.distance = P_u;

                mtx_595[u].lock();
                L_temp_595<hop_weight_type>[u].push_back(node);
                // std::cout << "u: " << u << " vk: " << v_k << " dis: " << node.distance << std::endl;
                mtx_595[u].unlock();
                new_label_num++;

                for (const auto &xx: input_graph.ADJs[u]) {
                    int adj_v = xx.first, ec = xx.second;
                    if (P_dij[adj_v] == std::numeric_limits<hop_weight_type>::max()) {
                        node.vertex = adj_v;
                        node.distance = P_u + ec;

                        Q_handles[adj_v] = Q.push(node);
                        P_dij[adj_v] = node.distance;
                        P_changed_vertices.push_back(adj_v);
                    } else if (P_dij[adj_v] > P_u + ec) {
                        node.vertex = adj_v;
                        node.distance = P_u + ec;
                        Q.update(Q_handles[adj_v], node);
                        P_dij[adj_v] = node.distance;
                    }
                }
            } else {
                if (common_hub_for_query_v_k_u != v_k) {
                    ppr_595[u].lock();
                    PPR_TYPE::PPR_insert(&PPR_595, u, common_hub_for_query_v_k_u, v_k);
                    ppr_595[u].unlock();
                }
                if (common_hub_for_query_v_k_u != u) {
                    ppr_595[v_k].lock();
                    PPR_TYPE::PPR_insert(&PPR_595, v_k, common_hub_for_query_v_k_u, u);
                    ppr_595[v_k].unlock();
                }
            }
        }

        for (const auto i: P_changed_vertices) {
            P_dij[i] = std::numeric_limits<int>::max(); // reverse-allocate P values
        }
        for (const auto i: T_changed_vertices) {
            T_dij[i] = std::numeric_limits<int>::max(); // reverse-allocate T values
        }

        mtx_595[max_N_ID_for_mtx_595 - 1].lock();
        std::cout << "calculate pll v: " << v_k << std::endl;
        Qid_595.push(used_id);
        mtx_595[max_N_ID_for_mtx_595 - 1].unlock();
    };

    template<typename weight_type, typename hop_weight_type>
    std::vector<std::vector<two_hop_label<hop_weight_type> > > GeneratorNonHop<weight_type,
        hop_weight_type>::sortL(int num_of_threads) const {
        /*time complexity: O(V*L*logL), where L is average number of labels per vertex*/

        int N = L_temp_595<hop_weight_type>.size();
        std::vector<std::vector<two_hop_label<hop_weight_type> > > output_L(N);

        /*time complexity: O(V*L*logL), where L is average number of labels per vertex*/
        ThreadPool pool(num_of_threads);
        std::vector<std::future<int> > results; // return typename: xxx
        for (int v_k = 0; v_k < N; v_k++) {
            results.emplace_back(
                pool.enqueue([this, &output_L, v_k] {
                    // pass const type value j to thread; [] can be empty
                    std::sort(L_temp_595<hop_weight_type>[v_k].begin(), L_temp_595<hop_weight_type>[v_k].end(),
                              experiment::nonhop::compare_two_hop_label_small_to_large<hop_weight_type>);
                    std::vector<two_hop_label<hop_weight_type> >(L_temp_595<hop_weight_type>[v_k]).swap(
                        L_temp_595<hop_weight_type>[v_k]);
                    output_L[v_k] = L_temp_595<hop_weight_type>[v_k];
                    std::vector<two_hop_label<hop_weight_type> >().swap(L_temp_595<hop_weight_type>[v_k]);
                    // clear new labels for RAM efficiency

                    return 1; // return to results; the return type must be the same with results
                }));
        }
        for (auto &&result: results)
            result.get(); // all threads finish here
        results.clear();

        return output_L;
    };

    template<typename weight_type, typename hop_weight_type>
    void GeneratorNonHop<weight_type, hop_weight_type>::clean_L(
        two_hop_case_info<hop_weight_type> &case_info, int thread_num) const {
        auto &L = case_info.L;
        const int N = L.size();

        ThreadPool pool(thread_num);
        std::vector<std::future<int> > results;

        for (int v = 0; v < N; v++) {
            results.emplace_back(
                pool.enqueue([v, &L] {
                    // pass const type value j to thread; [] can be empty
                    mtx_595[max_N_ID_for_mtx_595 - 1].lock();
                    int used_id = Qid_595.front();
                    Qid_595.pop();
                    mtx_595[max_N_ID_for_mtx_595 - 1].unlock();

                    std::vector<two_hop_label<hop_weight_type> > Lv_final_inner;

                    std::vector<two_hop_label<hop_weight_type> > &Lv = L[v];

                    auto &T = T_dij_595<hop_weight_type>[used_id];

                    for (const auto &Lvi: Lv) {
                        int u = Lvi.vertex;
                        if (v == u) {
                            Lv_final_inner.push_back(two_hop_label(Lvi));
                            T[v] = Lvi.distance;
                            continue;
                        }
                        mtx_595[u].lock();
                        const auto &Lu = L[u];

                        int min_dis = std::numeric_limits<hop_weight_type>::max();
                        for (const auto &label: Lu) {
                            if (T[label.vertex] == std::numeric_limits<hop_weight_type>::max()) {
                                continue;
                            }
                            hop_weight_type query_dis = label.distance + T[label.vertex];
                            if (query_dis < 0) {
                                std::cout << "overflow happen in clean with nonhop" << std::endl;
                            }
                            if (query_dis < min_dis) {
                                min_dis = query_dis;
                            }
                        }
                        mtx_595[u].unlock();
                        if (min_dis > Lvi.distance) {
                            Lv_final_inner.push_back(two_hop_label(Lvi));
                            T[u] = Lvi.distance;
                        }
                    }

                    for (const auto &label: Lv_final_inner) {
                        T[label.vertex] = std::numeric_limits<hop_weight_type>::max();
                    }

                    mtx_595[v].lock();
                    Lv = std::move(Lv_final_inner);
                    mtx_595[v].unlock();

                    mtx_595[max_N_ID_for_mtx_595 - 1].lock();
                    Qid_595.push(used_id);
                    std::cout << "print pll v: " << v << std::endl;
                    mtx_595[max_N_ID_for_mtx_595 - 1].unlock();

                    return 1; // return to results; the return type must be the same with results
                }));
        }

        for (auto &&result: results)
            result.get(); // all threads finish here
        results.clear();
    }

    template<typename weight_type, typename hop_weight_type>
    void GeneratorNonHop<weight_type, hop_weight_type>::PLL_clear_global_values() const {
        std::vector<std::vector<two_hop_label<hop_weight_type> > >().swap(L_temp_595<hop_weight_type>);
        PPR_TYPE::PPR_type().swap(PPR_595);
        std::queue<int>().swap(Qid_595);
        std::vector<std::vector<hop_weight_type> >().swap(P_dij_595<hop_weight_type>);
        std::vector<std::vector<hop_weight_type> >().swap(T_dij_595<hop_weight_type>);
        std::vector<std::vector<PLL_handle_t_for_sp<hop_weight_type> > >().swap(Q_handles_595<hop_weight_type>);
    }
}
