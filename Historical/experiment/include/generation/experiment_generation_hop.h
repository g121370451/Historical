#pragma once
#include <boost/heap/fibonacci_heap.hpp>
#include <algorithm>
#include "utils/ExecutionTimer.h"
#include "entity/graph.h"
#include "entity/two_hop_label.h"
#include "utils/ThreadPool.h"
#include "entity/hop_global_params.h"
namespace experiment
{
    namespace hop
    {
        inline ExecutionTimer generator_timer;
        template <typename weight_type, typename hop_weight_type>
        class GeneratorHop
        {
        public:
            void operator()(graph<weight_type> &input_graph, two_hop_case_info<hop_weight_type> &case_info);

        private:
            int global_upper_k = 0;
            std::vector<std::vector<two_hop_label<hop_weight_type>>> L_temp_599;
            std::vector<std::vector<two_hop_label<hop_weight_type>>> Lv_final_599;
            void HSDL_thread_function(int v_k, graph<weight_type> &graph);

            std::vector<std::vector<two_hop_label<hop_weight_type>>> hop_constrained_sortL(int num_of_threads);

            /*canonical_repair*/
            void hop_constrained_clean_L(two_hop_case_info<hop_weight_type> &case_info, int thread_num);

            void hop_constrained_clear_global_values();
        };
        template <typename weight_type, typename hop_weight_type>
        void GeneratorHop<weight_type, hop_weight_type>::operator()(graph<weight_type> &graph, two_hop_case_info<hop_weight_type> &case_info)
        {
            //----------------------------------- step 1: initialization -----------------------------------
            generator_timer.startSubtask("step 1: initialization");
            int N = graph.size();
            /* store the L Label and PPR*/
            this->L_temp_599.resize(N);
            this->Lv_final_599.resize(N);
            PPR_599.resize(N);

            int num_of_threads = case_info.thread_num;
            ThreadPool pool(num_of_threads);
            std::vector<std::future<int>> results;
            generator_timer.endSubtask();
            //----------------------------------------------- step 2: generate labels ---------------------------------------------------------------
            /**
             * Use Temp_L_vk_599 to mark the positional relationship between the iterated nodes
             * and the nodes with already generated labels. This is done to reduce the process of
             * traversing the L labels of the iterated nodes. Additionally, register multithreaded
             * tasks and retrieve the results
             */
            generator_timer.startSubtask("step 2: generate labels");
            global_upper_k = case_info.upper_k;
            Temp_L_vk_599<hop_weight_type>.resize(num_of_threads);
            dist_hop_599<hop_weight_type>.resize(num_of_threads);
            Q_handle_priorities_599<hop_weight_type>.resize(num_of_threads);
            for (int i = 0; i < num_of_threads; i++)
            {
                Temp_L_vk_599<hop_weight_type>[i].resize(N);
                dist_hop_599<hop_weight_type>[i].resize(N, {std::numeric_limits<hop_weight_type>::max(), 0});
                Q_handle_priorities_599<hop_weight_type>[i].resize(N);
                for (int j = 0; j < N; j++) {
                    hop_constrained_node_handle<hop_weight_type> handle_x;
                    Q_handle_priorities_599<hop_weight_type>[i][j].resize(global_upper_k + 1, {handle_x, std::numeric_limits<hop_weight_type>::max()});
                }
                Qid_599.push(i);
            }

            int last_check_vID = N - 1;

            for (int v_k = 0; v_k <= last_check_vID; v_k++)
            {
                results.emplace_back(
                    pool.enqueue([this, v_k, &graph]
                                 {
							this->HSDL_thread_function(v_k,graph);
							return 1; }));
            }

            for (auto &&result : results)
                result.get();
            generator_timer.endSubtask();
            //----------------------------------------------- step 3: sortL ---------------------------------------------------------------
            generator_timer.startSubtask("step 3: sortL");
            case_info.L = L_temp_599;
            case_info.L = hop_constrained_sortL(num_of_threads);
            case_info.PPR = PPR_599;
            generator_timer.endSubtask();
            //----------------------------------------------- step 4: canonical_repair---------------------------------------------------------------
            generator_timer.startSubtask("step 4: canonical_repair");
//            hop_constrained_clean_L(case_info, num_of_threads);
            generator_timer.endSubtask();
            //---------------------------------------------------------------------------------------------------------------------------------------
            hop_constrained_clear_global_values();
            std::cout << "end" << std::endl;
        }
        template <typename weight_type, typename hop_weight_type>
        inline void GeneratorHop<weight_type, hop_weight_type>::hop_constrained_clean_L(two_hop_case_info<hop_weight_type> &case_info, int thread_num)
        {

            std::vector<std::vector<experiment::hop::two_hop_label<hop_weight_type>>> &L = case_info.L;
            const int N = L.size();

            ThreadPool pool(thread_num);
            std::vector<std::future<int>> results;
            for (int v = 0; v < N; v++)
            {
                results.emplace_back(
                    pool.enqueue([this, v, &L] { // pass const type value j to thread; [] can be empty
                        mtx_599[max_N_ID_for_mtx_599 - 1].lock();
                        std::cout << "clean pll v is " << v << std::endl;
                        int used_id = Qid_599.front();
                        Qid_599.pop();
                        mtx_599[max_N_ID_for_mtx_599 - 1].unlock();

                        std::vector<two_hop_label<hop_weight_type>> &Lv_final = this->Lv_final_599[v];

                        /**
                         * get the L result of the current vertex
                         */

                        std::vector<two_hop_label<hop_weight_type>> &Lv = L[v];

                        /**
                         * the temp_L in this thread
                         */
                        auto &T = Temp_L_vk_599<hop_weight_type>[used_id];

                        /**
                         * Traverse the L-list of the current vertex
                         */
                        for (const auto &Lvi : Lv)
                        {
                            int u = Lvi.hub_vertex;
                            int u_hop = Lvi.hop;

                            /**
                             * Traverse the L on the opposite vertex of the current label.
                             */
                            const auto &Lu = L[u];

                            /**
                             * traverse downward from the perfectly correct first vertex
                             */
                            int min_dis = std::numeric_limits<hop_weight_type>::max();
                            for (auto &label1 : Lu)
                            {
                                for (auto &label2 : T[label1.hub_vertex])
                                {
                                    if (label1.hop + label2.second <= u_hop)
                                    {
                                        hop_weight_type query_dis = label1.distance + label2.first;
                                        if (query_dis < 0) {
                                            std::cout << "overflow happen in clean with hop" << std::endl;
                                        }
                                        if (query_dis < min_dis)
                                        {
                                            min_dis = query_dis;
                                        }
                                    }
                                }
                            }

                            if (min_dis > Lvi.distance)
                            {
                                Lv_final.push_back(two_hop_label(Lvi));
                                T[u].push_back({Lvi.distance, Lvi.hop});
                            }
                        }

                        for (const auto &label : Lv_final)
                        {
                            std::vector<std::pair<hop_weight_type, int>>().swap(T[label.hub_vertex]);
                        }

                        mtx_599[max_N_ID_for_mtx_599 - 1].lock();
                        Qid_599.push(used_id);
                        mtx_599[max_N_ID_for_mtx_599 - 1].unlock();

                        return 1; // return to results; the return type must be the same with results
                    }));
            }

            for (auto &&result : results)
                result.get(); // all threads finish here
            std::cout << "start move" << std::endl;
            case_info.L = std::move(Lv_final_599);
            std::cout << "end move" << std::endl;
            results.clear();
        };
        template <typename weight_type, typename hop_weight_type>
        void GeneratorHop<weight_type, hop_weight_type>::hop_constrained_clear_global_values()
        {
            std::vector<std::vector<two_hop_label<hop_weight_type>>>().swap(L_temp_599);
            std::vector<std::vector<two_hop_label<hop_weight_type>>>().swap(Lv_final_599);
            std::vector<std::vector<std::vector<std::pair<hop_weight_type, int>>>>().swap(Temp_L_vk_599<hop_weight_type>);
            std::vector<std::vector<std::pair<hop_weight_type, int>>>().swap(dist_hop_599<hop_weight_type>);
            std::vector<std::vector<std::vector<std::pair<hop_constrained_node_handle<hop_weight_type>, hop_weight_type>>>>().swap(Q_handle_priorities_599<hop_weight_type>);
            std::queue<int>().swap(Qid_599);
            PPR_TYPE::PPR_type().swap(PPR_599);
        }
        template <typename weight_type, typename hop_weight_type>
        void GeneratorHop<weight_type, hop_weight_type>::HSDL_thread_function(int v_k, graph<weight_type> &graph)
        {
            /* get unique thread id */
            /* critical section obtain array index  */
            mtx_599[max_N_ID_for_mtx_599 - 1].lock();
            std::cout << "calculate pll vk is " << v_k << std::endl;
            int used_id = Qid_599.front();
            Qid_599.pop();
            mtx_599[max_N_ID_for_mtx_599 - 1].unlock();

            /* store the Temp L and dist_hop in current thread*/
            std::vector<int> Temp_L_vk_changes, dist_hop_changes;
            /* Temp_L_vk stores the dest_vertex_id and distance and hop */
            auto &Temp_L_vk = Temp_L_vk_599<hop_weight_type>[used_id];
            auto &dist_hop = dist_hop_599<hop_weight_type>[used_id]; // record the minimum distance (and the corresponding hop) of a searched vertex in Q
            std::vector<std::pair<int, int>> Q_handle_priorities_changes;
            /* get the label list in current thread*/
            auto &Q_handle_priorities = Q_handle_priorities_599<hop_weight_type>[used_id];

            /* a class contains information about destination vertex, hop count, and cost in the priority queue */
            boost::heap::fibonacci_heap<two_hop_label<hop_weight_type>> Q;

            /* generate a label for the starting vertex itself */
            two_hop_label<hop_weight_type> node;
            node.hub_vertex = v_k;
            node.hop = 0;
            node.distance = 0;
            Q_handle_priorities[v_k][0] = {Q.push(node), node.distance};
            Q_handle_priorities_changes.emplace_back(v_k, 0);
            // size_t size = 0;
            /* Temp_L_vk_599 stores the label (dist and hop) of vertex v_k */
            mtx_599[v_k].lock();
            /* root is vk-> vk->obj info -> vector<obj> -> index-> vertexId obj-><distance,hop> */
            for (auto &xx : L_temp_599[v_k])
            {
                int L_vk_vertex = xx.hub_vertex;
                Temp_L_vk[L_vk_vertex].push_back({xx.distance, xx.hop});
                Temp_L_vk_changes.push_back(L_vk_vertex);
            }
            mtx_599[v_k].unlock();
            /*  dist_hop_599 stores the shortest distance from vk to any other vertices with its hop_cst,
                note that the hop_cst is determined by the shortest distance */
            dist_hop[v_k] = {0, 0};
            dist_hop_changes.push_back(v_k);
            while (!Q.empty())
            {
                // size = std::max(size, Q.size());
                /* poll the vertex from heap.In other words, poll the vertex with the minimal cost */
                node = Q.top();
                Q.pop();

                /* current node, u, which is the node generating the labels*/
                int u = node.hub_vertex;

                if (v_k > u)
                {
                    continue;
                }

                int u_hop = node.hop;
                hop_weight_type P_u = node.distance;
                int common_hub_for_query_v_k_u = -1;
                hop_weight_type query_v_k_u = std::numeric_limits<hop_weight_type>::max();
                mtx_599[u].lock();
                for (auto &xx : L_temp_599[u])
                {
                    int common_v = xx.hub_vertex;
                    for (auto &yy : Temp_L_vk[common_v])
                    {
                        if (xx.hop + yy.second <= u_hop)
                        {
                            hop_weight_type dis = xx.distance + yy.first;
                            if (dis < 0) {
                                std::cout <<"overflow happen in generation with hop" << std::endl;
                            }
                            if (query_v_k_u > dis)
                            {
                                query_v_k_u = dis;
                                common_hub_for_query_v_k_u = xx.hub_vertex;
                            }
                        }
                    }
                }
                mtx_599[u].unlock();
                if (P_u < query_v_k_u)
                {
                    node.hub_vertex = v_k;
                    node.hop = u_hop;
                    node.distance = P_u;
                    mtx_599[u].lock();
                    L_temp_599[u].push_back(node);
                    mtx_599[u].unlock();

                    if (u_hop + 1 > global_upper_k)
                    {
                        continue;
                    }

                    /* update adj */
                    /* Traverse neighboring nodes */
                    for (auto &xx : graph[u])
                    {
                        /* adh_v is the neighborhood and the ec is the distance from u to ajd_v*/
                        int adj_v = xx.first, ec = xx.second;

                        /* update node info */
                        node.hub_vertex = adj_v;
                        node.distance = P_u + ec;
                        node.hop = u_hop + 1;

                        auto &yy = Q_handle_priorities[adj_v][node.hop];

                        if (yy.second <= node.distance)
                        { // adj_v has been reached with a smaller distance and the same hop
                            continue;
                        }
                        /* the vertex has not been visited*/
                        if (dist_hop[adj_v].first == std::numeric_limits<hop_weight_type>::max())
                        { // adj_v has not been reached
                            yy = {Q.push(node), node.distance};
                            Q_handle_priorities_changes.emplace_back(adj_v, node.hop);
                            dist_hop[adj_v].first = node.distance;
                            dist_hop[adj_v].second = node.hop;
                            dist_hop_changes.push_back(adj_v);
                        }
                        else
                        {
                            if (node.distance < dist_hop[adj_v].first)
                            { // adj_v has been reached with a less distance
                                if (yy.second != std::numeric_limits<hop_weight_type>::max())
                                {
                                    Q.update(yy.first, node);
                                    yy.second = node.distance;
                                }
                                else
                                {
                                    yy = {Q.push(node), node.distance};
                                    Q_handle_priorities_changes.push_back({adj_v, node.hop});
                                }
                                dist_hop[adj_v].first = node.distance;
                                dist_hop[adj_v].second = node.hop;
                            }
                            else if (node.hop < dist_hop[adj_v].second)
                            { // adj_v has been reached with a less hop
                                if (yy.second != std::numeric_limits<hop_weight_type>::max())
                                {
                                    Q.update(yy.first, node);
                                    yy.second = node.distance;
                                }
                                else
                                {
                                    yy = {Q.push(node), node.distance};
                                    Q_handle_priorities_changes.push_back({adj_v, node.hop});
                                }
                            }
                        }
                    }
                }
                else
                {
                    /* add v_k into PPR(u,common_hub_for_query_v_k_u), and add u into PPR(v_k,common_hub_for_query_v_k_u)*/
                    if (common_hub_for_query_v_k_u != v_k)
                    {
                        ppr_599[u].lock();
                        PPR_TYPE::PPR_insert(&PPR_599, u, common_hub_for_query_v_k_u, v_k);
                        ppr_599[u].unlock();
                    }
                    if (common_hub_for_query_v_k_u != u)
                    {
                        ppr_599[v_k].lock();
                        PPR_TYPE::PPR_insert(&PPR_599, v_k, common_hub_for_query_v_k_u, u);
                        ppr_599[v_k].unlock();
                    }
                }
            }
            for (const auto &xx : Temp_L_vk_changes)
            {
                std::vector<std::pair<hop_weight_type, int>>().swap(Temp_L_vk[xx]);
            }
            for (const auto &xx : dist_hop_changes)
            {
                dist_hop[xx] = {std::numeric_limits<hop_weight_type>::max(), 0};
            }
            for (auto &[fst, snd] : Q_handle_priorities_changes) {
                hop_constrained_node_handle<hop_weight_type> handle_x;
                Q_handle_priorities[fst][snd] = {handle_x, std::numeric_limits<hop_weight_type>::max()};
            }

            // mtx_599[v_k].lock();
            // std::vector<two_hop_label>(L_temp_599[v_k]).swap(L_temp_599[v_k]);
            // mtx_599[v_k].unlock();

            mtx_599[max_N_ID_for_mtx_599 - 1].lock();
            Qid_599.push(used_id);
            mtx_599[max_N_ID_for_mtx_599 - 1].unlock();
        }
        template <typename weight_type, typename hop_weight_type>
        std::vector<std::vector<two_hop_label<hop_weight_type>>> GeneratorHop<weight_type, hop_weight_type>::hop_constrained_sortL(int num_of_threads)
        {

            /*time complexity: O(V*L*logL), where L is average number of labels per vertex*/

            int N = L_temp_599.size();
            std::vector<std::vector<two_hop_label<hop_weight_type>>> output_L(N);

            /*time complexity: O(V*L*logL), where L is average number of labels per vertex*/
            ThreadPool pool(num_of_threads);
            std::vector<std::future<int>> results; // return typename: xxx
            for (int v_k = 0; v_k < N; v_k++)
            {
                results.emplace_back(
                    pool.enqueue([this, &output_L, v_k] { // pass const type value j to thread; [] can be empty
                        std::sort(this->L_temp_599[v_k].begin(), this->L_temp_599[v_k].end(), experiment::hop::compare_hop_constrained_two_hop_label<hop_weight_type>);
                        std::vector<two_hop_label<hop_weight_type>>(this->L_temp_599[v_k]).swap(this->L_temp_599[v_k]); // ʹ��vector��swap�Ż��ڴ�ռ�ã��ͷŶ���Ŀռ�
                        output_L[v_k] = this->L_temp_599[v_k];
                        std::vector<two_hop_label<hop_weight_type>>().swap(this->L_temp_599[v_k]); // clear new labels for RAM efficiency

                        return 1; // return to results; the return type must be the same with results
                    }));
            }
            for (auto &&result : results)
                result.get(); // all threads finish here

            return output_L;
        }
    };
}