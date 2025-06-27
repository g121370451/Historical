#pragma once

#include <limits>
#include <vector>
#include <iostream>
#include "utils/BinaryPersistence.h"
#include "utils/global.h"
#include "utils/vector_operations.h"


namespace experiment {
    namespace PPR_TYPE {
        using PPR_type = std::vector<std::vector<std::pair<int, std::vector<int>>>>;

        long long int getSize(PPR_type &PPR);

        int PPR_binary_operations_insert(std::vector<int> &input_vector, int key);

        void PPR_insert(PPR_type *PPR, int v1, int v2, int v3);

        void PPR_insert_with_csv(PPR_type *PPR, int v1, int v2, int v3, result::MaintainShard &shard);

        std::vector<int> PPR_retrieve(PPR_type &PPR, int v1, int v2);

        void PPR_replace(PPR_type &PPR, int v1, int v2, std::vector<int> &loads);

        void PPR_erase(PPR_type &PPR, int v1, int v2, int v3);
    }

    namespace nonhop {
        template<typename weight_type>
        class two_hop_label {
        public:
            int vertex;
            weight_type distance;
            int t_s, t_e;

            two_hop_label();

            two_hop_label(const two_hop_label<weight_type> &other);

            explicit two_hop_label(int start_time);

            two_hop_label &operator=(const two_hop_label &) = default;

            bool operator==(const two_hop_label<weight_type> &other) const {
                return vertex == other.vertex && distance == other.distance && t_s == other.t_s && t_e == other.t_e;
            };

            bool operator<(const two_hop_label<weight_type> &other) const {
                return distance > other.distance; // < is the max-heap; > is the min heap
            }

            void serialize(std::ofstream &out) const {
                experiment::saveBinary(out, vertex);
                experiment::saveBinary(out, distance);
                experiment::saveBinary(out, t_s);
                experiment::saveBinary(out, t_e);
            };

            void deserialize(std::ifstream &in) {
                experiment::loadBinary(in, vertex);
                experiment::loadBinary(in, distance);
                experiment::loadBinary(in, t_s);
                experiment::loadBinary(in, t_e);
            };
        };

        template<typename weight_type>
        bool compare_two_hop_label_small_to_large(two_hop_label<weight_type> &i, two_hop_label<weight_type> &j) {
            if (i.t_e != j.t_e)
                return i.t_e > j.t_e;    // t_e降序
            return i.vertex < j.vertex; // < is from small to big; > is from big to small
        };

        // method of 2hop label
        /*result as std::pair<int,int> means <distance,hub>*/
        template<typename hop_weight_type>
        std::pair<hop_weight_type, int> graph_weighted_two_hop_extract_distance_and_hub_in_current_with_csv(
                std::vector<two_hop_label<hop_weight_type>> &L_s,
                std::vector<two_hop_label<hop_weight_type>> &L_t, int source, int terminal,
                experiment::result::MaintainShard &maintain_shard) {
            /*return std::numeric_limits<double>::max() is not connected*/
            if (experiment::status::currentTimeMode == experiment::status::MaintainTimeMode::SLOT1) {
                ++maintain_shard.cover_count_slot1;
            } else {
                ++maintain_shard.cover_count_slot2;
            }
            if (source == terminal) {
                return {0, source};
            }

            hop_weight_type distance = std::numeric_limits<hop_weight_type>::max(); // if disconnected, return this large value
            int common_hub = -1;

            auto vector1_check_pointer = L_s.begin();
            auto vector2_check_pointer = L_t.begin();
            auto pointer_L_s_end = L_s.end(), pointer_L_t_end = L_t.end();
            while (vector1_check_pointer != pointer_L_s_end && vector2_check_pointer != pointer_L_t_end &&
                   vector1_check_pointer->t_e == std::numeric_limits<int>::max() &&
                   vector2_check_pointer->t_e == std::numeric_limits<int>::max()) {
                if (vector1_check_pointer->vertex == vector2_check_pointer->vertex) {
                    hop_weight_type dis = vector1_check_pointer->distance + vector2_check_pointer->distance;
                    if (distance > dis) {
                        distance = dis;
                        common_hub = vector1_check_pointer->vertex;
                    }
                    ++vector1_check_pointer;
                } else if (vector1_check_pointer->vertex > vector2_check_pointer->vertex) {
                    ++vector2_check_pointer;
                } else {
                    ++vector1_check_pointer;
                }
            }

            return {distance, common_hub};
        };

        /*result as std::pair<int,int> means <distance,vertex>*/
        template<typename hop_weight_type>
        hop_weight_type
        search_sorted_two_hop_label_in_current_with_csv(std::vector<two_hop_label<hop_weight_type>> &input_vector,
                                                        int key, experiment::result::MaintainShard &maintain_shard) {
            if (experiment::status::currentTimeMode == experiment::status::MaintainTimeMode::SLOT1) {
                ++maintain_shard.label_count_slot1;
            } else {
                ++maintain_shard.label_count_slot2;
            }
            int left = 0, right = input_vector.size() - 1;

            while (left <= right) {
                int mid = left + ((right - left) / 2);

                if (input_vector[mid].t_e == std::numeric_limits<int>::max()) {
                    if (input_vector[mid].vertex == key) {
                        return input_vector[mid].distance;
                    } else if (input_vector[mid].vertex < key) {
                        left = mid + 1;
                    } else {
                        right = mid - 1;
                    }
                } else {
                    right = mid - 1;
                }
            }

            return std::numeric_limits<hop_weight_type>::max();
        }

        template<typename hop_weight_type>
        void
        insert_sorted_two_hop_label_with_csv(std::vector<two_hop_label<hop_weight_type>> &input_vector, int key,
                                             int value,
                                             int time, experiment::result::MaintainShard &maintain_shard) {
            if (experiment::status::currentTimeMode == experiment::status::MaintainTimeMode::SLOT1) {
                if (experiment::status::currentMaintainMode == experiment::status::MaintainMode::DECREASE) {
                    ++maintain_shard.label_decrease_insert_slot1;
                } else if (experiment::status::currentMaintainMode == experiment::status::MaintainMode::INCREASE) {
                    ++maintain_shard.label_increase_insert_slot1;
                }
            } else {
                if (experiment::status::currentMaintainMode == experiment::status::MaintainMode::DECREASE) {
                    ++maintain_shard.label_decrease_insert_slot2;
                } else if (experiment::status::currentMaintainMode == experiment::status::MaintainMode::INCREASE) {
                    ++maintain_shard.label_increase_insert_slot2;
                }
            }
            int left = 0, right = input_vector.size() - 1;

            while (left <= right) {
                int mid = left + ((right - left) / 2);

                if (input_vector[mid].t_e == std::numeric_limits<int>::max()) {
                    if (input_vector[mid].vertex == key) {
                        two_hop_label<hop_weight_type> old_label = input_vector[mid];
                        old_label.t_e = time - 1;
                        if (input_vector[mid].distance == value) {
                            std::cout << "nihop has same label,so dont insert new label to it" << std::endl;
                            return;
                        }
                        input_vector[mid].distance = value;
                        input_vector[mid].t_s = time;
                        if (old_label.distance != std::numeric_limits<hop_weight_type>::max()) {
                            int insert_left = mid + 1, insert_right = input_vector.size() - 1;

                            while (insert_left <= insert_right) {
                                int insert_mid = insert_left + ((insert_right - insert_left) / 2);

                                if (input_vector[insert_mid].t_e < time) {
                                    insert_right = insert_mid - 1;
                                } else if (input_vector[insert_mid].t_e == time) {
                                    if (input_vector[insert_mid].vertex > key) {
                                        insert_right = insert_mid - 1;
                                    } else if (input_vector[insert_mid].vertex < key) {
                                        insert_left = insert_mid + 1;
                                    } else {
                                        insert_left = insert_mid;
                                        break;
                                    }
                                } else {
                                    insert_left = insert_mid + 1;
                                }
                            }
                            input_vector.insert(input_vector.begin() + insert_left, old_label);
                        }
                        return;
                    } else if (input_vector[mid].vertex < key) {
                        left = mid + 1;
                    } else {
                        right = mid - 1;
                    }
                } else {
                    right = mid - 1;
                }
            }
            if (value != std::numeric_limits<hop_weight_type>::max()) {
                two_hop_label<hop_weight_type> new_label(time);
                new_label.vertex = key;
                new_label.distance = value;

                input_vector.insert(input_vector.begin() + left, new_label);
            }
        }

        template<typename hop_weight_type>
        class two_hop_case_info {
        public:
            int thread_num = 1;

            /*labels*/
            std::vector<std::vector<two_hop_label<hop_weight_type>>> L;
            PPR_TYPE::PPR_type PPR;

            bool operator==(const two_hop_case_info<hop_weight_type> &other) const {
                if (thread_num != other.thread_num || L != other.L || PPR != other.PPR) {
                    return false;
                }
                return true;
            }

            void serialize(std::ofstream &out) const {
                experiment::saveBinary(out, thread_num);
                size_t size = L.size();
                BinarySerializer<size_t>::saveBinary(out, size);
                for (auto &item: L) {
                    BinarySerializer<std::vector<experiment::nonhop::two_hop_label<hop_weight_type>>>::saveBinary(out,
                                                                                                                  item);
                    out.flush();
                }
                experiment::saveBinary(out, PPR);
            }

            void deserialize(std::ifstream &in) {
                experiment::loadBinary(in, thread_num);
                size_t size;
                BinarySerializer<size_t>::loadBinary(in, size);
                L.resize(size);
                for (auto &item: L) {
                    BinarySerializer<std::vector<experiment::nonhop::two_hop_label<hop_weight_type>>>::loadBinary(in,
                                                                                                                  item);
                }
                experiment::loadBinary(in, PPR);
            }

            long long int compute_L_size() {
                long long int res = 0;
                for (const std::vector<experiment::nonhop::two_hop_label<hop_weight_type>> &L_info: L) {
                    for (const experiment::nonhop::two_hop_label<hop_weight_type> &inner: L_info) {
                        ++res;
                    }
                }
                return res;
            }

            /*clear labels*/
            void clear_labels() {
                std::vector<std::vector<two_hop_label<hop_weight_type>>>().swap(L);
                PPR_TYPE::PPR_type().swap(PPR);
            }

            /*compute label size; this should equal label_size_after_canonical_repair when use_canonical_repair==true*/
            long long int compute_L_byte_size() {
                long long int size = 0;
                for (auto it = L.begin(); it != L.end(); it++) {
                    size = size + (*it).size() * sizeof(two_hop_label<hop_weight_type>); // 12 byte per two_hop_label
                }
                return size;
            }

            long long int compute_PPR_byte_size() {
                long long int size = 0;
                for (const auto &i: PPR) {
                    for (int j = 0; j < i.size(); j++) {
                        // TODO-GPY 如果正确要修改这里
                        // size = size + (PPR[i][j].second.size() + 1) * sizeof(int);
                    }
                }
                return size;
            }

            /*printing*/
            void print_L() {
                std::cout << "print_L:" << std::endl;
                for (int i = 0; i < L.size(); i++) {
                    std::cout << "L[" << i << "]=";
                    for (int j = 0; j < L[i].size(); j++) {
                        std::cout << "{" << L[i][j].vertex << "," << L[i][j].distance << "," << L[i][j].t_s << ","
                                  << L[i][j].t_e << "}";
                    }
                    std::cout << std::endl;
                }
            }

            /*record_all_details*/
            void record_all_details(const std::string &save_name) {
                std::ofstream outputFile;
                outputFile.precision(6);
                outputFile.setf(std::ios::fixed);
                outputFile.setf(std::ios::showpoint);
                outputFile.open(save_name + ".txt");

                outputFile << "PLL info:" << std::endl;
                outputFile << "thread_num=" << thread_num << std::endl;
                outputFile << "compute_label_byte_size()=" << compute_L_byte_size() << std::endl;

                outputFile.close();
            }

            void record_all_details_stream(std::ofstream &outputFile) {
                outputFile << "PLL info:" << std::endl;
                outputFile << "thread_num=" << thread_num << std::endl;
                outputFile << "compute_label_byte_size()=" << compute_L_byte_size() << std::endl;
            }

            long long int query(int source, int terminal, int t_s, int t_e) {
                if (source == terminal) {
                    return 0;
                }

                int distance = std::numeric_limits<int>::max();
                auto vector1_check_pointer = L[source].begin();
                auto vector2_check_pointer = L[terminal].begin();
                auto pointer_L_s_end = L[source].end(), pointer_L_t_end = L[terminal].end();

                for (auto vector1_begin = vector1_check_pointer; vector1_begin != pointer_L_s_end; ++vector1_begin) {
                    // cout << "x (" << vector1_begin->hub_vertex << "," << vector1_begin->hop << "," << vector1_begin->distance << "," << vector1_begin->parent_vertex << ") " << endl;
                    for (auto vector2_begin = vector2_check_pointer;
                         vector2_begin != pointer_L_t_end; ++vector2_begin) {
                        // cout << "y (" << vector2_begin->hub_vertex << "," << vector2_begin->hop << "," << vector2_begin->distance << "," << vector2_begin->parent_vertex << ") " << endl;
                        if (vector1_begin->vertex == vector2_begin->vertex &&
                            std::max(vector1_begin->t_s, std::max(vector2_begin->t_s, t_s)) <=
                            std::min(vector1_begin->t_e, std::min(vector2_begin->t_e, t_e))) {
                            long long int dis = (long long int) vector1_begin->distance + vector2_begin->distance;
                            if (distance > dis) {
                                distance = dis;
                                // cout << "x (" << vector1_begin->hub_vertex << "," << vector1_begin->hop << "," << vector1_begin->distance <<  ") " << endl;
                                // cout << "y (" << vector2_begin->hub_vertex << "," << vector2_begin->hop << "," << vector2_begin->distance <<  ") " << endl;
                            }
                        }
                    }
                }
                return distance;
            }
        };

        template<typename weight_type>
        inline two_hop_label<weight_type>::two_hop_label() {
            t_s = 0;
            t_e = INT_MAX;
            vertex = std::numeric_limits<int>::max();
            distance = std::numeric_limits<weight_type>::max();
        };

        template<typename weight_type>
        inline experiment::nonhop::two_hop_label<weight_type>::two_hop_label(const two_hop_label<weight_type> &other) {
            t_s = other.t_s;
            t_e = other.t_e;
            vertex = other.vertex;
            distance = other.distance;
        };

        template<typename weight_type>
        inline experiment::nonhop::two_hop_label<weight_type>::two_hop_label(int start_time) {
            t_s = start_time;
            t_e = INT_MAX;
            vertex = std::numeric_limits<int>::max();
            distance = std::numeric_limits<weight_type>::max();
        };
    };

    namespace hop {
        template<typename hop_weight_type>
        class two_hop_label {
        public:
            int hub_vertex, hop;
            hop_weight_type distance;
            int t_s, t_e;

            two_hop_label(const two_hop_label &other);

            two_hop_label();

            two_hop_label &operator=(const two_hop_label &) = default;

            bool operator<(const two_hop_label &other) const {
                if (distance != other.distance) {
                    return distance > other.distance; // < is the max-heap; > is the min heap
                }
                return hop > other.hop; // < is the max-heap; > is the min heap
            }

            void serialize(std::ofstream &out) const {
                experiment::saveBinary(out, hub_vertex);
                experiment::saveBinary(out, hop);
                experiment::saveBinary(out, distance);
                experiment::saveBinary(out, t_s);
                experiment::saveBinary(out, t_e);
            }

            void deserialize(std::ifstream &in) {
                experiment::loadBinary(in, hub_vertex);
                experiment::loadBinary(in, hop);
                experiment::loadBinary(in, distance);
                experiment::loadBinary(in, t_s);
                experiment::loadBinary(in, t_e);
            }
        };

        template<typename hop_weight_type>
        two_hop_label<hop_weight_type>::two_hop_label(const two_hop_label &other) {
            t_s = other.t_s;
            t_e = other.t_e;
            hub_vertex = other.hub_vertex;
            hop = other.hop;
            distance = other.distance;
        }

        template<typename hop_weight_type>
        two_hop_label<hop_weight_type>::two_hop_label() : t_s(0), t_e(std::numeric_limits<int>::max()) {
            hub_vertex = 0;
            hop = 0;
            distance = std::numeric_limits<hop_weight_type>::max();
        }

        template<typename weight_type>
        bool compare_hop_constrained_two_hop_label(two_hop_label<weight_type> &i, two_hop_label<weight_type> &j) {
            if (i.t_e != j.t_e) {
                return i.t_e > j.t_e;
            } else if (i.hub_vertex != j.hub_vertex) {
                return i.hub_vertex < j.hub_vertex;
            } else if (i.hop != j.hop) {
                return i.hop < j.hop;
            } else if (i.t_s != j.t_s) {
                return i.t_s < j.t_s;
            } else {
                return i.distance < j.distance;
            }
        };

        /** result is as \<dis,hop,hub\> .if dont find effective result, method will return the \<max,max,-1\> */
        template<typename hop_weight_type>
        std::tuple<hop_weight_type, int, int>
        graph_weighted_two_hop_extract_distance_and_hop_and_hub_in_current_with_csv(
                std::vector<two_hop_label<hop_weight_type>> &L_s, std::vector<two_hop_label<hop_weight_type>> &L_t,
                int source, int terminal,
                int hop_cst,
                result::MaintainShard &maintain_shard) {
            if (status::currentTimeMode == status::MaintainTimeMode::SLOT1) {
                ++maintain_shard.cover_count_slot1;
            } else {
                ++maintain_shard.cover_count_slot2;
            }
            if (source == terminal) {
                return {0, 0, source};
            }
            /*return std::numeric_limits<int>::max() is not connected*/
            int distance = std::numeric_limits<hop_weight_type>::max();
            int hop = std::numeric_limits<int>::max();
            int hub = -1;
            auto vector1_check_pointer = L_s.begin();
            auto vector2_check_pointer = L_t.begin();
            auto pointer_L_s_end = vector1_check_pointer, pointer_L_t_end = vector2_check_pointer;
            while (pointer_L_s_end != L_s.end() && pointer_L_s_end->t_e == std::numeric_limits<int>::max()) {
                ++pointer_L_s_end;
            }
            while (pointer_L_t_end != L_t.end() && pointer_L_t_end->t_e == std::numeric_limits<int>::max()) {
                ++pointer_L_t_end;
            }

            while (vector1_check_pointer != pointer_L_s_end && vector2_check_pointer != pointer_L_t_end) {
                if (vector1_check_pointer->hub_vertex == vector2_check_pointer->hub_vertex) {
                    auto vector1_end = vector1_check_pointer;
                    while (vector1_end != pointer_L_s_end &&
                           vector1_check_pointer->hub_vertex == vector1_end->hub_vertex &&
                           vector1_end->t_e == std::numeric_limits<int>::max()) {
                        ++vector1_end;
                    }
                    auto vector2_end = vector2_check_pointer;
                    while (vector2_end != pointer_L_t_end &&
                           vector2_check_pointer->hub_vertex == vector2_end->hub_vertex &&
                           vector2_end->t_e == std::numeric_limits<int>::max()) {
                        ++vector2_end;
                    }

                    for (auto vector1_begin = vector1_check_pointer; vector1_begin != vector1_end; ++vector1_begin) {
                        // cout << "x (" << vector1_begin->hub_vertex << "," << vector1_begin->hop << "," << vector1_begin->distance << "," << vector1_begin->parent_vertex << ") " << endl;
                        for (auto vector2_begin = vector2_check_pointer;
                             vector2_begin != vector2_end; ++vector2_begin) {
                            // cout << "y (" << vector2_begin->hub_vertex << "," << vector2_begin->hop << "," << vector2_begin->distance << "," << vector2_begin->parent_vertex << ") " << endl;
                            if (vector1_begin->hop + vector2_begin->hop <= hop_cst) {
                                hop_weight_type dis = vector1_begin->distance + vector2_begin->distance;
                                if (distance > dis) {
                                    distance = dis;
                                    hop = vector1_begin->hop + vector2_begin->hop;
                                    hub = vector1_begin->hub_vertex;
                                }
                            } else {
                                break;
                            }
                        }
                    }

                    vector1_check_pointer = vector1_end;
                    vector2_check_pointer = vector2_end;
                } else if (vector1_check_pointer->hub_vertex > vector2_check_pointer->hub_vertex) {
                    ++vector2_check_pointer;
                } else {
                    ++vector1_check_pointer;
                }
            }

            return {distance, hop, hub};
        };

        /**result is <dis,hop> .if dont find effective result, method will return the <max,max>*/
        template<typename hop_weight_type>
        std::pair<hop_weight_type, int>
        search_sorted_two_hop_label_in_current_with_csv(std::vector<two_hop_label<hop_weight_type>> &input_vector,
                                                        int key, result::MaintainShard &maintain_shard) {
            if (status::currentTimeMode == status::MaintainTimeMode::SLOT1) {
                ++maintain_shard.label_count_slot1;
            } else {
                ++maintain_shard.label_count_slot2;
            }
            int left = 0, right = static_cast<int>(input_vector.size()) - 1;
            hop_weight_type mindis = std::numeric_limits<hop_weight_type>::max();
            int hop_val = std::numeric_limits<int>::max();
            if (input_vector.empty()) {
                return {mindis, hop_val};
            }
//            std::cout << "diffuse 1.1.1" <<input_vector.size()<< std::endl;
            while (left < right) {
                int mid = (right - left) / 2 + left;
//                std::cout << "left is " << left <<" right is " << right <<" mid is "<< mid << std::endl;
                if (input_vector[mid].t_e != std::numeric_limits<int>::max()) {
                    right = mid - 1;
                } else {
                    if (input_vector[mid].hub_vertex < key) {
                        left = mid + 1;
                    } else {
                        // mindis = input_vector[mid].distance;
                        // hop_val = input_vector[mid].hop;
                        // left = mid + 1;
                        // if (hop_val == hop_k) {
                        //     return {mindis, hop_val};
                        // }
                        right = mid - 1;
                    }
                }
            }
//            std::cout << "diffuse 1.1.2" << std::endl;
            while (left < input_vector.size() &&
                   input_vector[left].t_e == std::numeric_limits<int>::max() &&
                   input_vector[left].hub_vertex == key) {
                if (input_vector[left].distance < mindis) {
                    mindis = input_vector[left].distance;
                    hop_val = input_vector[left].hop;
                }
                left++;
            }
            return {mindis, hop_val};
        }

        /**result is <dis,hop> .if dont find effective result, method will return the <max,max>*/
        template<typename hop_weight_type>
        std::pair<hop_weight_type, int>
        search_sorted_two_hop_label_in_current_with_less_than_k_limit_with_csv(
                std::vector<two_hop_label<hop_weight_type>> &input_vector, int key,
                int hop_k, result::MaintainShard &maintain_shard) {
            if (status::currentTimeMode == status::MaintainTimeMode::SLOT1) {
                ++maintain_shard.label_count_slot1;
            } else {
                ++maintain_shard.label_count_slot2;
            }
            int left = 0, right = static_cast<int>(input_vector.size()) - 1;
            hop_weight_type mindis = std::numeric_limits<hop_weight_type>::max();
            int hop_val = std::numeric_limits<int>::max();
            if (input_vector.empty()) {
                return {mindis, hop_val};
            }

            while (left < right) {
                int mid = (right - left) / 2 + left;
                if (input_vector[mid].t_e != std::numeric_limits<int>::max()) {
                    right = mid - 1;
                } else {
                    if (input_vector[mid].hub_vertex < key) {
                        left = mid + 1;
                    } else {
                        // mindis = input_vector[mid].distance;
                        // hop_val = input_vector[mid].hop;
                        // left = mid + 1;
                        // if (hop_val == hop_k) {
                        //     return {mindis, hop_val};
                        // }
                        right = mid - 1;
                    }
                }
            }
            while (left < input_vector.size() &&
                   input_vector[left].t_e == std::numeric_limits<int>::max() &&
                   input_vector[left].hub_vertex == key && input_vector[left].hop <= hop_k) {
                if (input_vector[left].distance < mindis) {
                    mindis = input_vector[left].distance;
                    hop_val = input_vector[left].hop;
                }
                left++;
            }
            return {mindis, hop_val};
        }

        /**result is <dis,hop> .if dont find effective result, method will return the <max,max>*/
        template<typename hop_weight_type>
        std::pair<hop_weight_type, int>
        search_sorted_two_hop_label_in_current_with_equal_k_limit_with_csv(
                std::vector<two_hop_label<hop_weight_type>> &input_vector, int key,
                int hop_k, result::MaintainShard &maintain_shard) {
            if (status::currentTimeMode == status::MaintainTimeMode::SLOT1) {
                ++maintain_shard.label_count_slot1;
            } else {
                ++maintain_shard.label_count_slot2;
            }
            int left = 0, right = input_vector.size() - 1;
            if (input_vector.empty()) {
                return {std::numeric_limits<hop_weight_type>::max(), std::numeric_limits<int>::max()};
            }
            while (left <= right) {
                int mid = (right - left) / 2 + left;
                if (input_vector[mid].t_e != std::numeric_limits<int>::max()) {
                    right = mid - 1;
                } else {
                    if (input_vector[mid].hub_vertex < key) {
                        left = mid + 1;
                    } else if (input_vector[mid].hub_vertex > key) {
                        right = mid - 1;
                    } else {
                        if (input_vector[mid].hop < hop_k) {
                            left = mid + 1;
                        } else if (input_vector[mid].hop > hop_k) {
                            right = mid - 1;
                        } else {
                            return std::make_pair(input_vector[mid].distance, input_vector[mid].hop);
                        }
                    }
                }
            }
            return {std::numeric_limits<hop_weight_type>::max(), std::numeric_limits<int>::max()};
        }

        template<typename hop_weight_type>
        void
        insert_sorted_hop_constrained_two_hop_label_with_csv(std::vector<two_hop_label<hop_weight_type>> &input_vector,
                                                             int key,
                                                             int hop, hop_weight_type new_distance, int t,
                                                             result::MaintainShard &maintain_shard) {
            if (experiment::status::currentTimeMode == experiment::status::MaintainTimeMode::SLOT1) {
                if (experiment::status::currentMaintainMode == experiment::status::MaintainMode::DECREASE) {
                    ++maintain_shard.label_decrease_insert_slot1;
                } else if (experiment::status::currentMaintainMode == experiment::status::MaintainMode::INCREASE) {
                    ++maintain_shard.label_increase_insert_slot1;
                }
            } else {
                if (experiment::status::currentMaintainMode == experiment::status::MaintainMode::DECREASE) {
                    ++maintain_shard.label_decrease_insert_slot2;
                } else if (experiment::status::currentMaintainMode == experiment::status::MaintainMode::INCREASE) {
                    ++maintain_shard.label_increase_insert_slot2;
                }
            }
            int left = 0, right = input_vector.size() - 1;
            if (input_vector.empty()) {
                if (new_distance != std::numeric_limits<hop_weight_type>::max()) {
                    two_hop_label<hop_weight_type> new_label;
                    new_label.hub_vertex = key;
                    new_label.hop = hop;
                    new_label.distance = new_distance;
                    new_label.t_s = t;
                    new_label.t_e = std::numeric_limits<int>::max();

                    input_vector.insert(input_vector.begin(), new_label);
                }
                return;
            }
            while (left <= right) {
                int mid = left + ((right - left) / 2);

                if (input_vector[mid].t_e == std::numeric_limits<int>::max()) {
                    if (input_vector[mid].hub_vertex == key) {
                        if (input_vector[mid].hop == hop) {
                            two_hop_label old_label = input_vector[mid];
                            old_label.t_e = t - 1;
                            if (input_vector[mid].distance < new_distance &&
                                new_distance != std::numeric_limits<hop_weight_type>::max()) {
//                                std::cout << "error input " << std::endl;
                            }
                            input_vector[mid].distance = new_distance;
                            input_vector[mid].t_s = t;
                            if (old_label.t_s != t &&
                                old_label.distance != std::numeric_limits<hop_weight_type>::max()) {
                                int insert_left = mid + 1, insert_right = input_vector.size() - 1;

                                while (insert_left <= insert_right) {
                                    int insert_mid = insert_left + ((insert_right - insert_left) / 2);

                                    if (input_vector[insert_mid].t_e < t - 1) {
                                        insert_right = insert_mid - 1;
                                    } else if (input_vector[insert_mid].t_e == t - 1) {
                                        if (input_vector[insert_mid].hub_vertex > key ||
                                            (input_vector[insert_mid].hub_vertex == key &&
                                             input_vector[insert_mid].hop > hop)) {
                                            insert_right = insert_mid - 1;
                                        } else {
                                            insert_left = insert_mid + 1;
                                        }
                                    } else {
                                        insert_left = insert_mid + 1;
                                    }
                                }
                                if (insert_left > input_vector.size()) {
                                    std::cerr << "insert_left out of range: " << insert_left << " > "
                                              << input_vector.size() << std::endl;
                                    return;
                                }
                                input_vector.insert(input_vector.begin() + insert_left, old_label);
                            }
                            return;
                        } else if (input_vector[mid].hop < hop) {
                            left = mid + 1;
                        } else {
                            right = mid - 1;
                        }
                    } else if (input_vector[mid].hub_vertex < key) {
                        left = mid + 1;
                    } else {
                        right = mid - 1;
                    }
                } else {
                    right = mid - 1;
                }
            }
            if (new_distance != std::numeric_limits<hop_weight_type>::max()) {
                two_hop_label<hop_weight_type> new_label;
                new_label.hub_vertex = key;
                new_label.hop = hop;
                new_label.distance = new_distance;
                new_label.t_s = t;
                new_label.t_e = std::numeric_limits<int>::max();

                input_vector.insert(input_vector.begin() + left, new_label);
            }
        }

        template<typename weight_type>
        class two_hop_case_info {
        public:
            /*hop bounded*/
            int thread_num = 1;
            int upper_k = 0;

            /*labels*/
            std::vector<std::vector<two_hop_label<weight_type>>> L;
            PPR_TYPE::PPR_type PPR;

            void serialize(std::ofstream &out) const {
                experiment::saveBinary(out, thread_num);
                experiment::saveBinary(out, upper_k);
                size_t size = L.size();
                BinarySerializer<size_t>::saveBinary(out, size);
                for (const std::vector<experiment::hop::two_hop_label<weight_type>> &item: L) {
                    BinarySerializer<std::vector<experiment::hop::two_hop_label<weight_type>>>::saveBinary(out, item);
                    out.flush();
                }
                experiment::saveBinary(out, PPR);
            }

            void deserialize(std::ifstream &in) {
                experiment::loadBinary(in, thread_num);
                experiment::loadBinary(in, upper_k);
                size_t size;
                BinarySerializer<size_t>::loadBinary(in, size);
                L.resize(size);
                for (std::vector<experiment::hop::two_hop_label<weight_type>> &item: L) {
                    BinarySerializer<std::vector<experiment::hop::two_hop_label<weight_type>>>::loadBinary(in, item);
                }
                experiment::loadBinary(in, PPR);
            }

            long long int compute_label_bit_size() {
                long long int size = 0;
                for (auto &xx: L) {
                    size = size + xx.size() * sizeof(two_hop_label<weight_type>);
                }
                return size;
            }

            /*clear labels*/
            void clear_labels() {
                std::vector<std::vector<two_hop_label<weight_type>>>().swap(L);
                PPR_TYPE::PPR_type().swap(PPR);
            }

            void print_L() {
                int index = 0;
                std::cout << "print_L: (hub_vertex, hop, distance)" << std::endl;
                for (auto &xx: L) {
                    std::cout << "vertex " << index++ << ": ";
                    for (auto &yy: xx) {
                        std::cout << "(" << yy.hub_vertex << "," << yy.hop << "," << yy.distance << "," << yy.t_s << ","
                                  << yy.t_e << ")";
                    }
                    std::cout << std::endl;
                }
            }

            void print_PPR() {
                std::cout << "print_PPR:" << std::endl;
                for (int i = 0; i < PPR.size(); i++) {
                    for (int j = 0; j < PPR[i].size(); j++) {
                        std::cout << "PPR(" << i << "," << PPR[i][j].first << "): ";
                        for (int k = 0; k < PPR[i][j].second.size(); k++) {
                            std::cout << PPR[i][j].second[k] << " ";
                        }
                        std::cout << std::endl;
                    }
                }
            }

            void print_L_vk(int v_k) {
                for (auto it = L[v_k].begin(); it != L[v_k].end(); it++) {
                    std::cout << "(" << it->hub_vertex << "," << it->hop << "," << it->distance << "," << it->t_s << ","
                              << it->t_e << ")";
                }
                std::cout << std::endl;
            }

            /*record_all_details*/
            void record_all_details(std::string save_name) {
                std::ofstream outputFile;
                outputFile.precision(6);
                outputFile.setf(std::ios::fixed);
                outputFile.setf(std::ios::showpoint);
                outputFile.open(save_name + ".txt");

                outputFile << "hop_constrained_case_info:" << std::endl;
                outputFile << "thread_num=" << thread_num << std::endl;
                outputFile << "upper_k=" << upper_k << std::endl;

                outputFile << "compute_label_bit_size()=" << compute_label_bit_size() << std::endl;

                outputFile.close();
            }

            void record_all_details_stream(std::ofstream &outputFile) {
                outputFile << "hop_constrained_case_info:" << std::endl;
                outputFile << "thread_num=" << thread_num << std::endl;
                outputFile << "upper_k=" << upper_k << std::endl;

                outputFile << "compute_label_bit_size()=" << compute_label_bit_size() << std::endl;
            }

            long long int query(int source, int terminal, int t_s, int t_e, int hop_cst) {
                /*return std::numeric_limits<int>::max() is not connected*/

                if (hop_cst < 0) {
                    return std::numeric_limits<int>::max();
                }
                if (source == terminal) {
                    return 0;
                } else if (hop_cst == 0) {
                    return std::numeric_limits<int>::max();
                }

                int distance = std::numeric_limits<int>::max();
                auto vector1_check_pointer = L[source].begin();
                auto vector2_check_pointer = L[terminal].begin();
                auto pointer_L_s_end = L[source].end(), pointer_L_t_end = L[terminal].end();

                for (auto vector1_begin = vector1_check_pointer; vector1_begin != pointer_L_s_end; vector1_begin++) {
                    // cout << "x (" << vector1_begin->hub_vertex << "," << vector1_begin->hop << "," << vector1_begin->distance << "," << vector1_begin->parent_vertex << ") " << endl;
                    for (auto vector2_begin = vector2_check_pointer;
                         vector2_begin != pointer_L_t_end; vector2_begin++) {
                        // cout << "y (" << vector2_begin->hub_vertex << "," << vector2_begin->hop << "," << vector2_begin->distance << "," << vector2_begin->parent_vertex << ") " << endl;
                        if (vector1_begin->hub_vertex == vector2_begin->hub_vertex &&
                            vector1_begin->hop + vector2_begin->hop <= hop_cst &&
                            std::max(vector1_begin->t_s, std::max(vector2_begin->t_s, t_s)) <=
                            std::min(vector1_begin->t_e, std::min(vector2_begin->t_e, t_e))) {
                            long long int dis = (long long int) vector1_begin->distance + vector2_begin->distance;
                            if (distance > dis) {
                                distance = dis;
                                // cout << "x (" << vector1_begin->hub_vertex << "," << vector1_begin->hop << "," << vector1_begin->distance <<  ") " << endl;
                                // cout << "y (" << vector2_begin->hub_vertex << "," << vector2_begin->hop << "," << vector2_begin->distance <<  ") " << endl;
                            }
                        }
                    }
                }
                return distance;
            }
        };

    }

    template<typename weight_type>
    class BinarySerializer<experiment::nonhop::two_hop_label<weight_type>> {
    public:
        static void saveBinary(std::ofstream &out, const experiment::nonhop::two_hop_label<weight_type> &label) {
            label.serialize(out);
        }

        static void loadBinary(std::ifstream &in, experiment::nonhop::two_hop_label<weight_type> &label) {
            label.deserialize(in);
        }
    };

    template<typename weight_type>
    class BinarySerializer<experiment::hop::two_hop_label<weight_type>> {
    public:
        static void saveBinary(std::ofstream &out, const experiment::hop::two_hop_label<weight_type> &label) {
            label.serialize(out);
        }

        static void loadBinary(std::ifstream &in, experiment::hop::two_hop_label<weight_type> &label) {
            label.deserialize(in);
        }
    };

    template<typename weight_type>
    class BinarySerializer<experiment::nonhop::two_hop_case_info<weight_type>> {
    public:
        static void saveBinary(std::ofstream &out, const experiment::nonhop::two_hop_case_info<weight_type> &vec) {
            vec.serialize(out);
        }

        static void loadBinary(std::ifstream &in, experiment::nonhop::two_hop_case_info<weight_type> &vec) {
            vec.deserialize(in);
        }
    };

    template<typename weight_type>
    class BinarySerializer<experiment::hop::two_hop_case_info<weight_type>> {
    public:
        static void saveBinary(std::ofstream &out, const experiment::hop::two_hop_case_info<weight_type> &vec) {
            vec.serialize(out);
        }

        static void loadBinary(std::ifstream &in, experiment::hop::two_hop_case_info<weight_type> &vec) {
            vec.deserialize(in);
        }
    };
}
