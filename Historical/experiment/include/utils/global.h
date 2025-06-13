#pragma once

#include <string>
#include <atomic>
#include <array>
#include <thread>
#include <functional>
#include <mutex>
#include <boost/random.hpp>

namespace experiment {
    namespace result {
        struct basicData {
            std::string dataset;
            std::string changeName;
            double ruc_time_slot1 = 0;
            double ruc_time_slot2 = 0;
            double a2021_time_slot1 = 0;
            double a2021_time_slot2 = 0;
        };

        // 对齐缓存行的分片数据单元
        struct alignas(128) MaintainShard {
            int label_count_slot1 = 0; //4
            int label_count_slot2 = 0; //4
            int cover_count_slot1 = 0; //4
            int cover_count_slot2 = 0; // 4
            size_t ppr_insert_slot1 = 0; // 8
            size_t ppr_insert_slot2 = 0; // 8
            size_t label_decrease_insert_slot1 = 0; // 8
            size_t label_increase_insert_slot1 = 0; // 8
            size_t label_decrease_insert_slot2 = 0; // 8
            size_t label_increase_insert_slot2 = 0; // 8
            size_t diffuse_count_slot1 = 0; // 8
            size_t diffuse_count_slot2 = 0; // 8
            char padding[128 - 80]{};
        };
        struct alignas(64) QueryShard {
            long long int compare_count = 0;
            char padding[64 - (1 * sizeof(long long int))]{}; // 显式填充
        };

        static_assert(sizeof(MaintainShard) == 128, "Maintain MethodShard size mismatch");
        static_assert(sizeof(QueryShard) == 64, "Query MethodShard size mismatch");

        // 分片计数器管理器（固定线程数优化版）
        class FixedShardedCounter {
        private:
            std::vector<MaintainShard> maintain_shards_; // 动态数量分片
            std::vector<QueryShard> query_shards_; // 动态数量分片
            std::once_flag init_flag_;

        public:
            // 初始化固定分片数（线程安全）
            void initialize(size_t thread_count) {
                std::call_once(init_flag_, [=, this] {
                    maintain_shards_.resize(thread_count);
                });
            }

            // 获取线程专属分片（需先初始化）
            MaintainShard &get_thread_maintain_shard(size_t thread_idx) noexcept {
                return maintain_shards_.at(thread_idx); // 带边界检查
            }

            QueryShard &get_thread_query_shard(size_t thread_idx) noexcept {
                return query_shards_.at(thread_idx); // 带边界检查
            }

            // 合并所有分片数据
            void merge_to(MaintainShard &target) const noexcept;
        };

        // 对外暴露的统计接口
        struct MethodData {
            int label_count_slot1 = 0;
            int label_count_slot2 = 0;
            int cover_count_slot1 = 0;
            int cover_count_slot2 = 0;
            size_t ppr_insert_slot1 = 0;
            size_t ppr_insert_slot2 = 0;
            size_t label_decrease_insert_slot1 = 0;
            size_t label_increase_insert_slot1 = 0;
            size_t label_decrease_insert_slot2 = 0;
            size_t label_increase_insert_slot2 = 0;
            // 记录diffuse中向外查询的数量
            size_t diffuse_count_slot1 = 0;
            size_t diffuse_count_slot2 = 0;

            double total_label_count() const noexcept;

            double total_cover_count() const noexcept;

            double total_ppr_insert_count() const noexcept;

            double total_increase_label_insert_count() const noexcept;

            double total_decrease_label_insert_count() const noexcept;

            double total_diffuse_count() const noexcept;

            explicit operator MaintainShard &();
        };

        // 全局配置
        struct CSVConfig {
            basicData basic_data;    // 不可变基础数据
            FixedShardedCounter ruc_counter; // RUC方法计数器
            FixedShardedCounter old_counter; // 旧方法计数器
            MethodData ruc_data;         // 合并后的RUC结果
            MethodData old_data;         // 合并后的旧方法结果
        };

        // Global configuration instance
        extern CSVConfig global_csv_config;

        // Initialization function
        void init_config(const std::string &dataset, int threadNum);

    } // namespace config
    namespace status {
        enum MaintainTimeMode {
            SLOT1 = 0,
            SLOT2,
        };
        enum class HopMode {
            WithHop,
            NoHop
        };
        enum class MaintainAlgorithmMode {
            Algorithm2021,
            Algorithm2024,
        };
        enum class MaintainMode {
            DECREASE,
            INCREASE
        };
        extern MaintainTimeMode currentTimeMode;
        extern HopMode currentHopMode;
        extern MaintainAlgorithmMode currentAlgorithmMode;
        extern MaintainMode currentMaintainMode;
        extern boost::random::mt19937 boost_random_time_seed;
    }
} // namespace experiment