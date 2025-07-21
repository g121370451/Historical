#include "utils/global.h"

experiment::result::CSVConfig experiment::result::global_csv_config;
experiment::status::MaintainTimeMode experiment::status::currentTimeMode = experiment::status::MaintainTimeMode::SLOT1;
experiment::status::HopMode experiment::status::currentHopMode = experiment::status::HopMode::NoHop;
experiment::status::MaintainAlgorithmMode experiment::status::currentAlgorithmMode = experiment::status::MaintainAlgorithmMode::Algorithm2024;
experiment::status::MaintainMode experiment::status::currentMaintainMode = experiment::status::MaintainMode::DECREASE;
boost::random::mt19937 experiment::status::boost_random_time_seed{static_cast<std::uint32_t>(std::time(0))};
void experiment::result::init_config(const std::string &dataset,int thread_count,int iteration_count,int change_count,int hop_limit)
{
    global_csv_config.basic_data = {
        dataset,
        thread_count,
        iteration_count,
        change_count,
        hop_limit,
        "",
        0, 0};
    global_csv_config.ruc_counter.initialize(thread_count);
    global_csv_config.old_counter.initialize(thread_count);
    // Initialize with default values (can be modified later)
    global_csv_config.ruc_data = {0, 0, 0, 0, 0, 0, 0, 0};
    global_csv_config.old_data = {0, 0, 0, 0, 0, 0, 0, 0};
}
void experiment::result::FixedShardedCounter::merge_to(MaintainShard &target) const noexcept
{
    // 内存屏障更新
    std::atomic_thread_fence(std::memory_order_acquire);
    for (const auto &shard : this->maintain_shards_)
    {
        target.label_count_slot1 += shard.label_count_slot1;
        target.label_count_slot2 += shard.label_count_slot2;
        target.cover_count_slot1 += shard.cover_count_slot1;
        target.cover_count_slot2 += shard.cover_count_slot2;
        target.label_decrease_insert_slot1 += shard.label_decrease_insert_slot1;
        target.label_decrease_insert_slot2 += shard.label_decrease_insert_slot2;
        target.label_increase_insert_slot1 += shard.label_increase_insert_slot1;
        target.label_increase_insert_slot2 += shard.label_increase_insert_slot2;
        target.ppr_insert_slot1 += shard.ppr_insert_slot1;
        target.ppr_insert_slot2 += shard.ppr_insert_slot2;
        target.diffuse_count_slot1 += shard.diffuse_count_slot1;
        target.diffuse_count_slot2 += shard.diffuse_count_slot2;
    }
}

double experiment::result::MethodData::total_label_count() const noexcept
{
    return this->label_count_slot1 + this->label_count_slot2;
}

double experiment::result::MethodData::total_cover_count() const noexcept
{
    std::cout << "cover is " << this->cover_count_slot1 << "+" << this->cover_count_slot2 << std::endl;
    return this->cover_count_slot1 + this->cover_count_slot2;
}

double experiment::result::MethodData::total_ppr_insert_count() const noexcept
{
    return this->ppr_insert_slot1 + this->ppr_insert_slot2;
}

double experiment::result::MethodData::total_increase_label_insert_count() const noexcept
{
    return this->label_increase_insert_slot1 + this->label_increase_insert_slot2;
}

double experiment::result::MethodData::total_decrease_label_insert_count() const noexcept
{
    return this->label_decrease_insert_slot1 + this->label_decrease_insert_slot2;
}

double experiment::result::MethodData::total_diffuse_count() const noexcept
{
    return this->diffuse_count_slot1 + this->diffuse_count_slot2; 
}

experiment::result::MethodData::operator MaintainShard &()
{
    return *reinterpret_cast<MaintainShard *>(this);
}