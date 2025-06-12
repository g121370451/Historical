#pragma once
#include <filesystem>
#include <string>

namespace experiment
{
    enum Mode
	{
        MODE_ERROR,
		GENERATE_LABEL,
		MAINTAIN_LABEL,
		QUERY_RESULT
	};

	struct ExperimentConfig
	{
		// 成员变量需要默认初始化
		Mode mode = MODE_ERROR;
		int threads = 0;
		std::filesystem::path data_source{};
		std::filesystem::path save_path{};
		std::string datasetName{};
		int hop_limit = 0;

		// maintain-label param
		int iterations = 0;
		int change_count = 0;
		int max_value = 0;
		int min_value = 0;
		// 构造函数需要初始化所有成员
		ExperimentConfig();
	};
    void parse_arguments(ExperimentConfig *config,int argc, char *argv[]);
} // namespace experiment
