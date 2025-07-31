#include "parse/experiment_argparse.h"
#include "third/argparse.h"
#include <iostream>

void experiment::parse_arguments(ExperimentConfig *config, int argc, char *argv[])
{
    std::cout << "argc: " << argc << std::endl;
    // 如果没有参数，设置默认值并返回
    if (argc <= 1) {
        std::cout << "No arguments provided. Using default values." << std::endl;
        config->mode = Mode::MODE_ERROR; // 或者其他默认模式;
        return;
    }

    argparse::ArgumentParser program("experiment");

    // // generate-label
    argparse::ArgumentParser generate_label("generate-label");
    generate_label.add_argument("-t", "--threads").required().scan<'i', int>();
    generate_label.add_argument("-f", "--data_source").required();
    generate_label.add_argument("-p", "--save_path").required();
    generate_label.add_argument("-max", "--max_value").required().scan<'i', int>();
    generate_label.add_argument("-min", "--min_value").required().scan<'i', int>();
    generate_label.add_argument("-k", "--hop_limit").required().default_value(0).scan<'i', int>();

    // maintain-label
    argparse::ArgumentParser maintain_label("maintain-label");
    maintain_label.add_argument("-t", "--threads").required().scan<'i', int>().help("开启的线程数");
    maintain_label.add_argument("-f", "--data_source").required().help("数据来源");
    maintain_label.add_argument("-p", "--save_path").required().help("结果集的保存路径，包括统计数据、图、索引数据");
    maintain_label.add_argument("-k", "--hop_limit").required().scan<'i', int>().help("k-hop限制，如果k等于0使用non-hop算法，反之使用with-hop算法");
    maintain_label.add_argument("-m", "--iterations").required().scan<'i', int>().help("演化图的演化次数");
    maintain_label.add_argument("-c", "--change_count").required().scan<'i', int>().help("每次演化的边变化次数");
    maintain_label.add_argument("-max", "--max_value").required().scan<'i', int>().help("变变化的最大值");
    maintain_label.add_argument("-min", "--min_value").required().scan<'i', int>().help("变变化的最小值 不为负数");
    maintain_label.add_argument("-n", "--dataset").required().help("csv统计中的数据集记录值");
    maintain_label.add_argument("--check")
            .default_value(false)
            .implicit_value(true)
            .nargs(0)
            .help("开启正确性检查");

    maintain_label.add_argument("--read_old")
            .default_value("")
            .help("读取旧数据路径，若为空则不读取");

    maintain_label.add_argument("--strategy")
            .default_value("")
            .help("边变化策略，如：high_high_increase, low_low_mixed 等");

//  query-result
    argparse::ArgumentParser query_label("query-result");
    query_label.add_argument("-t", "--threads").required().scan<'i', int>().help("开启的线程数");
    query_label.add_argument("-f", "--data_source").required().help("数据来源");
    query_label.add_argument("-p", "--save_path").required().help("结果csv的保存路径 追加写");
    query_label.add_argument("-k", "--hop_limit").required().scan<'i', int>().help("k-hop限制，查找不同hop目录下的二进制文件");
    query_label.add_argument("-c", "--query_count").required().scan<'i', int>().help("每次演化的边变化次数");

    program.add_subparser(generate_label);
    program.add_subparser(maintain_label);
    program.add_subparser(query_label);

    try
    {
        program.parse_args(argc, argv);
    }
    catch (const std::runtime_error &err)
    {
        std::cerr << "Error: " << err.what() << std::endl;
        std::cerr << program;
        return;
    }
    if (program.is_subcommand_used("generate-label"))
    {
        config->mode = GENERATE_LABEL;
        config->threads = generate_label.get<int>("-t");
        config->data_source = std::filesystem::path(generate_label.get<std::string>("-f"));
        config->save_path = std::filesystem::path(generate_label.get<std::string>("-p"));
        config->hop_limit = generate_label.get<int>("-k");
        config->max_value = generate_label.get<int>("-max");
        config->min_value = generate_label.get<int>("-min");
        if (config->hop_limit < 0)
        {
            throw std::invalid_argument("Error: hop_constrained (-k) must be >= 0.");
        }
    }
    else if (program.is_subcommand_used("maintain-label"))
    {
        config->mode = MAINTAIN_LABEL;
        config->threads = maintain_label.get<int>("-t");
        config->data_source = maintain_label.get<std::string>("-f");
        config->save_path = maintain_label.get<std::string>("-p");
        config->hop_limit = maintain_label.get<int>("-k");
        config->iterations = maintain_label.get<int>("-m");
        config->change_count = maintain_label.get<int>("-c");
        config->max_value = maintain_label.get<int>("-max");
        config->min_value = maintain_label.get<int>("-min");
        config->datasetName = maintain_label.get<std::string>("-n");
        config->enableCorrectnessCheck = maintain_label.get<bool>("--check");
        auto readOldPath = maintain_label.get<std::string>("--read_old");
        if (!readOldPath.empty()) {
            config->generatedFilePath = readOldPath;
        }
        auto strategyInput = maintain_label.get<std::string>("--strategy");
        if (!strategyInput.empty()) {
            config->changeStrategy = strategyInput;
        }
        if (config->max_value <= 0 || config->min_value <= 0)
        {
            throw std::invalid_argument("Error: max_value and min_value must be greater than 0.");
        }
    }
    else if (program.is_subcommand_used("query-result"))
    {
        config->mode = QUERY_RESULT;
        config->threads = query_label.get<int>("-t");
        config->data_source = query_label.get<std::string>("-f");
        config->save_path = query_label.get<std::string>("-p");
        config->hop_limit = query_label.get<int>("-k");
        config->change_count = query_label.get<int>("-c");
    }
    else
    {
        std::cerr << "Error: Unknown subcommand." <<std::endl;
        std::cerr << program;
    }
}

experiment::ExperimentConfig::ExperimentConfig()
{
    std::cout << "ExperimentConfig constructor called" << std::endl;
}
