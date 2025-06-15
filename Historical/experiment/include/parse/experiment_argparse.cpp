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
    maintain_label.add_argument("-t", "--threads").required().scan<'i', int>();
    maintain_label.add_argument("-f", "--data_source").required();
    maintain_label.add_argument("-p", "--save_path").required();
    maintain_label.add_argument("-k", "--hop_limit").required().scan<'i', int>();
    maintain_label.add_argument("-m", "--iterations").required().scan<'i', int>();
    maintain_label.add_argument("-c", "--change_count").required().scan<'i', int>();
    maintain_label.add_argument("-max", "--max_value").required().scan<'i', int>();
    maintain_label.add_argument("-min", "--min_value").required().scan<'i', int>();
    maintain_label.add_argument("-n", "--dataset").required();

    // query-result
    argparse::ArgumentParser query_label("query-result");
    query_label.add_argument("-f", "--data_source").required();
    query_label.add_argument("-t", "--threads").required().scan<'i', int>();
    query_label.add_argument("-c", "--search_count").required().scan<'i', int>();
    query_label.add_argument("-k", "--hop_limit").required().scan<'i', int>();

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
        if (config->max_value <= 0 || config->min_value <= 0)
        {
            throw std::invalid_argument("Error: max_value and min_value must be greater than 0.");
        }
    }
    else if (program.is_subcommand_used("query-result"))
    {
        config->mode = QUERY_RESULT;
        config->data_source = query_label.get<std::string>("-f");
        config->change_count = query_label.get<int>("-c");
        config->hop_limit = query_label.get<int>("-k");
        config->threads = query_label.get<int>("-t");
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
