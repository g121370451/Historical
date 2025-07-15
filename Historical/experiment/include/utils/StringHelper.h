#pragma once

#include <vector>
#include <string>

namespace experiment{
	std::vector<std::string> parse_string(std::string parse_target, std::string delimiter);
	std::string get_current_time_string();
    std::string toLower(const std::string& s);
}