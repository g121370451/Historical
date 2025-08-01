#include "StringHelper.h"
#include <cstdio>
#include <ctime>
#include <algorithm>

std::vector<std::string> experiment::parse_string(std::string parse_target, std::string delimiter) {
    std::vector<std::string> Parsed_content;
    size_t pos = 0;
    std::string token;
    while ((pos = parse_target.find(delimiter)) != std::string::npos) {
        // find(const string& str, size_t pos = 0) function returns the position of the first occurrence of str in the string, or npos if the string is not found.
        token = parse_target.substr(0, pos);
        // The substr(size_t pos = 0, size_t n = npos) function returns a substring of the object, starting at position pos and of length npos
        Parsed_content.push_back(token);                 // store the subtr to the list
        parse_target.erase(0, pos + delimiter.length()); // remove the front substr and the first delimiter
    }
    Parsed_content.push_back(parse_target); // store the subtr to the list

    return Parsed_content;
}

std::string experiment::get_current_time_string() {
    time_t now = time(nullptr);
    tm *curr_tm = localtime(&now);
    char time[80] = {0};
    strftime(time, 80, "%Y-%m-%d-%H-%M-%S", curr_tm);
    return time;
}

std::string experiment::toLower(const std::string &s) {
    std::string out = s;
    std::transform(out.begin(), out.end(), out.begin(), ::tolower);
    return out;
}
