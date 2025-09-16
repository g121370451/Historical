#pragma once
#include <string>
#include "third/spdlog/spdlog.h"
class LogTransaction {
public:
    void log(const std::string& msg);

    void commit();

    void rollback();

private:
    std::vector<std::string> buffer;
};