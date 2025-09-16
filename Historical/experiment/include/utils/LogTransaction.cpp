//
// Created by gpy on 2025/9/15.
//

#include "LogTransaction.h"

void LogTransaction::rollback() {
    buffer.clear();
}

void LogTransaction::log(const std::string &msg) {
    buffer.emplace_back(msg);
}

void LogTransaction::commit() {
    for (auto& msg : buffer) {
        spdlog::info("{}", msg);
    }
    buffer.clear();
}
