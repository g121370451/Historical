#include "ExecutionTimer.h"

void experiment::ExecutionTimer::writeStatsToFileRecursively(std::ofstream &outFile,  const std::vector<std::pair<std::string, std::shared_ptr<Subtask>>> &tasks, int level)
{
    for (const auto &task : tasks)
    {
        outFile << std::string(level * 4, ' ') << "-> " << task.first << ", Time: "
                << task.second->getTotalDuration() << " seconds" << std::endl;
        writeStatsToFileRecursively(outFile, task.second->subtasks, level + 1);
    }
}

void experiment::ExecutionTimer::printTaskStats(const std::vector<std::pair<std::string, std::shared_ptr<Subtask>>> &tasks, int level)
{
    for (const auto &task : tasks)
    {
        std::cout << std::string(level * 4, ' ') << "-> " << task.first << ": "
                  << std::fixed << std::setprecision(3) << task.second->getTotalDuration() << "s" << std::endl;
        printTaskStats(task.second->subtasks, level + 1);
    }
}