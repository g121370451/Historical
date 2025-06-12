#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <chrono>
#include <iomanip>
#include <sstream>
#include <memory>

namespace experiment
{

	struct Subtask
	{
		std::string name;
		std::chrono::steady_clock::time_point startTime;
		std::chrono::steady_clock::time_point endTime;
		std::vector<std::pair<std::string, std::shared_ptr<Subtask>>> subtasks;
		std::weak_ptr<Subtask> parent;
		bool isEnded = false;

		inline Subtask() : startTime(std::chrono::steady_clock::now())
		{
			name = "task";
		}

		inline Subtask(const std::string &taskName, std::weak_ptr<Subtask> parentTask = std::weak_ptr<Subtask>())
			: name(taskName), startTime(std::chrono::steady_clock::now()), parent(parentTask)
		{
		}

		inline void end()
		{
			endTime = std::chrono::steady_clock::now();
			isEnded = true;
		}

		inline double getDuration() const
		{
			if (!isEnded)
				return 0.0;
			return std::chrono::duration_cast<std::chrono::duration<double>>(endTime - startTime).count();
		}

		inline double getTotalDuration() const
		{
			if (!subtasks.empty())
			{
				double totalDuration = 0.0;
				for (const auto &subtask : subtasks)
				{
					totalDuration += subtask.second->getTotalDuration();
				}
				return totalDuration;
			}
			return getDuration();
		}
	};

	class ExecutionTimer
	{
	public:
		inline void startTask(const std::string &taskName)
		{
			currentMainTaskName = taskName;
			currentSubtask = nullptr;
			allTasks.clear();
		}

		inline void startSubtask(const std::string &subtaskName)
		{
			std::cout << subtaskName << std::endl;
			auto newSubtask = std::make_shared<Subtask>(subtaskName, currentSubtask);

			if (currentSubtask)
			{
				currentSubtask->subtasks.emplace_back(subtaskName, newSubtask);
			}
			else
			{
				allTasks.emplace_back(subtaskName, newSubtask);
			}
			currentSubtask = newSubtask;
		}

		inline double endSubtask()
		{
			if (currentSubtask)
			{
				currentSubtask->end();
				double res = currentSubtask->getTotalDuration();
				auto parentTask = currentSubtask->parent.lock();
				currentSubtask = parentTask ? parentTask : nullptr;
				return res;
			}
			else
			{
				std::cerr << "Error: No subtask is currently running.\n";
				return 0;
			}
		}

		inline double getTaskDuration() const
		{
			double totalDuration = 0;
			for (const auto &task : allTasks)
			{
				totalDuration += task.second->getTotalDuration();
			}
			return totalDuration;
		}

		inline void printStats()
		{
			std::cout << "Task Timing Statistics: " << getTaskDuration() << " seconds\n";
			printTaskStats(allTasks, 0);
		}

		inline void writeStatsToFile(std::ofstream &outFile)
		{
			if (!outFile.is_open())
			{
				std::cerr << "Error opening file for writing." << std::endl;
				return;
			}

			outFile.precision(6);
			outFile.setf(std::ios::fixed);
			outFile.setf(std::ios::showpoint);
			outFile << "Task Timing Statistics:" << std::endl;
			writeStatsToFileRecursively(outFile, allTasks, 0);
		}

	private:
		void writeStatsToFileRecursively(std::ofstream &outFile, const std::vector<std::pair<std::string, std::shared_ptr<Subtask>>> &tasks, int level);

		void printTaskStats( const std::vector<std::pair<std::string, std::shared_ptr<Subtask>>> &tasks, int level);
	private:
		std::string currentMainTaskName;
		std::shared_ptr<Subtask> currentSubtask;
		std::vector<std::pair<std::string, std::shared_ptr<Subtask>>> allTasks;
	};

}
