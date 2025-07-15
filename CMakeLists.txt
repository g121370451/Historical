cmake_minimum_required(VERSION 3.12)
project(HistoricalHopHL)

set(CMAKE_CXX_STANDARD 20)
# 在顶部添加调试模式配置
# 修改编译选项部分
if(DEBUG_MODE)
    message(STATUS "Debug mode enabled")
    set(CMAKE_BUILD_TYPE Debug)
    add_compile_options(-g3 -O0 -ggdb -pthread)  # 更详细的调试信息
    add_compile_definitions(_DEBUG)
else()
    message(STATUS "Release mode enabled")
    set(CMAKE_BUILD_TYPE Release)
    add_compile_options(-O2 -pthread)
endif()
find_package(Boost REQUIRED COMPONENTS random)

if(Boost_FOUND)
    message(STATUS "Boost found: ${Boost_VERSION}")
    message(STATUS "Include dirs: ${Boost_INCLUDE_DIRS}")
    message(STATUS "Libraries: ${Boost_LIBRARIES}")
    message(STATUS "Boost库目录: ${Boost_LIBRARY_DIRS}")           # 库文件所在目录
else()
    message(FATAL_ERROR "Boost NOT found!")
endif()

include_directories(
        ${PROJECT_SOURCE_DIR}/Historical/experiment/include
        ${Boost_INCLUDE_DIRS}
)

# 打印包含目录
message(STATUS "Include directories for ${TARGET_NAME}:")
if(${INCLUDE_DIRECTORIES})
    foreach(dir ${INCLUDE_DIRECTORIES})
        message(STATUS "  ${dir}")
    endforeach()
else()
    message(STATUS "  No include directories found")
endif()
#查找Python 3
add_executable(experiment
        ${PROJECT_SOURCE_DIR}/Historical/experiment/main.cpp  # 主入口文件
        ${PROJECT_SOURCE_DIR}/Historical/experiment/include/entity/two_hop_label.cpp
        ${PROJECT_SOURCE_DIR}/Historical/experiment/include/parse/experiment_argparse.cpp  # 头文件的实现
        ${PROJECT_SOURCE_DIR}/Historical/experiment/include/utils/StringHelper.cpp
        ${PROJECT_SOURCE_DIR}/Historical/experiment/include/utils/ExecutionTimer.cpp
        ${PROJECT_SOURCE_DIR}/Historical/experiment/include/utils/global.cpp
        ${PROJECT_SOURCE_DIR}/Historical/experiment/include/utils/CSVWriter.cpp
)
target_compile_options(experiment PRIVATE
        $<$<CONFIG:Debug>:-DDEBUG>  # 定义DEBUG宏
        -Wall -Wextra -Wpedantic    # 更严格的警告
)
# 链接库
target_link_libraries(experiment
        ${Boost_LIBRARIES}
        pthread
)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/experiment)
