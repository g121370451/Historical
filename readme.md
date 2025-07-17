# Historical 2-hop label maintain

## 项目结构说明

源代码集中在 `Historical/` 中，构建文件与可执行程序分别输出至 `cmake-build-*` 和 `bin/` 目录。

其中：

- `data/` 用于组织不同数据集，按数据集名称（如 `enron-email/`）建立子目录，每个子目录下包含：
    - `raw/`：原始图数据、边变化文件、初始 2-hop 索引、分析图等元数据；
    - `processed/`：由程序维护后生成的 2-hop 索引数据，仅保存一次运行结果（默认覆盖）；

  数据结构统一便于程序识别与批处理。

- `output/` 用于存储程序运行生成的结果，包括：
    - 结果统计（如变化影响分析的 CSV 文件）；
    - Debug 模式下的检查输出（如标签正确性验证日志）；
    - 其他辅助输出信息（如中间状态、运行摘要等）。

- `build-*.bat/build-*.sh` 用于调用cmake构建项目
- `run-generator.bat/run-generator.sh` 用于调用generator生成原始图的初始2hop标签
- 

## 项目构建与运行
