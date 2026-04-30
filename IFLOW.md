# 🫚 ginger-cpp 项目上下文

## 项目概述

`ginger-cpp` 是一个用 C++ 实现的多项式求根算法库，支持并行化处理。该项目基于现代 CMake 实践构建，提供了多种多项式求根算法的实现，包括：

- **Aberth-Ehrlich 方法**：用于查找多项式根的数值方法，支持单线程和多线程版本
- **Bairstow 方法**：用于查找多项式根的另一种数值方法
- **自相关函数的特殊处理**：针对自相关函数的优化版本

## 构建与运行

### 依赖项
- CMake (版本 3.14...3.22)
- C++17 编译器
- 通过 CPM.cmake 管理的依赖项，如 fmt 库 (版本 12.1.0)

### 构建命令

#### 构建独立可执行文件
```bash
cmake -S standalone -B build/standalone
cmake --build build/standalone
./build/standalone/Ginger --help
```

#### 构建并运行测试套件
```bash
cmake -S test -B build/test
cmake --build build/test
CTEST_OUTPUT_ON_FAILURE=1 cmake --build build/test --target test
# 或直接运行可执行文件
./build/test/GingerTests
```

#### 构建基准测试
```bash
cmake -S bench -B build/bench
cmake --build build/bench
./build/bench/BenchmarkTarget
```

#### 一键构建所有内容
```bash
cmake -S all -B build
cmake --build build
```

### 代码格式化
```bash
# 检查格式
cmake --build build/test --target format
# 应用格式
cmake --build build/test --target fix-format
```

## 项目结构

- `include/ginger/` - 头文件目录，包含库的公共接口
  - `aberth.hpp` - Aberth-Ehrlich 方法的声明
  - `bairstow.hpp` - Bairstow 方法的声明
  - `config.hpp` - 配置类（Options）的定义
  - `matrix2.hpp`, `vector2.hpp` - 矩阵和向量的实现
- `source/` - 源文件目录，包含库的实现
  - `aberth.cpp` - Aberth-Ehrlich 方法的实现
  - `autocorr.cpp` - 自相关函数相关实现
  - `rootfinding.cpp` - 根查找算法的通用实现
- `test/source/` - 测试文件目录
  - `test_aberth.cpp` - Aberth 方法的测试
  - `test_autocorr.cpp` - 自相关函数的测试
  - `test_matrix2.cpp` - 矩阵相关测试
  - `test_rootfinding.cpp` - 根查找算法测试
  - `test_stress.cpp` - 压力测试
  - `test_vector2_ref.cpp`, `test_vector2.cpp` - 向量相关测试
- `bench/` - 基准测试目录
  - `BM_aberth.cpp` - Aberth 方法的基准测试
  - `BM_autocorr.cpp` - 自相关函数的基准测试
  - `BM_fir.cpp` - FIR 滤波器相关的基准测试
- `standalone/` - 独立可执行文件目录
- `documentation/` - 文档目录，使用 Doxygen 生成文档
- `cmake/` - CMake 模块，包含 CPM.cmake、doctest 配置等

## 开发约定

- **C++17**：使用 C++17 标准
- **现代 CMake**：遵循现代 CMake 实践
- **代码格式化**：使用 clang-format 和 cmake-format 进行代码格式化
- **测试**：集成 doctest 测试框架
- **并行化**：使用 ThreadPool.h 实现多线程处理
- **依赖管理**：使用 CPM.cmake 进行依赖管理

## 核心功能

- **Options 类**：包含 `max_iters`（最大迭代次数，默认 2000）和 `tolerance`（容差，默认 1e-14）
- **Aberth 方法**：包括单线程和多线程版本，以及针对自相关函数的特殊版本
- **Bairstow 方法**：实现 Bairstow 方法用于求多项式的根
- **Horner 方法**：实现 Horner 规则用于多项式求值

## 特殊工具

- **静态分析器**：支持 clang-tidy、iwyu、cppcheck
- **Sanitizers**：支持 AddressSanitizer、MemorySanitizer、UndefinedBehaviorSanitizer 等
- **Ccache**：支持使用 ccache 加速编译

## 注意事项

- 项目使用 glob 方式自动包含源文件，这在 CMake 中被认为不是最佳实践，但提供了简单性
- 项目模板设计为易于扩展，可以轻松添加新的算法或修改现有算法
- 项目支持安装，并可通过 find_package() 找到