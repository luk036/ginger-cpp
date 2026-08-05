# AGENTS.md - Agentic Coding Guidelines for ginger-cpp

## Project Overview

ginger-cpp is a C++ library for polynomial root-finding algorithms (parallelizable). It uses modern CMake practices, Google Test-style (doctest) for unit testing, and RapidCheck for property-based testing.

---

## Build Commands

### Build the Library
```bash
cmake -B build
cmake --build build
```

### Build and Run All Tests
```bash
cmake -B build
cmake --build build
ctest --test-dir build --output-on-failure

# Or run the executable directly:
./build/GingerTests
```

### Run a Single Test (doctest filter)
```bash
# Run specific test case
./build/GingerTests -tc="test_case_name"

# Run specific test suite
./build/GingerTests -ts="test_suite_name"

# List all test cases
./build/GingerTests --list-test-cases
```

### Build with Coverage (Linux/macOS)
```bash
cmake -B build -DGINGER_ENABLE_COVERAGE=1 -DCMAKE_BUILD_TYPE=Debug
cmake --build build
./build/GingerTests
```

### Build Standalone
```bash
cmake -B build
cmake --build build
./build/Ginger --help
```

---

## Code Formatting & Linting

### Check Format
```bash
cmake -B build
cmake --build build --target format
```

### Auto-Fix Format
```bash
cmake -B build
cmake --build build --target fix-format
```

Requirements: `clang-format==18.1.2`, `cmake-format==0.6.13`, `pyyaml`

---

## Code Style Guidelines

### File Organization
- **Headers**: `include/ginger/*.hpp` (or `.h` for C interop)
- **Sources**: `source/*.cpp`
- **Tests**: `test/source/*.cpp`

### Naming Conventions
- **Files**: lowercase with underscores: `rootfinding.hpp`, `vector2_ref.hpp`
- **Classes/Types**: PascalCase: `class Options`, `struct Vector2`
- **Functions**: snake_case: `initial_guess()`, `pbairstow_even()`
- **Variables**: snake_case: `max_iters`, `tolerance`
- **Constants**: UPPER_SNAKE_CASE for enum values, camelCase otherwise

### Header Includes
```cpp
#pragma once

#include <vector>        // STL
#include <utility>       // STL

#include "matrix2.hpp"  // Local (quoted)
#include "vector2.hpp"  // Local (quoted)

// Use include order: 1. Standard library, 2. External libs, 3. Local
```

### C++ Standard
- **Minimum**: C++17
- **Target**: Modern C++ (prefer `<filesystem>`, structured bindings, `std::optional`, etc.)

### Type Usage
- Prefer STL containers: `std::vector<T>`, `std::optional<T>`
- Use aliases for common types:
  ```cpp
  using Vec2 = ginger::Vector2<double>;
  using Mat2 = ginger::Matrix2<Vec2>;
  ```
- Use `auto` for return type deduction where clear

### Documentation (Doxygen)
```cpp
/**
 * @brief Brief description
 *
 * Longer description if needed.
 *
 * @param[in] param_name Description
 * @param[out] param_name Description
 * @return Return description
 */
extern auto function_name(const std::vector<double> &coeffs) -> std::vector<Vec2>;
```

### Formatting (per .clang-format)
- **Base Style**: Google
- **Indent Width**: 4 spaces
- **Column Limit**: 100
- **Brace Style**: Attach
- **Namespace Indentation**: All

### Error Handling
- No exceptions in this codebase (numerical algorithms)
- Return status via `std::pair<T, bool>` or output parameters
- Document convergence via boolean return value

### Testing Patterns
- Use **doctest** for unit tests
- Use **RapidCheck** for property-based tests (macro: `RAPIDCHECK`)
- Test files: `test/source/test_*.cpp`
- Test naming: `TEST_CASE("description")` or `TEST_SUITE("suite_name")`

---

## Project Structure

```
ginger-cpp/
├── include/ginger/       # Public headers
│   ├── aberth.hpp
│   ├── rootfinding.hpp
│   ├── matrix2.hpp
│   └── vector2.hpp
├── source/               # Implementation (.cpp)
├── test/source/          # Test suite (doctest + RapidCheck)
├── standalone/source/    # Example executable
├── test_installed/       # Consumer install test (CI only)
├── CMakeLists.txt       # Single-root build configuration
└── .clang-format        # Code formatter config
```

---

## Dependencies (via CPM.cmake)

- **fmt**@12.1.0 - Formatting (system-installed first)
- **spdlog**@v1.17.0 - Logging (system-installed first, `SPDLOG_FMT_EXTERNAL` always set)
- **doctest**@2.5.2 - Testing framework
- **rapidcheck** (master) - Property-based testing
- **cxxopts**@3.2.1 - CLI parsing (standalone only)
- **nanobench**@4.3.11 - Benchmarking (header-only download)

---

## CI/CD (GitHub Actions)

- **ubuntu.yml**: Build + test + coverage on Ubuntu
- **windows.yml**: Build + test on Windows
- **macos.yml**: Build + test on macOS
- **install.yml**: Build/install library, then test consumer via `find_package(Ginger)`
- **benchmark.yml**: Performance benchmarking
- **documentation.yaml**: Build + publish Doxygen docs on release

All CI runs with `CTEST_OUTPUT_ON_FAILURE=1`.

---

## Additional Tools

### Static Analyzers

```bash
cmake -B build -DGINGER_ENABLE_CLANG_TIDY=ON
cmake --build build --target clang-tidy
```

### Coverage

```bash
cmake -B build -DGINGER_ENABLE_COVERAGE=1
cmake --build build
cmake --build build --target coverage
```