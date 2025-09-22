## Project Overview
This is a C++ project that provides polynomial root-finding algorithms, specifically the Aberth-Ehrlich method. It is designed as a template for other projects and follows modern CMake practices. The project is well-structured, with a clear separation of the library code, tests, and a standalone executable. It also includes support for multi-threading to accelerate the root-finding process.

## Building and Running
The project uses CMake for building and provides several targets for different purposes.

### Build and run the standalone target
```bash
cmake -S standalone -B build/standalone
cmake --build build/standalone
./build/standalone/Ginger --help
```

### Build and run test suite
```bash
cmake -S test -B build/test
cmake --build build/test
CTEST_OUTPUT_ON_FAILURE=1 cmake --build build/test --target test
```

### Run clang-format
```bash
cmake -S test -B build/test
cmake --build build/test --target format
cmake --build build/test --target fix-format
```

### Build the documentation
```bash
cmake -S documentation -B build/doc
cmake --build build/doc --target GenerateDocs
open build/doc/doxygen/html/index.html
```

### Build everything at once
```bash
cmake -S all -B build
cmake --build build
```

## Development Conventions
*   **Coding Style:** The project uses `clang-format` to enforce a consistent coding style.
*   **Testing:** The project uses `doctest` for unit testing. Tests are located in the `test` directory.
*   **Dependencies:** Dependencies are managed using `CPM.cmake`.
*   **Documentation:** The project uses Doxygen for generating documentation from the source code.
*   **Continuous Integration:** The project uses GitHub Actions for continuous integration, with workflows for different operating systems.
*   **Alternative Build System:** The project also includes a `xmake.lua` file, which suggests that it can also be built using the `xmake` build system.
