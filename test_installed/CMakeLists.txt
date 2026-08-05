cmake_minimum_required(VERSION 3.14)
project(GingerInstalledTest LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Find dependencies first - the installed config references fmt.
find_package(fmt REQUIRED)
find_package(spdlog REQUIRED)
find_package(Ginger REQUIRED)

add_executable(ginger_installed_test main.cpp)
target_link_libraries(ginger_installed_test PRIVATE Ginger::Ginger)

enable_testing()
add_test(NAME ginger_installed_test COMMAND ginger_installed_test)
