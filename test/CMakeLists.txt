cmake_minimum_required(VERSION 3.14...3.22)

project(GingerTests LANGUAGES CXX)

# ---- Options ----

option(ENABLE_TEST_COVERAGE "Enable test coverage" OFF)
option(TEST_INSTALLED_VERSION "Test the version found by find_package" OFF)

# --- Import tools ----

include(../cmake/tools.cmake)

# ---- Dependencies ----

include(../cmake/CPM.cmake)
include(../specific.cmake)

CPMAddPackage("gh:doctest/doctest@2.4.11")
CPMAddPackage("gh:TheLartians/Format.cmake@1.7.3")
CPMAddPackage(
  NAME rapidcheck
  GITHUB_REPOSITORY emil-e/rapidcheck
  GIT_TAG master
  OPTIONS "RC_BUILD_TESTS OFF" "RC_BUILD_EXAMPLES OFF" "RC_INSTALL OFF"
)

if(TEST_INSTALLED_VERSION)
  find_package(Ginger REQUIRED)
else()
  CPMAddPackage(NAME Ginger SOURCE_DIR ${CMAKE_CURRENT_LIST_DIR}/..)
endif()

# ---- Create binary ----

file(GLOB sources CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/source/*.cpp)
add_executable(${PROJECT_NAME} ${sources})
target_link_libraries(${PROJECT_NAME} doctest::doctest Ginger::Ginger ${SPECIFIC_LIBS})
set_target_properties(${PROJECT_NAME} PROPERTIES CXX_STANDARD 17)

# ---- RapidCheck integration ----
if(rapidcheck_ADDED)
  target_compile_definitions(${PROJECT_NAME} PRIVATE RAPIDCHECK_H)
  if(TARGET rapidcheck::rapidcheck)
    target_link_libraries(${PROJECT_NAME} rapidcheck::rapidcheck)
  else()
    target_include_directories(${PROJECT_NAME} PRIVATE ${rapidcheck_SOURCE_DIR}/include)
    if(TARGET rapidcheck)
      target_link_libraries(${PROJECT_NAME} rapidcheck)
    endif()
  endif()
endif()

# enable compiler warnings
if(NOT TEST_INSTALLED_VERSION)
  if(CMAKE_CXX_COMPILER_ID MATCHES "Clang" OR CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    target_compile_options(Ginger PUBLIC -Wall -Wpedantic -Wextra -Werror)
  elseif(MSVC)
    target_compile_options(Ginger PUBLIC /W4 /WX /wd4819 /wd4996)
    target_compile_definitions(${PROJECT_NAME} PUBLIC DOCTEST_CONFIG_USE_STD_HEADERS)
  endif()
endif()

# ---- Add GingerTests ----

enable_testing()

# Note: doctest and similar testing frameworks can automatically configure CMake tests. For other
# testing frameworks add the tests target instead: add_test(NAME ${PROJECT_NAME} COMMAND
# ${PROJECT_NAME})

# Add a single test that runs all tests (both doctest and RapidCheck)
add_test(NAME GingerAllTests COMMAND ${PROJECT_NAME})

# ---- code coverage ----

if(ENABLE_TEST_COVERAGE)
  target_compile_options(Ginger PUBLIC -O0 -g -fprofile-arcs -ftest-coverage)
  target_link_options(Ginger PUBLIC -fprofile-arcs -ftest-coverage)
endif()
