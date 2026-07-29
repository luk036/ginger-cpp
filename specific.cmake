# Try system-installed fmt first (Ubuntu: libfmt-dev, macOS: brew install fmt, Termux: fmt)
find_package(fmt CONFIG QUIET)

if(fmt_FOUND)
  message(STATUS "Found system fmt: ${fmt_DIR}")
else()
  CPMAddPackage(
    NAME fmt
    GIT_TAG 12.1.0
    GITHUB_REPOSITORY fmtlib/fmt
    OPTIONS "FMT_INSTALL YES" # create an installable target
  )
endif()

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

# CPMAddPackage( NAME LdsGen GIT_TAG 1.2.1 GITHUB_REPOSITORY luk036/lds-gen-cpp OPTIONS
# "INSTALL_ONLY YES" # create an installable target )

# Try system-installed spdlog first (Ubuntu: libspdlog-dev, macOS: brew install spdlog, Termux: spdlog)
find_package(spdlog CONFIG QUIET)

if(spdlog_FOUND)
  message(STATUS "Found system spdlog: ${spdlog_DIR}")
else()
  # Add spdlog for logging functionality - use external fmt to avoid duplicate symbols
  CPMAddPackage(
    NAME spdlog
    GIT_TAG v1.17.0
    GITHUB_REPOSITORY gabime/spdlog
    OPTIONS "SPDLOG_INSTALL YES" "SPDLOG_FMT_EXTERNAL YES" # create an installable target
  )
endif()

set(SPECIFIC_LIBS Threads::Threads fmt::fmt spdlog::spdlog)

# cpmaddpackage( NAME GSL GITHUB_REPOSITORY "microsoft/GSL" GIT_TAG "v4.0.0" GIT_SHALLOW ON )

# include(FetchContent)
#
# FetchContent_Declare(GSL GIT_REPOSITORY "https://github.com/microsoft/GSL" GIT_TAG "v4.0.0"
# GIT_SHALLOW ON )

# FetchContent_MakeAvailable(GSL)

# find_package(Microsoft.GSL CONFIG REQUIRED)

# set(SPECIFIC_LIBS fmt::fmt Microsoft.GSL::GSL)
