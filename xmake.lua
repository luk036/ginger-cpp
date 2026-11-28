add_rules("mode.debug", "mode.release", "mode.coverage")
add_requires("doctest", "fmt", "nanobench", "spdlog")

set_languages("c++20")
add_includedirs("include")

if is_plat("linux", "macosx") then
    add_syslinks("pthread")
end

if is_mode("coverage") then
    add_cxflags("-ftest-coverage", "-fprofile-arcs", {force = true})
end

if is_plat("linux") then
    set_warnings("all", "error")
    -- Termux/Android ships its libc++ headers under $PREFIX
    local prefix = os.getenv("PREFIX")
    if prefix then
        add_cxflags("-Wno-unused-command-line-argument", {force = true})
        add_sysincludedirs(prefix .. "/include/c++/v1", prefix .. "/include", {public = true})
    end
    if is_mode("release") then
        add_cxflags("-march=native", "-ffast-math", {force = true})
    end
elseif is_plat("macosx") then
    if is_mode("release") then
        add_cxflags("-march=native", "-ffast-math", {force = true})
    end
elseif is_plat("windows") then
    add_cxflags("/EHsc /utf-8 /W4 /WX /wd4702", {force = true})
    add_ldflags("/FORCE:MULTIPLE", {force = true})
    if is_mode("release") then
        add_cxflags("/arch:AVX2", "/fp:fast", {force = true})
    end
end

target("Ginger")
    set_kind("static")
    -- aberth.hpp publicly includes <lds/lds.hpp>, so the Lds include dir must be public
    add_includedirs("../lds-cpp/include", {public = true})
    add_files("source/*.cpp")
    add_packages("fmt", "spdlog")

target("test_ginger")
    add_deps("Ginger")
    add_files("test/source/*.cpp")
    add_packages("doctest", "fmt", "spdlog")
    add_tests("default")

    -- rapidcheck is optionally provided by a CMake build
    local rapidcheck_dir = path.join(os.projectdir(), "build", "_deps", "rapidcheck-src")
    local rapidcheck_lib_dir = path.join(os.projectdir(), "build", "_deps", "rapidcheck-build")
    local rapidcheck_lib = nil

    if is_plat("windows") then
        rapidcheck_lib_dir = path.join(rapidcheck_lib_dir, "Release")
        rapidcheck_lib = path.join(rapidcheck_lib_dir, "rapidcheck.lib")
    else
        rapidcheck_lib = path.join(rapidcheck_lib_dir, "librapidcheck.a")
    end

    if os.isdir(rapidcheck_dir) and os.isfile(rapidcheck_lib) then
        add_includedirs(path.join(rapidcheck_dir, "include"))
        add_linkdirs(rapidcheck_lib_dir)
        add_links("rapidcheck")
        add_defines("RAPIDCHECK_H")
    end

target("test_fir")
    add_deps("Ginger")
    add_files("bench/BM_fir.cpp")
    add_packages("nanobench", "fmt", "spdlog")

target("test_autocorr")
    add_deps("Ginger")
    add_files("bench/BM_autocorr.cpp")
    add_packages("nanobench", "fmt", "spdlog")

target("test_aberth")
    add_deps("Ginger")
    add_files("bench/BM_aberth.cpp")
    add_packages("nanobench", "fmt", "spdlog")

target("test_scaling")
    add_deps("Ginger")
    add_files("bench/BM_scaling.cpp")
    add_packages("fmt", "spdlog")
