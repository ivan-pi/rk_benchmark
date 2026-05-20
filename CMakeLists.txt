cmake_minimum_required(VERSION 3.30.0 FATAL_ERROR)

project(rk_benchmark VERSION 0.1.0 LANGUAGES Fortran)

add_library(rk_lib rk_kinds.f90 rk_types.f90 rk_solvers.f90)

add_executable(rk_benchmark rk_benchmark.f90 robertson_models.f90)
target_link_libraries(rk_benchmark PRIVATE rk_lib)

#add_library(robertson_models PLUGIN robertson_models.f90)