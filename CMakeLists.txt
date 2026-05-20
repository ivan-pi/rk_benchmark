cmake_minimum_required(VERSION 3.30.0 FATAL_ERROR)

project(rk_benchmark VERSION 0.1.0 LANGUAGES C Fortran)

option(ENABLE_C_SOLVER "Build the external C rk23_cptr replacement" OFF)
option(ENABLE_LTO "Enable IPO/LTO" OFF)

add_library(rk_lib rk_kinds.f90 rk_types.f90 rk_solvers.f90)

if(ENABLE_C_SOLVER)
  target_sources(rk_lib PRIVATE rk_c_solver.f90 rk23_cptr_external.c)
endif()

add_executable(rk_benchmark rk_benchmark.f90 robertson_models.f90)
target_link_libraries(rk_benchmark PRIVATE rk_lib)
set_source_files_properties(rk_benchmark.f90 PROPERTIES Fortran_PREPROCESS ON)

if(ENABLE_C_SOLVER)
  target_compile_definitions(rk_benchmark PRIVATE USE_EXTERNAL_C_RK23)
endif()

#add_library(robertson_models PLUGIN robertson_models.f90)
if(ENABLE_LTO)
  include(CheckIPOSupported)
  check_ipo_supported(RESULT result OUTPUT output)
  if(result)
    set_property(TARGET rk_lib PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
    set_property(TARGET rk_benchmark PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
  else()
    message(WARNING "IPO is not supported: ${output}")
  endif()
endif()

message(STATUS "Current configuration: ${CMAKE_BUILD_TYPE}")
message(STATUS "ENABLE_C_SOLVER: ${ENABLE_C_SOLVER}")
message(STATUS "ENABLE_LTO: ${ENABLE_LTO}")
message(STATUS "CMAKE_Fortran_FLAGS: ${CMAKE_Fortran_FLAGS}")
