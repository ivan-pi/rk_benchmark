cmake_minimum_required(VERSION 3.30.0 FATAL_ERROR)

project(rk_benchmark VERSION 0.1.0 LANGUAGES Fortran)

if(NOT CMAKE_BUILD_TYPE)
  set(BUILD_TYPE "<not set>")
else()
  set(BUILD_TYPE "${CMAKE_BUILD_TYPE}")
endif()

add_library(rk_lib rk_kinds.f90 rk_types.f90 rk_solvers.f90)

add_executable(rk_benchmark rk_benchmark.f90 robertson_models.f90)
target_link_libraries(rk_benchmark PRIVATE rk_lib)

#add_library(robertson_models PLUGIN robertson_models.f90)

option(ENABLE_LTO "Enable IPO/LTO" OFF)
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

message(STATUS "Build settings:")
message(STATUS "  BUILD_TYPE: ${BUILD_TYPE}")
message(STATUS "  ENABLE_LTO: ${ENABLE_LTO}")
message(STATUS "  CMAKE_Fortran_FLAGS: ${CMAKE_Fortran_FLAGS}")
