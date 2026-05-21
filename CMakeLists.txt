cmake_minimum_required(VERSION 3.30.0 FATAL_ERROR)

project(rk_benchmark VERSION 0.1.0 LANGUAGES C Fortran)

option(ENABLE_C_SOLVER "Use the external C rk23_cptr replacement in the benchmark driver" OFF)
option(ENABLE_LTO "Enable IPO/LTO" OFF)

add_library(rk_lib rk_kinds.f90 rk_types.f90 rk_solvers.f90 rk23_cptr_external.c)

include(CheckLibraryExists)
check_library_exists(m nextafter "" HAVE_LIBM_NEXTAFTER)
if(HAVE_LIBM_NEXTAFTER)
  target_link_libraries(rk_lib PRIVATE m)
endif()

add_library(robertson_lib robertson_models.f90)
target_link_libraries(robertson_lib PRIVATE rk_lib)

add_executable(rk_benchmark rk_benchmark.F90)
target_link_libraries(rk_benchmark PRIVATE rk_lib robertson_lib)

add_executable(euler_benchmark rk_forward_euler.f90 euler_benchmark.f90)
target_link_libraries(euler_benchmark PRIVATE rk_lib robertson_lib)

if(ENABLE_C_SOLVER)
  target_compile_definitions(rk_benchmark PRIVATE USE_EXTERNAL_C_RK23)
endif()

#add_library(robertson_models PLUGIN robertson_models.f90)
if(ENABLE_LTO)
  include(CheckIPOSupported)
  check_ipo_supported(RESULT result OUTPUT output)
  if(result)
    set_property(TARGET rk_lib PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
    set_property(TARGET robertson_lib PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
    set_property(TARGET rk_benchmark PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
    set_property(TARGET euler_benchmark PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
  else()
    message(WARNING "IPO is not supported: ${output}")
  endif()
endif()

message(STATUS "Current configuration: ${CMAKE_BUILD_TYPE}")
message(STATUS "ENABLE_C_SOLVER: ${ENABLE_C_SOLVER}")
message(STATUS "ENABLE_LTO: ${ENABLE_LTO}")
message(STATUS "CMAKE_Fortran_FLAGS: ${CMAKE_Fortran_FLAGS}")
