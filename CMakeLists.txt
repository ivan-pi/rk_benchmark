cmake_minimum_required(VERSION 3.30.0 FATAL_ERROR)

project(rk_benchmark VERSION 0.1.0 LANGUAGES Fortran C)

add_library(rk_lib rk_kinds.f90 rk_types.f90 rk_solvers.f90)

add_executable(rk_benchmark rk_benchmark.f90 robertson_models.f90)
target_link_libraries(rk_benchmark PRIVATE rk_lib)

# ---------------------------------------------------------------------------
# dlopen benchmark support
# ---------------------------------------------------------------------------

# MODULE library: compiled to a shared object that may only be loaded at
# runtime via dlopen; it must not be linked against by other targets.
add_library(robertson_plugin MODULE robertson_plugin.f90)

# C shim providing rkb_dlopen / rkb_dlsym / rkb_dlclose to Fortran.
add_library(dlopen_shim STATIC dlopen_iface.c)

# Generate a small Fortran module that bakes in the absolute path of the
# plugin .so so the benchmark executable can find it without LD_LIBRARY_PATH.
file(GENERATE
  OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/plugin_path.f90"
  CONTENT
"module plugin_path_mod
  implicit none
  character(len=*), parameter :: robertson_plugin_path = &
      '$<TARGET_FILE:robertson_plugin>'
end module plugin_path_mod
")

target_sources(rk_benchmark PRIVATE
  "${CMAKE_CURRENT_BINARY_DIR}/plugin_path.f90")
target_link_libraries(rk_benchmark PRIVATE dlopen_shim ${CMAKE_DL_LIBS})