include_directories("${CMAKE_CURRENT_SOURCE_DIR}")

set(LUISA_COMPUTE_ENABLE_GUI OFF CACHE BOOL "" FORCE)

add_subdirectory(compute)
add_subdirectory(apps)
add_subdirectory(mdrun)



