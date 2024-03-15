function(luisa_md_add_application name)
    cmake_parse_arguments(APP "" "" "SOURCES" ${ARGN})
    add_executable(${name} ${APP_SOURCES})
    target_link_libraries(${name} PRIVATE luisa::compute)
endfunction()

luisa_md_add_application(test SOURCES test.cpp)