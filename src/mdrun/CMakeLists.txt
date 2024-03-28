set(LUISA_MD_SOURCES
        atom.cpp
        force_lj.cpp
        integrate.cpp
        input.cpp
        neighbor.cpp
        thermo.cpp
        setup.cpp
        output.cpp
        comm.cpp
        atom.h
        force.h
        force_lj.h
        integrate.h
        neighbor.h
        thermo.h
        types.h
        variant.h
        lj.h
        comm.h
        )

add_library(luisa-md SHARED ${LUISA_MD_SOURCES})
target_link_libraries(luisa-md PUBLIC luisa::compute) 
set_target_properties(luisa-md PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS ON)

function(luisa_md_add_executable name)
add_executable(${name} ${name}.cpp)
target_link_libraries(${name} PRIVATE luisa-md)
endfunction()

luisa_md_add_executable(lj)