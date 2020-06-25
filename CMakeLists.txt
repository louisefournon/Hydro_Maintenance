cmake_minimum_required(VERSION 3.8...3.15)
project(HydroUnitMaintenance VERSION 0.2.0 LANGUAGES CXX)

# ----- Settings ------------------------------------------------------------ #
# Since we are using the block factory, objects from statically linked libraries
# may wrongly appear unused, and by default the linker does not include them,
# so we have to force the linking.

if (NOT BUILD_SHARED_LIBS)
    if (CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
        add_link_options("-Wl,-all_load")
    else ()
        add_link_options("-Wl,--whole-archive,--allow-multiple-definition")
    endif ()
endif ()

# ----- Requirements -------------------------------------------------------- #
find_package(UCBlock REQUIRED)
find_package(MILPSolver REQUIRED)

# --------------------------------------------------------------------------- #
add_executable(hydro_maintenance HydroUnitMaintenance.cpp main.cpp)
target_compile_features(hydro_maintenance PUBLIC cxx_std_17)

target_link_libraries(hydro_maintenance PRIVATE
                      SMS++::UCBlock
                      SMS++::MILPSolver)

# ----- Install instructions -------------------------------------------- #
include(GNUInstallDirs)
install(TARGETS hydro_maintenance
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
