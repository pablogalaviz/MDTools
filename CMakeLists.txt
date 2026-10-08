cmake_minimum_required(VERSION 3.30)

project(MDTools
        VERSION 1.0.0
        DESCRIPTION "Molecular dynamics analysis tools"
        LANGUAGES C CXX)

#-------------------------------------------------------------------------------
# Language / standard settings
#-------------------------------------------------------------------------------
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Only meaningful for single-config generators; harmless otherwise.
if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
            Debug Release RelWithDebInfo MinSizeRel)
endif()

#-------------------------------------------------------------------------------
# Dependencies
#-------------------------------------------------------------------------------
find_package(OpenMP REQUIRED)

set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)

find_package(Boost 1.70 REQUIRED COMPONENTS
        program_options
        iostreams
        filesystem
        date_time
)

find_package(GSL REQUIRED)

find_package(PkgConfig REQUIRED)
pkg_check_modules(FFTW REQUIRED IMPORTED_TARGET fftw3)

#-------------------------------------------------------------------------------
# Target
#-------------------------------------------------------------------------------
add_executable(MDTools
        src/main.cpp
        src/io.cpp
        src/trajectoryReader.cpp
        # xdrfile (C sources)
        src/xdrfile/xdrfile.c
        src/xdrfile/xdrfile_trr.c
        src/xdrfile/xdrfile_xtc.c
        # Modules
        src/Modules/PhononDOS/mainPhononDOS.cpp
        src/Modules/DynamicStructureFactor/mainDynamicStructureFactor.cpp
        src/Modules/AxialDistributionHistogram/mainAxialDistributionHistogram.cpp
        src/Modules/PairDistributionHistogram/mainPairDistributionHistogram.cpp
        src/Modules/RadiusOfGyration/mainRadiusOfGyration.cpp
        src/Modules/RadialDistributionHistogram/mainRadialDistributionHistogram.cpp
)

# Public headers live here (propagated to consumers, not just globally).
target_include_directories(MDTools
        PUBLIC
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
)

#-------------------------------------------------------------------------------
# Link libraries (imported targets propagate includes, flags, and libs)
#-------------------------------------------------------------------------------
target_link_libraries(MDTools
        PRIVATE
        Boost::program_options
        Boost::iostreams
        Boost::filesystem
        Boost::date_time
        GSL::gsl
        PkgConfig::FFTW
        OpenMP::OpenMP_CXX
)

#-------------------------------------------------------------------------------
# Compiler flags (target-scoped, per-configuration where possible)
#-------------------------------------------------------------------------------
# Release optimization (generator-expression form keeps it per-config-safe).
target_compile_options(MDTools PRIVATE
        $<$<AND:$<CONFIG:Release>,$<COMPILE_LANGUAGE:C,CXX>>:-O3>
)

# Debug: no optimization.
target_compile_options(MDTools PRIVATE
        $<$<AND:$<CONFIG:Debug>,$<COMPILE_LANGUAGE:C,CXX>>:-O0>
)

install(TARGETS MDTools DESTINATION $ENV{HOME}/bin CONFIGURATIONS Release)
