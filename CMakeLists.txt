cmake_minimum_required(VERSION 3.30)
project(MDTools)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -O0")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0")

include_directories("${PROJECT_SOURCE_DIR}/include")

#-- Boost library setup
set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)
find_package(Boost COMPONENTS program_options iostreams filesystem system date_time REQUIRED)

#-- GSL library setup
find_package(GSL REQUIRED)
include_directories("${GSL_INCLUDE_DIRS}")


find_package(PkgConfig REQUIRED)
pkg_search_module(FFTW REQUIRED fftw3 IMPORTED_TARGET)
include_directories(PkgConfig::FFTW)
link_libraries(PkgConfig::FFTW)

add_executable(MDTools src/main.cpp
        include/logger.h
        include/parameters.h
        include/io.h
        src/io.cpp
        src/Modules/PhononDOS/mainPhononDOS.cpp
        src/Modules/PhononDOS/mainPhononDOS.h
        src/Modules/DynamicStructureFactor/mainDynamicStructureFactor.cpp
        src/Modules/DynamicStructureFactor/mainDynamicStructureFactor.h
        src/trajectoryReader.cpp
        include/trajectoryReader.h
)

target_link_libraries(MDTools
        Boost::boost
        Boost::program_options
        Boost::iostreams
        Boost::filesystem
        Boost::date_time
        GSL::gsl
)
