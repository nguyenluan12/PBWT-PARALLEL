cmake_minimum_required(VERSION 3.28)
project(pbwt)

set(CMAKE_CXX_STANDARD 11)

add_executable(randomizer
        randomizer.cpp
)

add_executable(pbwt
        pbwt.cpp pbwt.h
)

add_executable(pbwt_1
        pbwt_1.cpp
)

find_package(OpenMP REQUIRED)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

set(Boost_USE_STATIC_RUNTIME    OFF)
set(Boost_REALPATH    OFF)
find_package(Boost REQUIRED COMPONENTS iostreams)
if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
    target_link_libraries(pbwt_1 ${Boost_LIBRARIES})
    target_link_libraries(pbwt ${Boost_LIBRARIES})
endif()