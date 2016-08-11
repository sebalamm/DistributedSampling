cmake_minimum_required(VERSION 2.8.9)
project(sampling)

set(CMAKE_CXX_COMPILER g++)
 
add_subdirectory(extern)

include_directories(
    app
    lib
    lib/sampling
    lib/tools
)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DNDEBUG -Wall -march=native -O3 -g -std=c++0x -DHAVE_SSE2 -mavx2 -ffast-math")

add_executable(run_methodD "app/run_methodD.cpp")
add_executable(run_methodH "app/run_methodH.cpp")
add_executable(run_methodR "app/run_methodR.cpp")
add_executable(run_methodSH "app/run_methodSH.cpp")
add_executable(run_methodSR "app/run_methodSR.cpp")
# add_executable(run_experiments "app/run_experiments.cpp")

set(LIBS 
    stocc
    dsfmt
    spooky
    city
)
 
target_link_libraries(run_methodD LINK_PUBLIC ${LIBS})
target_link_libraries(run_methodH LINK_PUBLIC ${LIBS})
target_link_libraries(run_methodR LINK_PUBLIC ${LIBS})
target_link_libraries(run_methodSH LINK_PUBLIC ${LIBS})
target_link_libraries(run_methodSR LINK_PUBLIC ${LIBS})
# target_link_libraries(run_experiments LINK_PUBLIC ${LIBS})
