set(mt_sources "mt19937-64.c")
set(mt_headers 
    "mt64.h"
    "mt_wrapper.h"
)

add_library(mt ${mt_sources} ${mt_headers})

target_include_directories(mt PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
