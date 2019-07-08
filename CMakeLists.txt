cmake_minimum_required(VERSION 3.10)
project(metacache_mpi)

find_package(MPI REQUIRED)
find_package(Threads REQUIRED)

set(CMAKE_CXX_STANDARD 14)

set(THREADS_PREFER_PTHREAD_FLAG ON)

#SET(MEDIUM_TARGETS False CACHE BOOLEAN "Support for up to 4,294,967,295 reference sequences")
#SET(HIGH_TARGETS False CACHE BOOLEAN "Support for more than 4,294,967,295 reference sequences")

option(MEDIUM_TARGETS "Support for up to 4,294,967,295 reference sequences" ON)
option(HIGH_TARGETS "Support for more than 4,294,967,295 reference sequences" OFF)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake/modules")

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

file(GLOB sources  "src/*.h"  "src/*.cpp")

if(${MEDIUM_TARGETS})
    set(CMAKE_CXX_FLAGS "-DMC_TARGET_ID_TYPE=uint32_t -g -Wall -Wno-unused-function -DUSE_MALLOC_WRAPPERS")
    set(CMAKE_CXX_FLAGS_DEBUG "-DMC_TARGET_ID_TYPE=uint32_t -g -O3")
    set(CMAKE_CXX_FLAGS_RELEASE "-DMC_TARGET_ID_TYPE=uint32_t -O3")
elseif(${HIGH_TARGETS})
    set(CMAKE_CXX_FLAGS "-DMC_TARGET_ID_TYPE=uint64_t -g -Wall -Wno-unused-function -DUSE_MALLOC_WRAPPERS")
    set(CMAKE_CXX_FLAGS_DEBUG "-DMC_TARGET_ID_TYPE=uint64_t -g -O3")
    set(CMAKE_CXX_FLAGS_RELEASE "-DMC_TARGET_ID_TYPE=uint64_t -O3")
else(${MEDIUM_TARGETS})
    set(CMAKE_CXX_FLAGS "-g -Wall -Wno-unused-function -DUSE_MALLOC_WRAPPERS")
    set(CMAKE_CXX_FLAGS_DEBUG "-g -O3")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3")
endif(${MEDIUM_TARGETS})

include_directories(${MPI_CXX_INCLUDE_DIRS})

add_executable(metacache_mpi ${sources})

target_link_libraries(metacache_mpi ${MPI_LIBRARIES})
#target_link_libraries(metacache_mpi Threads::Threads)
target_link_libraries(metacache_mpi ${CMAKE_THREAD_LIBS_INIT})

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${CMAKE_SOURCE_DIR}/docs $<TARGET_FILE_DIR:${PROJECT_NAME}>/docs)