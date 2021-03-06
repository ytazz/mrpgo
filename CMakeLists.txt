cmake_minimum_required(VERSION 2.8.12) # Default version on Ubuntu 14.04

project(mrpgo)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)

set(USE_CHOLMOD false CACHE BOOL "use cholmod or not")

# let cmake find FindMKL.cmake
# note: links with lp64 version (not ilp64)
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR})

if(WIN32)
  if(NOT CMAKE_DEBUG_POSTFIX)
    set(CMAKE_DEBUG_POSTFIX d)
  endif()
endif()

if(WIN32)
    set(CMAKE_CXX_FLAGS "/openmp /wd4530")
endif()
if(UNIX)
	set(CMAKE_CXX_FLAGS "-std=c++11 -Wall")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-invalid-offsetof -Wno-sign-compare -Wno-delete-non-virtual-dtor -Wno-reorder")
	set(CMAKE_CXX_FLAGS_DEBUG "-g")
	set(CMAKE_CXX_FLAGS_RELEASE "-O2")
endif()

find_package(Boost  REQUIRED)
find_package(MKL    REQUIRED)
find_package(Eigen3 REQUIRED)

if(USE_CHOLMOD)
    find_package(SuiteSparse)
	add_compile_definitions(USE_CHOLMOD)
endif()

include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${Boost_INCLUDE_DIR})
include_directories(${MKL_INCLUDE_DIR})
include_directories(${EIGEN3_INCLUDE_DIR})

set(MRPGO_DATASET_FILES
  "dataset/intel.g2o"
  "dataset/city10000.g2o"
  "dataset/city10000b.g2o"
  "dataset/sphere2500.g2o"
  "dataset/cubicle.g2o"
  "dataset/cubicle_original.g2o"
  "dataset/grid3D.g2o"
  "dataset/m3500.g2o"
  "dataset/parking-garage.g2o"
  "dataset/rim.g2o"
  "dataset/rim_original.g2o"
  "dataset/torus3D.g2o"
  )

add_library(mrpgo
  include/loader.h
  include/node.h
  include/link.h
  include/posegraph.h
  include/variable.h
  include/solver.h
  include/corrector.h
  include/blas.h
  include/csv.h
  include/tokenizer.h
  include/timer.h
  include/mrpgo.h
  src/loader.cpp
  src/node.cpp
  src/link.cpp
  src/posegraph.cpp
  src/variable.cpp
  src/solver.cpp
  src/corrector.cpp
  src/blas.cpp
  src/timer.cpp
 )
 
 add_executable(mrpgo_app
   example/main.cpp
 )

target_link_libraries(mrpgo_app mrpgo)
target_link_libraries(mrpgo
    ${MKL_LIBRARIES}
    ${Boost_LIBRARIES}
)
if(USE_CHOLMOD)
	target_link_libraries(mrpgo SuiteSparse::cholmod)
endif()

install(TARGETS mrpgo mrpgo_app
  EXPORT        mrpgoExport
  RUNTIME       DESTINATION bin
  LIBRARY       DESTINATION lib
  INCLUDES      DESTINATION include/mrpgo
  PUBLIC_HEADER DESTINATION include/mrpgo
  )

install(EXPORT mrpgoExport
  FILE         mrpgoConfig.cmake
  DESTINATION  cmake
  )

install(FILES ${MRPGO_DATASET_FILES} DESTINATION bin/dataset)

## print all variables
get_cmake_property(_variableNames VARIABLES)
foreach(_variableName ${_variableNames})
	message(STATUS "${_variableName}=${${_variableName}}")
endforeach()
