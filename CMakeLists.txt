cmake_minimum_required(VERSION 3.3)

project(hair-rendering-thesis LANGUAGES CXX)
set (CMAKE_CXX_STANDARD 11) # use c++11 standard

# Setting source include directory
include_directories(src)

# Output appears in bin directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Boost libraries
set(Boost_USE_STATIC_LIBS ON) 
set(Boost_USE_MULTITHREADED ON)  
set(Boost_USE_STATIC_RUNTIME OFF) 

find_package(Boost 1.58.0 COMPONENTS program_options filesystem date_time  REQUIRED)

if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
  add_executable(replace_args src/scripts/replace_args.cpp
                              src/scripts/propertybag.h src/scripts/propertybag.cpp)
  target_link_libraries(replace_args ${Boost_LIBRARIES})

  add_executable(preprocess src/scripts/preprocess.cpp
                            src/scripts/propertybag.h src/scripts/propertybag.cpp
                            src/scripts/hairsimplify.h src/scripts/hairsimplify.cpp)
  target_link_libraries(preprocess ${Boost_LIBRARIES})
endif()

add_library(sayhello STATIC src/scripts/sayhello.cpp)
target_include_directories (sayhello PUBLIC src)

add_executable(simplify src/scripts/simplify.cpp
                        src/scripts/hairsimplify.h src/scripts/hairsimplify.cpp)
target_link_libraries(simplify LINK_PUBLIC sayhello)

#add_executable(renderframe src/scripts/renderframe.cpp)
#target_link_libraries(renderframe LINK_PUBLIC sayhello)


