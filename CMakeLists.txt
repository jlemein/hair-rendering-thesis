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

file(GLOB SRC_FILES
    "*.h"
    "*.cpp"
)

if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})

    add_executable(prepare src/scripts/preparescene.cpp
        src/scripts/property.h src/scripts/property.cpp
        src/scripts/propertybag.h src/scripts/propertybag.cpp
        src/scripts/constantproperty.h src/scripts/constantproperty.cpp
        src/scripts/animationpropertybag.h src/scripts/animationpropertybag.cpp
        src/scripts/animatedproperty.h src/scripts/animatedproperty.cpp)
    target_link_libraries(prepare ${Boost_LIBRARIES})

    add_executable(preprocess src/scripts/preprocess.cpp
        src/scripts/make_relative.h
        src/scripts/property.h src/scripts/property.cpp
        src/scripts/propertybag.h src/scripts/propertybag.cpp
        src/scripts/constantproperty.h src/scripts/constantproperty.cpp
        src/scripts/animationpropertybag.h src/scripts/animationpropertybag.cpp
        src/scripts/animatedproperty.h src/scripts/animatedproperty.cpp
        src/scripts/propertybag.h src/scripts/propertybag.cpp
        src/scripts/hairsimplify.h src/scripts/hairsimplify.cpp)
    target_link_libraries(preprocess ${Boost_LIBRARIES})

    add_executable(analyzehair src/scripts/analyze-hair.cpp
        src/scripts/hairstruct.h src/scripts/hairstruct.cpp
        src/scripts/pbrthairparser.h src/scripts/pbrthairparser.cpp)
    target_link_libraries(analyzehair ${Boost_LIBRARIES})

#    add_executable(processhair src/scripts/processhair.cpp
#        src/scripts/hairstruct.h src/scripts/hairstruct.cpp
#        src/scripts/pbrthairparser.h src/scripts/pbrthairparser.cpp)
#    target_link_libraries(processhair ${Boost_LIBRARIES})

endif()

add_library(sayhello STATIC src/scripts/sayhello.cpp)
target_include_directories (sayhello PUBLIC src)

add_executable(simplify src/scripts/simplify.cpp
                        src/scripts/hairsimplify.h src/scripts/hairsimplify.cpp)
target_link_libraries(simplify LINK_PUBLIC sayhello)

#add_executable(renderframe src/scripts/renderframe.cpp)
#target_link_libraries(renderframe LINK_PUBLIC sayhello)


