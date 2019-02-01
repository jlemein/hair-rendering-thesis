cmake_minimum_required(VERSION 3.3)

project(hair-rendering-thesis LANGUAGES CXX)
set (CMAKE_CXX_STANDARD 11) # use c++11 standard

# Setting source include directory
include_directories(src)

# Output appears in bin directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

# Boost libraries
set(Boost_USE_STATIC_LIBS ON) 
set(Boost_USE_MULTITHREADED ON)  
set(Boost_USE_STATIC_RUNTIME OFF) 

#OpenVDB libraries
find_package(Boost 1.58.0 COMPONENTS program_options filesystem date_time REQUIRED)
find_package(OpenVDB REQUIRED)
include_directories(${OpenVDB_INCLUDE_DIRS})
set(LIBS ${LIBS} ${OpenVDB_LIBRARIES})

#set(GLFW_BUILD_DOCS OFF CACHE BOOL "" FORCE)
#set(GLFW_BUILD_TESTS OFF CACHE BOOL "" FORCE)
#set(GLFW_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
#add_subdirectory("../ext/glfw-3.2.1")

find_package(GLEW REQUIRED)
include_directories(${GLEW_INCLUDE_DIRS})

find_package(glfw3 3.2 REQUIRED)

find_package(OpenGL REQUIRED)
include_directories(${OPENGL_INCLUDE_DIR})

find_package(glm REQUIRED)

##target_sources (gbuffer) PUBLIC plyloader.cpp ShaderSource.cpp SimpleGlUtil.cpp)
#target_include_directories(gbuffer ${OPENGL_INCLUDE_DIR})


#add_library (glfw STATIC IMPORTED)



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
        src/scripts/bezier.h src/scripts/bezier.cpp
        src/scripts/hairstruct.h src/scripts/hairstruct.cpp
        src/scripts/pbrthairparser.h src/scripts/pbrthairparser.cpp)
    target_link_libraries(analyzehair ${Boost_LIBRARIES})

    add_executable(hairtransform src/scripts/hairtransform.cpp
        src/scripts/bezier.h src/scripts/bezier.cpp
        src/scripts/hairstruct.h src/scripts/hairstruct.cpp
        src/scripts/pbrthairparser.h src/scripts/pbrthairparser.cpp)
    target_link_libraries(hairtransform ${Boost_LIBRARIES})

    add_executable(testvdb src/scripts/vdb/testvdb.cpp
        src/scripts/vdb/OpenVdbReader.h src/scripts/vdb/OpenVdbReader.cpp
        src/scripts/vdb/InputOutputUtil.h src/scripts/vdb/InputOutputUtil.cpp
        src/scripts/bezier.h src/scripts/bezier.cpp
        src/scripts/hairstruct.h src/scripts/hairstruct.cpp
        src/scripts/pbrthairparser.h src/scripts/pbrthairparser.cpp)
    target_link_libraries(testvdb ${OpenVDB_LIBRARIES} -ltbb -lHalf)

    add_executable(writevdb src/scripts/vdb/writevdb.cpp
        src/scripts/vdb/OpenVdbReader.h src/scripts/vdb/OpenVdbReader.cpp
        src/scripts/vdb/InputOutputUtil.h src/scripts/vdb/InputOutputUtil.cpp
        src/scripts/bezier.h src/scripts/bezier.cpp
        src/scripts/hairstruct.h src/scripts/hairstruct.cpp
        src/scripts/pbrthairparser.h src/scripts/pbrthairparser.cpp)
    target_link_libraries(writevdb ${OpenVDB_LIBRARIES} -ltbb -lHalf)

    add_executable(readvdb src/scripts/vdb/readvdb.cpp
        src/scripts/vdb/OpenVdbReader.h src/scripts/vdb/OpenVdbReader.cpp
        src/scripts/vdb/InputOutputUtil.h src/scripts/vdb/InputOutputUtil.cpp
        src/scripts/bezier.h src/scripts/bezier.cpp
        src/scripts/hairstruct.h src/scripts/hairstruct.cpp
        src/scripts/pbrthairparser.h src/scripts/pbrthairparser.cpp)
    target_link_libraries(readvdb ${OpenVDB_LIBRARIES} -ltbb -lHalf)

    add_executable(bezier src/scripts/bezier/main.cpp
        src/scripts/bezier.h src/scripts/bezier.cpp
        src/scripts/bezier/CurveRenderer.h src/scripts/bezier/CurveRenderer.cpp
        src/scripts/hairstruct.h src/scripts/hairstruct.cpp
        src/scripts/pbrthairparser.h src/scripts/pbrthairparser.cpp
        src/scripts/gl_util/ShaderSource.h src/scripts/gl_util/ShaderSource.cpp
        #src/scripts/gl_util/Renderable.h src/scripts/gl_util/Renderable.cpp
        src/scripts/gl_util/SimpleGlUtil.h src/scripts/gl_util/SimpleGlUtil.cpp)
    target_link_libraries(bezier ${OPENGL_gl_LIBRARY} glfw ${GLFW_LIBRARIES} GLEW glm)

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


