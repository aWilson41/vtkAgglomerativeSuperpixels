
cmake_minimum_required(VERSION 2.8)
 
PROJECT(VTKSuperpixels)
 
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

file(GLOB HEADER_FILES *.h *.hxx)
file(GLOB MXHEADER_FILES Mx/*.h Mx/*hxx)
file(GLOB SOURCE_FILES *.cpp *cxx)
file(GLOB MXSOURCE_FILES  Mx/*.cpp Mx/*cxx)

source_group("Mx" FILES ${MXHEADER_FILES} ${MXSOURCE_FILES})

add_executable(VTKSuperpixels MACOSX_BUNDLE ${HEADER_FILES} ${SOURCE_FILES} ${MXHEADER_FILES} ${MXSOURCE_FILES})
 
target_link_libraries(VTKSuperpixels ${VTK_LIBRARIES})
