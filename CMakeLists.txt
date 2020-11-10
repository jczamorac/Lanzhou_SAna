cmake_minimum_required(VERSION 2.6)
project(C11_Lanzhou)


set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O3 -Wall -Wextra")
include_directories(${PROJECT_SOURCE_DIR}/include)

#Find ROOT
find_package(ROOT REQUIRED)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ROOT_CXX_FLAGS}  -pthread -Wno-deprecated-declarations ")
include_directories(${ROOT_INCLUDE_DIR})
link_directories(${ROOT_LIBRARY_DIR})




set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/build")

## User code
#file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
#file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.h Tandar.h)




#root_generate_dictionary(RootLibDict ${PROJECT_SOURCE_DIR}/include/MSimpleRansac.h ${PROJECT_SOURCE_DIR}/include/MTrack.h LINKDEF ${PROJECT_SOURCE_DIR}/include/LinkDef.h)

#add_library(RootLib SHARED ${sources} RootLibDict)


add_executable(C11_Lanzhou   AnaClass.C  AnaClass.h)
#add_executable(AnaHits AnaHits.C )

target_link_libraries(C11_Lanzhou ${ROOT_LIBRARIES} -lSpectrum )

