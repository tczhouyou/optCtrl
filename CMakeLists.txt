cmake_minimum_required(VERSION 2.6)
if (NOT CMAKE_BUILD_TYPE)
    message(STATUS "No build type selected, default to Debug")
    set(CMAKE_BUILD_TYPE Debug)
endif()

project(OPTCTRL)

SET(OLD_POLICIES  CMP0043 CMP0062 CMP0012)
foreach(OLD_POLICY ${OLD_POLICIES})
#    message(STATUS "Checking ${OLD_POLICY} policy")
    if(POLICY ${OLD_POLICY})
        cmake_policy(SET ${OLD_POLICY} OLD)
#        message(STATUS "Setting ${OLD_POLICY} policy to OLD")
    endif()
endforeach()

set(LIBNAME optctrl)
set(OPTCTRL_LIBS_BUILD ${LIBNAME})



SET(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")

find_package(Boost REQUIRED)
find_package(GSL 1.15 REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(OpenMP QUIET)
find_package(Simox QUIET)
find_package(Armadillo REQUIRED)
find_package(SFML QUIET)

SET(CMAKE_SKIP_BUILD_RPATH  FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

SET(CMAKE_INSTALL_RPATH "${INSTALL_LIB_DIR}")

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH

set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

#set(GCOV_FLAG "--coverage")

add_definitions(-std=gnu++0x -Wno-unused-parameter -Wno-unused-function -Wno-unused-variable ${OpenMP_CXX_FLAGS} ${GCOV_FLAG})

set(OPTCTRL_BIN_DIR "${CMAKE_BINARY_DIR}/bin")
set(CONF_LIB_DIRS "${CMAKE_BINARY_DIR}")
set(CONF_INCLUDE_DIRS "${PROJECT_SOURCE_DIR}/src" "${PROJECT_BINARY_DIR}")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${OPTCTRL_BIN_DIR})

include("${CMAKE_SOURCE_DIR}/cmake/Testing.cmake")

include_directories("${CMAKE_SOURCE_DIR}/src/methods/")
include_directories(${Boost_INCLUDE_DIR})
include_directories(${LWPR_INCLUDE_DIRS})
include_directories(${GSL_INCLUDE_DIRS})
include_directories(SYSTEM ${Eigen3_INCLUDE_DIR})
include_directories(${ARMADILLO_INCLUDE_DIRS})


if(Simox_FOUND)
    setupSimoxExternalLibraries()
    message(STATUS "Found Simox!")
    message("Virtual Robot Libraries: " ${VirtualRobot_LIBRARIES})
    message("Simox Libraries: " ${Simox_LIBRARIES})
    list(APPEND LIBS ${Simox_LIBRARIES})
endif()

if(SFML_FOUND)
    message(STATUS "Found SFML")
    ADD_DEFINITIONS(-DSFML_FOUND)
    INCLUDE_DIRECTORIES(${SFML_INCLUDE_DIR})
    list(APPEND LIBS ${SFML_LIBRARIES})
endif()

# use, i.e. don't skip the full RPATH for the build tree
SET(CMAKE_SKIP_BUILD_RPATH  FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)


# the RPATH to be used when installing, but only if it's not a system directory
LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
IF("${isSystemDir}" STREQUAL "-1")
   SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
ENDIF("${isSystemDir}" STREQUAL "-1")

SET(DATA_DIR "${PROJECT_SOURCE_DIR}/data")
ADD_DEFINITIONS(-DDATA_DIR="${DATA_DIR}")

file(RELATIVE_PATH RELATIVE_BASE_DIR
     "${CMAKE_INSTALL_PREFIX}/share/cmake/${PROJECT_NAME}"
     "${CMAKE_INSTALL_PREFIX}")


configure_file(cmake/OPTCTRLConfig.cmake.in
  "${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config.cmake" @ONLY)
configure_file(cmake/OPTCTRLConfigVersion.cmake.in
  "${PROJECT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake" @ONLY)

set(LIBS ${LIBS} ${Boost_LIBRARIES} ${LWPR_LIBRARIES}  ${GSL_LIBRARIES} ${VirtualRobot_LIBRARIES} gomp blas pthread dl boost_serialization ${ARMADILLO_LIBRARIES})


list(APPEND HEADERS
	src/methods/ilqr.h
	)

list(APPEND SOURCES
	src/methods/ilqr.cpp
  ) 


include_directories(${PROJECT_SOURCE_DIR}/src)
set(CMAKE_CXX_FLAGS "-fPIC")
set(CMAKE_C_FLAGS "-fPIC")

add_definitions(${OPTCTRL_COMPILE_DEFINITIONS})

add_library(${LIBNAME} SHARED ${SOURCES} ${HEADERS})
#SET_TARGET_PROPERTIES(${LIBNAME} PROPERTIES LINK_FLAGS --coverage)
target_link_libraries(${LIBNAME} ${LIBS})

#exporting and installing
export(TARGETS ${OPTCTRL_LIBS_BUILD} FILE "${PROJECT_BINARY_DIR}/OPTCTRLTargets.cmake")

install(TARGETS ${OPTCTRL_LIBS_BUILD} EXPORT OPTCTRLTargets DESTINATION lib)
install(DIRECTORY src/ DESTINATION include FILES_MATCHING PATTERN *.h)
install(DIRECTORY src/ DESTINATION include FILES_MATCHING PATTERN *.hpp)

install(FILES
  "${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
  "${PROJECT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
  "${PROJECT_BINARY_DIR}/${PROJECT_NAME}Targets.cmake"
  DESTINATION "share/cmake/${PROJECT_NAME}"
  COMPONENT cmake)

export(PACKAGE OPTCTRL)
include("${CMAKE_SOURCE_DIR}/cmake/Packaging.cmake")

add_subdirectory(examples)

