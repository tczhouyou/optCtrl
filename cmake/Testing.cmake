# macros for testing

# Switch for test cases
option(DMP_BUILD_TESTS "Build tests" ON)

if( DMP_BUILD_TESTS )
    find_package(Boost 1.44 REQUIRED COMPONENTS unit_test_framework)
    link_directories(${Boost_LIBRARY_DIRS})

    enable_testing()
    
endif( DMP_BUILD_TESTS )

##
## This macro adds a testcase to the testsuite which is run by the command 'make test'.
## All test executables are stored in ${DMP_BIN_DIR}.
## The Output path of the log files is stored in Test.h which gets generated above.
##
## PARAM TEST_NAME name of the test and the executable which gets created
## PARAM TEST_FILE name of the cpp file containing the test code
## PARAM DEPENDENT_LIBRARIES the libraries which must be linked to the testcase executable
##
macro(DMP_add_test TEST_NAME TEST_FILE DEPENDENT_LIBRARIES)
    if (DMP_BUILD_TESTS)
        message(STATUS "    Building test ${TEST_NAME}")
#        message("${DEPENDENT_LIBRARIES}")
        add_executable(${TEST_NAME} ${TEST_FILE})
        include_directories(${PROJECT_SOURCE_DIR}/src)
        if(GCOV_FLAG)
            SET_TARGET_PROPERTIES(${TEST_NAME} PROPERTIES LINK_FLAGS ${GCOV_FLAG})
        endif()
        if (NOT DMP_OS_WIN)
            target_link_libraries(${TEST_NAME} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY} ${DEPENDENT_LIBRARIES})
        endif()
        message(STATUS "DMP_BIN_DIR: ${DMP_BIN_DIR}")
        add_test(NAME ${TEST_NAME}
                 COMMAND "${DMP_BIN_DIR}/${TEST_NAME}" --output_format=XML --log_level=all --report_level=detailed
                 WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
    endif()
endmacro()
