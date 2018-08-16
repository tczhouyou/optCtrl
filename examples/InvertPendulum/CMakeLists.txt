


link_directories(${CMAKE_BINARY_DIR})



add_executable(InvertPendulum main.cpp )


target_link_libraries(InvertPendulum optctrl)
