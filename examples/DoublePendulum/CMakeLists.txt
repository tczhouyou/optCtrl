


link_directories(${CMAKE_BINARY_DIR})



add_executable(DoublePendulum main.cpp )


target_link_libraries(DoublePendulum optctrl)
