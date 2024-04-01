add_executable(main main.cpp)

target_link_libraries(main PRIVATE ${PROJECT_NAME}::${PROJECT_NAME})
