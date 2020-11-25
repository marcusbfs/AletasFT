add_executable(AletasFT ${CMAKE_CURRENT_SOURCE_DIR}/Main.cpp)
target_link_libraries(AletasFT PUBLIC AletasTasklib)
target_link_libraries(AletasFT PUBLIC GeratrizesTasklib)
