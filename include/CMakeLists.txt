add_library(FunDerlib INTERFACE)
target_include_directories(FunDerlib INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})

add_library(quadlib INTERFACE)
target_include_directories(quadlib INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})

add_library(Aletaslib INTERFACE)
target_include_directories(Aletaslib INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(Aletaslib INTERFACE FunDerlib)
target_link_libraries(Aletaslib INTERFACE eigen)

add_library(AletasTasklib INTERFACE)
target_include_directories(AletasTasklib INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(AletasTasklib INTERFACE FunDerlib)
target_link_libraries(AletasTasklib INTERFACE Aletaslib)
target_link_libraries(AletasTasklib INTERFACE eigen)
target_link_libraries(AletasTasklib INTERFACE quadlib)

add_library(GeratrizesTasklib INTERFACE)
target_include_directories(GeratrizesTasklib
                           INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(GeratrizesTasklib INTERFACE FunDerlib)
target_link_libraries(GeratrizesTasklib INTERFACE Aletaslib)
target_link_libraries(GeratrizesTasklib INTERFACE AletasTasklib)
