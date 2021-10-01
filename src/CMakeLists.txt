project (huygens)

SET(CMAKE_CXX_FLAGS "-Wall -g -O2 -std=c++17")

set (LIBS ${LIBS} portaudio rtmidi) 
include_directories(/usr/local/include)
SET(CMAKE_EXE_LINKER_FLAGS "-L /usr/local/lib")

add_executable(huygens huygens.cpp)
target_link_libraries (huygens dl ${LIBS})


