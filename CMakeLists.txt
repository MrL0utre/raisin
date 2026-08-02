cmake_minimum_required(VERSION 3.16)

project(raisin VERSION 1.2 LANGUAGES C)

option(RAISIN_ENABLE_SANITIZERS "Enable AddressSanitizer and UndefinedBehaviorSanitizer" OFF)

file(GLOB RAISIN_SOURCES CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/src/*.c")

add_executable(raisin ${RAISIN_SOURCES})
target_include_directories(raisin PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_compile_features(raisin PRIVATE c_std_11)

if(MSVC)
  target_compile_options(raisin PRIVATE /W4)
else()
  target_compile_options(raisin PRIVATE -Wall -Wextra -Wpedantic)
  target_link_libraries(raisin PRIVATE m)
endif()

if(RAISIN_ENABLE_SANITIZERS)
  if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(raisin PRIVATE -fsanitize=address,undefined -fno-omit-frame-pointer)
    target_link_options(raisin PRIVATE -fsanitize=address,undefined)
  else()
    message(FATAL_ERROR "RAISIN_ENABLE_SANITIZERS requires GCC or Clang")
  endif()
endif()

set_target_properties(raisin PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/exe"
)
