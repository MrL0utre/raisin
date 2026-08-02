# SPDX-License-Identifier: GPL-3.0-only

cmake_minimum_required(VERSION 3.16)

project(raisin VERSION 1.3 LANGUAGES C)

option(RAISIN_ENABLE_SANITIZERS "Enable AddressSanitizer and UndefinedBehaviorSanitizer" OFF)
option(RAISIN_WARNINGS_AS_ERRORS "Treat compiler warnings as errors" ON)

file(GLOB RAISIN_CORE_SOURCES CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/src/*.c")
list(REMOVE_ITEM RAISIN_CORE_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/raisin.c")

add_library(raisin_core STATIC ${RAISIN_CORE_SOURCES})
target_include_directories(raisin_core PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_compile_features(raisin_core PUBLIC c_std_11)

add_executable(raisin "${CMAKE_CURRENT_SOURCE_DIR}/src/raisin.c")
target_link_libraries(raisin PRIVATE raisin_core)

if(MSVC)
  target_compile_options(raisin_core PRIVATE /W4)
  target_compile_options(raisin PRIVATE /W4)
  if(RAISIN_WARNINGS_AS_ERRORS)
    target_compile_options(raisin_core PRIVATE /WX)
    target_compile_options(raisin PRIVATE /WX)
  endif()
else()
  target_compile_options(raisin_core PRIVATE -Wall -Wextra -Wpedantic)
  target_compile_options(raisin PRIVATE -Wall -Wextra -Wpedantic)
  target_link_libraries(raisin_core PUBLIC m)
  if(RAISIN_WARNINGS_AS_ERRORS)
    target_compile_options(raisin_core PRIVATE -Werror)
    target_compile_options(raisin PRIVATE -Werror)
  endif()
endif()

if(RAISIN_ENABLE_SANITIZERS)
  if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(raisin_core PRIVATE -fsanitize=address,undefined -fno-omit-frame-pointer)
    target_link_options(raisin_core PRIVATE -fsanitize=address,undefined)
    target_compile_options(raisin PRIVATE -fsanitize=address,undefined -fno-omit-frame-pointer)
    target_link_options(raisin PRIVATE -fsanitize=address,undefined)
  else()
    message(FATAL_ERROR "RAISIN_ENABLE_SANITIZERS requires GCC or Clang")
  endif()
endif()

set_target_properties(raisin PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/exe"
)

include(CTest)
if(BUILD_TESTING)
  add_executable(test_invariants "${CMAKE_CURRENT_SOURCE_DIR}/test/test_invariants.c")
  target_link_libraries(test_invariants PRIVATE raisin_core)
  if(RAISIN_ENABLE_SANITIZERS AND CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(test_invariants PRIVATE -fsanitize=address,undefined -fno-omit-frame-pointer)
    target_link_options(test_invariants PRIVATE -fsanitize=address,undefined)
  endif()
  set_target_properties(test_invariants PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/test"
  )

  add_test(
    NAME rbh_invariants
    COMMAND test_invariants "${CMAKE_CURRENT_SOURCE_DIR}/hypergraphs/b14.rzn2"
  )

  add_executable(test_arch "${CMAKE_CURRENT_SOURCE_DIR}/test/test_arch.c")
  target_link_libraries(test_arch PRIVATE raisin_core)
  if(RAISIN_ENABLE_SANITIZERS AND CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(test_arch PRIVATE -fsanitize=address,undefined -fno-omit-frame-pointer)
    target_link_options(test_arch PRIVATE -fsanitize=address,undefined)
  endif()
  set_target_properties(test_arch PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/test"
  )
  add_test(
    NAME architecture_topology
    COMMAND test_arch
            "${CMAKE_CURRENT_SOURCE_DIR}/test/fixtures/arch_path.arch"
            "${CMAKE_CURRENT_SOURCE_DIR}/test/fixtures/arch_disconnected.arch"
            "${CMAKE_CURRENT_SOURCE_DIR}/test/fixtures/arch_duplicate.arch"
            "${CMAKE_CURRENT_SOURCE_DIR}/test/fixtures/arch_resources.arch"
            "${CMAKE_CURRENT_SOURCE_DIR}/test/fixtures/arch_bad_resources.arch"
  )

  add_executable(test_capacity "${CMAKE_CURRENT_SOURCE_DIR}/test/test_capacity.c")
  target_link_libraries(test_capacity PRIVATE raisin_core)
  if(RAISIN_ENABLE_SANITIZERS AND CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(test_capacity PRIVATE -fsanitize=address,undefined -fno-omit-frame-pointer)
    target_link_options(test_capacity PRIVATE -fsanitize=address,undefined)
  endif()
  set_target_properties(test_capacity PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/test"
  )
  add_test(
    NAME communication_capacity
    COMMAND test_capacity
            "${CMAKE_CURRENT_SOURCE_DIR}/test/fixtures/arch_path.arch"
            "${CMAKE_CURRENT_SOURCE_DIR}/test/fixtures/arch_resources.arch"
  )

  add_executable(test_scientific "${CMAKE_CURRENT_SOURCE_DIR}/test/test_scientific.c")
  target_link_libraries(test_scientific PRIVATE raisin_core)
  if(RAISIN_ENABLE_SANITIZERS AND CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(test_scientific PRIVATE -fsanitize=address,undefined -fno-omit-frame-pointer)
    target_link_options(test_scientific PRIVATE -fsanitize=address,undefined)
  endif()
  set_target_properties(test_scientific PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/test"
  )
  add_test(
    NAME scientific_reference
    COMMAND test_scientific
            "${CMAKE_CURRENT_SOURCE_DIR}/test/fixtures/chain.rzn2"
            "${CMAKE_CURRENT_SOURCE_DIR}/test/fixtures/arch_path.arch"
  )

  find_package(Python3 COMPONENTS Interpreter)
  if(Python3_Interpreter_FOUND)
    add_test(
      NAME cli_regression
      COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/test/test_cli.py"
              "$<TARGET_FILE:raisin>" "${CMAKE_CURRENT_SOURCE_DIR}"
    )
    add_test(
      NAME algorithm_regression
      COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/test/test_algorithms.py"
              "$<TARGET_FILE:raisin>" "${CMAKE_CURRENT_SOURCE_DIR}"
              "${CMAKE_CURRENT_SOURCE_DIR}/test/baselines/v1.2-b14-arch0.json"
    )
  endif()
endif()
