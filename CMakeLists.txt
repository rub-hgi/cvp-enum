cmake_minimum_required (VERSION 3.0)
project (cvp-enum)

# The version number.
set (VERSION_MAJOR 0)
set (VERSION_MINOR 0)
set (VERSION_JUNIOR 1)

set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)

# Set standard compile and link flags
set (CMAKE_C_FLAGS           "-Wall -Wpedantic -std=gnu11")
set (CMAKE_C_FLAGS_RELEASE   "-O3 -march=native -DNDEBUG")
set (CMAKE_C_FLAGS_DEBUG     "-O0 -g")

set (CMAKE_CXX_FLAGS         "-Wall -Wpedantic -pthread -std=gnu++11")
set (CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -DNDEBUG")
set (CMAKE_CXX_FLAGS_DEBUG   "-O0 -g")

set (CMAKE_LDFLAGS           "-O3 -flto -pthread -lstdc++")

# Configure files
# the cmdline header's version needs to be substituted
configure_file ("${PROJECT_SOURCE_DIR}/src.in/cmdline_decoding.h" "${PROJECT_SOURCE_DIR}/src/bindecoding/cmdline_decoding.h")
configure_file ("${PROJECT_SOURCE_DIR}/src.in/cmdline_decoding.c" "${PROJECT_SOURCE_DIR}/src/bindecoding/cmdline_decoding.c")
configure_file ("${PROJECT_SOURCE_DIR}/src.in/cmdline_lwesampler.h" "${PROJECT_SOURCE_DIR}/src/binlwesampler/cmdline_lwesampler.h")
configure_file ("${PROJECT_SOURCE_DIR}/src.in/cmdline_lwesampler.c" "${PROJECT_SOURCE_DIR}/src/binlwesampler/cmdline_lwesampler.c")
configure_file ("${PROJECT_SOURCE_DIR}/src.in/cmdline_reduction.h" "${PROJECT_SOURCE_DIR}/src/binreduction/cmdline_reduction.h")
configure_file ("${PROJECT_SOURCE_DIR}/src.in/cmdline_reduction.c" "${PROJECT_SOURCE_DIR}/src/binreduction/cmdline_reduction.c")
configure_file ("${PROJECT_SOURCE_DIR}/src.in/cmdline_enumeration.h" "${PROJECT_SOURCE_DIR}/src/binenumeration/cmdline_enumeration.h")
configure_file ("${PROJECT_SOURCE_DIR}/src.in/cmdline_enumeration.c" "${PROJECT_SOURCE_DIR}/src/binenumeration/cmdline_enumeration.c")

LIST(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/build/CMakeModules")
# Look for needed libraries
find_package(GMP REQUIRED)
find_package(NTL REQUIRED)
find_package(FPLLL REQUIRED)

# Helper Library
add_library(toolbox STATIC
	src/libtoolbox/conversions.h
	src/libtoolbox/io.h
	src/libtoolbox/matrix_helper.cpp
	src/libtoolbox/matrix_helper.h
	src/libtoolbox/matrix_operations.cpp
	src/libtoolbox/matrix_operations.h
	src/libtoolbox/signals.h
	src/libtoolbox/vector_templates.h
)
target_include_directories (toolbox PRIVATE
	src/libtoolbox
)

# Executables
add_executable (decoding
	src/bindecoding/main_decoding.cpp
	src/bindecoding/cmdline_decoding.c
	src/bindecoding/cmdline_decoding.h
	src/binlwesampler/lwesampler.cpp
	src/binlwesampler/lwesampler.h
	src/binlwesampler/ziggurat.cpp
	src/binlwesampler/ziggurat.h
	src/binreduction/reduction.cpp
	src/binreduction/reduction.h
	src/binenumeration/enumeration_helper.cpp
	src/binenumeration/enumeration_ntl.cpp
	src/binenumeration/enumeration_babai.cpp
	src/binenumeration/enumeration_lp.cpp
	src/binenumeration/enumeration_ln.cpp
	src/binenumeration/enumeration.h
)
target_link_libraries (decoding
	-pthread
	toolbox
	${NTL_LIBRARIES}
	${GMP_LIBRARIES}
)
target_include_directories (decoding PRIVATE
	src/libtoolbox
	src/bindecoding
	src/binlwesampler
	src/binreduction
	src/binenumeration
)

add_executable (lwesampler
	src/binlwesampler/main_lwesampler.cpp
	src/binlwesampler/cmdline_lwesampler.c
	src/binlwesampler/cmdline_lwesampler.h
	src/binlwesampler/lwesampler.cpp
	src/binlwesampler/lwesampler.h
	src/binlwesampler/ziggurat.cpp
	src/binlwesampler/ziggurat.h
)
target_link_libraries (lwesampler
	toolbox
	${NTL_LIBRARIES}
	${GMP_LIBRARIES}
)
target_include_directories (lwesampler PRIVATE
	src/libtoolbox
	src/bindecoding
	src/binlwesampler
	src/binreduction
	src/binenumeration
)

add_executable (reduction
	src/binreduction/main_reduction.cpp
	src/binreduction/cmdline_reduction.c
	src/binreduction/cmdline_reduction.h
	src/binreduction/reduction.cpp
	src/binreduction/reduction.h
)
target_link_libraries (reduction
	toolbox
	${NTL_LIBRARIES}
	${GMP_LIBRARIES}
)
target_include_directories (reduction PRIVATE
	src/libtoolbox
	src/bindecoding
	src/binlwesampler
	src/binreduction
	src/binenumeration
)

add_executable (enumeration
	src/binenumeration/main_enumeration.cpp
	src/binenumeration/cmdline_enumeration.c
	src/binenumeration/cmdline_enumeration.h
	src/binenumeration/enumeration_helper.cpp
	src/binenumeration/enumeration_ntl.cpp
	src/binenumeration/enumeration_babai.cpp
	src/binenumeration/enumeration_lp.cpp
	src/binenumeration/enumeration_ln.cpp
	src/binenumeration/enumeration.h
)
target_link_libraries (enumeration
	-pthread
	toolbox
	${NTL_LIBRARIES}
	${GMP_LIBRARIES}
)
target_include_directories (enumeration PRIVATE
	src/libtoolbox
	src/bindecoding
	src/binlwesampler
	src/binreduction
	src/binenumeration
)
