# Try to find the NTL library
#  NTL_FOUND - system has NTL lib
#  NTL_INCLUDE_DIR - the NTL include directory
#  NTL_LIBRARIES - Libraries needed to use NTL

# TODO: check for GMP, as NTL needs this...

if (NTL_INCLUDE_DIR AND NTL_LIBRARIES)
  # Already in cache, be silent
  set (NTL_FIND_QUIETLY TRUE)
endif (NTL_INCLUDE_DIR AND NTL_LIBRARIES)

find_path (NTL_INCLUDE_DIR NAMES NTL/ZZ.h )
find_library (NTL_LIBRARIES NAMES ntl libntl )
MESSAGE (STATUS "NTL libs: " ${NTL_LIBRARIES})

include (FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS (NTL DEFAULT_MSG NTL_INCLUDE_DIR NTL_LIBRARIES)

mark_as_advanced (NTL_INCLUDE_DIR NTL_LIBRARIES)
