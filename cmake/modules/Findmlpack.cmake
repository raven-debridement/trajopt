# - Try to find mlpack 
# Once done this will define
#  MLPACK_FOUND - System has mlpack 
#  MLPACK_INCLUDE_DIRS - The mlpack include directories
#  MLPACK_LIBRARIES - The libraries needed to use mlpack


find_path(MLPACK_INCLUDE_DIRS mlpack/core.hpp paths /usr/include /usr/local/include)

find_library(MLPACK_LIBRARIES mlpack paths /usr/lib /usr/local/lib)

if (MLPACK_INCLUDE_DIRS)
  if (MLPACK_LIBRARIES)
    set(MLPACK_FOUND "YES")
  endif()
endif()
