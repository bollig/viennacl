# - Find the OpenCL headers and library
#
# Defines the following if found:
#  OPENCL_FOUND        : TRUE if found, FALSE otherwise
#  OPENCL_INCLUDE_DIRS : Include directories for OpenCL
#  OPENCL_LIBRARIES    : The libraries to link against
#
# The user can set the OPENCLROOT environment variable to help finding OpenCL
# if it is installed in a non-standard place.

set(ENV_ATISTREAMSDKROOT $ENV{ATISTREAMSDKROOT})
if(ENV_ATISTREAMSDKROOT)
 set(ENV_OPENCLROOT $ENV{ATISTREAMSDKROOT})
endif(ENV_ATISTREAMSDKROOT)

set(ENV_AMDAPPSDKROOT $ENV{AMDAPPSDKROOT})
if(ENV_AMDAPPSDKROOT)
 set(ENV_OPENCLROOT $ENV{AMDAPPSDKROOT})
endif(ENV_AMDAPPSDKROOT)

set(ENV_OPENCLROOT2 $ENV{OPENCLROOT})
if(ENV_OPENCLROOT2)
 set(ENV_OPENCLROOT $ENV{OPENCLROOT})
endif(ENV_OPENCLROOT2)

set(ENV_OPENCLROOT3 $ENV{OPENCL_ROOT})
if(ENV_OPENCLROOT3)
  MESSAGE(STATUS "SET OPENCLROOT to ENV_OPENCLROOT3")
 set(ENV_OPENCLROOT $ENV{OPENCL_ROOT})
endif(ENV_OPENCLROOT3)


if(ENV_OPENCLROOT)
  find_path(
    OPENCL_INCLUDE_DIR
    NAMES CL/cl.h OpenCL/cl.h
    PATHS ${ENV_OPENCLROOT}/include ${ENV_OPENCLROOT}
    #NO_DEFAULT_PATH  #uncomment this is you wish to surpress the use of default paths for OpenCL
    )

  if (("${CMAKE_SYSTEM_NAME}" MATCHES "Linux") OR (${CMAKE_SYSTEM_NAME} MATCHES "Windows"))
    if(CMAKE_SIZEOF_VOID_P EQUAL 4)
      set(OPENCL_LIB_SEARCH_PATH
          ${OPENCL_LIB_SEARCH_PATH}
          ${ENV_OPENCLROOT}/lib/x86)
    else(CMAKE_SIZEOF_VOID_P EQUAL 4)
      set(OPENCL_LIB_SEARCH_PATH
          ${OPENCL_LIB_SEARCH_PATH}
          ${ENV_OPENCLROOT}/lib/x86_64)
    endif(CMAKE_SIZEOF_VOID_P EQUAL 4)
  endif(("${CMAKE_SYSTEM_NAME}" MATCHES "Linux") OR (${CMAKE_SYSTEM_NAME} MATCHES "Windows"))
  find_library(
    OPENCL_LIBRARY
    NAMES OpenCL
    PATHS ${OPENCL_LIB_SEARCH_PATH}
    #NO_DEFAULT_PATH  #uncomment this is you wish to surpress the use of default paths for OpenCL
    )
else(ENV_OPENCLROOT)
  find_path(
    OPENCL_INCLUDE_DIR
    NAMES CL/cl.h OpenCL/cl.h
    )

  find_library(
    OPENCL_LIBRARY
    NAMES OpenCL
    )
endif(ENV_OPENCLROOT)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  OPENCL
  DEFAULT_MSG
  OPENCL_LIBRARY OPENCL_INCLUDE_DIR
  )

if(OPENCL_FOUND)
  set(OPENCL_INCLUDE_DIRS ${OPENCL_INCLUDE_DIR})
  set(OPENCL_LIBRARIES ${OPENCL_LIBRARY})
else(OPENCL_FOUND)
  set(OPENCL_INCLUDE_DIRS)
  set(OPENCL_LIBRARIES)
endif(OPENCL_FOUND)

mark_as_advanced(
  OPENCL_INCLUDE_DIR
  OPENCL_LIBRARY
  )

