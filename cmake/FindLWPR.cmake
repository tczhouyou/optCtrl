# - Try to find LWPR
# Once done this will define
#  LWPR_FOUND - System has LibLWPR
#  LWPR_INCLUDE_DIRS - The LibLWPR include directories
#  LWPR_LIBRARIES - The libraries needed to use LibLWPR

find_path(LWPR_INCLUDE_DIR NAMES lwpr.hh
          HINTS ${LWPR_DIR}/include/
          lwpr/include
          )

find_library(LWPR_LIBRARY NAMES lwpr
          HINTS ${LWPR_DIR}/lib/
          lwpr/lib
          )

set(LWPR_LIBRARIES ${LWPR_LIBRARY} )
set(LWPR_INCLUDE_DIRS ${LWPR_INCLUDE_DIR} )

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set LWPR_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(LWPR  DEFAULT_MSG
                                  LWPR_LIBRARY LWPR_INCLUDE_DIR)

mark_as_advanced(LWPR_INCLUDE_DIR LWPR_LIBRARY )
