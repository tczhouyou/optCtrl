# - Try to find SFML
# Once done this will define
#  SFML_FOUND - System has LibSFML
#  SFML_INCLUDE_DIRS - The LibSFML include directories
#  SFML_LIBRARIES - The libraries needed to use LibSFML
set(SFML_DIR "${CMAKE_SOURCE_DIR}/SFML")
find_path(SFML_INCLUDE_DIR NAMES config.hh
          HINTS ${SFML_DIR}/include/
	  SFML/include
          )
  message(${SFML_INCLUDE_DIR})

  find_library(SFML_LIBRARY_GRAPHICS NAMES sfml-graphics  
          HINTS ${SFML_DIR}/lib/
	  SFML/lib
          )

  find_library(SFML_LIBRARY_SYSTEM NAMES sfml-system 
          HINTS ${SFML_DIR}/lib/
	  SFML/lib
          )

  find_library(SFML_LIBRARY_WINDOW NAMES sfml-window 
          HINTS ${SFML_DIR}/lib/
	  SFML/lib
          ) 
  
  set(SFML_LIBRARIES ${SFML_LIBRARY_GRAPHICS} ${SFML_LIBRARY_SYSTEM} ${SFML_LIBRARY_WINDOW})
set(SFML_INCLUDE_DIRS ${SFML_INCLUDE_DIR} )

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set SFML_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(SFML  DEFAULT_MSG
                                  SFML_LIBRARY SFML_INCLUDE_DIR)

mark_as_advanced(SFML_INCLUDE_DIR SFML_LIBRARY )
