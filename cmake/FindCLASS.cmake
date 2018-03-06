#  CLASS_FOUND - system has the CLASS Core library
#  CLASS_INCLUDE_DIR - the CLASS include directory
#  CLASS_LIBRARIES - The libraries needed to use the CLASS Core Library
#  CLASS_AGENT_TEST_LIBRARIES - A test library for agents
#  CLASS_TEST_LIBRARIES - All testing libraries
#  CLASS_DEFAULT_TEST_DRIVER - The default CLASS unit test driver

# Check if we have an environment variable to CLASS root
IF(DEFINED ENV{CLASS_PATH})
    IF(NOT DEFINED CLASS_ROOT_DIR)
        SET(CLASS_ROOT_DIR "$ENV{CLASS_PATH}")
    ELSE(NOT DEFINED CLASS_ROOT_DIR)
        # Yell if both exist
        MESSAGE(STATUS "\tTwo CLASS_ROOT_DIRs have been found:")
        MESSAGE(STATUS "\t\tThe defined cmake variable CLASS_ROOT_DIR: ${CLASS_ROOT_DIR}")
        MESSAGE(STATUS "\t\tThe environment variable CLASS_ROOT_DIR: $ENV{CLASS_PATH}")
    ENDIF(NOT DEFINED CLASS_ROOT_DIR)
ELSE(DEFINED ENV{CLASS_PATH})
MESSAGE(STATUS "\tTwo CLASS_ROOT_DIRs have not been found:")
MESSAGE(STATUS "\t\tThe  cmake variable CLASS_ROOT_DIR")
MESSAGE(STATUS "\t\tThe environment variable CLASS_ROOT_DIR")
MESSAGE(STATUS "\t\tOne at least need to be defined...")
ENDIF(DEFINED ENV{CLASS_PATH})

# Let the user know if we're using a hint
MESSAGE(STATUS "Using ${CLASS_ROOT_DIR} as CLASS_ROOT_DIR.")

# Set the include dir, this will be the future basis for other
# defined dirs
FIND_PATH(CLASS_INCLUDE_DIR CLASSHeaders.hxx
    HINTS "${CLASS_ROOT_DIR}" "${CLASS_ROOT_DIR}/source/include"
    "${CLASS_ROOT_DIR}/include"
    "${CLASS_ROOT_DIR}/include/CLASS"
    /usr/local/CLASS /opt/local/CLASS
    PATH_SUFFIXES CLASS/include include include/CLASS source/include CLASS/source/include)

FIND_PATH(CLASS_EQ_INCLUDE_DIR EQ_OneParameter.hxx
    HINTS "${CLASS_ROOT_DIR}" "${CLASS_ROOT_DIR}/source/Model/Equivalence"
    "${CLASS_ROOT_DIR}/Equivalence"
    "${CLASS_ROOT_DIR}/Equivalence/CLASS"
    /usr/local/CLASS /opt/local/CLASS
    PATH_SUFFIXES CLASS/Equivalence Equivalence Equivalence/CLASS source/Equivalence CLASS/source/Equivalence)

FIND_PATH(CLASS_XS_INCLUDE_DIR XSM_MLP.hxx
    HINTS "${CLASS_ROOT_DIR}" "${CLASS_ROOT_DIR}/source/Model/XS"
    "${CLASS_ROOT_DIR}/XS"
    "${CLASS_ROOT_DIR}/XS/CLASS"
    /usr/local/CLASS /opt/local/CLASS
    PATH_SUFFIXES CLASS/XS XS XS/CLASS source/XS CLASS/source/XS)

FIND_PATH(CLASS_IM_INCLUDE_DIR IM_RK4.hxx
    HINTS "${CLASS_ROOT_DIR}" "${CLASS_ROOT_DIR}/source/Model/Irradiation"
    "${CLASS_ROOT_DIR}/Irradiation"
    "${CLASS_ROOT_DIR}/Irradiation/CLASS"
    /usr/local/CLASS /opt/local/CLASS
    PATH_SUFFIXES CLASS/Irradiation Irradiation Irradiation/CLASS source/Irradiation CLASS/source/Irradiation)

# Add the root dir to the hints
SET(CLASS_ROOT_DIR "${CLASS_INCLUDE_DIR}/../..")

# Look for the library
FIND_LIBRARY(CLASS_LIBRARY 
    NAMES CLASSpkg
    HINTS "${CLASS_ROOT_DIR}" "${CLASS_ROOT_DIR}/lib" "$ENV{CLASS_PATH}/lib"
    "$ENV{CLASS_lib}"
    /usr/local/CLASS/lib /usr/local/CLASS
    /opt/CLASS /opt/local/CLASS
    PATH_SUFFIXES CLASS/lib lib)

# Copy the results to the output variables.
IF(CLASS_INCLUDE_DIR AND CLASS_LIBRARY)
    SET(CLASS_FOUND 1)
    SET(CLASS_LIBRARIES "${CLASS_LIBRARY}")
    SET(CLASS_INCLUDE_DIRS "${CLASS_INCLUDE_DIR}")
    SET(CLASS_XS_INCLUDE_DIRS "${CLASS_XS_INCLUDE_DIR}")
    SET(CLASS_EQ_INCLUDE_DIRS "${CLASS_EQ_INCLUDE_DIR}")
    SET(CLASS_IM_INCLUDE_DIRS "${CLASS_IM_INCLUDE_DIR}")
ELSE()
    SET(CLASS_FOUND 0)
    SET(CLASS_LIBRARIES)
    SET(CLASS_INCLUDE_DIRS)
ENDIF()

SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CLASS_INCLUDE_DIR}/../../share/CLASS/cmake)

# Report the results.
IF(CLASS_FOUND)
    SET(CLASS_DIR_MESSAGE "Found CLASS Core Headers : "
        ${CLASS_INCLUDE_DIRS} )
    SET(CLASS_LIB_MESSAGE "Found CLASS Core Library : "
        ${CLASS_LIBRARIES} )
    MESSAGE(STATUS ${CLASS_DIR_MESSAGE})
    MESSAGE(STATUS ${CLASS_LIB_MESSAGE})
ELSE(CLASS_FOUND)
    SET(CLASS_DIR_MESSAGE
        "CLASS was not found. Make sure CLASS_LIBRARY and CLASS_INCLUDE_DIR are set.")
    IF(NOT CLASS_FIND_QUIETLY)
        MESSAGE(STATUS "${CLASS_DIR_MESSAGE}")
        MESSAGE(STATUS "CLASS_INCLUDE_DIR was set to : ${CLASS_INCLUDE_DIR}")
        MESSAGE(STATUS "CLASS_LIBRARY was set to : ${CLASS_LIBRARY}")
        IF(CLASS_FIND_REQUIRED)
            MESSAGE(FATAL_ERROR "${CLASS_DIR_MESSAGE}")
        ENDIF(CLASS_FIND_REQUIRED)
    ENDIF(NOT CLASS_FIND_QUIETLY)

ENDIF(CLASS_FOUND)


MARK_AS_ADVANCED(
    CLASS_INCLUDE_DIRS
    CLASS_LIBRARY
    )
