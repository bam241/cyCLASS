# ------------------- Add all Concrete Agents ----------------------------

EXECUTE_PROCESS(COMMAND git describe --tags OUTPUT_VARIABLE cycbam_version OUTPUT_STRIP_TRAILING_WHITESPACE)
CONFIGURE_FILE(cycbam_version.h.in "${CMAKE_CURRENT_SOURCE_DIR}/cycbam_version.h" @ONLY)

SET(CYCLUS_CUSTOM_HEADERS "cycbam_version.h")

USE_CYCLUS("cycbam" "reactor")

USE_CYCLUS("cycbam" "fuel_fab")

INSTALL_CYCLUS_MODULE("cycbam" "" "NONE")

SET(TestSource ${cycbam_TEST_CC} PARENT_SCOPE)

# install header files
FILE(GLOB h_files "${CMAKE_CURRENT_SOURCE_DIR}/*.h")
INSTALL(FILES ${h_files} DESTINATION include/cycbam COMPONENT cycbam)
