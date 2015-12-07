### DO NOT DELETE THIS COMMENT: INSERT_ARCHETYPES_HERE ###
CONFIGURE_FILE(cycbam_version.h.in "${CMAKE_CURRENT_SOURCE_DIR}/cybam_version.h" @ONLY)

SET(CYCLUS_CUSTOM_HEADERS "cybam_version.h")

USE_CYCLUS("cybam" "bu_solver_mlp")
USE_CYCLUS("cybam" "fuel_fab")
USE_CYCLUS("cybam" "reactor")

INSTALL_CYCLUS_MODULE("cybam" "./")

# install header files
FILE(GLOB h_files "${CMAKE_CURRENT_SOURCE_DIR}/*.h")
INSTALL(FILES ${h_files} DESTINATION include/cybam COMPONENT cybam)
