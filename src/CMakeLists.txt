### DO NOT DELETE THIS COMMENT: INSERT_ARCHETYPES_HERE ###
USE_CYCLUS("cybam" "bu_solver_mlp")

INSTALL_CYCLUS_MODULE("cybam" "./")

# install header files
FILE(GLOB h_files "${CMAKE_CURRENT_SOURCE_DIR}/*.h")
INSTALL(FILES ${h_files} DESTINATION include/cybam COMPONENT cybam)
