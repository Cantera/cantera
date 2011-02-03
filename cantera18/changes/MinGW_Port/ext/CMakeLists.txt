SET (CMAKE_C_FLAGS  -DSkip_f2c_Undefs)

INCLUDE_DIRECTORIES (${PROJECT_SOURCE_DIR}/ext/f2c_libs )

add_subdirectory ( f2c_libs )
add_subdirectory ( f2c_blas )
add_subdirectory ( f2c_lapack )
add_subdirectory ( f2c_math )
add_subdirectory ( cvode )
add_subdirectory ( tpx )


