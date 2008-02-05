SET (CMAKE_C_FLAGS  -DSkip_f2c_Undefs)
add_subdirectory( f2c_libs) 
add_subdirectory( f2c_blas)
add_subdirectory (f2c_lapack)
add_subdirectory ( tpx )


