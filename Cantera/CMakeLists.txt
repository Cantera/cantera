add_subdirectory(src)

if (BUILD_F90_INTERFACE)
    add_subdirectory(fortran)
endif (BUILD_F90_INTERFACE) 

if (BUILD_MATLAB_TOOLBOX)
    add_subdirectory(matlab)
endif (BUILD_MATLAB_TOOLBOX) 
