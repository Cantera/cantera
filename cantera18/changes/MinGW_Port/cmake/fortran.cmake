#### Cantera Fortran configuration file

#if (NOT BUILD_WITH_F2C)

#### Fortran 90
message("Fortran cmake")

if (BUILD_F90_INTERFACE)
    if (F90 STREQUAL "default")
         FIND_LIBRARY(GFORTRAN_LIB  gfortran ${F90_LIB_DIR} /usr/local/lib)
         IF (GFORTRAN_LIB)
             MESSAGE("${GFORTRAN_LIB}")
         ENDIF (GFORTRAN_LIB)
    endif (F90 STREQUAL "default")
endif (BUILD_F90_INTERFACE)
