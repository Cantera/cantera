

SET (FRONT_ENDS  Cantera++   CanteraML)

FOREACH (front  ${FRONT_ENDS})
      ADD_SUBDIRECTORY (${front} )
ENDFOREACH (front)
