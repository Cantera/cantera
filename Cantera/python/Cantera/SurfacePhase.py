
from ThermoPhase import ThermoPhase
from exceptions import CanteraError
import Numeric
import _cantera

class SurfacePhase(ThermoPhase):
    """A class for surface phases."""

    def __init__(self, xml_phase=None, index=-1):
        ThermoPhase.__init__(self, xml_phase=xml_phase, index=index)

    def setSiteDensity(self, n0):
        _cantera.surf_setsitedensity(self._phase_id, n0)

    def siteDensity(self):
        return _cantera.surf_sitedensity(self._phase_id)
    
    def setCoverages(self, theta):
        nt = len(theta)
        if nt == self.nSpecies():
            _cantera.surf_setcoverages(self._phase_id,
                                     Numeric.asarray(theta,'d'))
        else:
            raise CanteraError('expected '+`self.nSpecies()`+
                               ' coverage values, but got '+`nt`)

    def coverages(self):
        return _cantera.surf_getcoverages(self._phase_id)

    def setConcentrations(self, theta):
        _cantera.surf_setconcentrations(self._phase_id, theta)

    def concentrations(self):
        return _cantera.surf_getconcentrations(self._phase_id)


    
