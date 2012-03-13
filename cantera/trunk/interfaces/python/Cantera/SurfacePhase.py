
from ThermoPhase import ThermoPhase
from exceptions import CanteraError
from Cantera.num import asarray
import _cantera

class SurfacePhase(ThermoPhase):
    """A class for surface phases."""

    def __init__(self, xml_phase=None, index=-1):
        ThermoPhase.__init__(self, xml_phase=xml_phase, index=index)

    def setSiteDensity(self, n0):
        """Set the site density."""
        _cantera.surf_setsitedensity(self._phase_id, n0)

    def siteDensity(self):
        """Site density [kmol/m2]"""
        return _cantera.surf_sitedensity(self._phase_id)

    def setCoverages(self, theta):
        """Set the surface coverages to the values in array *theta*."""
        nt = len(theta)
        if nt == self.nSpecies():
            _cantera.surf_setcoverages(self._phase_id,
                                     asarray(theta,'d'))
        else:
            raise CanteraError('expected '+`self.nSpecies()`+
                               ' coverage values, but got '+`nt`)

    def coverages(self):
        """Return the array of surface coverages."""
        return _cantera.surf_getcoverages(self._phase_id)

    def setConcentrations(self, conc):
        """Set the surface concentrations to the values in
        array *conc*."""
        _cantera.surf_setconcentrations(self._phase_id, conc)

    def concentrations(self):
        """Return the array of surface concentrations [kmol/m2]."""
        return _cantera.surf_getconcentrations(self._phase_id)


class EdgePhase(SurfacePhase):
    """A one-dimensonal edge."""

    def __init__(self, xml_phase=None, index=-1):
        SurfacePhase.__init__(self, xml_phase=xml_phase, index=index)
