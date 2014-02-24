""" Cantera provides a set of 'transport manager' classes that manage
the computation of transport properties. Every object representing a
phase of matter for which transport properties are needed has a
transport manager assigned to it. The transport manager has only one
job: to compute the values of the transport properties of its assigned
phase.

A transport manager may do things not apparent to the user in order to
improve the speed of transport property evaluation. For example, it
may cache intermediate results that depend only on temperature, so
that if it happens to be called again at the same temperature (a
common occurrence) it can skip over computing the stored
temperature-dependent intermediate properties.  This is why we use the
term 'manager' rather than 'calculator.'

In the Cantera kernel, each different transport model is implemented
by a different class derived from the base class Transport. A
highly simplified class structure is used in the Python interface --
there is only one class. """

import _cantera
from Cantera.num import asarray
import exceptions

class Transport:
    """Transport properties.

    This class provides the Python interface to the family of
    transport manager classes in the Cantera C++ kernel. A transport
    manager has one job: to compute transport properties of a phase of
    matter assigned to it. The phase is represented by an object
    belonging to a class derived from ThermoPhase.

    In the C++ kernel, a transport manager implements a single
    transport model, and is an instance of a subclass of the base
    class ``Transport``.  The structure in Python is a little
    different. A single class ``Transport`` represents any kernel-level
    transport manager. In addition, multiple kernel-kevel transport
    managers may be installed in one Python transport manager,
    although only one is active at any one time. This feature allows
    switching between transport models."""

    def __init__(self, xml_phase=None,
                 phase=None, model = "", loglevel=0):
        """Create a transport property manager.

        :param xml_phase:
            XML phase element
        :param phase:
            :class:`.ThermoPhase` instance representing the phase that the
            transport properties are for
        :param model:
            String specifying transport model. If omitted or set to ``Default``,
            the model will be read from the input file.
        :param loglevel:
            controls the amount of diagnostic output
        """

        # if the transport model is not specified, look for attribute
        # 'model' of the XML 'transport' element
        if model == "" or model == "Default" or model == "default":
            try:
                self.model = xml_phase.child('transport')['model']
            except:
                self.model = ""
        else:
            self.model = model

        self.__tr_id = 0
        self.__tr_id = _cantera.Transport(self.model,
                                             phase._phase_id, loglevel)
        self.trnsp = phase.nSpecies()
        self._phase_id = phase._phase_id

        # dictionary holding all installed transport managers
        self._models = {}
        self._models[self.model] = self.__tr_id

    def __del__(self):
        """Delete all installed transport models."""
        if hasattr(self,'_models'):
            for m in self._models.keys():
                try:
                    _cantera.tran_delete(self._models[m])
                except:
                    pass

    def addTransportModel(self, model, loglevel=1):
        """Add a new transport model.  Note that if *model* is the
        name of an already-installed transport model, the new
        transport manager will take the place of the old one, which
        will no longer be accessible. This method does not change the
        active model."""
        new_id = _cantera.Transport(model,
                                    self._phase_id, loglevel)
        self._models[model] = new_id


    def switchTransportModel(self, model):
        """Switch to a different transport model."""
        if self._models.has_key(model):
            self.__tr_id = self._models[model]
            self.model = model
        else:
            raise CanteraError("Transport model "+model+" not defined. Use "
                               +"method addTransportModel first.")

    def desc(self):
        """A short description of the active model."""
        if self.model == 'Multi':
            return 'Multicomponent'
        elif self.model == 'Mix':
            return 'Mixture-averaged'
        else:
            return self.model

    def transport_id(self):
        """For internal use."""
        return self.__tr_id

    def transport_hndl(self):
        """For internal use."""
        return self.__tr_id

    def viscosity(self):
        "Viscosity [Pa-s]."""
        return _cantera.tran_viscosity(self.__tr_id)

    def electricalConductivity(self):
        """electrical conductivity. [S/m]."""
        return _cantera.tran_electricalConductivity(self.__tr_id)

    def thermalConductivity(self):
        """Thermal conductivity. [W/m/K]."""
        return _cantera.tran_thermalConductivity(self.__tr_id)

    def thermalDiffCoeffs(self):
        """Return a one-dimensional array of the species thermal diffusion
        coefficients. Not implemented in all transport models."""
        return _cantera.tran_thermalDiffCoeffs(self.__tr_id,
                                             self.trnsp)

    def binaryDiffCoeffs(self):
        """Two-dimensional array of species binary diffusion coefficients."""
        return _cantera.tran_binaryDiffCoeffs(self.__tr_id,
                                             self.trnsp)

    def diffusionCoeffs(self):
        """Species diffusion coefficients. (m^2/s)."""
        return self.mixDiffCoeffs()


    def mixDiffCoeffs(self):
        """Mixture-averaged diffusion coefficients."""
        return _cantera.tran_mixDiffCoeffs(self.__tr_id,
                                         self.trnsp)

    def multiDiffCoeffs(self):
        """Two-dimensional array of species multicomponent diffusion
        coefficients. Not implemented in all transport managers."""
        return _cantera.tran_multiDiffCoeffs(self.__tr_id,
                                           self.trnsp)

    def setParameters(self, type, k, params):
        """Set model-specific parameters."""
        return _cantera.tran_setParameters(self.__tr_id,
                                           type, k, asarray(params))


    def molarFluxes(self, state1, state2, delta):
        return _cantera.tran_getMolarFluxes(self.__tr_id, self.trnsp,
                                            asarray(state1), asarray(state2),
                                            delta)

    def massFluxes(self, state1, state2, delta):
        return _cantera.tran_getMassFluxes(self.__tr_id, self.trnsp,
                                           asarray(state1), asarray(state2),
                                           delta)
