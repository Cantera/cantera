import _cantera
from Numeric import asarray
import exceptions

class Transport:

    """Transport property manager.
    A transport property manager is responsible for computing transport
    properties.
    """
    def __init__(self, xml_phase=None,
                 phase=None, model = "", loglevel=0):
        """Create a transport property manager.

        xml_phase --- XML phase element
        phase     --- ThermoPhase instance representing the phase that the
                      transport properties are for
        model     --- string specifying transport model. If omitted,
                      model will be taken from the input file.
        loglevel  --- controls amount of diagnostic output
        """
        if model == "" or model == "Default":
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
        self._models = {}
        self._models[self.model] = self.__tr_id

    def __del__(self):
        for m in self._models.keys():
            try:
                _cantera.tran_delete(self._models[m])
            except:
                pass

    def addTransportModel(self, model, loglevel=1):
        new_id = _cantera.Transport(model,
                                    self._phase_id, loglevel)
        self._models[model] = new_id
        

    def switchTransportModel(self, model):
        if self._models.has_key(model):
            self.__tr_id = self._models[model]
        else:
            raise CanteraError("Transport model "+model+" not defined. Use "
                               +"method addTransportModel first.")
            
    def desc(self):
        if self.model == 'Multi':
            return 'Multicomponent'
        elif self.model == 'Mix':
            return 'Mixture-averaged'
        else:
            return self.model
            
    def transport_id(self):
        return self.__tr_id

    def transport_hndl(self):
        return self.__tr_id
    
    def viscosity(self):
        return _cantera.tran_viscosity(self.__tr_id)

    def thermalConductivity(self):
        return _cantera.tran_thermalConductivity(self.__tr_id)

    def thermalDiffCoeffs(self):
        return _cantera.tran_thermalDiffCoeffs(self.__tr_id,
                                             self.trnsp)

    def binaryDiffCoeffs(self):
        return _cantera.tran_binaryDiffCoeffs(self.__tr_id,
                                             self.trnsp)

    def mixDiffCoeffs(self):
        return _cantera.tran_mixDiffCoeffs(self.__tr_id,
                                         self.trnsp)            

    def multiDiffCoeffs(self):
        return _cantera.tran_multiDiffCoeffs(self.__tr_id,
                                           self.trnsp)

    def setParameters(self, type, k, params):
        return _cantera.tran_setParameters(self.__tr_id,
                                           type, k, asarray(params))
    
    
