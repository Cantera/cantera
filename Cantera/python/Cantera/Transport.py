import _cantera

class Transport:
    """Transport property manager.

    A transport property manager is responsible for computing transport
    properties.
    """
    def __init__(self, xml_phase=None,
                 phase=None, model = "", loglevel=0):
        self._phase = phase
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
        self._models = {}
        self._models[self.model] = self.__tr_id

    def __del__(self):
        try:
            _cantera.tran_delete(self.__tr_id)
        except:
            pass

    def setTransportModel(self, model):
        if self._models.has_key(model):
            self.__tr_id = self._models[model]
        else:
            self.__tr_id = _cantera.Transport(model,
                                                 self._phase._phase_id, 0)
        self.model = model
        
    def desc(self):
        if self.model == 'Multi':
            return 'Multicomponent'
        elif self.model == 'Mix':
            return 'Mixture-averaged'
        else:
            return 'Unknown'
            
    def transport_id(self):
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
    
