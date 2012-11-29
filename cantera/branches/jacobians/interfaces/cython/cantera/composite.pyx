class Solution(ThermoPhase, Kinetics, Transport):
    pass

class Interface(InterfacePhase, InterfaceKinetics):
    pass

class DustyGas(ThermoPhase, Kinetics, DustyGasTransport):
    pass
