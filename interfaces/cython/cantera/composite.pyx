class Solution(ThermoPhase, Kinetics, Transport):
    pass

class Interface(InterfacePhase, Kinetics):
    pass

class DustyGas(ThermoPhase, Kinetics, DustyGasTransport):
    pass
