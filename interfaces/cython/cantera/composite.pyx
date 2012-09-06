class Solution(ThermoPhase, Kinetics, Transport):
    def __init__(self, *args, **kwars):
        pass

class Interface(InterfacePhase, Kinetics):
    pass
