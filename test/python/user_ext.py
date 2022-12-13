import cantera as ct

class SquareRateData(ct.ExtensibleRateData):
    __slots__ = ("Tsquared",)

    def update(self, gas):
        self.Tsquared = gas.T**2
        return True

@ct.extension(name="square-rate", data=SquareRateData)
class SquareRate(ct.ExtensibleRate):
    __slots__ = ("A",)

    def set_parameters(self, node, units):
        self.A = node["A"]

    def eval(self, data):
        return self.A * data.Tsquared
