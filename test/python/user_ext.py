import cantera as ct

class SquareRateData(ct.ExtensibleRateData):
    def replace_update(self, gas):
        self.Tsquared = gas.T**2
        return True

@ct.extension(name="square-rate", data=SquareRateData)
class SquareRate(ct.ExtensibleRate):
    def after_set_parameters(self, node, units):
        self.A = node["A"]

    def replace_eval(self, data):
        return self.A * data.Tsquared
