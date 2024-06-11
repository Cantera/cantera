import cantera as ct
import numpy as np

class SquareRateData(ct.ExtensibleRateData):
    __slots__ = ("Tsquared",)
    use_count = [0]  # used in test for memory leak

    def __init__(self):
        self.use_count[0] += 1

    def update(self, gas):
        self.Tsquared = gas.T**2
        return True

    def __del__(self):
        self.use_count[0] -= 1

@ct.extension(name="square-rate", data=SquareRateData)
class SquareRate(ct.ExtensibleRate):
    __slots__ = ("A",)
    use_count = [0]  # used in test for memory leak

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.use_count[0] += 1

    def set_parameters(self, node, units):
        self.A = node["A"]

    def get_parameters(self, node):
        node["A"] = self.A

    def eval(self, data):
        return self.A * data.Tsquared

    def __del__(self):
        self.use_count[0] -= 1

# A custom interface rate & rate data type

class FooRateData(ct.ExtensibleRateData):
    __slots__ = ("T",)

    def update(self, interface):
        gas = interface.adjacent['gas'] # assume gas phase is adjacent to us
        self.T = gas.T
        return True


@ct.extension(name="interface-foo-rate", data=FooRateData)
class FooRate(ct.ExtensibleRate):
    __slots__ = ("A", "E")

    def set_parameters(self, node, units):
        self.A = node["A"]
        self.E = node["E"]

    def get_parameters(self, node):
        node["A"] = self.A
        node["E"] = self.E

    def eval(self, data):
        return self.A * np.exp(-self.E / data.T)
