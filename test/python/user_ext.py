import cantera as ct

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
