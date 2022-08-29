import cantera as ct

@ct.extension(name="square-rate")
class SquareRate(ct.ExtensibleRate):
    def after_set_parameters(self, node, units):
        self.A = node["A"]

    def replace_eval(self, T):
        return self.A * T**2
