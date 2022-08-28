import cantera as ct

@ct.extension(name="square-rate")
class SquareRate(ct.ExtensibleRate):
    def replace_eval(self, T):
        return T**2
