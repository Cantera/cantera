from onedim import *

class BurnerFlame(Stack):
    
    def __init__(self, gas = None, type = 'burner'):
        self.type = type
        self.left = Inlet('burner')
        self.right = Outlet('outlet')
        self.flow = AxisymmetricFlow('flow',gas = gas)
        Stack.__init__(self, [self.left, self.flow, self.right])
        self.showSolution()

        
        
