from cantera import gas_constant
from numpy import zeros, ones
from .utilities import handleError

def spdict(phase, x):
    nm = phase.speciesNames()
    data = {}
    for k in range(len(nm)):
        data[nm[k]] = x[k]
    return data

class Species:
    def __init__(self,g,name):
        self.g = g
        t = g.T
        p = g.P
        x = g.X
        self.name = name
        self.symbol = name
        self.index = g.species_index(name)
        #self.minTemp = g.minTemp(self.index)
        #self.maxTemp = g.maxTemp(self.index)
        self.molecularWeight = g.molecular_weights[self.index]
        self.c = []
        self.e = g.element_names
        self.hf0 = self.enthalpy_RT(298.15)*gas_constant*298.15
        g.TPX = t,p,x
        for n in range(len(self.e)):
            na = g.n_atoms(self.index, n)
            if na > 0:
                self.c.append((self.e[n],na))

    def composition(self):
        return self.c

    def enthalpy_RT(self,t):
        self.g.TP = t, None
        return self.g.partial_molar_enthalpies[self.index] / (gas_constant*t)

    def cp_R(self,t):
        self.g.setTemperature(t)
        return self.g.cp_R()[self.index]

    def entropy_R(self,t):
        self.g.setTemperature(t)
        return self.g.entropies_R()[self.index]

class Mix:
    def __init__(self,g):
        self.g = g
        self._mech = g
        self.nsp = g.n_species
        self._moles = zeros(self.nsp,'d')
        self.wt = g.molecular_weights

    def setMoles(self, m):
        self._moles = m
        self.g.X = self._moles

    def moles(self):
        return self._moles

    def totalMoles(self):
        sum = 0.0
        for k in range(self.nsp):
            sum += self._moles[k]
        return sum

    def totalMass(self):
        sum = 0.0
        for k in range(self.nsp):
            sum += self._moles[k]*self.wt[k]
        return sum

    def moleDict(self):
        d = {}
        nm = self.g.species_names
        for e in range(self.nsp):
            d[nm[e]] = self._moles[e]
        return d

    def setMass(self, m):
        self.setMoles( m/self.wt)

    def mass(self):
        return self.wt*self._moles

    def speciesNames(self):
        return self.g.species_names

    def massDict(self):
        d = {}
        nm = self.g.species_names
        for e in range(self.nsp):
            d[nm[e]] = self._moles[e]*self.wt[e]
        return d

    def set(self, temperature = None, pressure = None,
            density = None, enthalpy = None,
            entropy = None, intEnergy = None, equil = 0):
        total_mass = self.totalMass()

        if temperature and pressure:
            self.g.TP = temperature, pressure
            if equil:
                self.g.equilibrate('TP',solver=0)

        elif temperature and density:
            self.g.TD = temperature, density
            if equil:
                self.g.equilibrate('TV',solver=0)

        elif pressure and enthalpy:
            self.g.HP = enthalpy, pressure
            if equil:
                self.g.equilibrate('HP',solver=0)

        elif pressure and entropy:
            self.g.SP = entropy, pressure
            if equil:
                self.g.equilibrate('SP',solver=0)

        elif density and entropy:
            self.g.SV = entropy, 1.0/density
            if equil:
                self.g.equilibrate('SV',solver=0)

        elif density and intEnergy:
            self.g.UV = intEnergy, 1.0/density
            if equil:
                self.g.equilibrate('UV',solver=0)

#       else:
#               handleError('unsupported property pair', warning=1)


        total_moles = total_mass/self.g.mean_molecular_weight
        self._moles = self.g.X*total_moles
