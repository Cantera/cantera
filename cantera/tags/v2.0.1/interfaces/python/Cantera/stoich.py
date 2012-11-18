from Cantera import exceptions
from Cantera.num import array
from Cantera.elements import elementMoles

def det3(A):
    """Determinant of a 3x3 matrix."""
    return (A[0,0]*(A[1,1]*A[2,2] - A[1,2]*A[2,1])
            - A[0,1]*(A[1,0]*A[2,2] - A[1,2]*A[2,0])
            + A[0,2]*(A[1,0]*A[2,1] - A[2,0]*A[1,1]))

def stoich_fuel_to_oxidizer(mix, fuel, oxidizer):
    """Fuel to oxidizer ratio for stoichiometric combustion.

    This function only works for fuels composed of carbon, hydrogen,
    and/or oxygen. The fuel to oxidizer ratio is returned that results in

    """

    # fuel
    mix.setMoleFractions(fuel)
    f_carbon = elementMoles(mix, 'C')
    f_oxygen = elementMoles(mix, 'O')
    f_hydrogen = elementMoles(mix, 'H')

    #oxidizer
    mix.setMoleFractions(oxidizer)
    o_carbon = elementMoles(mix, 'C')
    o_oxygen = elementMoles(mix, 'O')
    o_hydrogen = elementMoles(mix, 'H')

    B = array([f_carbon, f_hydrogen, f_oxygen],'d')
    A = array([[1.0, 0.0, -o_carbon],
               [0.0, 2.0, -o_hydrogen],
               [2.0, 1.0, -o_oxygen]], 'd')

    num = array(A,'d')
    num[:,2] = B
    r = det3(num)/det3(A)
    if r <= 0.0:
        raise CanteraError('negative or zero computed stoichiometric fuel/oxidizer ratio!')
    return 1.0/r

if __name__ == "__main__":
    g = GRI30()
    print stoich_fuel_to_oxidizer(g, 'CH4:1', 'O2:1')
