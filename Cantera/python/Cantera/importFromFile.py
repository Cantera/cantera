
import _cantera
from Thermo import Thermo
from Kinetics import Kinetics

def importFromFile(t, k, params):
    return _cantera.importFromFile(t.cthermo, k.ckin, params['import_file'],'',1)

