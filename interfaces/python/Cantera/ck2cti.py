import _cantera

"""
Convert a Chemkin-format input file to CTI format.
Parameters:

infile - name of the Chemkin-format input file.

thermodb - Thermodynamic database. This may be a standard
           Chemkin-format thermo database, or may be any
           Chemkin-format input file containing a THERMO section.

trandb - Transport database. File containing species transport
         parameters in Chemkin format. If this argument is omitted,
         the CTI file will not contain transport property information.

idtag - ID tag. Used to identify the ideal_gas entry in the CTI file. Optional.

debug - If set to 1, extra debugging output will be written. This
        should only be used if ck2cti fails, in order to view
        intermediate output of the parser. Default: off (0).

validate - If set to 1, the mechanism will be checked for errors. This
           is recommended, but for very large mechanisms may slow down
           the conversion process. Default: on (1).

The translated file is written to the standard output.
"""

def ck2cti(infile = "chem.inp", thermodb = "", trandb
           = "", idtag = "", debug = 0, validate = 1):
    _cantera.ct_ck2cti(infile,
                       thermodb, trandb, idtag, debug, validate)
