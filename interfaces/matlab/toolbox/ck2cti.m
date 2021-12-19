function f = ck2cti(infile, thermo, transport)
% CK2CTI  Convert a CHEMKIN input file to Cantera format.
% f = ck2cti(infile, thermo, transport)
% Examples::
%
%    f = ck2cti('chem.inp')
%    f = ck2cti('chem.inp', 'therm.dat')
%    f = ck2cti('chem.inp', 'therm.dat', 'tran.dat')
%
% These 3 statements all create a Cantera input file 'chem.cti.' In
% the first case, the CK-format file contains all required species
% thermo data, while in the second case some or all thermo data is
% read from file 'therm.dat.' In the third form, the input file
% created will also contain transport property parameters. The
% function return value is a string containing the output file
% name.
%
% :param infile:
%     Chemistry input file in CHEMKIN format. Required.
% :param thermo:
%     Thermodynamic input file in CHEMKIN format. Optional if
%     thermodynamic data is specified in the chemistry input file.
% :param transport:
%     Transport input file in CHEMKIN format. Optional.
% :return:
%     String with CTML output filename.
%

warning(['The CTI input file format is deprecated and will be removed in ' ...
         'Cantera 3.0. Use ''ck2yaml.py'' to convert Chemkin-format input files to ' ...
         'the YAML format. See https://cantera.org/tutorials/ck2yaml-tutorial.html ' ...
         'for more information.']);

if nargin == 0
    error('input file name must be supplied')
elseif nargin == 1
    thermo = '-';
    transport = '-';
elseif nargin == 2
    transport = '-';
end

dotloc = findstr(infile, '.');
if dotloc(end) > 1
    idtag = infile(1:dotloc(end)-1);
    outfile = [idtag '.cti'];
else
    idtag = infile;
    outfile = [infile '.cti'];
end

iok = ctmethods(0, 1, infile, thermo, transport, idtag, 0, 0);

f = outfile;
