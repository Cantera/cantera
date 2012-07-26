function f = ck2cti(infile, thermo, transport)
% CK2CTI - Convert a Chemkin-compatible reaction mechanism file to
%           Cantera format.
%
%       f = ck2cti('chem.inp')
%       f = ck2cti('chem.inp', 'therm.dat')
%       f = ck2cti('chem.inp', 'therm.dat', 'tran.dat')
%
%    These 3 statements all create a Cantera input file 'chem.cti.' In
%    the first case, the CK-format file contains all required species
%    thermo data, while in the second case some or all thermo data is
%    read from file 'therm.dat.' In the third form, the input file
%    created will also contain transport property parameters. The
%    function return value is a string containing the output file
%    name.

if nargin == 0
    error('input file name must be supplied')
elseif nargin == 1
    thermo = '-';
    transport = '-';
elseif nargin == 2
    transport = '-';
end

dotloc = findstr(infile,'.');
if dotloc(end) > 1
    idtag = infile(1:dotloc(end)-1);
    outfile = [idtag '.cti'];
else
    idtag = infile;
    outfile = [infile '.cti'];
end

iok = ctmethods(0,1, infile, thermo, transport, idtag, 0, 0);

f = outfile;
