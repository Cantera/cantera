function f = ck2ctml(infile, thermo, transport)
% CK2CTML - Convert a Chemkin-compatible reaction mechanism file to
%           CTML.
%
%    Cantera uses an XML-based file format (CTML) for specifying input
%    parameters of any type, including specifying reaction mechanism
%    attributes. This function reads a reaction mechanism file in "CK
%    format" - the file format used by the Chemkin software package -
%    and writes an equivalent description in CTML.
%
%       f = ck2ctml('chem.inp')
%       f = ck2ctml('chem.inp', 'therm.dat')
%       f = ck2ctml('chem.inp', 'therm.dat', 'tran.dat')
%
%    These 3 statements all create a CTML file 'chem.xml.' In the
%    first case, the CK-format file contains all required species
%    thermo data, while in the second case some or all thermo data is
%    read from file 'therm.dat.' In the third form, the CTML file
%    created will also contain transport property parameters. The
%    function return value is a string containing the output file
%    name.
%
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
    outfile = [idtag '.xml'];
else
    idtag = infile;
    outfile = [infile '.xml'];
end

iok = ctmethods(0, 1, infile, thermo, transport, outfile, idtag);
if iok
    error(geterr);
end
f = outfile;
