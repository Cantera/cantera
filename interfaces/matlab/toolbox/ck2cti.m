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
%
%prog = [ctbin,'/ck2cti'];

% set this to zero to turn off mechanism validation
validate = 1;

% set this to one to turn on debugging.  Use only if ck2cti
% fails, and you want to see how the parser is parsing the input file.
debug = 0;

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

iok = ctmethods(0,1, infile, thermo, transport, idtag, debug, validate);

%iok = system([prog,' -i ',infile,' -t ',thermo,' -tr ',transport, ...
%	      ' -id ',idtag,' > ',outfile]);
if iok
    %ierr2 = system([prog,' > log'])
    %if ierr2
    %  error(['Program ck2cti is not found at ',prog,['. Edit file' ...
    %		    [' ctbin.m to point to the Cantera bin directory.']]])
    % else
    error(['Error occurred while running ck2cti. Check file ck2cti.log' ...
        ' for error messages.']);
    %end
end
f = outfile;
