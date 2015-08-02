function s = IdealGasMix(infile, b, c)
% IDEALGASMIX  Create a mixture of ideal gases.
% s = IdealGasMix(infile, b, c)
% Create a :mat:func:`Solution` instance representing an ideal gas mixture. ::
%
%     gas1 = IdealGasMix('ctml_file'[,'transport_model'])
%     gas2 = IdealGasMix('ck_file'[,'thermo_db'[,'tran_db'[,'transport_model']]])
%
% creates an object that represents an ideal gas mixture. The
% species in the mixture, their properties, and the reactions among
% the species, if any, are specified in file 'ctml_file' and 'ck_file'.  Two
% input file formats are supported - CTML and CK
% (CHEMKIN-compatible).  Examples::
%
%     g1a = IdealGasMix('mech.xml')
%     g1b = IdealGasMix('mech.xml', 'Multi')
%     g2 = IdealGasMix('mech2.inp')
%     g3 = IdealGasMix('mech3.inp', 'therm.dat')
%     g4 = IdealGasMix('mech4.inp', 'therm.dat', 'tran.dat')
%
% Objects ``g1a`` and ``g1b`` are created from a CTML file. CTML files
% contain all data required to build the object, and do not require
% any additional database files. Objects ``g2`` - ``g4`` are created
% from CK-format input files. For ``g2``, 'mech2.inp' contains all required
% species thermo data. File 'mech3.inp' is missing some or all
% species thermo data, and requires database file 'therm.dat.'
% Object ``g4`` is created including transport data.
%
% Note that calling :mat:func:`IdealGasMix` with a CK-format input file
% also creates an equivalent CTML file that may be used in future
% calls. If the initial call includes a transport database, then
% the CTML file will contain transport data.
%
% See also: :mat:func:`ck2cti`, :mat:func:`Solution`
%
% :param infile:
%     Input file, either CTI, CTML, or CHEMKIN format
% :param b:
%     If a CTI or CTML file is specified with ``infile``, this can be
%     the transport modeling to be used. If a CHEMKIN format file is
%     specified with ``infile``, this is the filename of the
%     thermodynamic database, if required.
% :param c:
%     If a CHEMKIN format file is specified with ``infile``, this is the
%     filename of the transport database, if required.
% :return:
%     Instance of class :mat:func:`Solution`
%

dotloc = findstr(infile, '.');
if dotloc(end) > 1
    ext = infile(dotloc(end):end);
    if strcmp(ext, '.cti') || strcmp(ext, '.xml')
        if nargin == 1
            s = Solution(infile);
        elseif nargin == 2
            s = Solution(infile, b);
        end
        return
    end
end

if nargin == 1
    b = '-';
    c = '-';
elseif nargin == 2
    c = '-';
end
xml = ck2cti(infile, b, c);
s = Solution(xml);
set(s, 'P', oneatm);
