function set(tp, varargin)
% SET  Set properties of a phase.
% set(tp,varargin)
% The properties that may be set are
%
% * Temperature     (T)
% * Density         (Rho)
% * Volume          (V)
% * Pressure        (P)
% * Enthalpy        (H)
% * Entropy         (S)
% * Mole Fractions  (X)
% * Mass Fractions  (Y)
% * Vapor Fraction  (Vapor)
% * Liquid Fraction (Liquid)
%
% Either the full property name or the symbol may be
% specified. For the extensive properties (V,H,U,S), the values
% must be given per unit mass. H, U, and S must be set in
% conjunction with pressure (for H,S) or volume (for U,S). Either
% (specific) volume or density may be specified. Mole and mass
% fractions must be input as vectors (either row or column) with
% length equal to the number of species. Two properties may be
% specified in a single call to :mat:func:`set`, plus one of
% mass fractions or mole fractions.
%
% Examples::
%
%    >> set(gas,'Temperature',600.0);
%    >> set(gas,'T',600.0);
%    >> set(gas,'T',600.0,'P',2*oneatm,'Y',massfracs);
%    >> set(gas,'H',0.5*enthalpy_mass(gas),'P',pressure(gas));
%    >> set(gas,'S',entropy_mass(gas),'P',0.5*pressure(gas));
%    >> set(gas,'X',ones(nSpecies(gas),1));
%    >> set(gas,'T',500.0,'Vapor',0.8)
%
% Alternatively, individual methods to set properties may be
% called (setTemperature, setMoleFractions, etc.)
%
% See also: :mat:func:`setDensity`, :mat:func:`setMassFractions`,
% :mat:func:`setMoleFractions`, :mat:func:`setPressure`, :mat:func:`setState_HP`,
% :mat:func:`setState_Psat`, :mat:func:`setState_SP`, :mat:func:`setState_SV`,
% :mat:func:`setState_Tsat`, :mat:func:`setState_UV`, :mat:func:`setTemperature`
%
% :param tp:
%     Instance of class :mat:func:`ThermoPhase` (or another
%     object that derives from ThermoPhase)
% :param varargin:
%     Comma separated list of ``property, value`` pairs to be set
%

property_argin = varargin;
tval = -999;
pval = -999;
hval = -999;
uval = -999;
sval = -999;
vval = -999;
qval = -999;

np = 0;
nt = 0;
nv = 0;
nx = 0;
ny = 0;
ns = 0;
nh = 0;
nu = 0;
nq = 0;

while length(property_argin) >= 2,
    prop = property_argin{1};
    val = property_argin{2};
    if issparse(val)
        val = full(val);
    end
    property_argin = property_argin(3:end);
    switch prop
        case 'Temperature'
            nt = nt + 1;
            tval = val;
        case 'T'
            nt = nt + 1;
            tval = val;
        case 'Density'
            nv = nv + 1;
            vval = 1.0/val;
        case 'Rho'
            nv = nv + 1;
            vval = 1.0/val;
        case 'V'
            nv = nv + 1;
            vval = val;
        case 'MoleFractions'
            nx = nx + 1;
            setMoleFractions(tp, val);
        case 'X'
            nx = nx + 1;
            setMoleFractions(tp, val);
        case 'MassFractions'
            ny = ny + 1;
            setMassFractions(tp, val);
        case 'Y'
            ny = ny + 1;
            setMassFractions(tp, val);
        case 'Pressure'
            pval = val;
            np = np + 1;
        case 'P'
            pval = val;
            np = np + 1;
        case 'Enthalpy'
            hval = val;
            nh = nh + 1;
        case 'H'
            hval = val;
            nh = nh + 1;
        case 'IntEnergy'
            uval = val;
            nu = nu + 1;
        case 'U'
            uval = val;
            nu = nu + 1;
        case 'Entropy'
            sval = val;
            ns = ns + 1;
        case 'S'
            sval = val;
            ns = ns + 1;
        case 'Sat'
            qval = val;
            nq = nq + 1;
        case 'Vapor'
            qval = val;
            nq = nq + 1;
        case 'Liquid'
            qval = 1.0 - val;
            nq = nq + 1;
        otherwise
            error(['Unknown property ' char(prop)])
    end
end

if nx + ny > 1
    error('Composition specified multiple times.');
end

ntot = nt + np + nv + ns + nh + nu + nq;

if ntot == 1
    %
    % set T, v, or P individually
    %
    if nt == 1
        setTemperature(tp, tval);   % density held fixed
    elseif nv == 1
        setDensity(tp, 1.0/vval);   % temperature held fixed
    elseif np == 1
        setPressure(tp, pval);     % temperature held fixed
    else
        error('Pressure, Volume, or Density must also be specified.');
    end
elseif ntot == 2
    %
    % set property pairs
    %
    if nt == 1 && nv == 1
        setTemperature(tp, tval);
        setDensity(tp, 1.0/vval);
    elseif nt == 1 && np == 1
        setTemperature(tp, tval);
        setPressure(tp, pval);
    elseif np == 1 && nv == 1
        setState_RP(tp, [1.0/vval, pval])
    elseif nt == 1 && nq == 1
        setState_Tsat(tp, [tval,qval]);
    elseif np == 1 && nq == 1
        setState_Psat(tp, [pval,qval]);
    elseif np == 1 && nh == 1
        setState_HP(tp, [hval,pval]);
    elseif nu == 1 && nv == 1
        setState_UV(tp, [uval,vval]);
    elseif ns == 1 && np == 1
        setState_SP(tp, [sval,pval]);
    elseif ns == 1 && nv == 1
        setState_SV(tp, [sval,vval]);
    elseif ns == 1 && nt == 1
        setState_ST(tp, [sval,tval]);
    elseif nt == 1 && nv == 1
        setState_TV(tp, [tval,vval]);
    elseif np == 1 && nv == 1
        setState_PV(tp, [pval,vval]);
    elseif nu == 1 && np == 1
        setState_UP(tp, [uval,pval]);
    elseif nv == 1 && nh == 1
        setState_VH(tp, [vval,hval]);
    elseif nt == 1 && nh == 1
        setState_TH(tp, [tval,hval]);
    elseif ns == 1 && nh == 1
        setState_SH(tp, [sval,hval]);
    else
        error('Unimplemented property pair.');
    end
elseif ntot > 2
    error('Too many properties specified.');
end
