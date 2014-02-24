function a = set(a,varargin)
% SET -  Set properties.
%
%   The properties that may be set are
%
%   Temperature    (T)
%   Density        (Rho)
%   Volume         (V)
%   Pressure       (P)
%   Enthalpy       (H)
%   Entropy        (S)
%   MoleFractions  (X)
%   MassFractions  (Y)
%   Vapor Fraction (Vapor)
%   Liquid Fractio (Liquid)
%
%   Either the full property name or the symbol may be
%   specified. For the extensive properties (V,H,U,S), the values
%   must be given per unit mass. H, U, and S must be set in
%   conjunction with pressure (for H,S) or volume (for U,S). Either
%   (specific) volume or density may be specified. Mole and mass
%   fractions must be input as vectors (either row or column) with
%   length equal to the number of species.
%
%   Examples:
%
%      set(gas,'Temperature',600.0);
%      set(gas,'T',600.0);
%      set(gas,'T',600.0,'P',2*oneatm,'Y',massfracs);
%      set(gas,'H',0.5*enthalpy_mass(gas),'P',pressure(gas));
%      set(gas,'S',entropy_mass(gas),'P',0.5*pressure(gas));
%      set(gas,'X',ones(nSpecies(gas),1));
%      set(gas,'T',500.0,'Vapor',0.8)
%
%  Alternatively, individual methods to set properties may be
%  called (setTemperature, setMoleFractions, etc.)
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
            setMoleFractions(a,val);
        case 'X'
            nx = nx + 1;
            setMoleFractions(a,val);
        case 'MassFractions'
            ny = ny + 1;
            setMassFractions(a,val);
        case 'Y'
            ny = ny + 1;
            setMassFractions(a,val);
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
            error(['unknown property ' char(prop)])
    end
end

if nx + ny > 1
    error('composition specified multiple times');
end

ntot = nt + np + nv + ns + nh + nu + nq;

if ntot == 1
    %
    % set T, v, or P individually
    %
    if nt == 1
        setTemperature(a,tval);   % density held fixed
    elseif nv == 1
        setDensity(a,1.0/vval);   % temperature held fixed
    elseif np == 1
        setPressure(a, pval);     % temperature held fixed
    else
        error('pressure, volume, or density must also be specified');
    end
elseif ntot == 2
    %
    % set property pairs
    %
    if nt == 1 && nv == 1
        setTemperature(a,tval);
        setDensity(a,1.0/vval);
    elseif nt == 1 && np == 1
        setTemperature(a,tval);
        setPressure(a, pval);
    elseif nt == 1 && nq == 1
        setState_Tsat(a, [tval,qval]);
    elseif np == 1 && nq == 1
        setState_Psat(a, [pval,qval]);
    elseif np == 1 && nh == 1
        setState_HP(a,[hval,pval]);
    elseif nu == 1 && nv == 1
        setState_UV(a,[uval,vval]);
    elseif ns == 1 && np == 1
        setState_SP(a,[sval,pval]);
    elseif ns == 1 && nv == 1
        setState_SV(a,[sval,vval]);
    else
        error('unimplemented property pair');
    end
elseif ntot > 2
    error('too many properties specified');
end
