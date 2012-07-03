function a = set(a,varargin)
% SET -  Set properties.
%
%   The properties that may be set are
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
%
%  Alternatively, individual methods to set properties may be
%  called (setTemperature, setMoleFractions, etc.)
%

property_argin = varargin;

while length(property_argin) >= 2,
    prop = property_argin{1};
    val = property_argin{2};
    property_argin = property_argin(3:end);
    switch prop
        case 'Temperature'
            setTemperature(a,val);
        case 'T'
            setTemperature(a,val);
        case 'MassFractions'
            setMassFractions(a,val);
        case 'Y'
            setMassFractions(a,val);
        case 'mdot'
            setMdot(a,val);
        case 'MassFlux'
            setMdot(a,val);
        case 'P'
            setPressure(a,val);
        case 'Pressure'
            setPressure(a,val);
        case 'tol'
            sz = size(val);
            if sz == nComponents(a)
                setTolerances(a, val(1,:), val(2,:));
            elseif length(val) == 2
                setTolerances(a, 'default', val(1), val(2));
            else
                error('wrong array size for error tolerances');
            end
        case 'tol-time'
            sz = size(val);
            if sz == nComponents(a)
                setTolerances(a, val(1,:), val(2,:));
            elseif length(val) == 2
                rt = val(1);
                at = val(2);
                setTolerances(a, 'default', rt, at, 'ts');
            else
                error('wrong array size for error tolerances');
            end
        case 'grid'
            setupGrid(a, val);
        case 'bounds'
            setBounds(a, val(1,:), val(2,:));
        case 'X'
            setMoleFractions(a, val);
        case 'MoleFractions'
            setMoleFractions(a, val);
        case 'T_fixed'
            setFixedTempProfile(a, val);
        case 'ID'
            setID(a, val);
        otherwise
            error(['unknown property ' char(prop)]);
    end
end

