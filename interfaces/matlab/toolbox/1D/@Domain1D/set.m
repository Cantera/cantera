function a = set(a,varargin)
% SET  Set properties of a Domain1D.
% a = set(a,varargin)
% The properties that may be set are
%
% * Temperature     (T)
% * Pressure        (P)
% * Mole Fractions  (X)
% * Mass Flux       (mdot)
% * tol
% * tol-time
% * grid
% * bounds
% * T_fixed
% * ID
%
% Either the full property name or the symbol may be
% specified. Mole and mass
% fractions must be input as vectors (either row or column) with
% length equal to the number of species.
%
% Examples::
%
%     >> set(gas,'Temperature',600.0);
%     >> set(gas,'T',600.0);
%     >> set(gas,'T',600.0,'P',2*oneatm,'Y',massfracs);
%     >> set(gas,'X',ones(nSpecies(gas),1));
%
% Alternatively, individual methods to set properties may be
% called (setTemperature, setMoleFractions, etc.)
%
% See also: :mat:func:`setBounds`, :mat:func:`setFixedTempProfile` :mat:func:`setID`,
% :mat:func:`setMdot`, :mat:func:`setMoleFractions`, :mat:func:`setPressure`,
% :mat:func:`setProfile`, :mat:func:`setSteadyTolerances`, :mat:func:`setTemperature`,
% :mat:func:`setTransientTolerances`, :mat:func:`setupGrid`
%
% :param a:
%     Instance of class :mat:func:`Domain1D`
% :param varargin:
%     Comma separated list of ``property, value`` pairs to be set
%

property_argin = varargin;

while length(property_argin) >= 2,
    prop = property_argin{1};
    val = property_argin{2};
    property_argin = property_argin(3:end);
    switch prop
        case 'Temperature'
            setTemperature(a, val);
        case 'T'
            setTemperature(a, val);
        case 'mdot'
            setMdot(a, val);
        case 'MassFlux'
            setMdot(a, val);
        case 'P'
            setPressure(a, val);
        case 'Pressure'
            setPressure(a, val);
        case 'tol'
            sz = size(val);
            if sz == nComponents(a)
                setSteadyTolerances(a, val(1,:), val(2,:));
            elseif length(val) == 2
                setSteadyTolerances(a, 'default', val(1), val(2));
            else
                error('Wrong array size for error tolerances.');
            end
        case 'tol-time'
            sz = size(val);
            if sz == nComponents(a)
                setTransientTolerances(a, val(1,:), val(2,:));
            elseif length(val) == 2
                rt = val(1);
                at = val(2);
                setTransientTolerances(a, 'default', rt, at);
            else
                error('Wrong array size for error tolerances.');
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

