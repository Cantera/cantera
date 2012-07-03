function a = set(a,varargin)
% SET Set phase properties and return the updated object
property_argin = varargin;
while length(property_argin) >= 2,
    prop = property_argin{1};
    val = property_argin{2};
    property_argin = property_argin(3:end);
    switch prop
        case 'Temperature'
            phase_set(a.tp_id,1,val);
        case 'Density'
            phase_set(a.tp_id,2,val);
        otherwise
            error('Phase properties: Temperature, Density')
    end
end
