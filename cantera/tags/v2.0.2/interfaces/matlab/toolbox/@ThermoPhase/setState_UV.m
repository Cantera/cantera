function a = setState_UV(a,uv)
% SETSTATE_UV    Set the specific internal energy [J/kg] and
% specific volume [m^3/kg].
%
%    setState_UV(a, uv) sets the specific internal energy and
%    specific volume of object a, holding its composition
%    fixed. Argument 'uv' must be a vector of length 2 containing
%    the desired values for the specific internal energy (J/kg) and
%    specific volume (m^3/kg).
%
if uv(2) <= 0.0
    error('the specific volume must be positive');
end
thermo_set(a.tp_id,21,uv);
