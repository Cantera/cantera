function setThermalConductivity(tr, lam)
% SETTHERMALCONDUCTIVITY - Set the thermal conductivity.
%
%    This method can only be used with transport models that
%    support directly setting the value of the thermal
%    conductivity.
%
setParameters(tr, 1, 0, lam);
