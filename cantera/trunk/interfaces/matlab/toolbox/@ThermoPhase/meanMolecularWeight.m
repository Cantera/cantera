function mmw = meanMolecularWeight(tp)
% MEANMOLECULARWEIGHT - Mean molecular weight.
%
%    This method is a synonym for method meanMolarMass and is
%    provided for backward compatibility.
%
%    See also; meanMolarMass
%

mmw = phase_get(tp.tp_id,4);
