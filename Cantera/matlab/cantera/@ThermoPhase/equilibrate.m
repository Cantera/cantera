function a = equilibrate(a, xy)
% EQUILIBRATE  Set the phase to a state of chemical equilibrium.
%
%   The second argument must be one of the strings 'TP', 'TV',
%   'HP', 'SP', 'SV', 'UV', 'PT', 'VT', 'PH', 'PS', 'VS', or 'VU',
%   and specifies the two thermodynamic properties held fixed at
%   the values in the initial state. Note that if U, H, V, or S is
%   specified, it is the specific value (per unit mass), not the
%   molar value, that is held fixed.
%
if nargin ~= 2
  error('two arguments required')
end

iok = 0;
switch xy
 case 'TP'
  iok = thermo_set(a.tp_id, 50, 104);
 case 'TV'
  iok = thermo_set(a.tp_id, 50, 100);
 case 'HP'
  iok = thermo_set(a.tp_id, 50, 101);
 case 'SP'
  iok = thermo_set(a.tp_id, 50, 102);
 case 'SV'
  iok = thermo_set(a.tp_id, 50, 107); 
 case 'UV'
  iok = thermo_set(a.tp_id, 50, 105);
 case 'PT'
  iok = thermo_set(a.tp_id, 50, 104);
 case 'VT'
  iok = thermo_set(a.tp_id, 50, 100);
 case 'PH'
  iok = thermo_set(a.tp_id, 50, 101);
 case 'PS'
  iok = thermo_set(a.tp_id, 50, 102);
 case 'VS'
  iok = thermo_set(a.tp_id, 50, 107); 
 case 'VU'
  iok = thermo_set(a.tp_id, 50, 105); 
 otherwise
  error('unsupported option')
end
if iok < 0
  e = geterr;
  if e == 0
    e = 'unknown error';
  end
   error(e);
end

