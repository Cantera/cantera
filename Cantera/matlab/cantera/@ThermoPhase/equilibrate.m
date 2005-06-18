function a = equilibrate(a, xy, solver, rtol, maxsteps, loglevel)
% EQUILIBRATE  Set the phase to a state of chemical equilibrium.
%
%   The second argument must be one of the strings 'TP', 'TV',
%   'HP', 'SP', 'SV', 'UV', 'PT', 'VT', 'PH', 'PS', 'VS', or 'VU',
%   and specifies the two thermodynamic properties held fixed at
%   the values in the initial state. Note that if U, H, V, or S is
%   specified, it is the specific value (per unit mass), not the
%   molar value, that is held fixed.
%

% use the ChemEquil solver by default
if nargin < 3
  solver = 0;
end
if nargin < 4
  rtol = 1.0e-9;
end
if nargin < 5
  maxsteps = 1000;
end
if nargin < 6
  loglevel = 0;
end

iok = 0;
switch xy
 case 'TP'
  iok = thermo_set(a.tp_id, 50, 104, solver, rtol, maxsteps, loglevel);
 case 'TV'
  iok = thermo_set(a.tp_id, 50, 100, solver, rtol, maxsteps, loglevel);
 case 'HP'
  iok = thermo_set(a.tp_id, 50, 101, solver, rtol, maxsteps, loglevel);
 case 'SP'
  iok = thermo_set(a.tp_id, 50, 102, solver, rtol, maxsteps, loglevel);
 case 'SV'
  iok = thermo_set(a.tp_id, 50, 107, solver, rtol, maxsteps, loglevel); 
 case 'UV'
  iok = thermo_set(a.tp_id, 50, 105, solver, rtol, maxsteps, loglevel);
 case 'PT'
  iok = thermo_set(a.tp_id, 50, 104, solver, rtol, maxsteps, loglevel);
 case 'VT'
  iok = thermo_set(a.tp_id, 50, 100, solver, rtol, maxsteps, loglevel);
 case 'PH'
  iok = thermo_set(a.tp_id, 50, 101, solver, rtol, maxsteps, loglevel);
 case 'PS'
  iok = thermo_set(a.tp_id, 50, 102, solver, rtol, maxsteps, loglevel);
 case 'VS'
  iok = thermo_set(a.tp_id, 50, 107, solver, rtol, maxsteps, loglevel); 
 case 'VU'
  iok = thermo_set(a.tp_id, 50, 105, solver, rtol, maxsteps, loglevel); 
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

