function s = GRI30(tr)
% GRI30 - Create a Solution instance representing
%         reaction mechanism GRI-Mech 3.0.
%
%     GRI-Mech 3.0 is a widely-used reaction mechanism for natural gas
%     combustion. It contains 53 species composed of the elements H,
%     C, O, N, and/or Ar, and contains 325 reactions, most of which
%     are reversible. GRI-Mech 3.0, like most combustion mechanisms,
%     is designed for use at pressures where the ideal gas law holds.
%
%     Function GRI30 creates the solution according to the
%     specifications in file gri30.xml. The ideal gas equation of
%     state is used. Transport property evaluation is disabled by
%     default. To enable transport properties, supply the name of
%     the transport model to use.
%
%     g1 = GRI30           % no transport properties
%     g2 = GRI30('Mix')    % mixture-averaged transport properties
%     g3 = GRI30('Multi')  % miulticomponent transport properties
%
if nargin == 0
  s = Solution('gri30.xml');
elseif nargin == 1
  s = Solution('gri30.xml',tr); 
else
  error('wrong number of arguments')
end
