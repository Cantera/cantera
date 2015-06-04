function setFixedTempProfile(d, profile)
% SETFIXEDTEMPPROFILE  Set a fixed temperature profile.
% d = setFixedTempProfile(d, profile)
% Set the temperature profile to use when the
% energy equation is not being solved. The profile must be entered
% as an array of positions / temperatures, which may be in rows or
% columns.
%
% :param d:
%     Instance of class :mat:func:`Domain1D`
% :param profile:
%     n x 2 or 2 x n array of ``n`` points at which the temperature
%     is specified.
%

sz = size(profile);
if sz(1) == 2
  domain_methods(d.dom_id, 64, profile(1,:), profile(2,:));
elseif sz(2) == 2
  domain_methods(d.dom_id, 64, profile(:,1), profile(:,2));
else
  error('Wrong temperature profile array shape.');
end

