function d = setFixedTempProfile(d, profile)
% SETFIXEDTEMPPROFILE - set the temperature profile to use when the
% energy equation is not being solved. The profile must be entered
% as an array of positions / temperatures, which may be in rows or
% columns. 
%   
sz = size(profile);
if sz(1) == 2
  domain_methods(d.dom_id, 64, profile(1,:), profile(2,:));
elseif sz(2) == 2
  domain_methods(d.dom_id, 64, profile(:,1), profile(:,2));
else
  error('wrong temperature profile array shape');
end

