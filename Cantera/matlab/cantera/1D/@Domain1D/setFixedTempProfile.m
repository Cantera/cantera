function d = setFixedTempProfile(d, profile)
% SETFIXEDTEMPPROFILE - 
%   
sz = size(profile);
if sz(1) == 2
  domain_methods(d.dom_id, 64, profile(1,:), profile(2,:));
elseif sz(2) == 2
  domain_methods(d.dom_id, 64, profile(:,1), profile(:,2));
else
  error('wrong temperature profile array shape');
end

