function setProfile(s, name, comp, p)
% SETPROFILE - Specify a profile for one component.
%
%   name   --   domain name
%   comp   --   component number
%   zr     --   array of relative positions (0.0 to 1.0)
%   v      --   array of values
%
%   The solution vector values for this component will be linearly
%   interpolated from the discrete function defined by v vs. zr.
%   Note that zr = 0.0 corresponds to the leftmost grid point in
%   the specified domain, and zr = 1.0 corresponds to the rightmost
%   grid point. This method can be called at any time, but is
%   usually used to set the initial guess for the solution.
%
%   Example:
%
%      zr = [0 0.1 0.2 0.4 0.8 1];
%      v  = [500 650 700 730 800 900];
%      setProfile(1, 2, zr, v);
%
if isa(name,'double')
    n = name;
else
    n = domainIndex(s, name);
end

d = s.domains(n);

if isa(comp,'double') || isa(comp,'cell')
    c = comp;
elseif isa(comp,'char')
    c = {comp};
else
    error('wrong type');
end

np = length(c);
sz = size(p);
if sz(1) == np + 1;
    for j = 1:np
        ic = componentIndex(d,c{j});
        stack_methods(s.stack_id, 101, n, ic, p(1,:), p(j+1,:));
    end
elseif sz(2) == np + 1;
    ic = componentIndex(d,c{j});
    stack_methods(s.stack_id, 101, n, ic, p(:,1), p(:,j+1));
else
    error('wrong profile shape');
end
