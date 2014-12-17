function setProfile(s, name, comp, p)
% SETPROFILE  Specify a profile for one component.
% setProfile(s, name, comp, p)
% The solution vector values for this component will be linearly
% interpolated from the discrete function defined by p(:,1) vs. p(:,2).
% Note that ``p(1,1) = 0.0`` corresponds to the leftmost grid point in
% the specified domain, and ``p(1,n) = 1.0`` corresponds to the rightmost
% grid point. This method can be called at any time, but is
% usually used to set the initial guess for the solution.
%
% Example (assuming ``s`` is an instance of :mat:func:`Stack`)::
%
%     >> zr = [0 0.1 0.2 0.4 0.8 1];
%     >> v  = [500 650 700 730 800 900];
%     >> setProfile(s, 1, 2, [zr, v]);
%
% :param s:
%     Instance of class :mat:func:`Stack`
% :param name:
%     Domain name
% :param comp:
%     component number
% :param p:
%     n x 2 array, whose columns are the relative (normalized) positions
%     and the component values at those points. The number of positions
%     ``n`` is arbitrary.
%

if isa(name, 'double')
    n = name;
else
    n = domainIndex(s, name);
end

d = s.domains(n);

if isa(comp, 'double') || isa(comp, 'cell')
    c = comp;
elseif isa(comp, 'char')
    c = {comp};
else
    error('Wrong type.');
end

np = length(c);
sz = size(p);
if sz(1) == np + 1;
    for j = 1:np
        ic = componentIndex(d, c{j});
        stack_methods(s.stack_id, 101, n, ic, p(1,:), p(j+1,:));
    end
elseif sz(2) == np + 1;
    ic = componentIndex(d,c{j});
    stack_methods(s.stack_id, 101, n, ic, p(:,1), p(:,j+1));
else
    error('Wrong profile shape.');
end
