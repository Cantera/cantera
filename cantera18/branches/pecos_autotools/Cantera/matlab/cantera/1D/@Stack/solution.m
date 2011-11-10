function x = solution(s, domain, component)
% SOLUTION - get a solution component in one domain.
%
%    x = solution(s, 'flow', 'T') returns in vector x the values of
%    solution component 'T' in domain 'flow'.
%   
idom = domainIndex(s, domain);
d = s.domains(idom);

if nargin == 3
  icomp = componentIndex(d, component);
  for n = 1:nPoints(d)
    x(n) = stack_methods(s.stack_id, 30, idom, icomp, n);
  end
else
  nc = nComponents(d);
  np = nPoints(d);
  for m = 1:nc
    for n = 1:np
      x(m,n) = stack_methods(s.stack_id, 30, idom, m, n);
    end    
  end
end
