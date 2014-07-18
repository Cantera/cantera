function x = solution(s, domain, component)
% SOLUTION  Get a solution component in one domain.
% x = solution(s, domain, component)
% :param s:
%     Instance of class :mat:func:`Stack`
% :param domain:
%     String, name of the domain from which the solution is desired
% :param component:
%     String, component for which the solution is desired. If omitted,
%     solutions for all of the components will be returned in an
%     :mat:func:`nPoints` x :mat:func:`nComponents` array.
% :return:
%     Either an :mat:func:`nPoints` x 1 vector, or
%     :mat:func:`nPoints` x :mat:func:`nComponents` array.
%

idom = domainIndex(s, domain);
d = s.domains(idom);
np = nPoints(d);
if nargin == 3
    icomp = componentIndex(d, component);
    x = zeros(1, np);
    for n = 1:np
        x(n) = stack_methods(s.stack_id, 30, idom, icomp, n);
    end
else
    nc = nComponents(d);
    x = zeros(nc, np);
    for m = 1:nc
        for n = 1:np
            x(m,n) = stack_methods(s.stack_id, 30, idom, m, n);
        end
    end
end
