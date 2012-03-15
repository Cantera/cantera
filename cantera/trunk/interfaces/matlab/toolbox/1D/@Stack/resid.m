function r = resid(s, domain, rdt, count)
if nargin == 2
    rdt = 0.0;
    count = 0;
end

idom = domainIndex(s, domain);
d = s.domains(idom);

r = zeros(nComponents(d), nPoints(d));
stack_methods(s.stack_id, 113, rdt, count);
for m = 1:nComponents(d)
    for n = 1:nPoints(d)
        r(m,n) = stack_methods(s.stack_id, 31, idom, m, n);
    end
end
