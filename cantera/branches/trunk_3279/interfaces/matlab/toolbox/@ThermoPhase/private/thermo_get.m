function i = thermo_get(n, job, a, b)
if nargin == 2
    i = ctmethods(20, n, job);
elseif nargin == 3
    i = ctmethods(20, n, job, a);
else
    i = ctmethods(20, n, job, a, b);
end
