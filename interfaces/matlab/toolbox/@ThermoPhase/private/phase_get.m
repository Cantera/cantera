function i = phase_get(n, job, a, b)
if nargin == 2
    i = ctmethods(30, n, job);
elseif nargin == 3
    i = ctmethods(30, n, job, a);
else
    i = ctmethods(30, n, job, a, b);
end