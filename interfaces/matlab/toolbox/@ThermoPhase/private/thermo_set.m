function i = thermo_set(n, job, a, b, c, d, e, f)
if nargin == 2
    i = ctmethods(20, n, -job);
elseif nargin == 3
    i = ctmethods(20, n, -job,a);
elseif nargin == 4
    i = ctmethods(20, n, -job, a, b);
elseif nargin == 5
    i = ctmethods(20, n, -job, a, b, c);
elseif nargin == 6
    i = ctmethods(20, n, -job, a, b, c, d);
elseif nargin == 7
    i = ctmethods(20, n, -job, a, b, c, d, e);
elseif nargin == 8
    i = ctmethods(20, n, -job, a, b, c, d, e, f);
end
