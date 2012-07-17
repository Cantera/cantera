function v = mixturemethods(n, job, a, b, c, d, e, f)
%
if nargin == 2
    v = ctmethods(120, n, job);
elseif nargin == 3
    v = ctmethods(120, n, job, a);
elseif nargin == 4
    v = ctmethods(120, n, job, a, b);
elseif nargin == 5
    v = ctmethods(120, n, job, a, b, c);
elseif nargin == 6
    v = ctmethods(120, n, job, a, b, c, d);
elseif nargin == 7
    v = ctmethods(120, n, job, a, b, c, d, e);
elseif nargin == 8
    v = ctmethods(120, n, job, a, b, c, d, e, f);
end
