function v = kinetics_get(n, job, a, b, c, d, e, f, g)
% KINETICS_GET - get kinetics attributes
%
if nargin == 2
    v = ctmethods(40, n, job);
elseif nargin == 3
    v = ctmethods(40, n, job, a);
elseif nargin == 4
    v = ctmethods(40, n, job, a, b);
elseif nargin == 5
    v = ctmethods(40, n, job, a, b, c);
elseif nargin == 6
    v = ctmethods(40, n, job, a, b, c, d);
elseif nargin == 7
    v = ctmethods(40, n, job, a, b, c, d, e);
elseif nargin == 8
    v = ctmethods(40, n, job, a, b, c, d, e, f);
elseif nargin == 9
    v = ctmethods(40, n, job, a, b, c, d, e, f, g);
end
