function v = reactornetmethods(n, job, a, b, c, d)
%
if nargin == 2
    v = ctmethods(65, n, job);
elseif nargin == 3
    v = ctmethods(65, n, job, a);
elseif nargin == 4
    v = ctmethods(65, n, job, a, b);
elseif nargin == 5
    v = ctmethods(65, n, job, a, b, c);
elseif nargin == 6
    v = ctmethods(65, n, job, a, b, c, d);
end
