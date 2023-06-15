function phase_set(n, job, a, b)
if nargin == 2
    ctmethods(30, n, -job);
elseif nargin == 3
    ctmethods(30, n, -job, a);
else
    ctmethods(30, n,-job, a, b);
end
