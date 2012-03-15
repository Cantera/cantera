function v = stack_methods(n, job, a, b, c, d, e, f)
% STACK_METHODS - converter function for methods of class Stack
%
%    All Cantera functions and methods are handled by the single MEX
%    file 'ctmethods.' This function is provided only for convenience,
%    and simply calls ctmethods with a flag associated with this class
%    as the first parameter, followed by the input arguments.
if nargin == 2
    v = ctmethods(90, n, job);
elseif nargin == 3
    v = ctmethods(90, n, job, a);
elseif nargin == 4
    v = ctmethods(90, n, job, a, b);
elseif nargin == 5
    v = ctmethods(90, n, job, a, b, c);
elseif nargin == 6
    v = ctmethods(90, n, job, a, b, c, d);
elseif nargin == 7
    v = ctmethods(90, n, job, a, b, c, d, e);
elseif nargin == 8
    v = ctmethods(90, n, job, a, b, c, d, e, f);
else
    error('too many arguments');
end
