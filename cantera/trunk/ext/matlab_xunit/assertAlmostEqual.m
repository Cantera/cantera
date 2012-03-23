function assertAlmostEqual(A, B, reltol, message)
%assertEqual Assert that inputs are equal within relative tolerance
%   assertEqual(A, B, RELTOL) throws an exception of any of the values in A and
%   B are not equal within the specified tolerance.  NaN values are considered
%   to be equal.  A and B have to have the same class and sparsity to be
%   considered equal.
%
%   assertEqual(A, B) uses the following relative tolerance value:
%
%       100 * eps(class(A))
%
%   assertEqual(A, B, RELTOL, MESSAGE) uses the specified message string when
%   throwing the exception.  With this syntax, use RELTOL = [] to specify the
%   default relative tolerance.
%
%   Note that if either A or B are not floating-point arrays, then A and B are
%   compared using ISEQUALWITHEQUALNANS and the relative tolerance value is not
%   used. 
%
%   Examples
%   --------
%   % This call returns silently.
%   assertAlmostEqual(1.0, 1.0 + eps);
%
%   % This call throws an error.
%   assertAlmostEqual(1.0, 1.1);
%
%   See also assertEqual, mtest.utils.isAlmostEqual

%   Steven L. Eddins
%   Copyright 2008-2009 The MathWorks, Inc.

if ~(issparse(A) == issparse(B))
   throw(MException('assertAlmostEqual:sparsityNotEqual', message));
end

if ~strcmp(class(A), class(B))
   throw(MException('assertAlmostEqual:classNotEqual', message));
end

if nargin < 3 || isempty(reltol)
    reltol = 100 * eps(class(A));
end

if nargin < 4
    message = sprintf('Inputs are not equal within relative tolerance: %g', ...
        reltol);
end

if ~xunit.utils.isAlmostEqual(A, B, reltol)
   throw(MException('assertAlmostEqual:tolExceeded', message));
end
