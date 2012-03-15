function a = equilibrate(a, xy, solver, rtol, maxsteps, maxiter, loglevel)
% EQUILIBRATE  Set the phase to a state of chemical equilibrium.
%
%        XY -- A two-letter string, which must be one of the set
%        ['TP','TV','HP','SP','SV','UV','PT','VT','PH','PS','VS','VU'].
%        If H, U, S, or V is specified, the value must be the specific
%        value (per unit mass).
%
%        solver -- specifies the equilibrium solver to use. If solver
%        = 0, a fast solver using the element potential method will be
%        used. If solver > 0, a slower but more robust Gibbs
%        minimization solver will be used. If solver < 0 or is
%        unspecified, the fast solver will be tried first, then if it
%        fails the other will be tried.
%
%        rtol -- the relative error tolerance.
%
%        maxsteps -- maximum number of steps in composition to take to
%        find a converged solution.
%
%        maxiter -- for the Gibbs minimization solver only, this
%        specifies the number of 'outer' iterations on T or P when
%        some property pair other than TP is specified.
%
%        loglevel -- set to a value > 0 to write diagnostic output to
%        a file in HTML format. Larger values generate more detailed
%        information. The file will be named 'equilibrate_log.html.'
%        Subsequent files will be named 'equilibrate_log1.html',
%        'equilibrate_log2.html', etc., so that log files are not
%        overwritten.
%
%

% use the ChemEquil solver by default
if nargin < 3
    solver = -1;
end
if nargin < 4
    rtol = 1.0e-9;
end
if nargin < 5
    maxsteps = 1000;
end
if nargin < 6
    maxiter = 100;
end
if nargin < 7
    loglevel = 0;
end

iok = thermo_set(a.tp_id, 50, xy, solver, rtol, maxsteps, maxiter, loglevel);
if iok < 0
    e = geterr;
    if e == 0
        e = 'unknown error';
    end
    error(e);
end

