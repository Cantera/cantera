function tp = equilibrate(tp, xy, solver, rtol, maxsteps, maxiter, loglevel)
% EQUILIBRATE  Set the phase to a state of chemical equilibrium.
% tp = equilibrate(tp, xy, solver, rtol, maxsteps, maxiter, loglevel)
% :param XY:
%     A two-letter string, which must be one of the set
%     ``['TP','TV','HP','SP','SV','UV','UP']``,
%     indicating which pair of properties should be held constant.
%     Not all of the properties to be held constant are available with
%     all of the solvers.
% :param solver:
%     Specifies the equilibrium solver to use. If solver = 0, a fast
%     solver using the element potential method will be used. If
%     solver = 1, a slower but more robust Gibbs minimization solver
%     will be used. If solver >= 2, a version of the VCS algorithm will
%     be used. If solver < 0 or is unspecified, the fast solver
%     will be tried first, then if it fails the Gibbs minimization solver
%     will be tried.
% :param rtol:
%     The relative error tolerance.
% :param maxsteps:
%     Maximum number of steps in composition to take to find a
%     converged solution.
% :param maxiter:
%     For the Gibbs minimization solver only, this specifies the number
%     of 'outer' iterations on T or P when some property pair other than
%     TP is specified.
% :param loglevel:
%     Set to a value > 0 to write diagnostic output. Larger values
%     generate more detailed information.
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

iok = thermo_set(tp.tp_id, 50, xy, solver, rtol, maxsteps, maxiter, loglevel);
if iok < 0
    e = geterr;
    if e == 0
        e = 'unknown error';
    end
    error(e);
end

