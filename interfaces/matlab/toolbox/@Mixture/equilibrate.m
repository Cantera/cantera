function r = equilibrate(self, XY, err, maxsteps, maxiter, ...
    loglevel)
%
% EQUILIBRATE - Set the mixture to a state of chemical equilibrium.
%
%       This method uses a version of the VCS algorithm to find the
%       composition that minimizes the total Gibbs free energy of the
%       mixture, subject to element conservation constraints. For a
%       description of the theory, see Smith and Missen, "Chemical
%       Reaction Equilibrium."  The VCS algorithm is implemented in
%       Cantera kernel class MultiPhaseEquil.
%
%       The VCS algorithm solves for the equilibrium composition for
%       specified temperature and pressure. If any other property pair
%       other than "TP" is specified, then an outer iteration loop is
%       used to adjust T and/or P so that the specified property
%       values are obtained.
%
%       XY - Two-letter string specifying the two properties to hold
%       fixed.  Currently, 'TP', 'HP', and 'SP' are
%       implemented. Default: 'TP'.
%
%       err - Error tolerance. Iteration will continue until (Delta
%       mu)/RT is less than this value for each reaction. Default:
%       1.0e-9. Note that this default is very conservative, and good
%       equilibrium solutions may be obtained with larger error
%       tolerances.
%
%       maxsteps - Maximum number of steps to take while solving the
%       equilibrium problem for specified T and P. Default: 1000.
%
%       maxiter - Maximum number of temperature and/or pressure
%       iterations.  This is only relevant if a property pair other
%       than (T,P) is specified. Default: 200.
%
%       loglevel -- set to a value > 0 to write diagnostic output.
%       Larger values generate more detailed information.
%
%          >> equilibrate(mix, 'TP')
%          >> equilibrate('TP', 1.0e-6, 500)
%
if nargin < 6
    loglevel = 0;
end
if nargin < 5
    maxiter = 200;
end
if nargin < 4
    maxsteps = 1000;
end
if nargin < 3
    err = 1.0e-9;
end
if nargin < 2
    XY = 'TP';
end
r = mixturemethods(31, mix_hndl(self), XY, err, maxsteps, maxiter, loglevel);
