function setInitialTime(r, t)
% SETINITIALTIME - Set the initial time at which to begin the time
% integration.
%
% If the time integration has already begun, this restarts the
% integrator using the current solution as the initial condition,
% setting the time to 't'.
%
reactornetmethods(5, reactornet_hndl(r), t);
