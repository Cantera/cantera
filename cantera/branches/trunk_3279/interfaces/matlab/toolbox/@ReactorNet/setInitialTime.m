function setInitialTime(r, t)
% SETINITIALTIME  Set the initial time of the integration.
% setInitialTime(r, t)
% If the time integration has already begun, this restarts the
% integrator using the current solution as the initial condition,
% setting the time to ``t``.
%
% :param r:
%     Instance of class :mat:func:`ReactorNet`
% :param t:
%     Time at which integration should be restarted, using the
%     current state as the initial condition. Units: s
%

reactornetmethods(5, reactornet_hndl(r), t);
