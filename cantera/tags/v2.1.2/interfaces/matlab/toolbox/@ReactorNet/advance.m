function advance(n, tout)
% ADVANCE - Advance the state of the reactor network in time.
%
%    Method advance integrates the system of ordinary differential
%    equations that determine the rate of change of the volume, the
%    mass of each species, and the total energy for each reactor. The
%    integration is carried out from the current time to time
%    'tout.' (Note 'tout' is an absolute time, not a time interval.)
%    The integrator may take many internal time steps before reaching
%    tout.
%
%       for i in 1:10
%          tout = 0.1*i
%          advance(n, tout)
%          ...
%          <add output commands here>
%          ...
%       end
%
%    See also: ReactorNet/step
%
reactornetmethods(8, reactornet_hndl(n), tout);
