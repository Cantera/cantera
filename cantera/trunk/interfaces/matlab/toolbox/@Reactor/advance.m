function advance(r, tout)
% ADVANCE - Advance the state of the reactor in time.
%
%    Method advance integrates the system of ordinary differential
%    equations that determine the rate of change of the reactor
%    volume, the mass of each species, and the total energy. The
%    integration is carried out from the current reactor time to time
%    'tout.' (Note 'tout' is an absolute time, not a time interval.)  The
%    integrator may take many internal time steps before reaching
%    tout.
%
%       for i in 1:10
%          tout = 0.1*i
%          advance(r, tout)
%          ...
%          <add output commands here>
%          ...
%       end
%
%    See also: Reactor/step
%
reactormethods(8, reactor_hndl(r), tout);
