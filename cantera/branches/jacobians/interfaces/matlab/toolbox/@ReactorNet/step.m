function t = step(r, tout)
% STEP - Take one internal time step toward tout.
%
%    The integrator used to integrate the ODEs (CVODE) takes
%    variable-size steps, chosen so that a specified error
%    tolerance is maintained. At times when the solution is rapidly
%    changing, the time step becomes smaller to resolve the
%    solution.
%
%    Method 'step' takes one internal time step and returns. This
%    can be useful when it is desired to resolve a rapidly-changing
%    solution.
%
%    This method can be used as follows:
%
%       t = 0.0
%       tout = 0.1
%       while t < tout
%          t = step(r, tout)
%          ,,,
%          <commands to save desired variables>
%          ...
%       end
%
%    See also: Reactor/advance
%
t = reactornetmethods(21, reactornet_hndl(r), tout);
