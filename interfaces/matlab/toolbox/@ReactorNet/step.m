function t = step(r)
% STEP  Take one internal time step.
% t = step(r)
% The integrator used to integrate the ODEs (CVODE) takes
% variable-size steps, chosen so that a specified error
% tolerance is maintained. At times when the solution is rapidly
% changing, the time step becomes smaller to resolve the
% solution.
%
% Method :mat:func:`step` takes one internal time step and returns
% the network time after taking that step. This
% can be useful when it is desired to resolve a rapidly-changing
% solution.
%
% This method can be used as follows:
%
% .. code-block:: matlab
%
%     t = 0.0
%     tout = 0.1
%     while t < tout
%        t = step(r)
%        ,,,
%        <commands to save desired variables>
%        ...
%     end
%
% See also: :mat:func:`advance`
%
% :param r:
%     Instance of class :mat:func:`ReactorNet`
% :return:
%     Network time after the internal time step. Units: s
%

t = reactornetmethods(21, reactornet_hndl(r));
