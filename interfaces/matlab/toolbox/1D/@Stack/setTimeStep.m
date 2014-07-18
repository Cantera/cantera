function setTimeStep(s, stepsize, steps)
% SETTIMESTEP  Specify a sequence of time steps.
% setTimeStep(s, stepsize, steps)
% :param stepsize:
%     Initial step size (s)
% :param steps:
%     Vector of number of steps to take before
%     re-attempting solution of steady-state problem. For
%     example, steps = [1, 2, 5, 10] would cause one time
%     step to be taken first the the steady-state
%     solution attempted. If this failed, two time steps
%     would be taken, etc.
%

stack_methods(s.stack_id, 112, stepsize, length(steps), steps)
