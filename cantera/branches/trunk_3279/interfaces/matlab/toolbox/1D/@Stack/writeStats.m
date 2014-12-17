function writeStats(s)
% WRITESTATS  Print statistics for the current solution.
% writeStats(s)
% Prints a summary of the number of function and
% Jacobian evaluations for each grid, and the CPU time spent on
% each one.
%
% :param s:
%     Instance of class :mat:func:`Stack`
%

stack_methods(s.stack_id, 108);
