function solve(s, loglevel, refine_grid)
% SOLVE  Solve the problem.
% solve(s, loglevel, refine_grid)
% :param s:
%     Instance of class :mat:func:`Stack`
% :param loglevel:
%     Integer flag controlling the amount of diagnostic
%     output. Zero supresses all output, and 5 produces
%     very verbose output.
% :param refine_grid:
%     Integer, 1 to allow grid refinement, 0 to disallow.
%

stack_methods(s.stack_id, 104, loglevel, refine_grid);
