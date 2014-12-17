function setupGrid(d, grid)
% SETUPGRID  Set up the solution grid.
% d = setupGrid(d, grid)
% :param d:
%     Instance of class :mat:func:`Domain1D`
% :param grid:
%

domain_methods(d.dom_id, 53, grid);
