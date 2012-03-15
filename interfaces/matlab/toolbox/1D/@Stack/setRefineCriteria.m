function d = setRefineCriteria(d, n, ratio, slope, curve, prune)
% SETREFINECRITERIA - Set the criteria used to refine the grid.
%
%     n  --  domain number beginning with domain 1 at the left
%     ratio  --   maximum size ratio between adjacent cells
%     slope  --   maximum relative difference in value between
%                 adjacent points
%     curve  --   maximum relative difference in slope between
%                 adjacent cells
%     prune  --   minimum value for slope or curve for which points
%                 will be retained in the grid. If the computed
%                 slope or curve value is below prune for all
%                 components, it will be deleted, unless either
%                 neighboring point is already marked for deletion.
%
if nargin < 3
    ratio = 10.0;
end
if nargin < 4
    slope = 0.8;
end
if nargin < 5
    curve = 0.8;
end
if nargin < 6
    prune = -0.1;
end

stack_methods(d.stack_id, 106, n, ratio, slope, curve, prune);
