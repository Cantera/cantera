function display(self, threshold)
if nargin < 2
    threshold = 1e-14
end
phase_get(thermo_hndl(self), 15, 1, threshold);
