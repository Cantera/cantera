function setKineticsMgr(r, k)
% SETKINETICSMGR - set the kinetics manager
%
if ~isa(k,'Kinetics')
  error('wrong object type');
end

reactormethods(7, reactor_hndl(r), kinetics_hndl(k));

