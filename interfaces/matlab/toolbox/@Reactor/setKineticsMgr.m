function setKineticsMgr(r, k)
% SETKINETICSMGR - set the kinetics manager. This method is used
% internally during Reactor initialization, but is usually not
% called by users.
%
if ~isa(k,'Kinetics')
    error('wrong object type');
end

reactormethods(7, reactor_hndl(r), kinetics_hndl(k));
