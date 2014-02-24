function insert(r, gas)
% INSERT - insert a mixture into the reactor
%
r.contents = gas;
setThermoMgr(r, gas);
setKineticsMgr(r, gas);
