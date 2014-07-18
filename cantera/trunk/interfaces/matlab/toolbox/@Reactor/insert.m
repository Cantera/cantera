function insert(r, gas)
% INSERT  Insert a solution or mixture into a reactor.
% insert(r, gas)
% :param r:
%     Instance of class :mat:func:`Reactor`
% :param gas:
%     Instance of class :mat:func:`Solution` to be inserted
%

r.contents = gas;
setThermoMgr(r, gas);
setKineticsMgr(r, gas);
