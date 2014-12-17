function setFunction(f, mf)
% SETFUNCTION  Set the mass flow rate with class :mat:func:`Func`.
% setFunction(f, mf)
%
% See also: :mat:func:`MassFlowController`, :mat:func:`Func`
%
% :param f:
%     Instance of class :mat:func:`MassFlowController`
% :param mf:
%     Instance of class :mat:func:`Func`
%

if f.type == 1
    k = flowdevicemethods(5, f.index, func_hndl(mf));
    if k < 0
        error(geterr);
    end
else
    error('Mass flow rate can only be set for mass flow controllers')
end
