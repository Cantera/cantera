function setFunction(f, mf)
% SETMASSFLOWRATE -
%
if f.type == 1
    k = flowdevicemethods(5, f.index, func_hndl(mf));
    if k < 0
        error(geterr);
    end
else
    error('Mass flow rate can only be set for mass flow controllers')
end
