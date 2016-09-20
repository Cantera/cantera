function setChemistry(r, flag)
% SETCHEMISTRY  Enable or disable changing reactor composition by reactions.
% setChemistry(r, flag)
% If the chemistry is disabled, then the reactor composition is
% constant. The parameter should be the string ``'on'`` to enable the
% species equations, or ``'off'`` to disable it.
%
% By default, Reactor objects are created with the species equations
% enabled if there are reactions present in the mechanism file, and
% disabled otherwise. ::
%
%     >> setChemistry(r, 'on');
%     >> setChemistry(r, 'off');
%
% :param r:
%     Instance of class :mat:func:`Reactor`
% :param flag:
%     String, either ``'on'`` or ``'off'`` to enable and disable
%     solving the energy equation, respectively
%

if strcmp(flag, {'on'})
    iflag = true;
elseif strcmp(flag, {'off'})
    iflag = false;
else
    error('Input to setChemistry not understood. It must be either "on" or "off".');
end
reactormethods(8, r.index, iflag);
