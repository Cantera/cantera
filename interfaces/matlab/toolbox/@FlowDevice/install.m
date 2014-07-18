function install(f, upstream, downstream)
% INSTALL  Install a flow device between reactors or reservoirs.
% install(f, upstream, downstream)
% :param f:
%     Instance of class :mat:func:`FlowDevice` to install
% :param upstream:
%     Upstream :mat:func:`Reactor` or :mat:func:`Reservoir`
% :param downstream:
%     Downstream :mat:func:`Reactor` or :mat:func:`Reservoir`
% :return:
%     Instance of class :mat:func:`FlowDevice`
%

if nargin == 3
    if ~isa(upstream, 'Reactor') || ~isa(downstream, 'Reactor')
        error(['Flow devices can only be installed between reactors or' ...
            ' reservoirs'])
    end
    i = reactor_hndl(upstream);
    j = reactor_hndl(downstream);
    ok = flowdevicemethods(2, f.index, i, j);
    if ok < 0
        error(geterr)
    end
else
    error('install requires 3 arguments')
end
