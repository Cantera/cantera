function install(f, upstream, downstream)

if nargin == 3
    if ~isa(upstream,'Reactor') || ~isa(downstream,'Reactor')
        error(['Flow devices can only be installed between reactors or' ...
            ' reservoirs'])
    end
    i = hndl(upstream);
    j = hndl(downstream);
    ok = flowdevicemethods(2, f.index, i, j);
    if ok < 0
        error(geterr)
    end
else
    error('install requires 3 arguments')
end
