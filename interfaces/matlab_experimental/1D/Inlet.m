function m = Inlet(id)
    % Create an inlet domain.
    m = Domain1D('Inlet1D');
    if nargin == 0
        m.setID('inlet');
    else
        m.setID(id);
    end
end
