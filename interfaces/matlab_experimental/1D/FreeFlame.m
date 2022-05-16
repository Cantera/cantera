function m = FreeFlame(gas, id)
    % Create a freely-propagating flat flame.
    m = Domain1D('StagnationFlow', gas, 2);
    if nargin == 1
        m.setID('flame');
    else
        m.setID(id);
    end
end
