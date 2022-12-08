function m = FreeFlame(gas, id)
    % Create a freely-propagating flat flame.
    % m = FreeFlame(gas, id)
    %
    % :param gas:
    %     Instance of class :mat:func:`Solution`
    % :param id:
    %     String, ID of the flow
    % :return:
    %     Domain1D instance representing a freely propagating,
    %     adiabatic flame
    %
    m = Domain1D('StagnationFlow', gas, 2);

    if nargin == 1
        m.setID('flame');
    else
        m.setID(id);
    end

end
