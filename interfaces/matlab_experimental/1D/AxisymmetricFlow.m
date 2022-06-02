function m = AxisymmetricFlow(gas, id)
    % Create an axisymmetric flow domain.
    % m = AxisymmetricFlow(gas, id)
    %
    % :param gas:
    %     Instance of class :mat:func:`Solution`
    % :param id:
    %     String, ID of the flow
    % :return:
    %     Domain1D instance representing an axisymmetric flow.
    %
    m = Domain1D('StagnationFlow', gas);
    if nargin == 1
        m.setID('flow');
    else
        m.setID(id);
    end
end
