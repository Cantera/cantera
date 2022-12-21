classdef AxisymmetricFlow < Domain1D
    % Create an axisymmetric flow domain.
    %
    % m = AxisymmetricFlow(gas, id)
    %
    % :param gas:
    %     Instance of class :mat:class:`Solution`
    % :param id:
    %     String, ID of the flow
    % :return:
    %     Instance of class :mat:class:`AxisymmetricFlow`.
    %

    methods

        % Constructor
        function m = AxisymmetricFlow(gas, id)

            m = m@Domain1D('StagnationFlow', gas);

            if nargin == 1
                m.setID('flow');
            else
                m.setID(id);
            end

        end

    end

end
