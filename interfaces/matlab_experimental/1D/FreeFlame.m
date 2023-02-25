classdef FreeFlame < Domain1D
    % Create a freely-propagating flat flame. ::
    %
    %     >> m = FreeFlame(gas, id)
    %
    % :param gas:
    %     Instance of class :mat:class:`Solution`.
    % :param id:
    %     String, ID of the flow.
    % :return:
    %     Instance of class :mat:class:`FreeFlame` representing
    %     a freely propagating, adiabatic flame.

    methods

        function m = FreeFlame(gas, id)
            % Constructor

            m@Domain1D('StagnationFlow', gas, 2);

            if nargin == 1
                m.setID('flame');
            else
                m.setID(id);
            end

        end

    end

end
