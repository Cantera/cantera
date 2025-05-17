function w = Water(backend)
    % Return an object representing water. ::
    %
    %     >> w = Water()
    %
    % The object returned by this method implements an accurate equation of
    % state for water that can be used in the liquid, vapor, saturated
    % liquid/vapor, and supercritical regions of the phase diagram. The
    % equation of state is taken from
    %
    % Reynolds, W. C. *Thermodynamic Properties in SI: graphs, tables, and
    % computational equations for forty substances.* Stanford: Stanford
    % University, 1979. Print.
    %
    % For more details, see classes Cantera::PureFluid and tpx::water in the
    % Cantera C++ source code documentation.
    %
    % :return:
    %     Instance of class :mat:class:`Solution`.
    if nargin == 0 | strcmp(backend, 'Reynolds')
        w = Solution('liquidvapor.yaml', 'water', 'water');
    elseif strcmp(backend, 'IAPWS95')
        w = Solution('liquidvapor.yaml', 'liquid-water-IAPWS95', 'water');
    else
        error(['Unknow backend: ', backend]);
    end
end
