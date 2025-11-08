function m = Methane()
    % Return an object representing methane. ::
    %
    %     >> m = ct.Methane()
    %
    % The object returned by this method implements an accurate equation of
    % state for methane that can be used in the liquid, vapor, saturated
    % liquid/vapor, and supercritical regions of the phase diagram.  The
    % equation of state is taken from
    %
    % Reynolds, W. C. *Thermodynamic Properties in SI: graphs, tables, and
    % computational equations for forty substances* Stanford: Stanford
    % University, 1979. Print.
    %
    % For more details, see classes :ct:`PureFluidPhase` and :ct:`methane` in the
    % Cantera C++ source code documentation.
    %
    % :return:
    %     Instance of class :mat:class:`ct.Solution`.

    m = ct.Solution('liquidvapor.yaml', 'methane');
end
