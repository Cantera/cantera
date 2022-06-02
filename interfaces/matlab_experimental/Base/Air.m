function gas = Air()
    % Create an object representing air.
    % gas = Air()
    % Air is modeled as an ideal gas mixture. The specification is taken
    % from file ``air.yaml``. Several reactions among oxygen and nitrogen are
    % defined. Mixture-averaged transport is specified by default.
    %
    % :return:
    %     Instance of class :mat:func:`Solution`
    %
    gas = Solution('air.yaml', 'air');
end
