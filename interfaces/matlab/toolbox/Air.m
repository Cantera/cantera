function gas = Air()
% AIR  Create an object representing air.
% gas = Air()
% Air is modeled as an ideal gas mixture. The specification is taken
% from file ``air.xml``. Several reactions among oxygen and nitrogen are
% defined.
%
% :return:
%     Instance of class :mat:func:`Solution`
%

gas = Solution('air.xml', 'air');
