classdef Tolerances
    % Container for relative tolerances used in thermodynamic validation tests.
    % Default values can be overridden as needed via name-value arguments.

    properties
        P
        U
        S
        dUdS
        dAdV
        dPdT
        hTs
    end

    methods

        function t = Tolerances(arg)
            arguments
                arg.P (1,1) double {mustBeNumeric} = 2e-5
                arg.U (1,1) double {mustBeNumeric} = 2e-6
                arg.S (1,1) double {mustBeNumeric} = 2e-6
                arg.dUdS (1,1) double {mustBeNumeric} = 2e-6
                arg.dAdV (1,1) double {mustBeNumeric} = 2e-6
                arg.dPdT (1,1) double {mustBeNumeric} = 2e-4
                arg.hTs (1,1) double {mustBeNumeric} = 2e-4
            end

            % Assign properties
            t.P = arg.P;
            t.U = arg.U;
            t.S = arg.S;
            t.dUdS = arg.dUdS;
            t.dAdV = arg.dAdV;
            t.dPdT = arg.dPdT;
            t.hTs = arg.hTs;
        end

    end
end
