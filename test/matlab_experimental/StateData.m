classdef StateData

    properties
        phase
        T
        P
        D
        U
        S
        tolMod
    end

    methods

        function statedata = StateData(phase, T, P, arg)
            arguments
                phase char
                T double
                P double
                arg.relax (1,1) logical = false
                arg.D (1,1) double {mustBeNumeric} = NaN
                arg.V (1,1) double {mustBeNumeric} = NaN
                arg.U (1,1) double {mustBeNumeric} = NaN
                arg.H (1,1) double {mustBeNumeric} = NaN
                arg.S (1,1) double {mustBeNumeric} = NaN
            end

            % Assign properties
            statedata.phase = phase;
            statedata.T = T;
            statedata.P = P * 1e6;

            if ~isnan(arg.D)
                statedata.D = arg.D;
            elseif ~isnan(arg.V)
                statedata.D = 1.0 / arg.V;
            else
                statedata.D = NaN;
            end

            if ~isnan(arg.U)
                statedata.U = 1e3 * arg.U;
            elseif ~isnan(arg.H)
                statedata.U = 1e3 * arg.H - statedata.P / statedata.D;
            else
                statedata.U = NaN;
            end

            statedata.S = 1e3 * arg.S;
            statedata.tolMod = arg.relax * 10.0 + ~arg.relax * 1.0;
        end

    end
end
