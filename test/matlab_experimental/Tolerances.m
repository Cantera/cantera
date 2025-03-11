classdef Tolerances

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

        function t = Tolerances(varargin)
            % Set up input parser
            param = inputParser;
            addParameter(param, 'P', 2e-5, @isnumeric);
            addParameter(param, 'U', 2e-6, @isnumeric);
            addParameter(param, 'S', 2e-6, @isnumeric);
            addParameter(param, 'dUdS', 2e-6, @isnumeric);
            addParameter(param, 'dAdV', 2e-6, @isnumeric);
            addParameter(param, 'dPdT', 2e-4, @isnumeric);
            addParameter(param, 'hTs', 2e-4, @isnumeric);

            % Parse inputs
            parse(param, varargin{:});
            args = param.Results;

            % Assign properties
            t.P = args.P;
            t.U = args.U;
            t.S = args.S;
            t.dUdS = args.dUdS;
            t.dAdV = args.dAdV;
            t.dPdT = args.dPdT;
            t.hTs = args.hTs;
        end

    end
end
