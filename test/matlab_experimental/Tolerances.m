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
            param.addParameter('P', 2e-5, @isnumeric);
            param.addParameter('U', 2e-6, @isnumeric);
            param.addParameter('S', 2e-6, @isnumeric);
            param.addParameter('dUdS', 2e-6, @isnumeric);
            param.addParameter('dAdV', 2e-6, @isnumeric);
            param.addParameter('dPdT', 2e-4, @isnumeric);
            param.addParameter('hTs', 2e-4, @isnumeric);

            % Parse inputs
            param.parse(varargin{:});
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
