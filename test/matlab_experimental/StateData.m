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

        function statedata = StateData(phase, T, P, varargin)
            param = inputParser;
            param.addRequired('phase', @ischar);
            param.addRequired('T', @isnumeric);
            param.addRequired('P', @isnumeric);
            param.addParameter('D', [], @isnumeric);
            param.addParameter('V', [], @isnumeric);
            param.addParameter('U', [], @isnumeric);
            param.addParameter('H', [], @isnumeric);
            param.addParameter('S', [], @isnumeric);
            param.addParameter('relax', false, @islogical);

            % Parse inputs
            param.parse(phase, T, P, varargin{:});
            args = param.Results;

            % Assign properties
            statedata.phase = args.phase;
            statedata.T = args.T;
            statedata.P = args.P * 1e6;

            if ~isempty(args.D)
                statedata.D = args.D;
            elseif ~isempty(args.V)
                statedata.D = 1.0 / args.V;
            else
                statedata.D = NaN;
            end

            if ~isempty(args.U)
                statedata.U = 1e3 * args.U;
            elseif ~isempty(args.H)
                statedata.U = 1e3 * args.H - statedata.P / statedata.D;
            else
                statedata.U = NaN;
            end

            statedata.S = 1e3 * args.S;
            statedata.tolMod = args.relax * 10.0 + ~args.relax * 1.0;
        end

    end
end
