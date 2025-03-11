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
            addRequired(param, 'phase', @ischar);
            addRequired(param, 'T', @isnumeric);
            addRequired(param, 'P', @isnumeric);
            addParameter(param, 'D', [], @isnumeric);
            addParameter(param, 'V', [], @isnumeric);
            addParameter(param, 'U', [], @isnumeric);
            addParameter(param, 'H', [], @isnumeric);
            addParameter(param, 'S', [], @isnumeric);
            addParameter(param, 'relax', false, @islogical);

            % Parse inputs
            parse(param, phase, T, P, varargin{:});
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
