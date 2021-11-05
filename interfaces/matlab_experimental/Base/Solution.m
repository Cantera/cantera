classdef Solution < handle
    properties
        thermo
        kinetics
        transport
    end
    methods
        % Solution class constructor
        function s = Solution(src, id, trans)
            if nargin == 1
                id = '-';
            end
            tp = ThermoPhase(src, id);
            kin = Kinetics(tp, src, id);
            s.kinetics = kin;
            s.thermo = tp;
            if nargin == 3
                if (strcmp(trans, 'default') || strcmp(trans, 'None')...
                    || strcmp(trans, 'Mix') || strcmp(trans, 'Multi'))
                    tr = Transport(tp, trans, 0);
                else
                    error('Unknown transport modelling specified.');
                end
            else
                tr = Transport(tp, 'default', 0);
            end
            s.transport = tr;
        end

        % Delete the kernel objects associated with a solution
        function clear(s)
            s.thermo.clear;
            s.kinetics.clear;
            s.transport.clear;
        end
    end
end
