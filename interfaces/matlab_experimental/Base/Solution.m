classdef Solution < handle & ThermoPhase & Kinetics & Transport
    properties
        tp
    end
    methods
        % Solution class constructor
        function s = Solution(src, id, trans)
            if nargin == 1
                id = '-';
            end
            tp = ThermoPhase(src, id);
            s@ThermoPhase(src, id);
            s@Kinetics(tp, src, id);
            if nargin == 3
                if ~(strcmp(trans, 'default') || strcmp(trans, 'None')...
                     || strcmp(trans, 'Mix') || strcmp(trans, 'Multi'))
                    error('Unknown transport modelling specified.');
                end
            else
                trans = 'default';
            end
            s@Transport(tp, trans, 0);
            s.tp_id = tp.tp_id;
        end

        % Delete the kernel objects associated with a solution
        function clear(s)
            s.tp_clear;
            s.kin_clear;
            s.tr_clear;
        end

    end
end
