function r = fplus(a, b)
    % Get a functor representing the sum of two input functors.
    %
    % r = fplus(a, b)
    %
    % :param a:
    %     Instance of class :mat:class:`Func`
    % :param b:
    %     Instance of class :mat:class:`Func`
    % :return:
    %     Instance of class :mat:class:`fplus`
    %
    methods

        % Constructor
        function f = fplus(a, b)
            f = f@Func('sum', a, b);
        end

    end

end
