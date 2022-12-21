function r = frdivide(a, b)
    % Get a functor that is the ratio of the input functors.
    %
    % r = frdivide(a,b)
    %
    % :param a:
    %     Instance of class :mat:class:`Func`
    % :param b:
    %     Instance of class :mat:class:`Func`
    % :return:
    %     Instance of class :mat:class:`frdivide`
    %
    methods

        % Constructor
        function f = frdivide(a, b)
            f = f@Func('ratio', a, b);
        end

    end

end
