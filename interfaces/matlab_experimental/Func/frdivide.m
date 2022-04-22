function r = frdivide(a, b)
    % Get a functor that is the ratio of the input functors 'a' and
    % 'b'.
    r = Func('ratio', a, b);
end
