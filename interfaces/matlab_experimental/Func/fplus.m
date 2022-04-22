function r = fplus(a, b)
    % Get a functor representing the sum two input functors 'a' and
    % 'b'.
    r = Func('sum', a, b);
end
