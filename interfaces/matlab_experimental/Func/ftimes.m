function r = ftimes(a, b)
     % Get a functor that multiplies two functors 'a' and 'b'
     r = Func('prod', a, b);
end
