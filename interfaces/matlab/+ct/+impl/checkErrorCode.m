function checkErrorCode(code)
    switch code
        case {-1, -2}
            error('Cantera:ctError', ct.impl.getError());
        case -999.999
            error('Cantera:ctError', ct.impl.getError());
        case intmax('uint64')
            error('Cantera:ctError', ct.impl.getError());
    end
end
