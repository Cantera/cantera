function s = char(p)
% CHAR  Get the formatted string to display the function.
% s = char(p)
% :param p:
%     Instance of class :mat:func:`Func`
% :return:
%     Formatted string displaying the function
%

if strcmp(p.typ,'sum')
    s = ['(' char(p.f1) ') + (' char(p.f2) ')'];
elseif strcmp(p.typ,'diff')
    s = ['(' char(p.f1) ') - (' char(p.f2) ')'];
elseif strcmp(p.typ,'prod')
    s = ['(' char(p.f1) ') * (' char(p.f2) ')'];
elseif strcmp(p.typ,'ratio')
    s = ['(' char(p.f1) ') / (' char(p.f2) ')'];
elseif all(p.coeffs == 0)
    s = '0';
else
    if strcmp(p.typ,'polynomial')
        d = length(p.coeffs) - 1;
        s = [];
        nn = 0;
        for b = p.coeffs;
            cc(d+1-nn) = b;
            nn = nn + 1;
        end
        for a = cc;
            if a ~= 0;
                if ~isempty(s)
                    if a > 0
                        s = [s ' + '];
                    else
                        s = [s ' - '];
                        a = -a;
                    end
                end
                if a ~= 1 || d == 0
                    s = [s num2str(a)];
                    if d > 0
                        s = [s '*'];
                    end
                end
                if d >= 2
                    s = [s 'x^' int2str(d)];
                elseif d == 1
                    s = [s 'x'];
                end
            end
            d = d - 1;
        end
    elseif strcmp(p.typ, 'gaussian')
        s = ['Gaussian(' num2str(p.coeffs(1)) ',' ...
            num2str(p.coeffs(2)) ',' ...
            num2str(p.coeffs(3)) ')'];
    elseif strcmp(p.typ, 'fourier')
        c = reshape(p.coeffs, [], 2);
        Ao = c(1, 1);
        w = c(1, 2);
        A = c(2:end, 1);
        B = c(2:end, 2);
        N = size(c, 1) - 1;
        if Ao ~= 0
            s = num2str(Ao/2);
        else
            s = '';
        end
        for n=1:N
            if A(n) ~= 0
                if A(n) < 0
                    prefix = ' - ';
                elseif s
                    prefix = ' + ';
                else
                    prefix = '';
                end

                s = [s prefix num2str(abs(A(n))) '*cos(' num2str(n*w) '*x)'];
            end

            if B(n) ~= 0
                if B(n) < 0
                    prefix = ' - ';
                elseif s
                    prefix = ' + ';
                else
                    prefix = '';
                end

                s = [s prefix num2str(abs(B(n))) '*sin(' num2str(n*w) '*x)'];
            end
        end
    else
        s = ['*** char not yet implemented for' p.typ ' ***'];
    end
end
