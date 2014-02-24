function rop = rop(k)
% ROP - Forward and reverse rates of progress.
%
%    ROP(K) returns an M x 2 array of reaction rates of
%    progress. The first column contains the forward rates of progress,
%    and the second column the reverse rates. If this function
%    is called with no output argument, a bar graph is produced.
%
f = rop_f(k);
r = rop_r(k);
rop = [f r];
if nargout == 0
    figure
    set(gcf,'Name','Rates of Progress');
    bar(rop);
    xlabel('Reaction Number');
    ylabel('Rate of Progress [kmol/m^3-s]');
    title('Rates of Progress');
    legend('Forward', 'Reverse');
end
