function rop = rop(a)
% ROP  Get the forward and reverse rates of progress.
% rop = rop(a)
%
% See also: :mat:func:`rop_f`, :mat:func:`rop_r`, :mat:func:`rop_net`
%
% :param a:
%     Instance of class :mat:func:`Kinetics` (or another
%     object deriving from Kinetics)
%     for which forward and reverse rates of progress are desired.
% :return:
%     Returns an I x 2 array of reaction rates of
%     progress, where I is the number of reactions. The first
%     column contains the forward rates of progress, and the
%     second column the reverse rates. If this function
%     is called with no output argument, a bar graph is produced.
%

f = rop_f(a);
r = rop_r(a);
rop = [f r];
if nargout == 0
    figure
    set(gcf, 'Name', 'Rates of Progress');
    bar(rop);
    xlabel('Reaction Number');
    ylabel('Rate of Progress [kmol/m^3-s]');
    title('Rates of Progress');
    legend('Forward', 'Reverse');
end
