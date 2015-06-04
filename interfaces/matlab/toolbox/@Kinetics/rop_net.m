function q = rop_net(a)
% ROP_NET  Net rates of progress for all reactions.
% q = rop_net(a)
%
% See also: :mat:func:`rop_f`, :mat:func:`rop_r`, :mat:func:`rop`
%
% :param a:
%     Instance of class :mat:func:`Kinetics` (or another
%     object deriving from Kinetics)
%     for which the net rates of progress are desired.
% :return:
%     Returns a column vector of the net rates of progress
%     for all reactions. If this function
%     is called with no output argument, a bar graph is produced.
%

q = kinetics_get(a.id, 13, 0);
if nargout == 0
    figure
    set(gcf, 'Name', 'Net Rates of Progress')
    bar(q)
    xlabel('Reaction Number')
    ylabel('Net Rate of Progress [kmol/m^3]')
    title('Net Rates of Progress')
end
