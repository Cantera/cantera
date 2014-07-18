function q = rop_f(a)
% ROP_F  Forward rates of progress for all reactions.
% q = rop_f(a)
%
%    See also: :mat:func:`rop_r`, :mat:func:`rop_net`, :mat:func:`rop`
%
% :param a:
%     Instance of class :mat:func:`Kinetics` (or another
%     object deriving from Kinetics)
%     for which forward rates of progress are desired.
% :return:
%     Returns a column vector of the forward rates of progress
%     for all reactions. If this function
%     is called with no output argument, a bar graph is produced.
%

q = kinetics_get(a.id, 11, 0);
if nargout == 0
    figure
    set(gcf, 'Name', 'Rates of Progress')
    bar(q)
    xlabel('Reaction Number')
    ylabel('Forward Rate of Progress [kmol/m^3]')
    title('Forward Rates of Progress')
end
