function q = rop_r(a)
% ROP_R  Get the reverse rates of progress for all reactions.
% q = rop_r(a)
%
%    See also: :mat:func:`rop_f`, :mat:func:`rop_net`, :mat:func:`rop`
%
% :param a:
%     Instance of class :mat:func:`Kinetics` (or another
%     object deriving from Kinetics)
%     for which reverse rates of progress are desired.
% :return:
%     Returns a column vector of the reverse rates of progress
%     for all reactions. If this function
%     is called with no output argument, a bar graph is produced.
%

q = kinetics_get(a.id, 12, 0);
if nargout == 0
    figure
    set(gcf, 'Name', 'Reverse Rates of Progress')
    bar(q)
    xlabel('Reaction Number')
    ylabel('Reverse Rate of Progress [kmol/m^3]')
    title('Reverse Rates of Progress')
end
