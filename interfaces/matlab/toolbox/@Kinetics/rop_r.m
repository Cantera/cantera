function q = rop_r(a)
% ROP_R  Reverse rates of progress for all reactions.
%
%    Q = ROP_R(K)
%
%        Returns a column vector of the reverse rates of progress
%        for all reactions. The value is zero for irreversible
%        reactions.
%
%    See also: rop_r, rop_net.
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
