function kc = equil_Kc(a)
% equil_Kc(a)  equilibrium constants for all reactions
%
%    q = equil_Kc(a)
%
%        Returns a column vector of the equilibrium constants
%        for all reactions. The vector has an entry for every
%        reaction, whether reversible or not, but non-zero values
%        occur only for the reversible reactions.
%
%
kc = kinetics_get(a.id,14,0);
if nargout == 0
    figure
    set(gcf,'Name','Equilibrium Constants')
    bar(log10(kc))
    xlabel('Reaction Number')
    ylabel('log_1_0 Kc [kmol, m, s]')
    title('Equilibrium Constants Kc')
end
