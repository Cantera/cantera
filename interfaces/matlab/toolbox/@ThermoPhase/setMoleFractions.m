function a = setMoleFractions(a,x,norm)
% SETMOLEFRACTIONS  Set the species mole fractions.
%
%   setMoleFractions(a,x)
%
%      If x is a vector of doubles, this call sets the species mole
%      fractions to the values in x and then scales them so that they
%      sum to 1.0. The length of x must equal the number of species.
%
%      If x is a string, then the mole fractions of the listed
%      species are set to the specified values, and all others are
%      set to zero. The syntax for the string is as shown in these
%      examples:
%         'CH4:4, O2:2, AR:3'
%         'SIH2:0.9,H2:0.1'
%
%   setMoleFractions(a, x, 'nonorm')
%
%      As above, but the normalization step is skipped. This only
%      works if x is an array, not a string. Since unnormalized mole
%      fractions can lead to unphysical results, 'nonorm' should be
%      used only in rare cases, such as computing partial
%      derivatives with respect to a species mole fraction.
%
%  Note that calling setMoleFractions leaves the temperature and
%  density unchanged, and therefore the pressure changes if the new
%  composition has a molar mass that is different than the old
%  composition. If it is desired to change the composition and hold
%  the pressure fixed, use method 'set' (defined in class Solution)
%  instead, or call setPressure after calling setMoleFractions.

if isa(x,'double')
    if nargin == 3
        if strcmp(norm,'nonorm')
            phase_set(a.tp_id,22,x);
        else
            phase_set(a.tp_id,20,x);
        end
    else
        phase_set(a.tp_id,20,x);
    end
    %
    % string input
    %
elseif isa(x,'char')
    phase_set(a.tp_id,30,x);
end
