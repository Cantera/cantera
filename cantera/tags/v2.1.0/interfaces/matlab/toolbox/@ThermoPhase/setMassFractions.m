function a = setMassFractions(a, y, norm)
% SETMASSFRACTIONS  Set the species mass fractions.
%
%   setMassFractions(a,y)
%
%      If y is a vector of doubles, this call sets the species mass
%      fractions to the values in y and then scales them so that they
%      sum to 1.0. The length of y must equal the number of species.
%
%      If y is a string, then the mass fractions of the listed
%      species are set to the specified values, and all others are
%      set to zero. The syntax for the string is as shown in these
%      examples:
%         'CH4:4, O2:2, AR:3'
%         'SIH2:0.9,H2:0.1'
%
%   setMassFractions(a, y, 'nonorm')
%
%      As above, but the normalization step is skipped. This only
%      works if y is an array, not a string. Since unnormalized mass
%      fractions can lead to unphysical results, 'nonorm' should be
%      used only in rare cases, such as computing partial
%      derivatives with respect to a species mass fraction.
%
%  Note that calling setMassFractions leaves the temperature and
%  density unchanged, and therefore the pressure changes if the new
%  composition has a molar mass that is different than the old
%  composition. If it is desired to change the composition and hold
%  the pressure fixed, use method 'set' (defined in class Solution)
%  instead, or call setPressure after calling setMassFractions.
%
if isa(y,'double')
    if nargin == 3
        if strcmp(norm,'nonorm')
            phase_set(a.tp_id,23,y);
        else
            phase_set(a.tp_id,21,y);
        end
    else
        phase_set(a.tp_id,21,y);
    end
    %
    % string input
    %
elseif isa(y,'char')
    phase_set(a.tp_id,31,y);
end
