function r = ConstPressureReactor(contents)
% CONSTPRESSUREREACTOR - Create a ConstPressureReactor object.
%
%    A ConstPressureReactor is an instance of class Reactor where the pressure
%    is held constant. The volume is not a state variable, but instead takes
%    on whatever value is consistent with holding the pressure constant.
%
%        r1 = ConstPressureReactor         % an empty reactor
%        r2 = ConstPressureReactor(gas)    % a reactor containing a gas
%
%   See also: Reactor
%
if nargin == 0
    contents = 0;
end
r = Reactor(contents, 4);
