% Tutorial 3:   Getting Help
%
% Keywords: tutorial

help tut3

% Suppose you have created a Cantera object and want to know what
% methods are available for it, and get help on using the methods.

g = GRI30

% The first thing you need to know is the MATLAB class object g
% belongs to. Type:

class(g)

% This tells you that g belongs to a class called 'Solution'. To find
% the methods for this class, type

methods Solution

% This command returns only a few method names. These are the ones
% directly defined in this class. But solution inherits many other
% methods from base classes. To see all of its methods, type

methods Solution -full

% Now a long list is printed, along with a specification of the class
% the method is inherited from. For example, 'setPressure' is
% inherited from a class 'ThermoPhase'. Don't be concerned at this
% point about what these base classes are - we'll come back to them
% later.

% Now that you see what methods are available, you can type 'help
% <method_name>' to print help text for any method. For example,

help setTemperature
help setMassFractions
help rop_net

% For help on how to construct objects of a given class, type 'help
% <classname>'

help Solution

% Now that you know how to get help when you need it, you can
% explore using the Cantera Toolbox on your own. But there are a
% few more useful things to know, which are described in the next
% few tutorials.

clear all
cleanup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   end of tutorial 3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
