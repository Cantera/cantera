% Cantera Toolbox.
% Version 1.0   17-April-2002
%
% Cantera is an open-source, object-oriented software package to aid
% the analysis and simulation of chemically-reacting flows. It
% consists of a computational kernel written largely in C++, along
% with interface libraries for several programming languages,
% including Fortran 90, Python, and MATLAB. More information about
% Cantera is available at http://www.cantera.org.
%
% This toolbox is the Cantera MATLAB interface.
%
% Cantera object construction.
%   Solution      - construct a Solution object.
%
% Constants.
%   oneatm      - One atmosphere.
%   gasconstant - Universal gas constant.
%
% Gas mixture models.
%   IdealGasMix - Ideal gas mixtures.
%   GRI30       - GRI-Mech 3.0.
%   air         - air.
%
% Zero-dimensional reactor models.
%   reactor     - A general, customizable reactor.
%   conhp       - An adiabatic, constant pressure reactor.
%   conuv       - An adiabatic, constant volume reactor.
%
% Utilities.
%   geterr      - Get Cantera error message.
%   adddir      - Add a directory to Cantera's search path
%   ctclear     - Clear all objects from memory.
%
%   Copyright 2002 California Institute of Technology
