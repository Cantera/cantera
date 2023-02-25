%% DIAMOND_CVD - A CVD example simulating growth of a diamond film
%
% This example computes the growth rate of a diamond film according to a
% simplified version of a particular published growth mechanism (see file
% diamond.yaml for details). Only the surface coverage equations are solved
% here; the gas composition is fixed. (For an example of coupled gas-phase
% and surface, see catalytic_combustion.py.) Atomic hydrogen plays an
% important role in diamond CVD, and this example computes the growth rate
% and surface coverages as a function of [H] at the surface for
% fixed temperature and [CH3].
%
% Requires: cantera >= 2.6.0, pandas >= 0.25.0, matplotlib >= 2.0
% Keywords: surface chemistry, kinetics

%% Initialization

clear all
close all

help diamond_cvd
t0 = cputime; % record the starting time

%% Operation Parameters

t = 1200.0; % surface temperature
p = 20.0 * OneAtm / 760.0; % pressure

%% Create the gas object
%
% The gas phase will be taken from the definition of phase 'gas' in
% input file 'diamond.yaml'.

gas = Solution('diamond.yaml', 'gas');
gas.TP = {t, p};
x = gas.X;
ih = gas.speciesIndex('H');
xh0 = x(ih);

%% Create the bulk diamond object
%
% The bulk phase will be taken from the definition of phase 'diamond' in
% input file 'diamond.yaml'.

dbulk = Solution('diamond.yaml', 'diamond');
mw = dbulk.molecularWeights;

%% Create the interface object
%
% This object will be used to evaluate all surface chemical production
% rates. It will be created from the interface definition 'diamond_100'
% in input file 'diamond.yaml'.

surf_phase = Interface('diamond.yaml', 'diamond_100', gas, dbulk);

%% Advance Coverages

iC = surf_phase.kineticsSpeciesIndex('C(d)', 'diamond');

xx = [];
rr = [];
cov = [];

for i = 1:20
    x(ih) = x(ih) / 1.4;
    gas.TPX = {t, p, x};
    surf_phase.advanceCoverages(10.0);
    r = surf_phase.netProdRates;
    carbon_dot = r(iC);
    mdot = mw * carbon_dot;
    rate = mdot / dbulk.D;
    xx = [xx; x(ih)];
    rr = [rr; rate * 1.0e6 * 3600.0];
    cov = [cov; surf_phase.coverages];
end

%% Make Plots

clf;

subplot(1, 2, 1);
plot(xx, rr);
xlabel('H Mole Fraction');
ylabel('Growth Rate (microns/hr)');
title('Growth Rate');

subplot(1, 2, 2);
plot(xx, cov);
xlabel('H Mole Fraction');
ylabel('Coverage');
title('Coverages');
