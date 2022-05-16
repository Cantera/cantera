% This example file calculates the cell voltage of a lithium-ion battery
% at given temperature, pressure, current, and range of state of charge (SOC).
%
% The thermodynamics are based on a graphite anode and a LiCoO2 cathode,
% modeled using the 'BinarySolutionTabulatedThermo' class.
% Further required cell parameters are the electrolyte ionic resistance, the
% stoichiometry ranges of the active materials (electrode balancing), and the
% surface area of the active materials.
%
% The functionality of this example is presented in greater detail in the
% reference (which also describes the derivation of the
% BinarySolutionTabulatedThermo class):
%
% Reference:
% M. Mayur, S. C. DeCaluwe, B. L. Kee, W. G. Bessler, “Modeling and simulation
% of the thermodynamics of lithium-ion battery intercalation materials in the
% open-source software Cantera,” Electrochim. Acta 323, 134797 (2019),
% https://doi.org/10.1016/j.electacta.2019.134797

% -----------------------------------------------------------------------------
% Input
% -----------------------------------------------------------------------------

clear all
close all
cleanup
clc

% Operation parameters
SOC = 0:0.02:1; % [-] Input state of charge (0...1) (can be a vector)
I_app = -1; % [A] Externally-applied current, negative for discharge
T = 293; % [K] Temperature
P = oneatm; % [Pa] Pressure

% Cell properties
inputFile = 'lithium_ion_battery.yaml'; % Cantera input file name
R_elyt = 0.0384; % [Ohm] Electrolyte resistance
S_ca = 1.1167; % [m^2] Cathode total active material surface area
S_an = 0.7824; % [m^2] Anode total active material surface area

% Electrode balancing: The "balancing" of the electrodes relates the chemical
% composition (lithium mole fraction in the active materials) to the macroscopic
% cell-level state of charge.
X_Li_an_0 = 0.01; % [-] anode Li mole fraction at SOC = 0 %
X_Li_an_1 = 0.75; % [-] anode Li mole fraction at SOC = 100 %
X_Li_ca_0 = 0.99; % [-] cathode Li mole fraction at SOC = 0 %
X_Li_ca_1 = 0.49; % [-] cathode Li mole fraction at SOC = 100 %

% -----------------------------------------------------------------------------
% Calculations
% -----------------------------------------------------------------------------

% Calculate mole fractions from SOC
X_Li_an = (X_Li_an_1-X_Li_an_0)*SOC+X_Li_an_0; % anode balancing
X_Li_ca = (X_Li_ca_0-X_Li_ca_1)*(1-SOC)+X_Li_ca_1; % cathode balancing

% Import all Cantera phases
anode = Solution(inputFile, 'anode');
cathode = Solution(inputFile, 'cathode');
elde = Solution(inputFile, 'electron');
elyt = Solution(inputFile, 'electrolyte');
anode_interface = Interface(inputFile, 'edge_anode_electrolyte', anode, elde, elyt);
cathode_interface = Interface(inputFile, 'edge_cathode_electrolyte', cathode, elde, elyt);

% Set the temperatures and pressures of all phases
anode.TP = {T, P};
cathode.TP = {T, P};
elde.TP = {T, P};
elyt.TP = {T, P};
anode_interface.TP = {T, P};
cathode_interface.TP = {T, P};

% Calculate cell voltage, separately for each entry of the input vectors
V_cell = zeros(length(SOC),1);
phi_l_an = 0;
phi_s_ca = 0;
for i = 1:length(SOC)
    % Set anode electrode potential to 0
    phi_s_an = 0;

    % Calculate anode electrolyte potential
    phi_l_an = fzero(@(E) anode_curr(phi_s_an, E, X_Li_an(i), anode, elde,...
                                     elyt, anode_interface, S_an) - I_app, phi_l_an);

    % Calculate cathode electrolyte potential
    phi_l_ca = phi_l_an + I_app*R_elyt;

    % Calculate cathode electrode potential
    phi_s_ca = fzero(@(E) cathode_curr(E, phi_l_ca, X_Li_ca(i), ...
                                       cathode, elde, elyt, cathode_interface,...
                                       S_ca) - I_app, phi_s_ca);

    % Calculate cell voltage
    V_cell(i) = phi_s_ca - phi_s_an;
end

% Let's plot the cell voltage, as a function of the state of charge:
figure(1);
plot(SOC*100, V_cell, 'linewidth', 2.5)
ylim([2.5, 4.3])
xlabel('State of charge / %')
ylabel('Cell voltage / V')
set(gca, 'fontsize', 14)


%--------------------------------------------------------------------------
% Helper functions
% -----------------------------------------------------------------------------

% This function returns the Cantera calculated anode current (in A)
function anCurr = anode_curr(phi_s, phi_l, X_Li_an, anode, elde, elyt, anode_interface, S_an)
    % Set the active material mole fraction
    anode.X = ['Li[anode]:' num2str(X_Li_an) ', V[anode]:' num2str(1-X_Li_an)];

    % Set the electrode and electrolyte potential
    elde.setElectricPotential(phi_s);
    elyt.setElectricPotential(phi_l);

    % Get the net reaction rate at the anode-side interface
    % Reaction according to cti file: Li+[elyt] + V[anode] + electron <=> Li[anode]
    r = anode_interface.ropNet; % [kmol/m2/s]

    % Calculate the current. Should be negative for cell discharge.
    anCurr = r*faradayconstant*S_an; %
end

% This function returns the Cantera calculated cathode current (in A)
function caCurr = cathode_curr(phi_s, phi_l, X_Li_ca, cathode, elde, elyt, cathode_interface, S_ca)
    % Set the active material mole fractions
    cathode.X = ['Li[cathode]:' num2str(X_Li_ca) ', V[cathode]:' num2str(1-X_Li_ca)];

    % Set the electrode and electrolyte potential
    elde.setElectricPotential(phi_s);
    elyt.setElectricPotential(phi_l);

    % Get the net reaction rate at the cathode-side interface
    % Reaction according to cti file: Li+[elyt] + V[cathode] + electron <=> Li[cathode]
    r = cathode_interface.ropNet; % [kmol/m2/s]

    % Calculate the current. Should be negative for cell discharge.
    caCurr = r*faradayconstant*S_ca*(-1); %
end
