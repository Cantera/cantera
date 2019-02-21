% This example file calculates the open-circuit voltage for a lithium-ion
% battery over a range of compositions.
%
% The thermodynamics are based on a graphite anode and a LiCoO2 cathode,
% modeled using the 'BinarySolutionTabulatedThermo' class.
%
% Note that the function 'E_cell' below has even greater capabilities than
% what we use, here. It calculates the steady state cell voltage, at a
% given composition and cell current, for a given electrolyte ionic
% resistance.  This functionality is presented in greater detail in the
% reference (which also describes the derivation of the
% BinarySolutionTabulatedThermo class):
%
% Reference:
% M. Mayur, S. DeCaluwe, B. L. Kee, W. G. Bessler, "Modeling
% thermodynamics and kinetics of intercalation phases for lithium-ion
% batteries in Cantera", under review at Electrochimica Acta.

% For the sake of simplicity, we're going to assume that the anode and
% cathode capacities are perfectly balanced (i.e. if the cathode lithium
% content is X percent of it's max possible (i.e. its capacity), then we
% will assume that the anode is at 1-X percent.  Without loss of
% generality, we will define the anode composition:

% The routine below returns the cell voltage (in Volt) of a lithium-ion
% cell for a given cell current and active material lithium stoichiometries.
%
% Input:
% - stoichiometries X_Li_ca and X_Li_an [-] (can be vectors)
% - temperature T [K]
% - pressure P [Pa]
% - externally-applied current I_app [A]
% - electrolyte resistance R_elyt [Ohm]

X_Li_an = [0.005:0.025:0.995];
X_Li_ca = 1 - X_Li_an;

I_app = 0;
R_elyt = 0;
T = 300;
P = oneatm;

global F
% Parameters
inputCTI = 'lithium_ion_battery.cti'; % cantera input file name
F = 96485; % Faraday's constant [C/mol]
S_ca = 1.1167; % [m^2] Cathode total active material surface area
S_an = 0.7824; % [m^2] Anode total active material surface area

% Import all Cantera phases
anode = Solution(inputCTI, 'anode');
cathode = Solution(inputCTI, 'cathode');
elde = Solution(inputCTI, 'electron');
elyt = Solution(inputCTI, 'electrolyte');
anode_interface = Interface(inputCTI, 'edge_anode_electrolyte', anode, elde, elyt);
cathode_interface = Interface(inputCTI, 'edge_cathode_electrolyte', cathode, elde, elyt);

% Set the temperatures and pressures of all phases
phases = [anode elde elyt cathode];
for ph = phases
    set(ph,'T',T,'P',P);
end

% Calculate cell voltage, separately for each entry of the input vectors
E_cell = zeros(length(X_Li_ca),1);
for i = 1:length(X_Li_ca)
    % Set anode electrode potential to 0
    phi_s_an = 0;

    % Calculate anode electrolyte potential
    phi_l_an = fzero(@(E) anode_curr(phi_s_an,E,X_Li_an(i),anode,elde,elyt,anode_interface,S_an)+I_app, 0);

    % Calculate cathode electrolyte potential
    phi_l_ca = phi_l_an + I_app*R_elyt;

    % Calculate cathode electrode potential
    phi_s_ca = fzero(@(E) cathode_curr(E,phi_l_ca,X_Li_ca(i),cathode,elde,elyt,cathode_interface,S_ca)+I_app, 0);

    % Calculate cell voltage
    E_cell(i) = phi_s_ca - phi_s_an;
end

% Let's plot the cell voltage, as a function of the cathode stoichiometry:
plot(X_Li_ca,E_cell,'linewidth',2.5)
ylim([2.5,4.3])
xlabel('Li Fraction in Cathode')
ylabel('Cell potential [V]')
set(gca,'fontsize',14)


%--------------------------------------------------------------------------
% Helper functions

% This function returns the Cantera calculated anode current (in A)
function anCurr = anode_curr(phi_s,phi_l,X_Li_an,anode,elde,elyt,anode_interface,S_an)

    global F

    % Set the active material mole fraction
    set(anode,'X',['Li[anode]:' num2str(X_Li_an) ', V[anode]:' num2str(1-X_Li_an)]);

    % Set the electrode and electrolyte potential
    setElectricPotential(elde,phi_s);
    setElectricPotential(elyt,phi_l);

    % Get the net reaction rate at the cathode-side interface
    r = rop_net(anode_interface).*1e3; % [mol/m2/s]

    % Calculate the current
    anCurr = r*F*S_an*1;
end

% This function returns the Cantera calculated cathode current (in A)
function caCurr = cathode_curr(phi_s,phi_l,X_Li_ca,cathode,elde,elyt,cathode_interface,S_ca)

    global F

    % Set the active material mole fractions
    set(cathode,'X',['Li[cathode]:' num2str(X_Li_ca) ', V[cathode]:' num2str(1-X_Li_ca)]);

    % Set the electrode and electrolyte potential
    setElectricPotential(elde,phi_s);
    setElectricPotential(elyt,phi_l);

    % Get the net reaction rate at the cathode-side interface
    r = rop_net(cathode_interface).*1e3; % [mol/m2/s]

    % Calculate the current
    caCurr = r*F*S_ca*(-1);
end
