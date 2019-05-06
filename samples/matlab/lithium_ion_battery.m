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
%
% The routine below returns the cell voltage (in Volt) of a lithium-ion
% cell for a given cell current and active material lithium stoichiometries.
%
% Input:
% - stoichiometries X_Li_ca and X_Li_an [-] (can be vectors)
% - temperature T [K]
% - pressure P [Pa]
% - externally-applied current I_app [A]
% - electrolyte resistance R_elyt [Ohm]


% Input parameters
SOC = 0:0.02:1; % [-] Input state of charge (0...1)
X_Li_an = (0.75-0.01)*SOC+0.01; % anode balancing
X_Li_ca = (0.99-0.49)*(1-SOC)+0.49; % cathode balancing
I_app = 0; % [A] Externally-applied current
R_elyt = 0; % [Ohm] Electrolyte resistance
T = 300; % [K] Temperature
P = oneatm; % [Pa] Pressure
inputCTI = 'lithium_ion_battery.cti'; % cantera input file name
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
set(anode,'T',T,'P',P);
set(cathode,'T',T,'P',P);
set(elde,'T',T,'P',P);
set(elyt,'T',T,'P',P);
set(anode_interface,'T',T,'P',P);
set(cathode_interface,'T',T,'P',P);

% Calculate cell voltage, separately for each entry of the input vectors
E_cell = zeros(length(SOC),1);
for i = 1:length(SOC)
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

% Let's plot the cell voltage, as a function of the state of charge:
figure(1);
plot(SOC*100,E_cell,'linewidth',2.5)
ylim([2.5,4.3])
xlabel('State of charge / %')
ylabel('Cell voltage / V')
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
    anCurr = r*96485*S_an*1; % F = 96485 C/mol Faraday's constant
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
    caCurr = r*96485*S_ca*(-1); % F = 96485 C/mol Faraday's constant
end
