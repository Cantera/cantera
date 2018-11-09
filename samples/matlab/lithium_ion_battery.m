function E_cell = lithium_ion_battery(X_Li_ca, X_Li_an, T, P, I_app, R_elyt)
% Returns the cell voltage (in Volt) of a lithium-ion cell for a given cell
% current and active material lithium stoichiometries.
%
% Input:
% - stoichiometries X_Li_ca and X_Li_an [-] (can be vectors)
% - temperature T [K]
% - pressure P [Pa]
% - externally-applied current I_app [A]
% - electrolyte resistance R_elyt [Ohm]
%
% Reference:
% M. Mayur, S. DeCaluwe, B. L. Kee, W. G. Bessler, "Modeling
% thermodynamics and kinetics of intercalation phases for lithium-ion
% batteries in Cantera", Computer Physics Communications


% Parameteres
inputCTI = 'lithium_ion_battery.cti'; % cantera input file name
F = 96485; % Faraday's constant [C/mol]
S_ca = 1.1167; % [m^2] Cathode total active material surface area
S_an = 0.7824; % [m^2] Anode total active material surface area

% Import all Cantera phases
anode = importThermoPhase(inputCTI, 'anode');
cathode = importThermoPhase(inputCTI,'cathode');
elde = importThermoPhase(inputCTI,'electron');
elyt = importThermoPhase(inputCTI,'electrolyte');
anode_interface = importEdge(inputCTI, 'edge_anode_electrolyte', anode, elde, elyt);
cathode_interface = importEdge(inputCTI, 'edge_cathode_electrolyte', cathode, elde, elyt);

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
    phi_l_an = fzero(@(E) anode_curr(phi_s_an,E,X_Li_an(i))+I_app, 0);

    % Calculate cathode electrolyte potential
    phi_l_ca = phi_l_an + I_app*R_elyt;

    % Calculate cathode electrode potential
    phi_s_ca = fzero(@(E) cathode_curr(E,phi_l_ca,X_Li_ca(i))+I_app, 0);

    % Calculate cell voltage
    E_cell(i) = phi_s_ca - phi_s_an;
end


%--------------------------------------------------------------------------
% Sub-functions

% This function returns the ThermoPhase class instance from CTI file
    function phase = importThermoPhase(inputCTI, name)
        doc = XML_Node('doc', inputCTI);
        node = findByID(doc, name);
        phase = ThermoPhase(node);
    end

% This function returns the Cantera calculated anode current (in A)
    function anCurr = anode_curr(phi_s,phi_l,X_Li_an)
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
    function caCurr = cathode_curr(phi_s,phi_l,X_Li_ca)
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
end