function v = calculatePsat(a, T)
% calcPsat - Calculate saturation pressure for temperature T (in K) for a
%  fluid with the MFTP thermo model.
v = thermo_get(a.tp_id,26,T);
