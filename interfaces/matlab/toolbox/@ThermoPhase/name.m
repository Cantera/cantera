function nm = name(tp)
% NAME - user-specified phase name.

nm = phase_get(thermo_hndl(tp), 42);
