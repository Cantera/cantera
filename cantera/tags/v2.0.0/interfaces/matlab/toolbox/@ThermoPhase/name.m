function nm = name(self)
% NAME - user-specified phase name.
nm = phase_get(thermo_hndl(self), 42);
