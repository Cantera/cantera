# Library-level utility functions and physical constants.

"One standard atmosphere [Pa]."
const one_atm = 101325.0

"Universal gas constant [J/kmol/K] (Cantera's molar convention)."
const gas_constant = 8314.46261815324

"Avogadro's number [1/kmol]."
const avogadro = 6.02214076e26

"""
    cantera_version() -> String

Version string of the linked Cantera library.
"""
cantera_version() = get_string((n, b) -> LibCantera.ct_version(n, b))

"""
    git_commit() -> String

Git commit of the linked Cantera library.
"""
git_commit() = get_string((n, b) -> LibCantera.ct_gitCommit(n, b))

"""
    add_data_directory(dir)

Add `dir` to the list of directories searched for Cantera input files.
"""
function add_data_directory(dir::AbstractString)
    check(LibCantera.ct_addDataDirectory(dir))
    return nothing
end

"""
    suppress_thermo_warnings(flag=true)

Suppress (or re-enable) Cantera's thermodynamic warnings.
"""
function suppress_thermo_warnings(flag::Bool=true)
    check(LibCantera.ct_suppressThermoWarnings(Int32(flag)))
    return nothing
end

"Reset all Cantera CLib storage, invalidating every existing handle."
reset_storage() = (check(LibCantera.ct_resetStorage()); nothing)
