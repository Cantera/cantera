from textwrap import dedent, indent
from string import Template


def get_property_definition_strings():
    """This function is here to be imported by the sdist builder."""
    UNITS = {
        "cp_mass": '"J/kg/K"', "cp_mole": '"J/kmol/K"', "cv_mass": '"J/kg/K"',
        "cv_mole": '"J/kmol/K"', "density_mass": '"kg/m**3"', "density_mole": '"kmol/m**3"',
        "enthalpy_mass": '"J/kg"', "enthalpy_mole": '"J/kmol"', "entropy_mass": '"J/kg/K"',
        "entropy_mole": '"J/kmol/K"', "gibbs_mass": '"J/kg"', "gibbs_mole": '"J/kmol"',
        "int_energy_mass": '"J/kg"', "int_energy_mole": '"J/kmol"',
        "volume_mass": '"m**3/kg"', "volume_mole": '"m**3/kmol"', "T": '"K"', "P": '"Pa"',
        "X": '"dimensionless"', "Y": '"dimensionless"', "Q": '"dimensionless"',
        "cp": '"J/K/" + self.basis_units', "cv": '"J/K/" + self.basis_units',
        "density": 'self.basis_units + "/m**3"', "h": '"J/" + self.basis_units',
        "s": '"J/K/" + self.basis_units', "g": '"J/" + self.basis_units',
        "u": '"J/" + self.basis_units', "v": '"m**3/" + self.basis_units',
        "H": '"J/" + self.basis_units', "V": '"m**3/" + self.basis_units',
        "S": '"J/K/" + self.basis_units', "D": 'self.basis_units + "/m**3"',
        "U": '"J/" + self.basis_units', "P_sat": '"Pa"', "T_sat": '"K"',
        "atomic_weight": '"kg/kmol"', "chemical_potentials": '"J/kmol"',
        "concentrations": '"kmol/m**3"', "critical_pressure": '"Pa"',
        "critical_temperature": '"K"', "critical_density": 'self.basis_units + "/m**3"',
        "electric_potential": '"V"', "electrochemical_potentials": '"J/kmol"',
        "isothermal_compressibility": '"1/Pa"', "sound_speed": '"m/s"', "max_temp": '"K"',
        "mean_molecular_weight": '"kg/kmol"', "min_temp": '"K"',
        "molecular_weights": '"kg/kmol"', "partial_molar_cp": '"J/kmol/K"',
        "partial_molar_enthalpies": '"J/kmol"', "partial_molar_entropies": '"J/kmol/K"',
        "partial_molar_int_energies": '"J/kmol"', "partial_molar_volumes": '"m**3/kmol"',
        "reference_pressure": '"Pa"', "thermal_expansion_coeff": '"1/K"'
    }

    SYMBOL = {
        "T": "T", "P": "P", "D": "density", "H": "h", "S": "s",
        "V": "v", "U": "u", "Q": "Q", "X": "X", "Y": "Y"
    }

    getter_properties = [
        "density_mass", "density_mole", "enthalpy_mass", "enthalpy_mole", "entropy_mass",
        "entropy_mole", "int_energy_mass", "int_energy_mole", "volume_mass", "volume_mole",
        "gibbs_mass", "gibbs_mole", "cp_mass", "cp_mole", "cv_mass", "cv_mole", "P",
        "P_sat", "T", "T_sat", "atomic_weight", "chemical_potentials", "concentrations",
        "critical_pressure", "critical_temperature", "critical_density",
        "electric_potential", "electrochemical_potentials", "isothermal_compressibility",
        "sound_speed", "max_temp", "mean_molecular_weight", "min_temp", "molecular_weights",
        "partial_molar_cp", "partial_molar_enthalpies", "partial_molar_entropies",
        "partial_molar_int_energies", "partial_molar_volumes", "reference_pressure",
        "thermal_expansion_coeff", "cp", "cv", "density", "h", "s", "g", "u", "v",
    ]

    getter_template = Template(dedent("""
        @property
        @copy_doc
        def ${name}(self):
            return Q_(self._phase.${name}, ${units})
    """))

    thermophase_getters = []
    for name in getter_properties:
        thermophase_getters.append(getter_template.substitute(name=name, units=UNITS[name]))

    setter_2_template = Template(dedent("""
        @property
        @copy_doc
        def ${name}(self):
            ${n0}, ${n1} = self._phase.${name}
            return Q_(${n0}, ${u0}), Q_(${n1}, ${u1})

        @${name}.setter
        def ${name}(self, value):
            ${n0} = value[0] if value[0] is not None else self.${n0}
            ${n1} = value[1] if value[1] is not None else self.${n1}
            for val, unit in ((${n0}, ${u0}), (${n1}, ${u1})):
                try:
                    val.ito(unit)
                except AttributeError as e:
                    if "'ito'" in str(e):
                        raise CanteraError(
                            f"Value {val!r} must be an instance of a pint.Quantity class"
                        ) from None
                    else:
                        raise
            self._phase.${name} = ${n0}.magnitude, ${n1}.magnitude
    """))

    tp_setter_2_properties = ["TP", "DP", "HP", "SP", "SV", "TD", "UV"]
    pf_setter_2_properties = ["PQ", "TQ", "PV", "SH", "ST", "TH", "TV", "UP", "VH"]

    thermophase_2_setters = []
    for name in tp_setter_2_properties:
        d = dict(name=name, n0=SYMBOL[name[0]], u0=UNITS[name[0]], n1=SYMBOL[name[1]],
                u1=UNITS[name[1]])
        thermophase_2_setters.append(setter_2_template.substitute(d))

    purefluid_2_setters = []
    for name in pf_setter_2_properties:
        d = dict(name=name, n0=SYMBOL[name[0]], u0=UNITS[name[0]], n1=SYMBOL[name[1]],
                u1=UNITS[name[1]])
        purefluid_2_setters.append(setter_2_template.substitute(d))

    setter_3_template = Template(dedent("""
        @property
        @copy_doc
        def ${name}(self):
            ${n0}, ${n1}, ${n2} = self._phase.${name}
            return Q_(${n0}, ${u0}), Q_(${n1}, ${u1}), Q_(${n2}, ${u2})

        @${name}.setter
        def ${name}(self, value):
            ${n0} = value[0] if value[0] is not None else self.${n0}
            ${n1} = value[1] if value[1] is not None else self.${n1}
            for val, unit in ((${n0}, ${u0}), (${n1}, ${u1})):
                try:
                    val.ito(unit)
                except AttributeError as e:
                    if "'ito'" in str(e):
                        raise CanteraError(
                            f"Value {val!r} must be an instance of a pint.Quantity class"
                        ) from None
                    else:
                        raise
            if value[2] is not None:
                try:
                    ${n2} = value[2].to(${u2}).magnitude
                except AttributeError:
                    ${n2} = value[2]
            else:
                ${n2} = self.${n2}.magnitude
            self._phase.${name} = ${n0}.magnitude, ${n1}.magnitude, ${n2}
    """))

    tp_setter_3_properties = [
        "TPX", "TPY", "DPX", "DPY", "HPX", "HPY", "SPX", "SPY", "SVX", "SVY", "TDX", "TDY",
        "UVX", "UVY"
    ]

    thermophase_3_setters = []
    for name in tp_setter_3_properties:
        d = dict(name=name, n0=SYMBOL[name[0]], u0=UNITS[name[0]], n1=SYMBOL[name[1]],
                u1=UNITS[name[1]], n2=SYMBOL[name[2]], u2=UNITS[name[2]])
        thermophase_3_setters.append(setter_3_template.substitute(d))

    getter_3_template = Template(dedent("""
        @property
        @copy_doc
        def ${name}(self):
            ${n0}, ${n1}, ${n2} = self._phase.${name}
            return Q_(${n0}, ${u0}), Q_(${n1}, ${u1}), Q_(${n2}, ${u2})
    """))

    pf_getter_3_properties = ["DPQ", "HPQ", "SPQ", "SVQ", "TDQ", "UVQ"]

    purefluid_3_getters = []
    for name in pf_getter_3_properties:
        d = dict(name=name, n0=name[0], u0=UNITS[name[0]], n1=name[1], u1=UNITS[name[1]],
                n2=name[2], u2=UNITS[name[2]])
        purefluid_3_getters.append(getter_3_template.substitute(d))


    def recursive_join(*args, joiner=""):
        result = ""
        for arg in args:
            result = result + joiner.join(arg)
        return result


    thermophase_properties = recursive_join(thermophase_getters, thermophase_2_setters,
                                            thermophase_3_setters)
    purefluid_properties = recursive_join(purefluid_2_setters, purefluid_3_getters)

    common_properties = dedent("""
        def __getattr__(self, name):
            return getattr(self._phase, name)

        def __setattr__(self, name, value):
            if name in dir(self):
                object.__setattr__(self, name, value)
            else:
                setattr(self._phase, name, value)

        @copy_doc
        def report(self, *args, **kwargs):
            return self._phase.report(*args, **kwargs)

        def __call__(self, *args, **kwargs):
            print(self.report(*args, **kwargs))

        @property
        def basis_units(self):
            \"\"\"The units associated with the mass/molar basis of this phase.\"\"\"
            if self._phase.basis == "mass":
                return "kg"
            else:
                return "kmol"

        @property
        @copy_doc
        def X(self):
            \"\"\"If an array is used for setting, the units must be dimensionless.\"\"\"
            X = self._phase.X
            return Q_(X, "dimensionless")

        @X.setter
        def X(self, value):
            if value is not None:
                try:
                    X = value.to("dimensionless").magnitude
                except AttributeError:
                    X = value
            else:
                X = self.X.magnitude
            self._phase.X = X

        @property
        @copy_doc
        def Y(self):
            \"\"\"If an array is used for setting, the units must be dimensionless.\"\"\"
            Y = self._phase.Y
            return Q_(Y, "dimensionless")

        @Y.setter
        def Y(self, value):
            if value is not None:
                try:
                    Y = value.to("dimensionless").magnitude
                except AttributeError:
                    Y = value
            else:
                Y = self.Y.magnitude
            self._phase.Y = Y
    """)

    return common_properties, thermophase_properties, purefluid_properties


def UnitsInterfaceBuilder(env, target, source):
    """This builder creates the cantera.with_units interface for the Python package.

    This builder is meant to be called like
    ``env.UnitsInterfaceBuilder(target_file, template_source_file)``

    The builder will create string templates for all the ``ThermoPhase`` methods and
    fill them in with the appropriate SI + kmol units to be converted with pint.

    The return value is an item of type ``env.SubstFile`` which can be put into a
    ``Depends`` call, so that this builder will be called for a particular module.
    """
    common_properties, thermophase_properties, purefluid_properties = get_property_definition_strings()
    env["common_properties"] = indent(common_properties.strip(), " "*4)
    env["thermophase_properties"] = indent(thermophase_properties.strip(), " "*4)
    env["purefluid_properties"] = indent(purefluid_properties.strip(), " "*4)
    units = env.SubstFile(target, source)
    return units


def generate(env, **kw):
    env.AddMethod(UnitsInterfaceBuilder)


def exists(env):
    return True
