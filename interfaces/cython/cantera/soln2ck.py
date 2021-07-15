from __future__ import print_function
from __future__ import division

import os
from string import Template

import cantera as ct
from cantera.test.test_mechanism_from_solution import test
import ck2cti


def write(solution):
    """Function to write cantera solution object to inp file.

    :param solution:
        Cantera solution object

    :return:
        Name of trimmed Mechanism file (.inp)

    >>> soln2ck.write(gas)
    """
    trimmed_solution = solution
    input_file_name_stripped = trimmed_solution.name
    cwd = os.getcwd()
    output_file_name = os.path.join(cwd,
                                    'pym_' +
                                    input_file_name_stripped +
                                    '.inp')
    with open(output_file_name, 'w+') as f:

        #Work functions

        calories_constant = 4184.0 #number of calories in 1000 Joules of energy
        def eliminate(input_string, char_to_replace, spaces='single'):
            """
            Eliminate characters from a string

            :param input_string
                string to be modified
            :param char_to_replace
                array of character strings to be removed
            """
            for char in char_to_replace:
                input_string = input_string.replace(char, "")
            if spaces == 'double':
                input_string = input_string.replace(" ", "  ")
            return input_string


        def replace_multiple(input_string, replace_list):
            """
            Replace multiple characters in a string

            :param input_string
                string to be modified
            :param replace list
                list containing items to be replaced (value replaces key)
            """
            for original_character, new_character in replace_list.items():
                input_string = input_string.replace(original_character,
                                                    new_character)
            return input_string

        def build_arrhenius(equation_object, equation_type):
            """
            Builds Arrhenius coefficient string

            :param equation_objects
                cantera equation object
            :param equation_type:
                string of equation type
            """
            coeff_sum = sum(equation_object.reactants.values())
            pre_exponential_factor = equation_object.rate.pre_exponential_factor
            temperature_exponent = '{:.3f}'.format(equation_object.rate.temperature_exponent)
            activation_energy = '{:.2f}'.format(equation_object.rate.activation_energy / calories_constant)
            if equation_type == 'ElementaryReaction':
                if coeff_sum == 1:
                    pre_exponential_factor = str(
                                    '{:.3E}'.format(pre_exponential_factor))
                if coeff_sum == 2:
                    pre_exponential_factor = str(
                                    '{:.3E}'.format(pre_exponential_factor*10**3))
                if coeff_sum == 3:
                    pre_exponential_factor = str(
                                    '{:.3E}'.format(pre_exponential_factor*10**6))
            if equation_type == 'ThreeBodyReaction':
                if coeff_sum == 1:
                    pre_exponential_factor = str(
                                    '{:.3E}'.format(pre_exponential_factor*10**3))
                if coeff_sum == 2:
                    pre_exponential_factor = str(
                                    '{:.3E}'.format(pre_exponential_factor*10**6))
            if (equation_type != 'ElementaryReaction'
                and equation_type != 'ThreeBodyReaction'):
                pre_exponential_factor = str(
                                    '{:.3E}'.format(pre_exponential_factor))
            arrhenius = [pre_exponential_factor,
                    temperature_exponent,
                    activation_energy]
            return arrhenius

        def build_modified_arrhenius(equation_object, t_range):
            """
            Builds Arrhenius coefficient strings for high and low temperature ranges

            :param equation_objects
                cantera equation object
            :param t_range:
                simple string ('high' or 'low') to designate temperature range
            """
            if t_range == 'high':
                pre_exponential_factor = equation_object.high_rate.pre_exponential_factor
                temperature_exponent = '{:.3f}'.format(equation_object.high_rate.temperature_exponent)
                activation_energy = '{:.2f}'.format(equation_object.high_rate.activation_energy/calories_constant)
                if len(equation_object.products) == 1:
                    pre_exponential_factor = str(
                                    '{:.5E}'.format(pre_exponential_factor*10**3))
                else:
                    pre_exponential_factor = str(
                                    '{:.5E}'.format(pre_exponential_factor))
                arrhenius_high = [pre_exponential_factor,
                            temperature_exponent,
                            activation_energy]
                return arrhenius_high
            if t_range == 'low':

                pre_exponential_factor = equation_object.low_rate.pre_exponential_factor
                temperature_exponent = '{:.3f}'.format(equation_object.low_rate.temperature_exponent)
                activation_energy = '{:.2f}'.format(equation_object.low_rate.activation_energy/calories_constant)
                if len(equation_object.products) == 1:
                    pre_exponential_factor = str(
                                    '{:.5E}'.format(pre_exponential_factor*10**6))
                else:
                    pre_exponential_factor = str(
                                    '{:.5E}'.format(pre_exponential_factor*10**3))

                arrhenius_low = [pre_exponential_factor,
                            temperature_exponent,
                            activation_energy]
                return arrhenius_low


        def build_nasa(nasa_coeffs, row):
            """
            Creates string of nasa polynomial coefficients

            :param nasa_coeffs
                cantera species thermo coefficients object
            :param row
                which row to write coefficients in
            """
            line_coeffs = ''
            lines = [[1, 2, 3, 4, 5], [6, 7, 8, 9, 10], [11, 12, 13, 14]]
            line_index = lines[row-2]
            for ix, c in enumerate(nasa_coeffs):
                if ix in line_index:
                    if c >= 0:
                        line_coeffs += ' '
                    line_coeffs += str('{:.8e}'.format(c))
            return line_coeffs

        def build_species_string():
            """
            formats species string for writing
            """
            species_list_string = ''
            line = 1
            for sp_index, sp_string in enumerate(trimmed_solution.species_names):
                sp = ' '
                #get length of string next species is added
                length_new = len(sp_string)
                length_string = len(species_list_string)
                total = length_new + length_string + 3
                #if string will go over width, wrap to new line
                if total >= 70*line:
                    species_list_string += '\n'
                    line += 1
                species_list_string += sp_string + ((16-len(sp_string))*sp)
            return species_list_string

        title = ''
        section_break = ('!'+ "-"*75 + '\n'
                         '!  ' + title +'\n'
                         '!'+ "-"*75 + '\n')

        #Write title block to file
        title = 'Chemkin File converted from Solution Object by pyMARS'
        f.write(section_break)

        #Write phase definition to file
        element_names = eliminate(str(trimmed_solution.element_names),
                                            ['[', ']', '\'', ','])
        element_string = Template(
                                'ELEMENTS\n' +
                                '$element_names\n' +
                                'END\n')
        f.write(element_string.substitute(element_names=element_names))
        species_names = build_species_string()
        species_string = Template(
                                'SPECIES\n' +
                                '$species_names\n'+
                                'END\n')
        f.write(species_string.substitute(species_names=species_names))

        #Write species to file
        title = 'Species data'
        f.write(section_break)

        f.write('THERMO ALL' + '\n' +
                '   300.000  1000.000  5000.000' +'\n')
        phase_unknown_list = []

        #write data for each species in the Solution object
        for sp_index in xrange(len(trimmed_solution.species_names)):
            d = 3.33564e-30 #1 debye = d coulomb-meters
            species = trimmed_solution.species(sp_index)
            name = str(trimmed_solution.species(sp_index).name)
            nasa_coeffs = trimmed_solution.species(sp_index).thermo.coeffs
            #Species attributes from trimmed solution object
            t_low = '{0:.3f}'.format(species.thermo.min_temp)
            t_max = '{0:.3f}'.format(species.thermo.max_temp)
            t_mid = '{0:.3f}'.format(species.thermo.coeffs[0])
            temp_range = str(t_low) + '  ' + str(t_max) + '  ' + t_mid
            species_comp = ''
            for atom in species.composition:
                species_comp += '{:<4}'.format(atom)
                species_comp += str(int(species.composition[atom]))
            if type(species.transport).__name__ == 'GasTransportData':
                species_phase = 'G'
            else:
                phase_unknown_list.append(name)
                species_phase = 'G'
            line_1 = (
                    '{:<18}'.format(name) +
                    '{:<6}'.format('    ') +
                    '{:<20}'.format(species_comp) +
                    '{:<4}'.format(species_phase) +
                    '{:<31}'.format(temp_range) +
                    '{:<1}'.format('1') +
                    '\n')
            f.write(line_1)
            line_2_coeffs = build_nasa(nasa_coeffs, 2)
            line_2 = line_2_coeffs  + '    2\n'
            f.write(line_2)
            line_3_coeffs = build_nasa(nasa_coeffs, 3)
            line_3 = line_3_coeffs + '    3\n'
            f.write(line_3)
            line_4_coeffs = build_nasa(nasa_coeffs, 4)
            line_4 = line_4_coeffs + '                   4\n'
            f.write(line_4)

        f.write('END\n')

        #Write reactions to file
        title = 'Reaction Data'
        f.write(section_break)
        f.write('REACTIONS\n')
        #write data for each reaction in the Solution Object
        for reac_index in xrange(len(trimmed_solution.reaction_equations())):
            equation_string = str(trimmed_solution.reaction_equation(reac_index))
            equation_string = eliminate(equation_string, ' ', 'single')
            equation_object = trimmed_solution.reaction(reac_index)
            equation_type = type(equation_object).__name__
            m = str(reac_index+1)
            if equation_type == 'ThreeBodyReaction':
                arrhenius = build_arrhenius(equation_object, equation_type)
                main_line = (
                            '{:<51}'.format(equation_string) +
                            '{:>9}'.format(arrhenius[0]) +
                            '{:>9}'.format(arrhenius[1]) +
                            '{:>11}'.format(arrhenius[2]) +
                            '\n')
                f.write(main_line)
                #trimms efficiencies list
                efficiencies = equation_object.efficiencies
                trimmed_efficiencies = equation_object.efficiencies
                for s in efficiencies:
                    if s not in trimmed_solution.species_names:
                        del trimmed_efficiencies[s]
                replace_list_2 = {
                                    '{':'',
                                    '}':'/',
                                    '\'':'',
                                    ':':'/',
                                    ',':'/'}
                efficiencies_string = replace_multiple(
                                            str(trimmed_efficiencies),
                                                replace_list_2)
                secondary_line = str(efficiencies_string) + '\n'
                if bool(efficiencies) is True:
                    f.write(secondary_line)
            if equation_type == 'ElementaryReaction':
                arrhenius = build_arrhenius(equation_object, equation_type)
                main_line = (
                            '{:<51}'.format(equation_string) +
                            '{:>9}'.format(arrhenius[0]) +
                            '{:>9}'.format(arrhenius[1]) +
                            '{:>11}'.format(arrhenius[2]) +
                            '\n')
                f.write(main_line)
            if equation_type == 'FalloffReaction':
                arr_high = build_modified_arrhenius(equation_object, 'high')
                main_line = (
                            '{:<51}'.format(equation_string) +
                            '{:>9}'.format(arr_high[0]) +
                            '{:>9}'.format(arr_high[1]) +
                            '{:>11}'.format(arr_high[2]) +
                            '\n')
                f.write(main_line)
                arr_low = build_modified_arrhenius(equation_object, 'low')
                second_line = (
                                '     LOW  /' +
                                '  ' + arr_low[0] +
                                '  ' + arr_low[1] +
                                '  ' + arr_low[2] + '/\n')
                f.write(second_line)
                j = equation_object.falloff.parameters
                #If optional Arrhenius data included:
                try:
                    third_line = (
                                '     TROE/' +
                                '   ' + str(j[0]) +
                                '  ' + str(j[1]) +
                                '  ' + str(j[2]) +
                                '  ' + str(j[3]) +' /\n')
                    f.write(third_line)
                except IndexError:
                    pass
                #trimms efficiencies list
                efficiencies = equation_object.efficiencies
                trimmed_efficiencies = equation_object.efficiencies
                for s in efficiencies:
                    if s not in trimmed_solution.species_names:
                        del trimmed_efficiencies[s]
                replace_list_2 = {
                                    '{':'',
                                    '}':'/',
                                    '\'':'',
                                    ':':'/',
                                    ',':'/'}
                efficiencies_string = replace_multiple(
                                            str(trimmed_efficiencies),
                                                replace_list_2)

                fourth_line = str(efficiencies_string) + '\n'
                if bool(efficiencies) is True:
                    f.write(fourth_line)
            #dupluicate option
            if equation_object.duplicate is True:
                duplicate_line = ' DUPLICATE' +'\n'
                f.write(duplicate_line)
        f.write('END')
    return output_file_name
