import cantera as ct

def test(original_solution, new_solution):
    """Test written cti file against original cti file.

    Arguments
    -------------
    Original cti file
    New cti file

    Returns
    -----------
    List of file properties that do not match


    """
    original = original_solution
    new = new_solution

    comparison_list=[]
    for i, species in enumerate(new.species_names):
        if species in original.species_names:
            comparison_list.append(species)



    def test_species_n():
        n_species_new=len(new.species())
        n_species_original=len(original.species())
        if n_species_new != n_species_original:
            print 'Not all species accounted for\n'
            print 'Number of species in original file: %s'  %n_species_original
            print 'Number of species in new file: %s'   %n_species_new





    def test_species_def():
        num =0
        #make sure comparing same species
        for i, name1 in enumerate(new.species_names):

            #start comparison with same species
            new_species= new.species(i)
            for j, name2 in enumerate(original.species_names):
                if name1.upper().lower() == name2.upper().lower():
                    original_species=original.species(j)
                    num += 1
            if original_species.name.upper().lower() != new_species.name.upper().lower():
                print (j, original_species, i, new_species,)


            assert original_species.composition == new_species.composition

            assert original_species.thermo.coeffs.all() == new_species.thermo.coeffs.all()

            try:
                assert original_species.transport.geometry == new_species.transport.geometry
            except AttributeError:
                pass

            try:
                assert original_species.transport.diameter == new_species.transport.diameter
            except AttributeError:
                pass

            try:
                assert original_species.transport.well_depth == new_species.transport.well_depth
            except AttributeError:
                pass

            try:
                assert original_species.transport.polarizability == new_species.transport.polarizability
            except AttributeError:
                pass
            try:
                assert original_species.transport.rotational_relaxation == new_species.transport.rotational_relaxation
            except AttributeError:
                pass

            try:
                assert original_species.transport.dipole == new_species.transport.dipole
            except AttributeError:
                pass

        print ('\nSpecies definition tests finished \n\n')

    def test_reactions_n():
        n_reactions_new=len(new.reactions())
        n_reactions_original=len(original.reactions())
        if n_reactions_new != n_reactions_original:
            print 'Not all reactions accounted for\n'
            print 'Number of reactions in original file: %s'  %n_reactions_original
            print 'Number of reactions in new file: %s'   %n_reactions_new


    def test_reactions():
        c=4184.0
        num = 0
        #iterate through all reactions in new solution
        for k, name1 in enumerate(new.reaction_equations()):
            num += 1
            new_reaction=new.reaction(k)
            new_eq_str=new_reaction.equation
            new_equation_type = type(new_reaction).__name__
            new_eq_ID=new_reaction.ID

            #handle different reaction types
            if new_equation_type == 'ElementaryReaction' or \
                new_equation_type == 'ThreeBodyReaction':

                new_A='{:.1E}'.format(new_reaction.rate.pre_exponential_factor)
                new_b='{:.1E}'.format(new_reaction.rate.temperature_exponent)
                new_E='{:.1E}'.format(new_reaction.rate.activation_energy)

            if new_equation_type == 'FalloffReaction':
                new_A_low='{:.1E}'.format(new_reaction.low_rate.pre_exponential_factor)
                new_b_low='{:.1E}'.format(new_reaction.low_rate.temperature_exponent)
                new_E_low='{:.1E}'.format(new_reaction.low_rate.activation_energy)

                new_A_high='{:.1E}'.format(new_reaction.high_rate.pre_exponential_factor)
                new_b_high='{:.1E}'.format(new_reaction.high_rate.temperature_exponent)
                new_E_high='{:.1E}'.format(new_reaction.high_rate.activation_energy)

            #iterate through original reactions
            for l, name2 in enumerate(original.reaction_equations()):
                rx_n=l+1
                original_reaction=original.reaction(l)
                original_eq_str=original_reaction.equation
                original_eq_ID=original_reaction.ID

                #check that same reaction is being compared
                if new_eq_ID == original_eq_ID:
                    assert new_eq_str == original_eq_str.upper()
                    original_equation_type=type(original_reaction).__name__

                    #check that reaction type is the same
                    if original_equation_type != new_equation_type:
                        print 'Reaction type does not match'
                        print 'Reaction # %s : %s\n\n' %(rx_n, new_eq_str)

                    if original_equation_type == 'ElementaryReaction' or \
                        original_equation_type == 'ThreeBodyReaction':

                        original_A='{:.1E}'.format(original_reaction.rate.pre_exponential_factor)
                        original_b='{:.1E}'.format(original_reaction.rate.temperature_exponent)
                        original_E='{:.1E}'.format(original_reaction.rate.activation_energy)

                        if new_A != original_A:
                            print 'Reaction Pre-Exponential does not match'
                            print 'Reaction # %s : %s   || %s ' %(rx_n, new_eq_str, original_eq_str)
                            print 'Original: %s  New: %s\n\n' %(original_A, new_A)
                        if new_b != original_b:
                            print 'Reaction Temperature Exponent does not match'
                            print 'Reaction # %s : %s   || %s' %(rx_n, new_eq_str, original_eq_str)
                            print 'Original: %s  New: %s\n\n' %(original_b, new_b)
                        if new_E != original_E:
                            print 'Reaction Activation Energy does not match'
                            print 'Reaction # %s : %s   || %s' %(rx_n, new_eq_str, original_eq_str)
                            print 'Original: %s  New: %s\n\n' %(original_E, new_E)

                    if original_equation_type == 'FalloffReaction':
                        original_A_low='{:.1E}'.format(original_reaction.low_rate.pre_exponential_factor)
                        original_b_low='{:.1E}'.format(original_reaction.low_rate.temperature_exponent)
                        original_E_low='{:.1E}'.format(original_reaction.low_rate.activation_energy)

                        original_A_high='{:.1E}'.format(original_reaction.high_rate.pre_exponential_factor)
                        original_b_high='{:.1E}'.format(original_reaction.high_rate.temperature_exponent)
                        original_E_high='{:.1E}'.format(original_reaction.high_rate.activation_energy)

                        if new_A_low != original_A_low:
                            print 'Reaction Pre-Exponential Low does not match'
                            print 'Reaction # %s : %s' %(rx_n, new_eq_str)
                            print 'Original: %s  New: %s\n\n' %(original_A_low, new_A_low)
                        if new_b_low != original_b_low:
                            print 'Reaction Temperature Exponent Low does not match'
                            print 'Reaction # %s : %s' %(rx_n, new_eq_str)
                            print 'Original: %s  New: %s\n\n' %(original_b_low, new_b_low)
                        if new_E_low != original_E_low:
                            print 'Reaction Activation Energy Low does not match'
                            print 'Reaction # %s : %s' %(rx_n, new_eq_str)
                            print 'Original: %s  New: %s\n\n' %(original_E_low, new_E_low)
                        if new_A_high != original_A_high:
                            print 'Reaction Pre-Exponential High does not match'
                            print 'Reaction # %s : %s' %(rx_n, new_eq_str)
                            print 'Original: %s  New: %s\n\n' %(original_A_high, new_A_high)
                        if new_b_high != original_b_high:
                            print 'Reaction Temperature Exponent High does not match'
                            print 'Reaction # %s : %s' %(rx_n, new_eq_str)
                            print 'Original: %s  New: %s\n\n' %(original_b_high, new_b_high)
                        if new_E_high != original_E_high:
                            print 'Reaction Activation Energy High does not match'
                            print 'Reaction # %s : %s' %(rx_n, new_eq_str)
                            print 'Original: %s  New: %s\n\n' %(original_E_high, new_E_high)


            #try:
                #if new_reaction.rate.pre_exponential_factor != original_reaction.rate.pre_exponential_factor:
                    #if new_reaction.rate.pre_exponential_factor/ original_reaction.rate.pre_exponential_factor > .004:
                    #print (k, (new_reaction.rate.pre_exponential_factor/ original_reaction.rate.pre_exponential_factor), new_reaction.reaction_type, new_reaction.rate.temperature_exponent, (new_reaction.rate.activation_energy/c) , new_reaction )
            #except AttributeError:
                #pass
            #assert new_reaction.efficiencies == original_reaction.efficiencies
        print ('Equation definition tests finished\n')

    test_species_def()
    test_reactions()
    test_species_n()
    test_reactions_n()

#test('gri301.cti', 'pym_gri30.cti')
