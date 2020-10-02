"""
PyVRNA: A custom Python 2.7 wrapper for ViennaRNA
Authors: Ayaan Hossain (ain.hoss07@gmail.com), Alex Reis
Trivia: This wrapper was engineered by Ayaan with lots of love and care during his
        third rotation in the BG program at Salis Lab between Feb-Mar 2017.
        ^Oh please we're still working on this wrapper.
"""

# Required for PyVRNA
import RNA
import math
from   itertools   import izip,  imap, chain
from   collections import deque, namedtuple, OrderedDict
from   operator    import itemgetter

# Required for ViennaRNA
import os
import re
import random
import string

"""
NOTES:

Allowed characters in constraint strings:
.    (no pressure to base pair)
x    (must not base pair) 
|    (must base pair, in either direction)
<, > (must base pair, in < or > direction)
(, ) (base marked ( must pair with its ) base)
i    (must pair with base in same molecule, intramolecular pairing)
e    (must pair with base in another molecule, intermolecular pairing)
+    (marked bases form a g-quadruplex)

Adding aptamers:
aptamers            = list of aptamer sequences in sequence (str)
aptamer_constraints = list of constraints (strs)
dG_ligands          = list of ligand binding free energies (double)

For a given aptamer, the sequence and constraint can be separated by "&" if
they are not a part of a contiguous sequence.
"""


class PyVRNA(object):
    """
    This class abstracts for RNAcentroid, RNAcofold, RNAeval, RNAfold and RNAsubopt programs in ViennaRNA suite, and is designed
    to be used by the various programs developed in Salis Lab.

    NOTE:  None of the (non-)sequence information is stored in PyVRNA object, because the parameters are extrinsic to the energy model
           that is abstracted by the PyVRNA object. All functions in this object operate on supplied sequence(s) and other parameters.
           The philosophy behind this design is to build a model once, at the beginning, and have it operate on all tge sequences and 
           other parameters supplied in future, to retain consistency throughout. Also, testing inputs add overhead, and is not always
           necessary, so it must be turned on only when operating on artifically generated sequences or parameters that need testing.

    Usage: object_name = PyVRNA(temperature=float, dangles=int, gquad=bool, parameter_file=string, test_inputs=bool)
           object_name.function_name(parameters)
    """
    def __init__(self,
        temperature=37.0,
        dangles=2,
        noGU=False,
        gquad=False,
        parameter_file="rna_andronescu2007.par",
        duplex_adjustment=True,
        test_inputs=False,
        enforce_constraints=False,
        pyindex=False):
        """
        Initializes an PyVRNA object, after validating all parameters.

        Usage: energy_model = PyVRNA(temperature=int, dangles=int, gquad=bool, parameter_file=string, test_inputs=bool)
        """
        # Setup variables for helping in assertions
        self.parameter_files     = ["dna_mathews1999.par",
                                    "rna_turner1999.par",
                                    "dna_mathews2004.par",
                                    "rna_turner2004.par",
                                    "rna_andronescu2007.par"]
        self.parameter_directory =  "/usr/local/share/ViennaRNA/"
        
        # Assert all supplied variables for correctness
        assert 0 < temperature,                                'temperature must be greater than 0.'
        assert 0 <= dangles <= 3 and isinstance(dangles, int), 'dangles must be one of the three integers: 0 (none), 1 (some), 2 (all), or 3 (all + coaxial stacking of multi-branch loops)'
        assert isinstance(noGU, bool),                         'noGU must be a boolean value.'
        assert isinstance(gquad, bool),                        'gquad must be a boolean value.'
        assert parameter_file in self.parameter_files,         ''.join(['parameter file must be ', " or ".join(self.parameter_files), '.'])
        if parameter_file + ".par" in self.parameter_files:
            parameter_file += ".par"
        assert isinstance(duplex_adjustment, bool),            'duplex_adjustment must be a boolean value.'
        assert isinstance(test_inputs, bool),                  'test_inputs must be a boolean value.'        
        assert isinstance(enforce_constraints, bool),          'enforce_constraints must be a boolean value.'
        assert isinstance(pyindex, bool),                      'pyindex must be a boolean value.'
        
        # Setup model default settings and variables
        self.settings             = RNA.md()
        self.settings.temperature = temperature   # Temperature in Celsius; default=37.0 (float)
        self.settings.dangles     = dangles       # Dangling end energies (0,1,2); see RNAlib documentation; default=2 (int)
        self.settings.gquad       = gquad         # Incorporate G-Quadruplex formation into structure prediction; default=off/0/False (bool/int)
        self.settings.noGU        = noGU          # Toggles GU wobble for RNA folding processes ... probably need to test this
        RNA.cvar.temperature      = temperature   # Global setting of temperature
        RNA.cvar.dangles          = dangles       # Global setting of dangles
        RNA.cvar.gquad            = gquad         # Global setting of gquad
        RNA.cvar.noGU             = noGU         # Global setting of gquad
        
        '''
        Other settings variables available for manual overriding
        self.settings.betaScale       = 1.0       # Set the scaling of the Boltzmann factors; default=1.0 (float)
        self.settings.special_hp      = True      # Include tabulated free energies for special hairpin loops (Tri-, Tetra-, or Hexa-loops); default=on/1/True (bool/int)
        self.settings.noLP            = False     # Produce structures without lonely pairs (helices of length 1); default=off/0/False (bool/int)
        self.settings.noGU            = False     # Do not allow GU wobble pairs; default=off/0/False (bool/int)
        self.settings.noGUclosure     = False     # Do not allow GU pairs at the end of helices; default=off/0/False (bool/int)
        self.settings.logML           = False     # Recompute free energies of multi-branch loops using a logarithmic model; default=off/0/False (bool/int)
        self.settings.circ            = False     # Assume a circular (instead of linear) RNA molecule; default=off/0/False (bool/int)
        self.settings.canonicalBPonly = False     # Remove non-canonical base pairs from the structure constraint; default=off/0/False (bool/int)
        self.settings.uniq_ML         = False     # Create additional matrix for unique multi-branch loop prediction; default=off/0/False (bool/int)
        self.settings.energy_set      = 0         # Energy set, rarely used, see --energyModel; default=0 (int)
        self.settings.backtrack       = 1         # Whether to backtrack secondary structures; default=1 (int)
        self.settings.backtrack_type  = 'F'       # Set backtrack type, i.e. which DP matrix is used; default='F' (str)
        self.settings.compute_bpp     = True      # Compute base pair probabilities after partition function computation; default=on/1/True (bool/int)
        self.settings.max_bp_span     = -1        # Maximum base pair span; default=-1 (int)
        self.settings.min_loop_size   = 3         # Minimal loop size; default=3 (int)
        self.settings.window_size     = -1        # Window size for sliding window structure prediction approaches; default=-1 (int)
        self.settings.oldAliEn        = False     # Use old energy model for comparative structure prediction; default=off/0/False  (bool/int)
        self.settings.ribo            = False     # Use Ribosum Scoring in comparative structure prediction; default=off/0/False (bool/int)
        self.settings.cv_fact         = 1.0       # Co-variance scaling factor used in comparative structure prediction; default=1.0 (float)
        self.settings.nc_fact         = 1.0       # Unknown (not well defined in docs); likely for comparative structure prediction; default=1.0 (float)
        self.settings.sfact           = 1.07      # Scaling factor used to avodi under-/overflows in partition function computation; default=1.07 (float)
        '''

        # Setup model constraint variables
        self.enforce_constraints        = enforce_constraints
        self.constraints_options        = 16760832 # Not a mysterious value, see the macro values defined in '.../ViennaRNA-2.3.4/src/ViennaRNA/constraints_hard.h'
        self.pyindex                    = pyindex  # Global boolean for Python based 0-indexing or 1-indexing

        # Setup model parameter file and test_input variables
        RNA.read_parameter_file(self.parameter_directory+parameter_file)
        self.test_inputs                = test_inputs

        # Calculate dG_init adjustment for associating strands with specified temperature
        self.duplex_adjustment = duplex_adjustment
        self.dG_init_adjustment = self._calc_dG_init_adjustment()

        # Setup result structures as namedtuples
        self.PyVRNA_fold_result         = namedtuple('PyVRNA_fold_result',     'structure energy')
        self.PyVRNA_inverse_result      = namedtuple('PyVRNA_inverse_result',  'sequence distance')
        self.PyVRNA_centroid_result     = namedtuple('PyVRNA_centroid_result', 'structure energy distance')
        self.PyVRNA_ensemble_result     = namedtuple('PyVRNA_ensemble_result', 'structure energy')
        self.PyVRNA_bp_result           = namedtuple('PyVRNA_bp_result',      'length bpx bpy pkx pky gquad')

    def _test_sequences(self, sequence_list, func_name):
        """
        Tests if the sequences supplied to a particular function is valid for that function.

        NOTE:  This function is explictly called from other functions in this class that operates on sequences, when test_inputs=True.
               No external calls are required in order to ensure that sequences are valid.

        Usage: energy_model._test_sequences(sequence_list=list, func_name=string) # Validates sequences for func_name function
        """
        # Check if there are proper number of items in sequence_list
        len_sequence_dict = {'RNAcentroid': [1, 1], 'RNAcofold': [2, 2], 'RNAeval': [1, 2], 'RNAfold': [1, 1], 'RNAsubopt': [1, 2]}
        assert len_sequence_dict[func_name][0] <= len(sequence_list) <= len_sequence_dict[func_name][1], ''.join(['in ', func_name, ': ', ['at most', 'exactly'][not len_sequence_dict[func_name][0] < len_sequence_dict[func_name][1]], ' ', str(max(len_sequence_dict[func_name])), ' sequence(s) are accepted.'])
        
        # Check if first item is a string
        assert isinstance(sequence_list[0], str), 'in '+func_name+': sequence[_1] must be a string.'
        
        # Check if second item (if applicable) is also a string
        second_seq_funcs_set = {'RNAcofold', 'RNAeval', 'RNAsubopt'}
        if func_name in second_seq_funcs_set and len(sequence_list) > 1:
            assert isinstance(sequence_list[1], str), 'in {}: {}sequence[_2] must be a string.'.format(func_name, ['(optional) ', ''][func_name == 'RNAcofold'])
        
        # Check if all sequences are valid
        valid_charset = set('ATGCUatgcu&')
        for sequence in sequence_list:
            if not sequence is None:
                assert set(sequence) <= valid_charset, ''.join(['in ', func_name, ': sequence string(s) must only contain ATGCUatgcu& characters only.'])

    def _test_non_sequences(self, non_sequence_list, sequence_list, func_name, test_for='structure'):
        """
        Tests if the non-sequences supplied to a particular function is valid for that function. Non-sequences must be 'structure' or 'constraint'

        NOTE:  This function is explictly called from other functions in this class that operates on non-sequences, when test_inputs=True.
               No external calls are required in order to ensure that non-sequences are valid.

        Usage: energy_model._test_non_sequences(non_sequence_list=list, sequence_list=list, func_name=string, test_for=string) # Validates non-sequences for func_name function
        """

        if not non_sequence_list:
            return

        # Check if test_for is structures or constraints
        assert test_for == 'structure' or test_for == 'constraint', ''.join(['in ', func_name, ": test_for must be 'constraint' or 'structure'"])
        
        # Check if there are proper number of items in non_sequence_list
        assert 1 <= len(non_sequence_list) <= 2, ''.join(['in ', func_name, ': only one or two structures(s) are accepted.'])
        
        # Check if structures types are compatible with sequence types
        assert isinstance(non_sequence_list[0], str) or isinstance(non_sequence_list[0], type(None)), ''.join(['in ', func_name, ': ', test_for, '[_1] must be a string or a None object.'])
        if isinstance(non_sequence_list[0], type(None)):
            assert isinstance(non_sequence_list[1], type(None)), ''.join(['in ', func_name, ': ', test_for, '[_2] must be a None object since ', test_for,'[_1] is also None.'])
        if isinstance(sequence_list[1], str):
            assert isinstance(non_sequence_list[1], str) or isinstance(non_sequence_list[1], type(None)), ''.join(['in ', func_name, ': ', test_for, '[_2] must be a string or a None object.'])
        elif isinstance(sequence_list[1], type(None)):
            assert isinstance(non_sequence_list[1], type(None)), ''.join(['in ', func_name, ': ', test_for, '[_2] must be a None object since sequence[_2] is also None.'])
        
        # Check if structures are valid, and their lengths are same as their sequences
        valid_chars   = '.([+])&' if test_for == 'structure' else '.x([+])&'
        valid_charset = set(valid_chars)
        for i in xrange(2):
            string_i = str(i+1)
            if not (non_sequence_list[i] is None or sequence_list[i] is None):
                assert len(non_sequence_list[i]) == len(sequence_list[i]), ''.join(['in ', func_name, ': ', test_for,'_', string_i, ' and sequence_', string_i,' must be of same length.'])
                assert set(non_sequence_list[i]) <= valid_charset,         ''.join(['in ', func_name, ': ', test_for,'_', string_i, ' must be a string of ', valid_chars,' characters only.'])

    def _test_aptamer_inputs(self, aptamers, aptamer_constraints, dG_ligands):
        """
        Tests if the aptamer constraints supplied to a particular function are valid for that function.

        NOTE:  This function is explictly called from other functions in this class that operates on non-sequences, when test_inputs=True.
               No external calls are required in order to ensure that non-sequences are valid.
        
        Usage: energy_model._test_non_sequences(non_sequence_list=list, sequence_list=list, func_name=string, test_for=string) # Validates non-sequences for func_name function
        """

        if not aptamers and not aptamer_constraints:
            return

        assert len(aptamers)==len(aptamer_constraints), "Make sure to provide equal number of aptamers and aptamer_constraint values"
        assert len(aptamers)==len(dG_ligands), "Make sure to provide equal number of aptamers and dG_ligand values"
        for aptamer in aptamers:
            self._test_sequences(aptamer)
        assert all(isinstance(dG,float) for dG in dG_ligands), "Make sure dG_ligand values are floats in kcal/mol"
        self._test_non_sequences(aptamer_constraints,aptamers,'RNAfold','constraint')

    def _test_bp_tuple(self, length, bpx, bpy, pkx, pky, gquad, func_name):
        """
        Tests if all inputs for a bp_tuple are valid.

        Usage: energy_model._test_bp_tuple(bp_tuple.bpx, bp_tuple.bpy, bp_tuple.pkx, bp_tuple.pky, bp_tuple.gquad, func_name=string) # Validates bpx/bpy, pkx/pky and gquad in func_name
        """
        assert isinstance(length,list),                         ''.join(['in ', func_name, ': length must be an integer.'])
        assert length >= max([len(bpx), len(pkx), len(gquad)]), ''.join(['in ', func_name, ': length must be greater than or equal to the length of bpx/bpy, pkx/pky and quad lists.'])
        assert len(bpx) == len(bpy),                            ''.join(['in ', func_name, ': bpx and bpy must have same length.'])
        assert len(pkx) == len(pky),                            ''.join(['in ', func_name, ': pkx and pky must have same length.'])
        for index, gquad_tuple in enumerate(gquad):
            assert len(gquad_tuple) == 4, ''.join(['in ', func_name, ': gquad[', str(index), '] must contain 4 indices.'])

    def _calc_dG_init_adjustment(self):
        """
        Adjustment according to Dirks et al., Thermodynamic Analysis of Interacting Nucleic Acid Strands
        See footnote 13: "Based on dimensional analysis, we define our concentrations as mole fractions rather than
        molarities. Therefore, the free energy of strand association for a complex of L strands is..."
        
        Based on a partition function analysis of dilute solutions of interacting strands
        in a fixed volume (the "box")

        NOTE:  This function is explictly called from other functions in this class that returns minimum free energy of a duplex.
               No external calls are required in order to ensure that energies are adjusted.
        """
        
        kB = 0.00198717 # Boltzmann constant in kcal/mol/K
        T = self.settings.temperature
        a = [-3.983035, 301.797, 522528.9, 69.34881, 999.974950]

        # Calculate the number of moles of water per liter (molarity) at temperature (T in deg C)
        # Density of water calculated using data from 
        # Tanaka M., Girard, G., Davis, R., Peuto A., Bignell, N.
        # Recommended table for the density of water..., Metrologia, 2001, 38, 301-309
        pH2O = a[4] * (1 - (T+a[0])**2.0*(T+a[1])/a[2]/(T+a[3])) / 18.0152

        return -kB*(T+273.15)*math.log(pH2O)

    def RNAcentroid(self, sequence, constraint=None, aptamers=[], aptamer_constraints=[], dG_ligands=[]):
        """
        Computes the centroid structure, its energy and distance for a sequence.

        Usage: centroid_result = energy_model.RNAcentroid(sequence, constraint)  # Executes RNAcentroid
               centroid_result.structure                                         # Retrieves centroid structure
               centroid_result.energy                                            # Retrieves centroid energy
               centroid_result.distance                                          # Retrieves centroid distance
               centroid_structure = energy_model.RNAcentroid(sequence).structure # Executes RNAcentroid and retrieves structure only
        """
        # Test if sequence is valid
        if self.test_inputs:
            self._test_sequences(sequence_list=[sequence], func_name='RNAcentroid')
            self._test_aptamer_inputs(aptamers,aptamer_constraints,dG_ligands)

        # Setup ViennaRNA library objects and call centroid function
        fc_obj = RNA.fold_compound(sequence, self.settings)
        if not constraint is None:
            if not self.enforce_constraints:
                fc_obj.hc_add_from_db(constraint)
            else:
                fc_obj.hc_add_from_db(constraint, self.constraints_options)
        if aptamers:
            for aptamer,fld,dG in zip(aptamers,aptamer_constraints,dG_ligands):
                fc_obj.sc_add_hi_motif(aptamer,fld,dG,0) # 0 is VRNA_OPTION_DEFAULT                
        fc_obj.pf()
        (structure, distance) = fc_obj.centroid()

        # Return the centroid structure, energy and distance
        return self.PyVRNA_centroid_result(structure=structure, energy=fc_obj.eval_structure(structure), distance=distance)

    def RNAcofold(self, sequences, constraints=[], aptamers=[], aptamer_constraints=[], dG_ligands=[]):
        """
        Computes the minimum free energy structure and its corresponding energy for a duplex.
        
        Usage: cofold_result = energy_model.RNAcofold(sequences=[sequence1,sequence2], constraints=[constraint1,constraint2]) # Executes RNAcofold
               cofold_result.structure                                                                                        # Retrieves mfe structure
               cofold_result.energy                                                                                           # Retrieves mfe
               cofold_energy = energy_model.RNAcofold(sequence_1, sequence_2, constraint_1, constraint_2).structure           # Executes RNAcofold and retrieves the mfe energy only
        """
        # Test if sequences and constraints are valid
        if self.test_inputs:
            self._test_sequences(sequence_list=sequences, func_name='RNAcofold')
            self._test_non_sequences(non_sequence_list=constraints, sequence_list=sequences, func_name='RNAcofold', test_for='constraint')
            self._test_aptamer_inputs(aptamers,aptamer_constraints,dG_ligands)

        # Setup ViennaRNA library objects and call mfe_dimer (cofold) function
        fc_obj = RNA.fold_compound("&".join(sequences), self.settings)
        if constraints:
            if not self.enforce_constraints:
                fc_obj.hc_add_from_db("".join(constraints)) # Ayaan I'm not following this logic - ACR 04.21.17
            else:
                fc_obj.hc_add_from_db("".join(constraints), self.constraints_options)
        if aptamers:
            for aptamer,fld,dG in zip(aptamers,aptamer_constraints,dG_ligands):
                fc_obj.sc_add_hi_motif(aptamer,fld,dG,0) # 0 is VRNA_OPTION_DEFAULT                
        structure, energy = fc_obj.mfe_dimer()
        
        if self.duplex_adjustment and set('.') < set(structure):
            energy += self.dG_init_adjustment

        structure = structure[:len(sequences[0])] + "&" + structure[len(sequences[0]):]

        # Return the cofold structure and energy
        return self.PyVRNA_fold_result(structure=structure, energy=energy)

    def RNAensemble(self, sequence, constraint=None, aptamers=[], aptamer_constraints=[], dG_ligands=[]):
        """
        Computes the ensemble structure, the Gibbs free energy based on the partition fxn and distance.

        Usage: ensemble_result = energy_model.RNAensemble(sequence, constraint)  # Executes RNAensemble
               ensemble_result.structure                                         # Retrieves ensemble structure
               ensemble_result.energy                                            # Retrieves ensemble energy
               ensemble_structure = energy_model.RNAcentroid(sequence).structure # Executes RNAensemble and retrieves structure only
        """
        # Test if sequence is valid
        if self.test_inputs:
            self._test_sequences(sequence_list=[sequence], func_name='RNAensemble')
            self._test_aptamer_inputs(aptamers,aptamer_constraints,dG_ligands)

        # Setup ViennaRNA library objects and call centroid function
        fc_obj = RNA.fold_compound(sequence, self.settings)
        if not constraint is None:
            if not self.enforce_constraints:
                fc_obj.hc_add_from_db(constraint)
            else:
                fc_obj.hc_add_from_db(constraint, self.constraints_options)
        if aptamers:
            for aptamer,fld,dG in zip(aptamers,aptamer_constraints,dG_ligands):
                fc_obj.sc_add_hi_motif(aptamer,fld,dG,0) # 0 is VRNA_OPTION_DEFAULT        
        (structure, energy) = fc_obj.pf()

        # Return the ensemble structure, energy and distance
        return self.PyVRNA_ensemble_result(structure=structure, energy=energy)

    def RNAeval(self, sequences, structures, aptamers=[], aptamer_constraints=[], dG_ligands=[]):
        """
        Computes the minimum free energy of a given structure and its corresponding sequence or duplex.

        Usage: energy = energy_model.RNAeval(sequences, structures)        # Executes RNAeval
        """
        # Test if sequences and structures are valid
        if self.test_inputs:
            self._test_sequences(sequence_list=sequences, func_name='RNAeval')
            self._test_non_sequences(non_sequence_list=structures, sequence_list=sequences, func_name='RNAeval', test_for='structure')

        energy = RNA.fold_compound("&".join(sequences), self.settings).eval_structure("".join(structures))

        if self.duplex_adjustment and len(sequences) > 1 and set('.') < set("".join(structures)):
            energy += self.dG_init_adjustment

        # Return the energy (chained setup of ViennaRNA library object and call to eval_structure function)
        return energy

    def RNAfold(self, sequence, constraint=None, aptamers=[], aptamer_constraints=[], dG_ligands=[]):
        """
        Computes the minimum free energy structure and its corresponding energy for a sequence.

        Usage: fold_result = energy_model.RNAfold(sequence, constraint)             # Executes RNAfold
               fold_result.structure                                                # Retrieves mfe structure
               fold_result.energy                                                   # Retrieves mfe
               mfe_structure = energy_model.RNAfold(sequence, constraint).structure # Executes RNAfold and retrieves the mfe structure only
        """

        # Test if sequence and constraint is valid
        if self.test_inputs:
            self._test_sequences(sequence_list=[sequence], func_name='RNAfold')
            self._test_non_sequences(non_sequence_list=[constraint, None], sequence_list=[sequence, None], func_name='RNAfold', test_for='constraint')
            self._test_aptamer_inputs(aptamers,aptamer_constraints,dG_ligands)
        
        # Setup ViennaRNA library objects and call mfe (fold) function
        fc_obj = RNA.fold_compound(sequence, self.settings)
        
        # Add constraint if provided
        if not constraint is None:
            if not self.enforce_constraints:
                fc_obj.hc_add_from_db(constraint)
            else:
                fc_obj.hc_add_from_db(constraint, self.constraints_options)

        # Add aptamers if provided
        if aptamers:
            for aptamer,fld,dG in zip(aptamers,aptamer_constraints,dG_ligands):
                fc_obj.sc_add_hi_motif(aptamer,fld,dG,0) # 0 is VRNA_OPTION_DEFAULT

        structure, energy = fc_obj.mfe()

        # Return the fold structure and energy
        return self.PyVRNA_fold_result(structure=structure, energy=energy)

    def RNAinverse(self, sequence, structure):
        """
        Computes a sequence conforming to given minimum free energy structure starting from given sequence and its distance to starting sequence

        Usage: inverse_result = energy_model.RNAinverse(sequence, structure)           # Executes RNAinverse
               inverse_result.structure                                                # Retrieves inverted sequence
               inverse_result.distance                                                 # Retrieves distance from startng sequence
               inverse_sequence = energy_model.RNAfold(sequence, constraint).sequence  # Executes RNAinverse and retrieves the sequence only
        """

        # Test if sequence and constraint is valid
        # Coming soon

        sequence, distance = RNA.inverse_fold(sequence, structure)

        # Return the sequence and the distance
        return self.PyVRNA_inverse_result(sequence=sequence, distance=distance)

    def RNAsubopt(self, sequences, constraints=[], delta_energy=5, aptamers=[], aptamer_constraints=[], dG_ligands=[]):
        """
        Computes the suboptimal structures for a sequence or duplex with optional constraints within a delta_energy range from its mfe.
        
        NOTE:  delta_energy must always be an integer. The ViennaRNA C library returns all suboptimal structures
               within a certain delta_energy range of the mfe, for example +5 of mfe etc.

        Usage: subopt_list = energy_model.RNAsubopt(sequences, constraints, delta_energy) # Executes RNAsubopt
               subopt_list[0].structure                                                                               # Retrieves 0th suboptimal structure
               subopt_list[0].energy                                                                                 # Retrieves mfe of 10th suboptimal structure
               subopt_energy_list = [solution.energy for solution in subopt_list]                                     # Filters out mfe of every suboptimal structure
        """
        # if set("".join(constraints)) == set('.'):
            # constraints = []
        # Test if sequences, structures and constraints are valid
        if self.test_inputs:
            self._test_sequences(sequence_list=sequences, func_name='RNAsubopt')
            self._test_non_sequences(non_sequence_list=constraints, sequence_list=sequences, func_name='RNAsubopt', test_for='constraint')
            self._test_aptamer_inputs(aptamers,aptamer_constraints,dG_ligands)
        delta_energy = int(round(delta_energy))
        
        # Setup ViennaRNA library objects, call subopt function and sort results
        fc_obj = RNA.fold_compound("&".join(sequences), self.settings)
        if constraints:
            if not self.enforce_constraints:
                fc_obj.hc_add_from_db("".join(constraints))
            else:
                fc_obj.hc_add_from_db("".join(constraints), self.constraints_options)

        # Add aptamers if provided
        if aptamers:
            for aptamer,fld,dG in zip(aptamers,aptamer_constraints,dG_ligands):
                fc_obj.sc_add_hi_motif(aptamer,fld,dG,0) # 0 is VRNA_OPTION_DEFAULT

        if self.duplex_adjustment and len(sequences) > 1:
            subopt_list = [self.PyVRNA_fold_result(structure=solution.structure, energy=solution.energy + self.dG_init_adjustment if set('&.') < set(solution.structure) else 0.0) for solution in fc_obj.subopt(delta_energy*100)]
        else:
            subopt_list = [self.PyVRNA_fold_result(structure=solution.structure, energy=solution.energy) for solution in fc_obj.subopt(delta_energy*100)]

        subopt_list.sort(key=itemgetter(1))
        
        # Return the fold structure and energy of all suboptimal structures
        return subopt_list

    def create_bp_tuple(self,length=[], bpx=[], bpy=[], pkx=[], pky=[], gquad=[]):
        """
        Returns a customized bp_tuple from input bpx, bpy, pkx, pky and gquad lists.

        NOTE:  bp_tuple is a tuple of (length, base_pair_x_list, base_pair_y_list, pseudo_pair_x_list, pseudo_pair_y_list, gquad_list).

        Usage: custom_bp_tuple = energy_model.create_bp_tuple(bpx=some_list, bpy=some_list) # Executes create_bp_tuple and returns the specified bp_tuple
        """
        # Test if the inputs are valid
        if self.test_inputs:
            self._test_bp_tuple(length=length, bpx=bpx, bpy=bpy, pkx=pkx, pky=pky, gquad=gquad, func_name='create_bp_tuple')
        
        # Create and return the customized bp_tuple
        return self.PyVRNA_bp_result (length=length, bpx=list(bpx), bpy=list(bpy), pkx=list(pkx), pky=list(pky), gquad=list(gquad))

    def vienna2bp(self,vienna_string, pyindex=None):
        """
        Converts the vienna_string (dot bracket string representing a mfe structure) to bp_tuple.

        NOTE:  bp_tuple is a tuple of (length, base_pair_x_list, base_pair_y_list, pseudo_pair_x_list, pseudo_pair_y_list, gquad_list).
               This function can be used for involution with bp2vienna, meaning bptuple2vienna(vienna2bp(vienna_string)) = vienna_string.

        Usage: bp_tuple = energy_model.vienna2bp(vienna_string=some_result.structure) # Executes vienna2bp
               bp_tuple.length                                                          # Retrieves the length of vienna_string
               bp_tuple.bpx                                                             # Retrieves the base_pair_x_list
               bp_tuple.bpy                                                             # Retrieves the base_pair_y_list
               bp_tuple.pkx                                                             # Retrieves the pseudo_pair_x_list
               bp_tuple.pky                                                             # Retrieves the pseudo_pair_y_list
               bp_tuple.gquad                                                           # Retrieves the gquad_list
        """
        # Setup indexing
        if pyindex is None:
            pyindex = self.pyindex
        
        # Test if vienna_string has all valid characters
        if self.test_inputs:
            assert set(vienna_string) <= set('&.([+])'), 'in vienna2bp: vienna_string must be a string of .([+]) characters only.'
            assert vienna_string.count('&') <= 1, 'Odd input. Should only have 1 & if a cofold, you provided: {}.'.format(vienna_string)

        if '&' in vienna_string:
            split_string = vienna_string.split('&')
            length = [len(x) for x in split_string]
            vienna_string = "".join(split_string)
        else:
            length = [len(vienna_string)]

        # Setup data structures and variables for computing bp 
        bp_tuple          = self.create_bp_tuple(length=length)
        bp_stack, pk_stack = [], []
        bp_dict, pk_dict   = OrderedDict(), OrderedDict()
        has_gquad          = False
        gquad_matrix       = {0:[], 1:[], 2:[], 3:[]}
        string_index       = 0
        gquad_index        = 0
        
        # Extract all parentheses and bracket indices in order
        for index, char in enumerate(vienna_string):
            if not char in ['x', '.']:
                if   char == '(':
                    bp_stack.append(index+1)
                    bp_tuple.bpx.append(index+1)
                    bp_dict[index+1] = -1
                elif char == '[':
                    pk_stack.append(index+1)                    
                    bp_tuple.pkx.append(index+1)
                    pk_dict[index+1] = -1
                elif char == ')':
                    bp_dict[bp_stack.pop()] = index+1
                elif char == ']':
                    pk_dict[pk_stack.pop()] = index+1
                elif char == '+':
                    has_gquad = True
        else:
            for val in bp_dict.itervalues():
                bp_tuple.bpy.append(val)
            for val in pk_dict.itervalues():
                bp_tuple.pky.append(val)
        
        # Extract all quadruplex indices in order
        while has_gquad and string_index < len(vienna_string):
            if vienna_string[string_index] == '+' and \
                (gquad_index == 0 or gquad_index == 2):
                if gquad_index == 0:
                    gquad_matrix[0].append(deque()); gquad_matrix[1].append(deque()); gquad_matrix[2].append(deque()); gquad_matrix[3].append(deque())
                while string_index < len(vienna_string) and \
                    vienna_string[string_index] == '+':
                    gquad_matrix[gquad_index][-1].append(string_index+1)
                    string_index += 1
                gquad_index = (gquad_index + 1) % 4
            elif vienna_string[string_index] == '+' and \
                (gquad_index == 1 or gquad_index == 3):
                #print 'Initial', string_index, len(vienna_string)
                while string_index < len(vienna_string) and \
                    vienna_string[string_index] == '+':
                    gquad_matrix[gquad_index][-1].appendleft(string_index+1)
                    string_index += 1
                    #print 'Updated', string_index, len(vienna_string)
                #print 'Exited', string_index, len(vienna_string)
                gquad_index = (gquad_index + 1) % 4
            else:
                string_index += 1
        else:
            if has_gquad:
                for i in xrange(len(gquad_matrix[0])):
                    bp_tuple.gquad.extend(zip(gquad_matrix[0][i], gquad_matrix[1][i], gquad_matrix[2][i], gquad_matrix[3][i]))
        
        if pyindex:
            bpx = [i-1 for i in bp_tuple.bpx]
            bpy = [j-1 for j in bp_tuple.bpy]
            pkx = [i-1 for i in bp_tuple.pkx]
            pky = [j-1 for j in bp_tuple.pky]
            gquad = [tuple(j-1 for j in i) for i in bp_tuple.gquad]

            bp_tuple = self.create_bp_tuple(length,bpx,bpy,pkx,pky,gquad)
        #print bp_tuple.gquad
        return bp_tuple

    def bp2vienna(self,length, bpx=[], bpy=[], pkx=[], pky=[], gquad=[], pyindex=None):
        """
        Creates a customized bp_tuple from input bpx, bpy, pkx, pky and gquad lists and returns the vienna_string encoded by it.

        NOTE:  bp_tuple is a tuple of (length, base_pair_x_list, base_pair_y_list, pseudo_pair_x_list, pseudo_pair_y_list, gquad_list).

        Usage: custom_vienna_string = energy_model.bp2vienna(bpx=some_list, bpy=some_list) # Executes bp2vienna and returns the vienna_string
        """
        # Setup indexing
        if pyindex is None:
            pyindex = self.pyindex

        return self.bptuple2vienna(self.create_bp_tuple(length, bpx, bpy, pkx, pky, gquad))

    def bptuple2vienna(self, bp_tuple, pyindex=None):
        """
        Converts an bp_tuple to vienna_string (dot bracket string representing a mfe structure).

        NOTE:  bp_tuple is a tuple of (length, base_pair_x_list, base_pair_y_list, pseudo_pair_x_list, pseudo_pair_y_list, gquad_list).
               This function can be used for involution with vienna2bp, meaning vienna2bp(bptuple2vienna(bp_tuple)) = bp_tuple.

        Usage: vienna_string = energy_model.bptuple2vienna(bp_tuple) # Executes bptuple2vienna and retrieves the mfe structure
        """
        # Setup indexing
        if pyindex is None:
            pyindex = self.pyindex
            
        # Test if vienna_string has all valid characters
        if self.test_inputs:
            self._test_bp_tuple(length=bp_tuple.length, bpx=bp_tuple.bpx, bpy=bp_tuple.bpy, pkx=bp_tuple.pkx, pky=bp_tuple.pky, gquad=bp_tuple.gquad, func_name='bptuple2vienna')

        if pyindex:
            bpx = [i+1 for i in bp_tuple.bpx]
            bpy = [j+1 for j in bp_tuple.bpy]
            pkx = [i+1 for i in bp_tuple.pkx]
            pky = [j+1 for j in bp_tuple.pky]
            gquad = [tuple(j+1 for j in i) for i in bp_tuple.gquad]
            bp_tuple = self.create_bp_tuple(bp_tuple.length,bpx,bpy,pkx,pky,gquad)        

        # Place all appropriate symbols in vienna_string_list
        vienna_string_list = ['.'] * sum(bp_tuple.length)
        for i, j in izip(bp_tuple.bpx, bp_tuple.bpy):
            if i >= 0 and j >= 0: vienna_string_list[i-1], vienna_string_list[j-1] = '(', ')'
            else                : break
        for i, j in izip(bp_tuple.pkx, bp_tuple.pky):
            if i >= 0 and j >= 0: vienna_string_list[i-1], vienna_string_list[j-1] = '[', ']'
            else                : break
        for i, j, k, l in chain(bp_tuple.gquad):
            if i >= 0 and j >= 0 and k >= 0 and l >= 0: vienna_string_list[i-1], vienna_string_list[j-1], vienna_string_list[k-1], vienna_string_list[l-1] = '+', '+', '+', '+'
            else                                      : break
        
        if len(bp_tuple.length) == 2:
            vienna_string_list.insert(bp_tuple.length[0],'&')

        # Return the vienna_string from its list
        return "".join(vienna_string_list)


# Legacy
class ViennaRNA(dict):
    # Legacy
    def __init__(self, Sequence_List, material, Gquad=False):
        # Legacy: Check if sequences in Sequence_List are valid
        exp = re.compile('[ATGCU]',re.IGNORECASE)
        for seq in Sequence_List:
            if exp.match(seq) is None:
                raise ValueError("Invalid letters found in inputted sequences. Only ATGCU allowed. \n Sequence is \"" + str(seq) + "\".")
        
        # Legacy: Setting object parameters
        self.ran                = 0
        self["sequences"]       = Sequence_List
        self["material"]        = material      
        parameter_files         = ["dna_mathews1999", "rna_turner1999", "dna_mathews2004", "rna_turner2004", "rna_andronescu2007"]
        material_par            = material+".par"
        self["RNA_model_param"] = material_par if material in parameter_files else parameter_files[-1]+".par"
        self["Gquad"]           = Gquad
        self["Gquad_param"]     = ["-g"] if Gquad else []

        # New: Refactoring dangle setup for other functions
        self.dangles_dict = {'all':2, 'some':1, 'none': 0}

    # Legacy
    def centroid(self, strands, constraint=None, Temp=37.0, dangles="all", outputPS=False):
        # Legacy: Checks and setup
        if Temp <= 0:        raise ValueError("The specified temperature must be greater than zero.")
        if len(strands) > 1: raise ValueError("Two RNA strands are inputted. ViennaRNA does NOT return Centroid for RNAcofold.")
        self["Centroid_composition"] = strands

        # New: PyVRNA execution
        energy_model = PyVRNA(temperature=Temp, dangles=self.dangles_dict[dangles], gquad=self["Gquad"], parameter_file=self["RNA_model_param"], test_inputs=False)
        structure, energy, distance = energy_model.RNAcentroid(sequence=self["sequences"][0])

        # Legacy: Parsing and storing output
        bp_tuple = energy_model.vienna2bp(structure, pyindex=False)
        self["program"]                 = "Centroid"
        self["totalnt"]                 = bp_tuple.length
        self["Centroid_energy"]         = [energy]
        self['Centroid_bracket_string'] = structure
        self["Centroid_basepairing_x"]  = [bp_tuple.bpx]
        self["Centroid_basepairing_y"]  = [bp_tuple.bpy]

    # Legacy
    def convert_bracket_to_numbered_pairs(self, bracket_string):
        # New: PyVRNA execution
        bp_tuple = PyVRNA(test_inputs=False).vienna2bp(bracket_string, pyindex=False)
        return [bp_tuple.length, bp_tuple.bpx, bp_tuple.bpy, bp_tuple.pkx, bp_tuple.pky]

    # Legacy
    def convert_numbered_pairs_to_bracket(self, strands, bp_x, bp_y, PK_bp_x=[], PK_bp_y=[], Gquad_bp=[]):
        # New: PyVRNA execution
        energy_model = PyVRNA(test_inputs=False)
        bp_tuple = energy_model.PyVRNA_bp_result(length=strands, bpx=bp_x, bpy=bp_y, pkx=PK_bp_x, pky=PK_bp_y, gquad=Gquad_bp)
        return energy_model.bptuple2vienna(bp_tuple,pyindex=False)

    # Legacy
    def energy(self, strands, base_pairing_x, base_pairing_y, Temp=37.0, dangles="all"):
        # Legacy: Checks and setup
        if Temp <= 0: raise ValueError("The specified temperature must be greater than zero.")
        self["energy_composition"] = strands
        strands                    = map(len, self["sequences"])
        sequences                  = self["sequences"] if len(strands) == 2 else [self["sequences"][0]]

        # New: PyVRNA execution
        energy_model = PyVRNA(temperature=Temp, dangles=self.dangles_dict[dangles], gquad=self["Gquad"], parameter_file=self["RNA_model_param"], test_inputs=False)
        bp_tuple     = energy_model.PyVRNA_bp_result(length=strands, bpx=base_pairing_x, bpy=base_pairing_y, pkx=[], pky=[], gquad=[])
        energy       = energy_model.RNAeval(sequences, energy_model.bptuple2vienna(bp_tuple, pyindex=False).split('&'))

        # Legacy: Parsing and storing output
        self["program"]              = "energy"
        self["energy_energy"]        = [energy]
        self["energy_basepairing_x"] = [base_pairing_x]
        self["energy_basepairing_y"] = [base_pairing_y]
        return energy

    # Legacy
    def mfe(self, strands, constraints=None, Temp=37.0, dangles="all", outputPS=False, duplex=False):
        # Legacy: Checks and setup
        if Temp <= 0: raise ValueError("The specified temperature must be greater than zero.")
        self["mfe_composition"] = strands
        
        # New: PyVRNA execution
        constraint = constraints.split('&') if (not constraints is None and '&' in constraints) else constraints
        energy_model = PyVRNA(temperature=Temp, dangles=self.dangles_dict[dangles], gquad=self["Gquad"], parameter_file=self["RNA_model_param"], test_inputs=False)
        if   len(strands) == 1:
            structure, energy = energy_model.RNAfold(sequence=self["sequences"][0], constraint=constraint)
        elif len(strands) == 2:
            structure, energy = energy_model.RNAcofold(self["sequences"], constraints=constraint)
        else:
            raise ValueError("Three RNA strands are inputted. ViennaRNA does NOT return structure and energy for three sequences in RNA(co)fold.")

        # Legacy: Parsing and storing output
        bp_tuple = energy_model.vienna2bp(structure, pyindex=False)
        self["program"]            = "mfe"
        self["totalnt"]            = bp_tuple.length
        self["mfe_energy"]         = [energy]
        self["mfe_bracket_string"] = structure
        self["mfe_basepairing_x"]  = [bp_tuple.bpx]
        self["mfe_basepairing_y"]  = [bp_tuple.bpy]

    # Legacy
    def subopt(self, strands, energy_gap, Temp=37.0, dangles="all", constraints=None, outputPS=False):
        # Legacy: Checks and setup
        if Temp <= 0: raise ValueError("The specified temperature must be greater than zero.")
        self["subopt_composition"]   = strands
        self["subopt_energy"]        = []
        self["subopt_basepairing_x"] = []
        self["subopt_basepairing_y"] = []

        # New: PyVRNA execution
        constraint = constraints.split('&') if (not constraints is None and '&' in constraints) else constraints
        energy_model = PyVRNA(temperature=Temp, dangles=self.dangles_dict[dangles], gquad=self["Gquad"], parameter_file=self["RNA_model_param"], test_inputs=False)
        if   len(strands) == 1:
            results = energy_model.RNAsubopt(sequences=[self["sequences"][0]], constraints=constraint, delta_energy=energy_gap)
        elif len(strands) == 2:
            results = energy_model.RNAsubopt(sequences=self["sequences"], constraints=constraint, delta_energy=energy_gap)
        else:
            raise ValueError("Three RNA strands are inputted. ViennaRNA does NOT return suboptimal structures and energies for three sequences in RNAsubopt.")
        
        # Legacy: Parsing and storing output
        for structure, energy in results:
            bp_tuple = energy_model.vienna2bp(structure,pyindex=False)
            self["subopt_energy"].append(energy)
            self["subopt_basepairing_x"].append(bp_tuple.bpx)
            self["subopt_basepairing_y"].append(bp_tuple.bpy)
        self["program"]           = "subopt"
        self["totalnt"]           = strands
        self["subopt_NumStructs"] = len(results)


def tests():
    """
    Tests each function in PyVRNA class
    """
    print ""
    print "### === Testing PyVRNA class === ###"

    # Tests for parameter files
    print 'Testing all parameter files...',
    sequence     = "CGCAGGGAUACCCGCG"
    parameter_files = ["dna_mathews1999.par", "rna_turner1999.par", "dna_mathews2004.par", "rna_turner2004.par", "rna_andronescu2007.par"]
    energy_model = PyVRNA(dangles=0, gquad=False, parameter_file=parameter_files[0])
    fold_result  = energy_model.RNAfold(sequence)
    assert [fold_result.structure, fold_result.energy] == ['(((.(((...))))))', -3.200000047683716]
    energy_model = None
    energy_model = PyVRNA(dangles=0, gquad=False, parameter_file=parameter_files[1])
    fold_result  = energy_model.RNAfold(sequence)
    assert [fold_result.structure, fold_result.energy] == ['(((.((....)).)))', -5.5]
    energy_model = None
    energy_model = PyVRNA(dangles=0, gquad=False, parameter_file=parameter_files[2])
    fold_result  = energy_model.RNAfold(sequence)
    assert [fold_result.structure, fold_result.energy] == ['(((.(((...))))))', -3.9000000953674316]
    energy_model = None
    energy_model = PyVRNA(dangles=0, gquad=False, parameter_file=parameter_files[3])
    fold_result  = energy_model.RNAfold(sequence)
    assert [fold_result.structure, fold_result.energy] == ['(((.(((...))))))', -5.599999904632568]
    energy_model = None
    energy_model = PyVRNA(dangles=0, gquad=False, parameter_file=parameter_files[4])
    fold_result  = energy_model.RNAfold(sequence)
    assert [fold_result.structure, fold_result.energy] == ['(((.(((...))))))', -4.619999885559082]
    energy_model = None
    print 'Good.'

    # Tests for RNAcentroid
    print 'Testing RNAcentroid...',
    sequence     = 'CGACGUAGAUGCUAGCUGACUCGAUGC'
    energy_model = PyVRNA(dangles=2, gquad=False, parameter_file='rna_turner2004.par')
    cntrd_result = energy_model.RNAcentroid(sequence)
    assert [cntrd_result.structure, cntrd_result.energy, cntrd_result.distance] == ['(((.(.((.......))..))))....', 1.399999976158142, 3.345900802497659]
    energy_model = None
    sequence     = 'CGCAGGGAUACCCGCG' + 'GCGCCCAUAGGGACGC'
    energy_model = PyVRNA(dangles=2, gquad=False, parameter_file='rna_turner2004.par')
    cntrd_result = energy_model.RNAcentroid(sequence)
    assert [cntrd_result.structure, cntrd_result.energy, cntrd_result.distance] == ['(((.(((...))))))((((((...))).)))', -13.5, 2.4947844957445344]
    energy_model = None
    print 'Good.'

    print 'Testing RNAcentroid with aptamers...',
    # example 1
    sequence = "AGACAUAGCGAUCAAGUGAUACCAGCAUCGUCUUGAUGCCCUUGGCAGCACUUCAUAGCUAGCUACG"
    theophylline_aptamer  = "GAUACCAG&CCCUUGGCAGC"
    theophylline_constraint = "(...((((&)...)))...)"
    dG_theophylline = -9.22
    energy_model = PyVRNA(dangles=2, gquad=False, parameter_file='rna_turner2004.par')
    fold_result  = energy_model.RNAcentroid(sequence,aptamers=[theophylline_aptamer],aptamer_constraints=[theophylline_constraint],dG_ligands=[dG_theophylline])
    assert [fold_result.structure, fold_result.energy] == ['.....((((..(.(((((...((((((((.....)))))...)))...))))).)..))))......',-17.0]
    fold_result  = energy_model.RNAcentroid(sequence)
    assert [fold_result.structure, fold_result.energy] == ['.....((((..(.(((((...((((((((.....)))))...)))...))))).)..))))......',-17.0]
    energy_model = None
    
    # example 2
    sequence = "CCAGAAUCAUGCAAGUGCGUAAGAUAGUCGCGGGCUUUAAUUACGCUAUAUUAUAUACCCAAUUCU"
    streptavidin_aptamer =              'CCAGAAUCAUGCAAGUGCGUAAGAUAGUCGCGGG'
    streptavidin_aptamer_structure =    '((xxxxxxxxxxxxx((((xxxxxxxxx))))))'
    dG_streptavidin = -10.15
    energy_model = PyVRNA(dangles=2, gquad=False, parameter_file='rna_turner2004.par')
    fold_result  = energy_model.RNAcentroid(sequence,aptamers=[streptavidin_aptamer],aptamer_constraints=[streptavidin_aptamer_structure],dG_ligands=[dG_streptavidin])
    assert [fold_result.structure, fold_result.energy] == ['((.............((((.........))))))................................',-0.4000000059604645]
    fold_result  = energy_model.RNAcentroid(sequence)
    assert [fold_result.structure, fold_result.energy] == ['................((((((...((((...))))....))))))....................',-9.600000381469727]
    
    # example 3 - add multiple apatmers
    sequence = "AGACAUAGCGAUCAAGUGAUACCAGCAUCGUCUUGAUGCCCUUGGCAGCACUUCAUAGCUAGCCAGAAUCAUGCAAGUGCGUAAGAUAGUCGCGGGCUUUAAUUACGCUAUAUUAUAUACCCAAUUCU"
    energy_model = PyVRNA(dangles=2, gquad=False, parameter_file='rna_turner2004.par')
    fold_result  = energy_model.RNAcentroid(sequence,aptamers=[theophylline_aptamer,streptavidin_aptamer], \
        aptamer_constraints=[theophylline_constraint,streptavidin_aptamer_structure],dG_ligands=[0.0,0.0]) # what if the aptamer provides no additional free energy?
    assert [fold_result.structure, fold_result.energy] == ['....((((((((((((((((......)))).))))))((((((((((((........))).)))))......((....)).............))))........)))))).................',-29.399999618530273]
    fold_result  = energy_model.RNAcentroid(sequence,aptamers=[theophylline_aptamer,streptavidin_aptamer], \
        aptamer_constraints=[theophylline_constraint,streptavidin_aptamer_structure],dG_ligands=[dG_theophylline,dG_streptavidin]) # now with the aptamer binding free energies
    assert [fold_result.structure, fold_result.energy] == ['....((((((((((((((((......)))).))))))(((...)))(((........)))((((.............((((.........)))))))).......)))))).................',-24.799999237060547]
    print 'Good.'

    # Tests for RNAcofold
    print 'Testing RNAcofold...',
    sequences    = ["CGCAGGGAUACCCGCG","GCGCCCAUAGGGACGC"]
    energy_model = PyVRNA()
    fold_result  = energy_model.RNAcofold(sequences)
    assert [fold_result.structure, fold_result.energy] == ['.((.(((...))))).&((((((...))).)))', -10.25 + energy_model.dG_init_adjustment]
    energy_model = None
    sequences    = ["GCGCACAUAGUGACGC","GCGCCCAUAGGGACGC"]
    constraints  = ["....xx....xx....","....xx....xx...."]
    energy_model = PyVRNA(dangles=1, gquad=False, parameter_file='rna_andronescu2007.par')
    fold_result  = energy_model.RNAcofold(sequences,constraints)
    assert [fold_result.structure, fold_result.energy] == ['((((............&))))............', -5.869999885559082 + energy_model.dG_init_adjustment]
    energy_model = None
    print 'Good.'

    # Tests for RNAensemble
    print 'Testing RNAensemble...',
    sequence     = 'CGCAGGGAUACCCGCGGCGCCCAUAGGGACGCCGCAGGGAUACCCGCGGCGCCCAUAGGGACGC'
    energy_model = PyVRNA(dangles=2, gquad=False, parameter_file='rna_turner2004.par')
    ensemble_result = energy_model.RNAensemble(sequence)
    assert [ensemble_result.structure, ensemble_result.energy] == ['(((.(((...||||||((((((.,.{||.||||||.|||...))))))}||}|}.,.))).)))', -39.47269058227539]
    energy_model = None
    print 'Good.'

    print 'Testing ensemble: RNAensemble with aptamers...',
    # example 1
    sequence = "AGACAUAGCGAUCAAGUGAUACCAGCAUCGUCUUGAUGCCCUUGGCAGCACUUCAUAGCUAGCUACG"
    theophylline_aptamer  = "GAUACCAG&CCCUUGGCAGC"
    theophylline_constraint = "(...((((&)...)))...)"
    dG_theophylline = -9.22
    energy_model = PyVRNA(dangles=2, gquad=False, parameter_file='rna_turner2004.par')
    fold_result  = energy_model.RNAensemble(sequence,aptamers=[theophylline_aptamer],aptamer_constraints=[theophylline_constraint],dG_ligands=[dG_theophylline])
    assert [fold_result.structure, fold_result.energy] == [',,.,.((((.,{.(((((...((((((((,...,)))))...)))...))))).},.)))),,,...',-28.491722106933594]
    fold_result  = energy_model.RNAensemble(sequence)
    assert [fold_result.structure, fold_result.energy] == [',,.,.((((.,{.(((((...((((((((,...,)))))...)))...))))).},.)))),,,...',-19.407447814941406]
    energy_model = None
    
    # example 2
    sequence = "CCAGAAUCAUGCAAGUGCGUAAGAUAGUCGCGGGCUUUAAUUACGCUAUAUUAUAUACCCAAUUCU"
    streptavidin_aptamer =              'CCAGAAUCAUGCAAGUGCGUAAGAUAGUCGCGGG'
    streptavidin_aptamer_structure =    '((xxxxxxxxxxxxx((((xxxxxxxxx))))))'
    dG_streptavidin = -10.15
    energy_model = PyVRNA(dangles=2, gquad=False, parameter_file='rna_turner2004.par')
    fold_result  = energy_model.RNAensemble(sequence,aptamers=[streptavidin_aptamer],aptamer_constraints=[streptavidin_aptamer_structure],dG_ligands=[dG_streptavidin])
    assert [fold_result.structure, fold_result.energy] == ['{{.............{(((,,,...,,,|}}})),,....},,,,,....................',-11.105342864990234]
    fold_result  = energy_model.RNAensemble(sequence)
    assert [fold_result.structure, fold_result.energy] == ['.........,,.....((((((...((((...))))....))))))..,.................',-10.555535316467285]
    
    # example 3 - add multiple apatmers
    sequence = "AGACAUAGCGAUCAAGUGAUACCAGCAUCGUCUUGAUGCCCUUGGCAGCACUUCAUAGCUAGCCAGAAUCAUGCAAGUGCGUAAGAUAGUCGCGGGCUUUAAUUACGCUAUAUUAUAUACCCAAUUCU"
    energy_model = PyVRNA(dangles=2, gquad=False, parameter_file='rna_turner2004.par')
    fold_result  = energy_model.RNAensemble(sequence,aptamers=[theophylline_aptamer,streptavidin_aptamer], \
        aptamer_constraints=[theophylline_constraint,streptavidin_aptamer_structure],dG_ligands=[0.0,0.0]) # what if the aptamer provides no additional free energy?
    assert [fold_result.structure, fold_result.energy] == ['..,.((((((((((((((((......)))).))))))(((({(((((((........))).))))}.,,,,{((....}|}}..}},...,..))))........)))))).,,..............',-33.18623352050781]
    fold_result  = energy_model.RNAensemble(sequence,aptamers=[theophylline_aptamer,streptavidin_aptamer], \
        aptamer_constraints=[theophylline_constraint,streptavidin_aptamer_structure],dG_ligands=[dG_theophylline,dG_streptavidin]) # now with the aptamer binding free energies
    assert [fold_result.structure, fold_result.energy] == ['....((((((((((((((((......}))).))))))(((...))){((........)))((((.............((((.........)))))))).......)))))).,...............',-35.90715408325195]
    print 'Good.'

    # Tests for RNAeval
    print 'Testing RNAeval...',
    sequence     = 'CGACGUAGAUGCUAGCUGACUCGAUGC'
    energy_model = PyVRNA(dangles=0, gquad=False, parameter_file='rna_andronescu2007.par')
    energy  = energy_model.RNAeval([sequence], ['(((.(.((.......))..))))....'])
    assert float("{0:0.2f}".format(energy)) == 2.15
    energy_model = None
    sequences    = ['CGCAGGGAUACCCGCG','GCGCCCAUAGGGACGC']
    energy_model = PyVRNA(dangles=2, gquad=False, parameter_file='rna_andronescu2007.par')
    energy  = energy_model.RNAeval(sequences,['.((.(((...))))).','((((((...))).)))'])
    assert float("{0:0.2f}".format(energy)) == -12.72
    energy_model = None
    print 'Good.'

    # Tests for RNAfold
    print 'Testing RNAfold...',
    sequence     = "CGCAGGGAUACCCGCG"
    energy_model = PyVRNA(dangles=0, gquad=False, parameter_file='rna_andronescu2007.par')
    fold_result  = energy_model.RNAfold(sequence)
    assert [fold_result.structure, fold_result.energy] == ['(((.(((...))))))', -4.619999885559082]
    energy_model = None
    sequence     = "CGCAGGGAUACCCGCGGCGCCCAUAGGGACGC"
    constraints  = '.....xxx........................'
    energy_model = PyVRNA(dangles=2, gquad=False, parameter_file='rna_turner2004.par')
    fold_result  = energy_model.RNAfold(sequence, constraints)
    assert [fold_result.structure, fold_result.energy] == ['(((.(......).)))((((((...))).)))', -9.300000190734863]
    energy_model = None
    print 'Good.'

    print 'Testing mfe: RNAfold with aptamers...',
    # example 1
    sequence = "AGACAUAGCGAUCAAGUGAUACCAGCAUCGUCUUGAUGCCCUUGGCAGCACUUCAUAGCUAGCUACG"
    theophylline_aptamer  = "GAUACCAG&CCCUUGGCAGC"
    theophylline_constraint = "(...((((&)...)))...)"
    dG_theophylline = -9.22
    energy_model = PyVRNA(dangles=2, gquad=False, parameter_file='rna_turner2004.par')
    fold_result  = energy_model.RNAfold(sequence,aptamers=[theophylline_aptamer],aptamer_constraints=[theophylline_constraint],dG_ligands=[dG_theophylline])
    assert [fold_result.structure, fold_result.energy] == ['.....((((.((.(((((...((((((((.....)))))...)))...))))).)).))))......',-26.920000076293945]
    fold_result  = energy_model.RNAfold(sequence)
    assert [fold_result.structure, fold_result.energy] == ['.....((((.((.(((((...((((((((.....)))))...)))...))))).)).))))......',-17.700000762939453]
    energy_model = None
    
    # example 2
    sequence = "CCAGAAUCAUGCAAGUGCGUAAGAUAGUCGCGGGCUUUAAUUACGCUAUAUUAUAUACCCAAUUCU"
    streptavidin_aptamer =              'CCAGAAUCAUGCAAGUGCGUAAGAUAGUCGCGGG'
    streptavidin_aptamer_structure =    '((xxxxxxxxxxxxx((((xxxxxxxxx))))))'
    dG_streptavidin = -10.15
    energy_model = PyVRNA(dangles=2, gquad=False, parameter_file='rna_turner2004.par')
    fold_result  = energy_model.RNAfold(sequence,aptamers=[streptavidin_aptamer],aptamer_constraints=[streptavidin_aptamer_structure],dG_ligands=[dG_streptavidin])
    assert [fold_result.structure, fold_result.energy] == ['((.............((((.........))))))................................',-10.550000190734863]
    fold_result  = energy_model.RNAfold(sequence)
    assert [fold_result.structure, fold_result.energy] == ['................((((((...((((...))))....))))))....................',-9.600000381469727]
    
    # example 3 - add multiple apatmers
    sequence = "AGACAUAGCGAUCAAGUGAUACCAGCAUCGUCUUGAUGCCCUUGGCAGCACUUCAUAGCUAGCCAGAAUCAUGCAAGUGCGUAAGAUAGUCGCGGGCUUUAAUUACGCUAUAUUAUAUACCCAAUUCU"
    energy_model = PyVRNA(dangles=2, gquad=False, parameter_file='rna_turner2004.par')
    fold_result  = energy_model.RNAfold(sequence,aptamers=[theophylline_aptamer,streptavidin_aptamer], \
        aptamer_constraints=[theophylline_constraint,streptavidin_aptamer_structure],dG_ligands=[0.0,0.0]) # what if the aptamer provides no additional free energy?
    assert [fold_result.structure, fold_result.energy] == ['....((((((((((((((((......)))).))))))((((((((((((........))).))))).(((((((....))))..)))......))))........)))))).................',-31.399999618530273]
    fold_result  = energy_model.RNAfold(sequence,aptamers=[theophylline_aptamer,streptavidin_aptamer], \
        aptamer_constraints=[theophylline_constraint,streptavidin_aptamer_structure],dG_ligands=[dG_theophylline,dG_streptavidin]) # now with the aptamer binding free energies
    assert [fold_result.structure, fold_result.energy] == ['....((((((((((((((((......)))).))))))(((...)))(((........)))((((.............((((.........)))))))).......)))))).................',-34.95000076293945]
    print 'Good.'

    # Tests for RNAsubopt
    print 'Testing RNAsubopt...',
    sequences    = ["CGCAGGGAUACCCGCG","GCGCCCAUAGGGACGC"]
    constraints = ["..xx........xx..","..xx........xx.."]
    delta_energy = 10
    energy_model = PyVRNA(dangles=0, gquad=False, parameter_file='rna_andronescu2007.par')
    subopt_list = energy_model.RNAsubopt(sequences, constraints, delta_energy)
    assert len(subopt_list) == 10196
    energy_model = None
    sequence   = "CGCAGGGAUACCCGCG"
    delta_energy = 3
    energy_model = PyVRNA(dangles=2, gquad=False, parameter_file='rna_turner2004.par')
    assert energy_model.RNAsubopt([sequence], delta_energy=3) == [('(((.(((...))))))', -5.599999904632568),
                                                                ('(((.((....)).)))', -5.5),
                                                                ('.((.(((...))))).', -4.699999809265137),
                                                                ('.((.((....)).)).', -4.599999904632568),
                                                                ('(((.((.....)))))', -3.5),
                                                                ('.((.((.....)))).', -2.5999999046325684),
                                                                ('....(((...)))...', -2.5999999046325684)]
    print 'Good.'

    print 'Testing RNAsubopt with aptamers...',
    sequence = "AGACAUAGCGAUCAAGUGAUACCAGCAUCGUCUUGAUGCCCUUGGCAGCACUUCAUAGCUAGCUACGAUGUAGCUCGGUAUUAUUU"
    theophylline_aptamer  = "GAUACCAG&CCCUUGGCAGC"
    theophylline_constraint = "(...((((&)...)))...)"
    dG_theophylline = -9.22
    energy_model = PyVRNA(dangles=0, gquad=False, parameter_file='rna_andronescu2007.par')
    subopt_result = energy_model.RNAsubopt([sequence], delta_energy=10, aptamers=[theophylline_aptamer], aptamer_constraints=[theophylline_constraint], dG_ligands=[dG_theophylline])
    assert subopt_result[0].structure == '......(((.((.(((((......(((((.....))))).........))))).)).)))(((((.....)))))...........'
    assert subopt_result[0].energy == -30.719999313354492
    # Test without aptamer free energy to confirm different predictions
    # subopt_result = energy_model.RNAsubopt([sequence], delta_energy=10)
    # assert subopt_result[0].structure == '......(((.((.(((((...((((((((.....)))))...)))...))))).)).)))(((((.....)))))...........'
    # assert subopt_result[0].energy == -21.5
    energy_model = None
    print 'Good.'


    # Tests for vienna2bp and bp2vienna
    print 'Testing vienna2bp and bp2vienna...',
    energy_model  = PyVRNA()
    vienna_string = '..++++....++++....++++....++++....++++....++++....++++....++++..'
    assert energy_model.bptuple2vienna(energy_model.vienna2bp(vienna_string)) == vienna_string
    energy_model  = None
    energy_model  = PyVRNA()
    vienna_string = '....(((...)))...'
    assert energy_model.bptuple2vienna(energy_model.vienna2bp(vienna_string)) == vienna_string
    energy_model  = None
    energy_model  = PyVRNA()
    vienna_string = '(((.((....)).)))'
    assert energy_model.bptuple2vienna(energy_model.vienna2bp(vienna_string)) == vienna_string
    energy_model  = None
    energy_model  = PyVRNA()
    vienna_string = '..[.((....)).]..+++..+++..+++..+++..'
    assert energy_model.bptuple2vienna(energy_model.vienna2bp(vienna_string)) == vienna_string
    energy_model  = None
    energy_model  = PyVRNA()
    vienna_string = '(((((((....++..++..++..++..)))))))'
    assert energy_model.bptuple2vienna(energy_model.vienna2bp(vienna_string)) == vienna_string
    energy_model  = None
    energy_model  = PyVRNA()
    vienna_string = '[[[[(((((((....)))))))..((((....))))]]]]'
    assert energy_model.bptuple2vienna(energy_model.vienna2bp(vienna_string)) == vienna_string
    energy_model  = None
    energy_model  = PyVRNA()
    vienna_string = '...((..((..))..((..))..))....((..((..))..(.)..))..'
    assert energy_model.bptuple2vienna(energy_model.vienna2bp(vienna_string)) == vienna_string
    energy_model  = None
    print 'Good.'

def legacy_tests():
    """
    Tests each function in ViennaRNA class
    """
    print ""
    print "### === Testing ViennaRNA legacy class === ###"

    # Tests for centroid
    print 'Testing centroid...',
    sequence      = 'CGACGUAGAUGCUAGCUGACUCGAUGC'
    energy_model  = PyVRNA(dangles=2, gquad=False, parameter_file='rna_turner2004.par')
    vienna_model  = ViennaRNA(Sequence_List=[sequence], material='rna_turner2004', Gquad=False)
    cntrd_result  = energy_model.RNAcentroid(sequence)
    vienna_model.centroid(strands=[1], constraint=None, Temp=37.0, dangles="all", outputPS=False)   
    assert [cntrd_result.structure, cntrd_result.energy] == [vienna_model['Centroid_bracket_string'], vienna_model["Centroid_energy"][0]] == ['(((.(.((.......))..))))....', 1.399999976158142]
    energy_model  = None
    vienna_model  = None
    sequence      = 'CGCAGGGAUACCCGCG' + 'GCGCCCAUAGGGACGC'
    energy_model  = PyVRNA(dangles=2, gquad=False, parameter_file='rna_turner2004.par')
    vienna_model  = ViennaRNA(Sequence_List=[sequence], material='rna_turner2004', Gquad=False)
    cntrd_result  = energy_model.RNAcentroid(sequence)
    vienna_model.centroid(strands=[1], constraint=None, Temp=37.0, dangles="all", outputPS=False)
    assert [cntrd_result.structure, cntrd_result.energy] == [vienna_model['Centroid_bracket_string'], vienna_model["Centroid_energy"][0]] == ['(((.(((...))))))((((((...))).)))', -13.5]
    energy_model  = None
    vienna_model  = None
    print 'Good.'

    # Tests for energy
    print 'Testing energy...',
    sequence     = 'CGACGUAGAUGCUAGCUGACUCGAUGC'
    structure    = '(((.(.((.......))..))))....'
    energy_model = PyVRNA(dangles=0, gquad=False, parameter_file='rna_andronescu2007.par')
    vienna_model = ViennaRNA(Sequence_List=[sequence], material='rna_andronescu2007', Gquad=False)
    bp_tuple    = energy_model.vienna2bp(structure)
    vienna_model.energy(strands=[1], base_pairing_x=bp_tuple.bpx, base_pairing_y=bp_tuple.bpy, Temp=37.0, dangles="none")
    assert float("{0:0.2f}".format(energy_model.RNAeval([sequence], [structure]))) == float("{0:0.2f}".format(vienna_model["energy_energy"][0])) == 2.15
    energy_model = None
    vienna_model = None
    sequences    = ['CGCAGGGAUACCCGCG','GCGCCCAUAGGGACGC']
    structures  = ['.((.(((...))))).','((((((...))).)))']
    energy_model = PyVRNA()
    vienna_model = ViennaRNA(Sequence_List=sequences, material='rna_andronescu2007', Gquad=False)
    bp_tuple    = energy_model.vienna2bp("".join(structures))
    vienna_model.energy(strands=[0, 1], base_pairing_x=bp_tuple.bpx, base_pairing_y=bp_tuple.bpy, Temp=37.0, dangles="all")
    assert float("{0:0.2f}".format(energy_model.RNAeval(sequences, structures))) == float("{0:0.2f}".format(vienna_model["energy_energy"][0])) == -12.72
    energy_model = None
    vienna_model = None
    print 'Good.'

    # Tests for mfe
    print 'Testing mfe: RNAfold...',
    sequence     = "CGCAGGGAUACCCGCG"
    energy_model = PyVRNA(dangles=0, gquad=False, parameter_file='rna_andronescu2007.par')
    vienna_model = ViennaRNA(Sequence_List=[sequence], material='rna_andronescu2007', Gquad=False)
    vienna_model.mfe(strands=[1], constraints=None, Temp=37.0, dangles="none", outputPS=False, duplex=False)
    assert list(energy_model.RNAfold(sequence)) == [vienna_model["mfe_bracket_string"], vienna_model["mfe_energy"][0]] == ['(((.(((...))))))', -4.619999885559082]
    energy_model = None
    vienna_model = None
    sequence     = "CGCAGGGAUACCCGCGGCGCCCAUAGGGACGC"
    constraints  = '.....xxx........................'
    energy_model = PyVRNA(dangles=2, gquad=False, parameter_file='rna_turner2004.par')
    vienna_model = ViennaRNA(Sequence_List=[sequence], material='rna_turner2004', Gquad=False)
    vienna_model.mfe(strands=[1], constraints=constraints, Temp=37.0, dangles="all", outputPS=False, duplex=False)
    assert list(energy_model.RNAfold(sequence, constraints)) == [vienna_model["mfe_bracket_string"], vienna_model["mfe_energy"][0]] == ['(((.(......).)))((((((...))).)))', -9.300000190734863]
    energy_model = None
    vienna_model = None
    print 'Good.'
    
    print 'Testing mfe: RNAcofold...',
    sequences    = ["CGCAGGGAUACCCGCG","GCGCCCAUAGGGACGC"]
    energy_model = PyVRNA()
    vienna_model = ViennaRNA(Sequence_List=sequences, material='rna_andronescu2007', Gquad=False)
    vienna_model.mfe(strands=[0, 1], constraints=None, Temp=37.0, dangles="all", outputPS=False, duplex=False)
    assert list(energy_model.RNAcofold(sequences)) == [vienna_model["mfe_bracket_string"], vienna_model["mfe_energy"][0]] == ['.((.(((...))))).&((((((...))).)))', -10.25 + energy_model.dG_init_adjustment]
    energy_model = None
    vienna_model = None
    sequences    = ["GCGCACAUAGUGACGC","GCGCCCAUAGGGACGC"]
    constraints  = ["....xx....xx....","....xx....xx...."]
    energy_model = PyVRNA(dangles=1, gquad=False, parameter_file='rna_andronescu2007.par')
    vienna_model = ViennaRNA(Sequence_List=sequences, material='rna_andronescu2007', Gquad=False)
    vienna_model.mfe(strands=[0, 1], constraints=constraints[0]+constraints[1], Temp=37.0, dangles="some", outputPS=False, duplex=False)
    assert list(energy_model.RNAcofold(sequences, constraints)) == [vienna_model["mfe_bracket_string"], vienna_model["mfe_energy"][0]] == ['((((............&))))............', -5.869999885559082 + energy_model.dG_init_adjustment]
    energy_model = None
    vienna_model = None
    print 'Good.'

    # Tests for RNAsubopt
    print 'Testing subopt...',
    sequences    = ["CGCAGGGAUACCCGCG", "GCGCCCAUAGGGACGC"]
    constraints = ["..xx........xx..", "..xx........xx.."]
    delta_energy = 10
    energy_model = PyVRNA(dangles=0, gquad=False, parameter_file='rna_andronescu2007.par')
    vienna_model = ViennaRNA(Sequence_List=sequences, material='rna_andronescu2007', Gquad=False)
    vienna_model.subopt(strands=[0, 1], energy_gap=delta_energy, Temp=37.0, dangles="none", constraints="".join(constraints), outputPS=False)
    assert len(energy_model.RNAsubopt(sequences, constraints, delta_energy=delta_energy)) == len(vienna_model["subopt_energy"]) == 10196
    energy_model = None
    vienna_model = None
    sequence     = "CGCAGGGAUACCCGCG"
    delta_energy = 3
    energy_model = PyVRNA(dangles=2, gquad=False, parameter_file='rna_turner2004.par')
    vienna_model = ViennaRNA(Sequence_List=[sequence], material='rna_turner2004', Gquad=False)
    vienna_model.subopt(strands=[1], energy_gap=delta_energy, Temp=37.0, dangles="all", constraints=None, outputPS=False)
    assert map(lambda x: x[1], energy_model.RNAsubopt([sequence], delta_energy=delta_energy)) == vienna_model["subopt_energy"] == [-5.599999904632568,
                                                                                                                                  -5.5,
                                                                                                                                  -4.699999809265137,
                                                                                                                                  -4.599999904632568,
                                                                                                                                  -3.5,
                                                                                                                                  -2.5999999046325684,
                                                                                                                                  -2.5999999046325684]
    energy_model = None
    vienna_model = None
    print 'Good.'

    # Tests for convert_bracket_to_numbered_pairs and convert_numbered_pairs_to_bracket
    print 'Testing convert_bracket_to_numbered_pairs ...'
    print 'Testing convert_bracket_to_numbered_pairs ...',
    vienna_model  = ViennaRNA(Sequence_List=[], material='rna_turner2004', Gquad=False)
    vienna_string = '................................................................'
    assert vienna_model.convert_numbered_pairs_to_bracket(*vienna_model.convert_bracket_to_numbered_pairs(vienna_string)) == vienna_string
    vienna_model  = None
    vienna_model  = ViennaRNA(Sequence_List=[], material='rna_turner2004', Gquad=False)
    vienna_string = '....(((...)))...'
    assert vienna_model.convert_numbered_pairs_to_bracket(*vienna_model.convert_bracket_to_numbered_pairs(vienna_string)) == vienna_string
    vienna_model  = None
    vienna_model  = ViennaRNA(Sequence_List=[], material='rna_turner2004', Gquad=False)
    vienna_string = '(((.((....)).)))'
    assert vienna_model.convert_numbered_pairs_to_bracket(*vienna_model.convert_bracket_to_numbered_pairs(vienna_string)) == vienna_string
    vienna_model  = None
    vienna_model  = ViennaRNA(Sequence_List=[], material='rna_turner2004', Gquad=False)
    vienna_string = '..[.((....)).]......................'
    assert vienna_model.convert_numbered_pairs_to_bracket(*vienna_model.convert_bracket_to_numbered_pairs(vienna_string)) == vienna_string
    vienna_model  = None
    vienna_model  = ViennaRNA(Sequence_List=[], material='rna_turner2004', Gquad=False)
    vienna_string = '(((((((....................)))))))'
    assert vienna_model.convert_numbered_pairs_to_bracket(*vienna_model.convert_bracket_to_numbered_pairs(vienna_string)) == vienna_string
    vienna_model  = None
    vienna_model  = ViennaRNA(Sequence_List=[], material='rna_turner2004', Gquad=False)
    vienna_string = '[[[[(((((((....)))))))..((((....))))]]]]'
    assert vienna_model.convert_numbered_pairs_to_bracket(*vienna_model.convert_bracket_to_numbered_pairs(vienna_string)) == vienna_string
    vienna_model  = None
    vienna_model  = ViennaRNA(Sequence_List=[], material='rna_turner2004', Gquad=False)
    vienna_string = '...((..((..))..((..))..))....((..((..))..(.)..))..'
    assert vienna_model.convert_numbered_pairs_to_bracket(*vienna_model.convert_bracket_to_numbered_pairs(vienna_string)) == vienna_string
    vienna_model  = None
    print 'Good.'

if __name__ == '__main__':
    tests()
    legacy_tests()
    # model = PyVRNA(temperature=37.0,dangles=0,gquad=False,parameter_file="rna_turner2004.par")
    # model.RNAsubopt(sequences=['AAGGGCGCGAGC', 'ACCUCCUUA'],constraints=['............','.........'],delta_energy=8.0)
    # model.RNAsubopt(sequences=['AAGGGCGCGAGC', 'ACCUCCUUA'],constraints=[],delta_energy=8.0)

