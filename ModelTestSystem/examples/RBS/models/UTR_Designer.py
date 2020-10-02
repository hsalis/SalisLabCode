#from NuPACK import NuPACK
from ViennaRNA import ViennaRNA
RNAEnergyModel = ViennaRNA

import re
import math

class Class_RNA(object):
    """For bacteria: contains ribosome binding sites, anti-sense RNAs, RNAse binding sites, and terminators. For eukaryotes: contains ribosomal pause sites, microRNAs, RNAse binding sites, exons, introns, and terminators"""
    def __repr__(self):
        # output="RNA Sequence Information:\n"
        # for RBS in self.RBS_list:
            # output = output + RBS.__repr__()
        # return output
        return "RBS object here (too long!)"
    pass

class Class_Ribosome_Binding_Site(object):
    """Information about the RBS surrounding a start codon in bacteria, including its interactions with the Ribosome or other RNAs."""
    def __repr__(self):
        # output = "Ribosome_Binding_Site Instance Information:\n"
        # for attr in vars(self):
            # value = getattr(self,attr)
            # output += "Attribute: " + attr + "\t Value: " + str(value) + "\n"
        # return output
        return "RBS object here (too long!)"
    pass

class CalcError(Exception):
    """Base class for exceptions in this module."""

    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)

class RBS_Calculator(RNAEnergyModel):

    #From experimental characterization of ALL Predicted RBSs as of 1/20/08
    RT_eff = 2.222
    logK =  7.824
    K = 2500.0

    #Global parameters -- constants
    infinity = 1e12 #For all practical purposes, here.
    RNA_model = "rna37"
    start_codon_energies = {"ATG":-1.194, "AUG": -1.194, "GTG": -0.0748, "GUG": -0.0748, "TTG":-0.0435, "UUG": -0.0435, "CTG": -0.03406, "CUG":-0.03406} #hybridization to CAT
    auto_dangles = True
    dangles_default = "all"
    temp = 37.0
    optimal_spacing = 5 #aligned spacing

    dG_spacing_constant_push = [12.2, 2.5, 2.0, 3.0] #[0.46, 0.0, -0.46]
    dG_spacing_constant_pull = [0.048, 0.24, 0.0] #[0.1, -0.08, 0.0]
    cutoff = 35 #number of nt +- start codon considering for folding
    standby_site_length = 4 #Number of nt before SD sequence that must be unpaired for ribosome binding
    energy_cutoff = 3.0
    start_codons = ["ATG", "AUG", "GTG", "GUG","TTG","UUG","CTG","CUG"] #substituted U for T in actual calcs
    rRNA = "acctcctta" #"acctcctta" #"ggatcacctcctta" #acctcctta is the typical anti-SD sequence

    footprint = 1000

    Warning_List = {'NONE': ('green','OK','No warning issued.'),
                        'SHORT_CDS': ('red', 'CDS','The protein coding sequence is too short.'),
                        'NON_EQUILIBRIUM': ('orange', 'NEQ','The thermodynamic model''s assumptions may not be valid.'),
                        'CLOSE_START_CODONS': ('orange', 'OLS','The start codons are too close together for an accurate prediction.'),
                        }

    def __init__(self, mRNA = None, start_range = None, rRNA_seq = None, name = "Unnamed mRNA", verbose = False):

        #RNAEnergyModel.__init__(self,sequences,self.RNA_model)

        exp = re.compile('[ATGCU]',re.IGNORECASE)
        if exp.match(mRNA) == None:
            raise ValueError("Invalid letters found in sequence mRNA ""%s"". Only ATGCU accepted." % mRNA)
            
        if rRNA_seq is not None and exp.match(rRNA_seq) == None:
            raise ValueError("Invalid letters found in sequence rRNA ""%s"". Only ATGCU accepted." % rRNA_seq)

        if start_range[0] < 0: start_range[0] = 0
        if start_range[1] > len(mRNA): start_range[1] = len(mRNA)

        self.name = name
        self.mRNA_input = mRNA.upper()
        if rRNA_seq is not None: self.rRNA = rRNA_seq
        
        self.rRNA_len = len(self.rRNA)
        self.mRNA_len = len(self.mRNA_input)
        self.total_sequence_length = len(mRNA) + len(self.rRNA)
        self.dG_rRNA = self.calc_dG_rRNA()
        self.has_run = 0
        self.start_range = start_range
        self.verbose = verbose
        

    def find_min(self,input_list):

        min_item = self.infinity
        min_index = 0

        for i, item in enumerate(input_list):
            if item < min_item:
                min_item = item
                min_index = i

        return (min_item,min_index)


    def find_start_codons(self,sequence):

        self.start_position_list = []
        self.start_codon_list = []

        seq_len = len(sequence)
        end = min(self.start_range[1],seq_len-2)
        begin = min(self.start_range[0],end)

        for i in range(begin,end+1):
            codon = sequence[i:i+3]
            if codon.upper() in self.start_codons:
                self.start_position_list.append(i)
                self.start_codon_list.append(codon)
                yield (i,codon)
            else:
                pass

    def calc_aligned_spacing(self,mRNA,start_pos,bp_x,bp_y):

        #rRNA is the concatenated at the end of the sequence in 5' to 3' direction
        #first: identify the farthest 3' nt in the rRNA that binds to the mRNA and return its mRNA base pairer

        Ok = False
        seq_len = len(mRNA) + self.rRNA_len
        for (rRNA_nt) in range(seq_len,seq_len - self.rRNA_len,-1):

            if rRNA_nt in bp_y:
                rRNA_pos = bp_y.index(rRNA_nt)
                if bp_x[rRNA_pos] < start_pos:
                    Ok = True
                    farthest_3_prime_rRNA = rRNA_nt - len(mRNA)

                    mRNA_nt = bp_x[rRNA_pos]
                    distance_to_start = start_pos - mRNA_nt + 1 #start_pos is counting starting from 0 (python)

                    #print "farthest 3' rRNA = ", farthest_3_prime_rRNA
                    #print "mRNA_nt = ", mRNA_nt
                    #print "start_pos = ", start_pos
                    #print "distance to start = ", distance_to_start


                    break
                else:
                    break

        #second: the aligned spacing is the difference between the "start distance", which is the farthest 5' mRNA bound by rRNA, and the farthest 3' prime rRNA bound.
        if Ok:
            aligned_spacing = distance_to_start - farthest_3_prime_rRNA
        else:
            aligned_spacing = self.infinity

        return aligned_spacing

    def calc_dG_spacing(self, aligned_spacing):
        #Relationship between aligned spacing and dG penalty for suboptimal spacing
        #Determined from experiments in terms of kT units

        if (aligned_spacing < self.optimal_spacing):
            ds = aligned_spacing - self.optimal_spacing

            dG_spacing_penalty = self.dG_spacing_constant_push[0] / (1.0 + math.exp(self.dG_spacing_constant_push[1]*(ds + self.dG_spacing_constant_push[2] )))**self.dG_spacing_constant_push[3]

        else:
            ds = aligned_spacing - self.optimal_spacing
            dG_spacing_penalty = self.dG_spacing_constant_pull[0] * ds * ds + self.dG_spacing_constant_pull[1] * ds + self.dG_spacing_constant_pull[2]

        return dG_spacing_penalty

    def calc_dG_mRNA_rRNA(self,start_pos):
        
        begin = max(0,start_pos-self.cutoff)
        mRNA_len = min(len(self.mRNA_input),start_pos+self.cutoff)
        start_pos_in_subsequence = min(start_pos, self.cutoff)
        startpos_to_end_len = mRNA_len - start_pos_in_subsequence - begin

        #First, identify the rRNA-binding site
        #Constraints: the entire rRNA-binding site must be upstream of the start codon
        mRNA = self.mRNA_input[begin:start_pos]
        
        if len(mRNA) == 0:
            raise CalcError("Warning: There is a leaderless start codon, which is being ignored.")
        
        fold = RNAEnergyModel([mRNA,self.rRNA],material = self.RNA_model)
        fold.subopt([1, 2],self.energy_cutoff,dangles = self.dangles, Temp = self.temp)
        
        if len(fold["subopt_basepairing_x"]) == 0:
            raise CalcError("Warning: The 16S rRNA has no predicted binding site. Start codon is considered as leaderless and ignored.")

        #Calculate dG_spacing

        #First, calculate aligned spacing of rRNA:mRNA
        aligned_spacing = []
        for (bp_x, bp_y) in zip(fold["subopt_basepairing_x"], fold["subopt_basepairing_y"]):
            aligned_spacing.append(self.calc_aligned_spacing(mRNA, start_pos_in_subsequence, bp_x,bp_y))

        dG_spacing_list = []
        dG_mRNA_rRNA = []
        dG_mRNA_rRNA_withspacing = []

        for (counter) in range(len(fold["subopt_basepairing_x"])):

            dG_mRNA_rRNA.append(fold["subopt_energy"][counter])
            val = self.calc_dG_spacing(aligned_spacing[counter])
            dG_spacing_list.append(val)
            dG_mRNA_rRNA_withspacing.append(val + fold["subopt_energy"][counter])
        
        [dG_mRNA_rRNA_folding, index] = self.find_min(dG_mRNA_rRNA_withspacing)
        dG_spacing_final = dG_spacing_list[index]

        dG_mRNA_rRNA_nospacing = dG_mRNA_rRNA[index]

        #Is the dG spacing large compared to the energy gap?
        if dG_spacing_final > self.energy_cutoff:
            if self.verbose: print "Warning: The spacing penalty is greater than the energy gap. dG (spacing) = ", dG_spacing_final

        #Now identify where the rRNA bound
        most_5p_mRNA = start_pos #start_pos_in_subsequence
        most_3p_mRNA = 0

        #List of rRNA-mRNA base paired nucleotides
        bp_x_target = []
        bp_y_target = []

        bp_x = fold["subopt_basepairing_x"][index]
        bp_y = fold["subopt_basepairing_y"][index]
        
        for (nt_x, nt_y) in zip(bp_x, bp_y):
            if nt_y > len(mRNA): #nt is rRNA
                most_5p_mRNA = min(most_5p_mRNA, bp_x[bp_y.index(nt_y)])
                most_3p_mRNA = max(most_3p_mRNA, bp_x[bp_y.index(nt_y)])
                bp_x_target.append(nt_x)
                bp_y_target.append(nt_y)

        #The rRNA-binding site is between the nucleotides at positions most_5p_mRNA and most_3p_mRNA

        #Now, fold the pre-sequence, rRNA-binding-sequence and post-sequence separately and put them together.

        #We postulate that not all of the post-sequence can form secondary structures.

        #We define a post-sequence window. The post-sequence mRNA inside the window is allowed to form secondary structures. Otherwise, it is not.
        #The length of the dsRNA-allowed window is RBS_Calculator.post_window_length
        #The window begins at most_5p_mRNA and extends in the 5' direction for post_window_length nucleotides

        #We also define a pre-sequence window. The pre-sequence mRNA inside the window is allowed to form secondary structures. Otherwise, it is not.
        #The length of the dsRNA-allowed window is RBS_Calculator.pre_window_length
        #The window begins at most_3p_mRNA and extends in the 3' direction for pre_window_length nucleotides

        #Alternate hypothesis: The ribosome footprint extends from the most_5p_mRNA to +19 nucleotides, where +1 is the beginning of the start codon. The ribosome prevents any secondary structure from forming. This is equivalent to having a folding window from +20 nt to the end of the mRNA.

        #print "begin = ", begin
        #print "most_5p_mRNA = ", most_5p_mRNA
        
        #Extrema case: sometimes the mRNA binds to itself and not the rRNA, when forming the mRNA:rRNA complex. This code below catches that error.
        if len(bp_x_target) == 0:
        
            total_energy_withspacing = 0.0  #Set the energy to zero if the rRNA does not bind to mRNA.
            if self.verbose: print "Error in RBS Calculator v1.0 code: minimization did not identify the 16S rRNA binding site on the mRNA. Setting energy to zero. v1.1 and later do not have this problem."
            structure = fold
            structure["program"] = "subopt"
            structure["mRNA"] = mRNA
            structure["MinStructureID"] = 0
            structure["dG_mRNA_rRNA"] = 0.0
            structure["dG_SD_aSD"] = 0.0
            structure["dG_mRNA_rRNA_withspacing"] = 0.0
            structure["dG_spacing"] = 0.0
            structure["subopt_energy"] = [0.0]
            structure["subopt_basepairing_x"] = []
            structure["subopt_basepairing_y"] = []
            structure["subopt_composition"] = [1, 2]
            structure["bp_x"] = []
            structure["bp_y"] = []

            return (total_energy_withspacing, structure)
        
        
        mRNA_pre = self.mRNA_input[begin:begin+most_5p_mRNA-1]
        post_window_end = mRNA_len + 1
        post_window_begin = min(start_pos + self.footprint,post_window_end) #Footprint
        post_window_end = mRNA_len + 1
        mRNA_post = self.mRNA_input[post_window_begin:post_window_end]

        mRNA_pre_len = len(mRNA_pre)
        mRNA_post_len = len(mRNA_post)
        mRNA_rRNA_binding_len = most_3p_mRNA - most_5p_mRNA + 1
        total_folded_len = mRNA_pre_len + mRNA_post_len + mRNA_rRNA_binding_len

        total_bp_x = []
        total_bp_y = []

        #Calculate pre-sequence folding -- using pre_window_length
        if len(mRNA_pre) > 0:
            fold_pre = RNAEnergyModel([mRNA_pre], material = self.RNA_model)
            fold_pre.mfe([1], dangles = self.dangles, Temp = self.temp)
            bp_x_pre = fold_pre["mfe_basepairing_x"][0]
            bp_y_pre = fold_pre["mfe_basepairing_y"][0]

        else:
            bp_x_pre = []
            bp_y_pre = []

        #Add pre-sequence base pairings to total base pairings
        offset = 0 #Begins at 0
        for (nt_x, nt_y) in zip(bp_x_pre, bp_y_pre):
            total_bp_x.append(nt_x + offset)
            total_bp_y.append(nt_y + offset)

        #print "pre bp_x = ", [x + offset for x in bp_x_pre]
        #print "pre bp_y = ", [x + offset for x in bp_y_pre]
        
        #Add rRNA-binding site base pairings to total base pairings
        
        # Calculate dG_SD:aSD - for UTR Deisgner
        # Add rRNA-binding site base pairings
        dG_SD_bp_x = []
        dG_SD_bp_y = []
        
        offset = 0 #Begins at zero
        if startpos_to_end_len < self.cutoff:
            rRNA_offset = startpos_to_end_len
        else:
            rRNA_offset = startpos_to_end_len

        for (nt_x, nt_y) in zip(bp_x_target, bp_y_target):
            total_bp_x.append(nt_x + offset)
            total_bp_y.append(nt_y + rRNA_offset)
            
            #Added for UTR Designer
            dG_SD_bp_x.append(nt_x + offset)
            dG_SD_bp_y.append(nt_y + rRNA_offset)
        
        #print dG_SD_bp_x
        #print dG_SD_bp_y
        
        #print "target bp_x = ", [x +offset for x in bp_x_target]
        #print "target bp_y = ", [x + rRNA_offset for x in bp_y_target]
        
        #Calculate post-sequence folding -- using post_window_length
        if len(mRNA_post) > 0:
            fold_post = RNAEnergyModel([mRNA_post], material = self.RNA_model)
            fold_post.mfe([1], dangles = self.dangles, Temp = self.temp)
            bp_x_post = fold_post["mfe_basepairing_x"][0]
            bp_y_post = fold_post["mfe_basepairing_y"][0]
        else:
            bp_x_post = []
            bp_y_post = []
            
        offset = post_window_begin - begin
        for (nt_x, nt_y) in zip(bp_x_post, bp_y_post):
            total_bp_x.append(nt_x + offset)
            total_bp_y.append(nt_y + offset)

        #print "post bp_x = ", [x + offset for x in bp_x_post]
        #print "post bp_y = ", [x + offset for x in bp_y_post]
            
        mRNA = self.mRNA_input[begin:mRNA_len]
        fold = RNAEnergyModel([mRNA, self.rRNA], material = self.RNA_model)

        #print "total_bp_x = ", total_bp_x
        #print "total_bp_y = ", total_bp_y
        
        total_energy = fold.energy([1, 2], total_bp_x, total_bp_y, Temp = self.temp, dangles = self.dangles)
        
        #Calculate dG_SD:aSD for UTR Designer
        dG_SD_aSD = fold.energy([1,2], dG_SD_bp_x, dG_SD_bp_y, Temp = self.temp, dangles = self.dangles)
        #print dG_SD_aSD
        
        #Compare with and without using windows
        #print "total dG (using windows) = ", total_energy
        #fold.export_PDF(0, "rRNA:mRNA (using windows)", filename =  "Test" + "_rRNA_windows_2" + ".pdf",program="energy")

        energy_nowindows = dG_mRNA_rRNA_nospacing
        #fold.export_PDF(index, "rRNA:mRNA (no windows)", filename =  "Test" + "_rRNA_nowindows" + ".pdf")

        #print "total dG (no windows) = ", energy_nowindows

        total_energy_withspacing = total_energy + dG_spacing_final

        structure = fold
        structure["program"] = "subopt"
        structure["mRNA"] = mRNA
        structure["MinStructureID"] = 0
        structure["dG_mRNA_rRNA"] = total_energy
        structure["dG_SD_aSD"] = dG_SD_aSD
        structure["dG_mRNA_rRNA_withspacing"] = total_energy_withspacing
        structure["dG_spacing"] = dG_spacing_final
        structure["subopt_energy"] = [total_energy_withspacing]
        structure["subopt_basepairing_x"] = [total_bp_x]
        structure["subopt_basepairing_y"] = [total_bp_y]
        structure["subopt_composition"] = [1, 2]
        structure["bp_x"] = total_bp_x
        structure["bp_y"] = total_bp_y

        return (total_energy_withspacing, structure)

    def calc_dG_standby_site(self,structure_old, rRNA_binding = True):
        #Set rRNA_binding to False if you want to calculate the structure of the unfolding of the standby site without the rRNA binding

        import copy
        structure = copy.deepcopy(structure_old)
        mRNA = structure["mRNA"]
        bp_x = structure["bp_x"]
        bp_y = structure["bp_y"]
        energy_before = structure["dG_mRNA_rRNA"] #without spacing effects

#        print "bp_x = ", bp_x
#        print "bp_y = ", bp_y
        
        #Identify the most 5p mRNA nt that is bound to rRNA
        most_5p_mRNA = 0
        for (nt_x, nt_y) in zip(bp_x, bp_y):
            if nt_x <= len(mRNA) and nt_y > len(mRNA): #nt_x is mRNA, nt_y is rRNA, they are bound.
                most_5p_mRNA = nt_x #starts counting from 0
                break
    
        #To calculate the mfe structure while disallowing base pairing at the standby site, we extract an mRNA subsequence from the beginning of the mRNA to the beginning of the standby site. We fold this subsequence and use these base pairings to replace the ones that occurred before the end of the standby site

        #Extract the base pairings that are 3' of the most_5p_mRNA base pairing
        bp_x_3p = []
        bp_y_3p = []
        for (nt_x, nt_y) in zip(bp_x, bp_y):
            if nt_x >= most_5p_mRNA:
                bp_x_3p.append(nt_x)
                bp_y_3p.append(nt_y)

        #Create the mRNA subsequence
        mRNA_subsequence = mRNA[0:max(0,most_5p_mRNA - self.standby_site_length - 1)]
        standby_site = mRNA[most_5p_mRNA - self.standby_site_length - 1:most_5p_mRNA]

        #Fold it and extract the base pairings
        if (len(mRNA_subsequence)) > 0:
            fold = RNAEnergyModel([mRNA_subsequence], material = self.RNA_model)
            fold.mfe([1], dangles = self.dangles, Temp = self.temp)
            energy_after_5p = fold["mfe_energy"][0]
            bp_x_5p = fold["mfe_basepairing_x"][0]   #[0] added 12/13/07
            bp_y_5p = fold["mfe_basepairing_y"][0]
        else:
             bp_x_5p = []
             bp_y_5p = []
             energy_after_5p = 0.0

        #Put the two sets of base pairings together
        bp_x_after = []
        bp_y_after = []
        for (nt_x, nt_y) in zip(bp_x_5p, bp_y_5p):
            bp_x_after.append(nt_x)
            bp_y_after.append(nt_y)

        for (nt_x, nt_y) in zip(bp_x_3p, bp_y_3p):
            bp_x_after.append(nt_x)
            bp_y_after.append(nt_y)

        #Calculate its energy
        fold = RNAEnergyModel([mRNA, self.rRNA], material = self.RNA_model)
        energy_after = fold.energy([1, 2], bp_x_after, bp_y_after, dangles = self.dangles, Temp = self.temp)

        dG_standby_site = energy_before - energy_after

        if (dG_standby_site > 0.0): dG_standby_site = 0.0

        index = structure["MinStructureID"]
        structure["bp_x"] = bp_x_after
        structure["bp_y"] = bp_y_after
        structure["subopt_basepairing_x"][index] = bp_x_after
        structure["subopt_basepairing_y"][index] = bp_y_after
        structure["subopt_energy"][index] = energy_after
        structure["dG_mRNA_rRNA_corrected"] = energy_after

        return (dG_standby_site, structure)

    def calc_dG_mRNA(self,start_pos):

        mRNA = self.mRNA_input[max(0,start_pos-self.cutoff):min(len(self.mRNA_input),start_pos+self.cutoff)]
        fold = RNAEnergyModel([mRNA],self.RNA_model)
        fold.mfe([1], Temp = self.temp, dangles = self.dangles)

        structure = fold
        structure["mRNA"] = mRNA
        structure["bp_x"] = fold["mfe_basepairing_x"][0]
        structure["bp_y"] = fold["mfe_basepairing_y"][0]
        structure["dG_mRNA"] = fold["mfe_energy"][0]
        structure["MinStructureID"] = 0

        dG_mRNA_folding = fold["mfe_energy"][0]

        return (dG_mRNA_folding, structure)

    # ================================================================================================
    # UTR Designer replacement free energy terms, dG_direct and dG_indirect
    def calc_dG_direct(self,start_pos):
        
        mRNA = self.mRNA_input[max(0,start_pos-25):min(len(self.mRNA_input),start_pos+18)]
        fold = RNAEnergyModel([mRNA],self.RNA_model)
        fold.mfe([1], Temp = self.temp, dangles = self.dangles)

        structure = fold
        structure["mRNA"] = mRNA
        structure["bp_x"] = fold["mfe_basepairing_x"][0]
        structure["bp_y"] = fold["mfe_basepairing_y"][0]
        structure["dG_mRNA"] = fold["mfe_energy"][0]
        structure["MinStructureID"] = 0

        dG_direct = fold["mfe_energy"][0]

        return (dG_direct, structure)

    def calc_dG_indirect(self,start_pos):
        
        mRNA = self.mRNA_input[max(0,start_pos-10):min(len(self.mRNA_input),start_pos+35)]
        fold = RNAEnergyModel([mRNA],self.RNA_model)
        fold.mfe([1], Temp = self.temp, dangles = self.dangles)

        structure = fold
        structure["mRNA"] = mRNA
        structure["bp_x"] = fold["mfe_basepairing_x"][0]
        structure["bp_y"] = fold["mfe_basepairing_y"][0]
        structure["dG_mRNA"] = fold["mfe_energy"][0]
        structure["MinStructureID"] = 0

        dG_indirect = fold["mfe_energy"][0]

        return (dG_indirect, structure)
    # ================================================================================================        
        
    def calc_dG_rRNA(self):
        fold = RNAEnergyModel([self.rRNA],self.RNA_model)
        fold.mfe([1], Temp = self.temp, dangles = "all")
        dG_rRNA_folding = fold["mfe_energy"][0]
        return dG_rRNA_folding

    def calc_dG_SDopen(self, mRNA_structure, mRNA_rRNA_structure):
        """Calculate the dG required to unfold any nucleotides in the mRNA folding that are bound to the rRNA in the mRNA:rRNA folding. Between the mRNA and mRNA:rRNA states, the base pairing of these nucleotides are mutually exlusive."""

        mRNA = mRNA_structure["mRNA"]
        program = mRNA_structure["program"]
        index = mRNA_structure["MinStructureID"]
        dG_mRNA = mRNA_structure[program + "_energy"][index]

        index = mRNA_rRNA_structure["MinStructureID"]
        bp_x_1 = mRNA_rRNA_structure["subopt_basepairing_x"][index][:]
        bp_y_1 = mRNA_rRNA_structure["subopt_basepairing_y"][index][:]

        most_5p_mRNA = 0
        most_3p_mRNA = 0

        for (nt_x, nt_y) in zip(bp_x_1, bp_y_1):
            if nt_y > len(mRNA): #nt is rRNA
                most_5p_mRNA = min(most_5p_mRNA, bp_x_1[bp_y_1.index(nt_y)])
                most_3p_mRNA = max(most_3p_mRNA, bp_x_1[bp_y_1.index(nt_y)])

        pre_mRNA = mRNA[0:most_5p_mRNA]
        post_mRNA = mRNA[most_3p_mRNA+1:len(mRNA)+1]

        if len(pre_mRNA) > 0:
            pre_fold = RNAEnergyModel([pre_mRNA],material = self.RNA_model)
            pre_fold.mfe([1],dangles = self.dangles, Temp = self.temp)
            dG_pre = pre_fold["mfe_energy"][0]
        else:
            dG_pre = 0.0

        if len(post_mRNA) > 0:
            post_fold = RNAEnergyModel([post_mRNA],material = self.RNA_model)
            post_fold.mfe([1],dangles = self.dangles, Temp = self.temp)
            dG_post = post_fold["mfe_energy"][0]
        else:
            dG_post = 0.0

        energy = dG_pre + dG_post

        ddG_mRNA = energy - dG_mRNA #positive if work is required to unfold SD sequence
        return ddG_mRNA

    def calc_kinetic_score(self, structure = None, mRNA_in = None, bp_x_in = None, bp_y_in = None):
        """Calculate a "kinetic score", a measure of the longest time required for the mRNA secondary structure to form. Determine the largest distance (in primary sequence space) between two base paired nucleotides (excluding mRNA:rRNA base pairing) and divide it by the length of the mRNA."""

        if not (structure is None):
            program = structure["program"]
            mRNA = structure["mRNA"]
            index = structure["MinStructureID"]
            bp_x = structure[program + "_basepairing_x"][index]
            bp_y = structure[program + "_basepairing_y"][index]

        if not (bp_x_in is None) and not (bp_y_in is None) and not (mRNA_in is None):
            mRNA = mRNA_in[:]
            bp_x = bp_x_in[:]
            bp_y = bp_y_in[:]

        largest_range_helix = 0
        for (nt_x, nt_y) in zip(bp_x, bp_y):
            if nt_x <= len(mRNA) and nt_y <= len(mRNA):
                val = nt_y - nt_x
                largest_range_helix = max(val, largest_range_helix)

        kinetic_score = float(largest_range_helix) / float(len(mRNA))
        if float(largest_range_helix) > 0:
            min_bp_prob = float(largest_range_helix)**(-1.44)
        else:
            min_bp_prob = 1.0

        return (kinetic_score, min_bp_prob)

    def calc_most_5p_mRNA(self,structure_old):

        import copy
        structure = copy.deepcopy(structure_old)
        mRNA = structure["mRNA"]
        bp_x = structure["bp_x"]
        bp_y = structure["bp_y"]

        #Identify the most 5p mRNA nt that is bound to rRNA
        most_5p_mRNA = 0
        for (nt_x, nt_y) in zip(bp_x, bp_y):
            if nt_x <= len(mRNA) and nt_y > len(mRNA): #nt_x is mRNA, nt_y is rRNA, they are bound.
                most_5p_mRNA = nt_x
                break

        return most_5p_mRNA

    def cpu_time(self):
        import resource
        return resource.getrusage(resource.RUSAGE_SELF)[0]

    def combsort(self, input_list):
        """Sorts the input_list in increasing order (minimum first) according to the comb sort algorithm from Wikipedia. Outputs the corresponding sorted index. """

        index = range(len(input_list))

        #Implementation of comb sort (from Wikipedia)
        shrink_factor = 1.24733095

        #Init
        gap = len(input_list)
        swaps = False

        while not (gap <= 1 and not swaps):

            #update the gap value for a next comb
            if gap > 1:
                gap = int(gap / shrink_factor)
                if gap == 10 or gap == 9:
                    gap = 11

            i = 0
            swaps = False #see bubblesort for an explanation

            #a single "comb" over the input list
            while gap + i < len(input_list):
                val = input_list[i]
                if val > input_list[i+gap]:
                    #Swap values -- verify that lists are being altered correctly
                    input_list[i] = input_list[gap + i]
                    input_list[gap + i] = val

                    #Swap indices of index list accordingly
                    temp = index[i]
                    index[i] = index[gap + i]
                    index[gap + i] = temp

                    swaps = True

                i += 1

        return index
       
    def run(self):

        start = self.cpu_time()

        #Initialization of data structures
        self.codon_list = []
        self.start_pos_list = []
        self.dG_total_list = []
        self.dG_mRNA_list = []
        self.dG_direct_list = []
        self.dG_indirect_list = []
        self.dG_mRNA_rRNA_list = []
        self.dG_SD_aSD_list = []
        self.fold_x_list = []
        self.fold_y_list = []
        self.dG_start_energy_list = []
        self.dG_spacing_list = []
        self.mRNA_structure_list = []
        self.mRNA_rRNA_uncorrected_structure_list = []
        self.mRNA_rRNA_corrected_structure_list = []
        self.dG_standby_site_list = []      
        self.kinetic_score_list = []
        self.min_bp_prob_list = []
        self.three_state_indicator_list = []

        self.most_5p_mRNA_list = []
        self.Expression_list = []

        for (start_pos, codon) in self.find_start_codons(self.mRNA_input):
            try:

                #print "Top of calc_dG here"

                #Set dangles based on length between 5' end of mRNA and start codon
                if self.auto_dangles:

                    if start_pos > self.cutoff:
                        self.dangles = "none"

                    else:
                        self.dangles = "all"

                else:
                    self.dangles = self.dangles_default
                    #print "Auto Dangles set to ", self.dangles

                #Start codon energy
                dG_start_codon = self.start_codon_energies[codon]

                #Energy of mRNA folding
                [dG_mRNA,mRNA_structure] = self.calc_dG_mRNA(start_pos)

                #Energy of mRNA:rRNA hybridization & folding
                [dG_mRNA_rRNA_withspacing,mRNA_rRNA_structure] = self.calc_dG_mRNA_rRNA(start_pos)

                # Energy of mRNA:rRNA hybridization
                dG_SD_aSD = mRNA_rRNA_structure["dG_SD_aSD"]
                dG_spacing = mRNA_rRNA_structure["dG_spacing"]

                dG_mRNA_rRNA_nospacing = mRNA_rRNA_structure["dG_mRNA_rRNA"]
                if dG_mRNA_rRNA_nospacing == 0: raise CalcError("Error: RBS Calculator v1.0 could not calculate the correct 16S rRNA binding site.")
                
                #Standby site correction:
                [dG_standby_site, corrected_structure] = self.calc_dG_standby_site(mRNA_rRNA_structure, rRNA_binding = True)
                
                # Energy of mRNA folding of structure from [-25,18] (dG_direct)
                [dG_direct,mRNA_structure] = self.calc_dG_direct(start_pos)

                # Energy of mRNA folding of structure from [-10,35] (dG_indirect)
                [dG_indirect,mRNA_structure] = self.calc_dG_indirect(start_pos)
                
                #Total energy is mRNA:rRNA + start - rRNA - mRNA - standby_site
                dG_total = dG_SD_aSD + dG_spacing + dG_start_codon - 0.5*dG_direct - 0.5*dG_indirect # Reis

                #Calculate 'kinetic score': directly related to probability of base pair formation
                (kinetic_score,min_bp_prob) = self.calc_kinetic_score(mRNA_structure)

                #Calculate dG to open SD sequence
                ddG_SD_open = self.calc_dG_SDopen(mRNA_structure, mRNA_rRNA_structure)


                self.codon_list.append(codon)
                
                self.mRNA_structure_list.append(mRNA_structure)
                self.mRNA_rRNA_uncorrected_structure_list.append(mRNA_rRNA_structure)
                self.mRNA_rRNA_corrected_structure_list.append(corrected_structure)
                self.most_5p_mRNA_list.append(self.calc_most_5p_mRNA(mRNA_rRNA_structure))

                self.dG_start_energy_list.append(dG_start_codon)
                self.dG_mRNA_list.append(dG_mRNA)
                self.dG_mRNA_rRNA_list.append(dG_mRNA_rRNA_nospacing)
                self.dG_SD_aSD_list.append(dG_SD_aSD)
                self.dG_spacing_list.append(mRNA_rRNA_structure["dG_spacing"])
                self.dG_standby_site_list.append(dG_standby_site)
                self.dG_total_list.append(dG_total)

                self.min_bp_prob_list.append(min_bp_prob)
                self.kinetic_score_list.append(kinetic_score)
                self.three_state_indicator_list.append(ddG_SD_open)

                #Start positions
                self.start_pos_list.append(start_pos)

                #Expression levels
                self.Expression_list.append(self.calc_expression_level(dG_total))
                
                # UTR Designer
                self.dG_direct_list.append(dG_direct)
                self.dG_indirect_list.append(dG_indirect)
                
                #For exporting the relevant structure to a PDF
                #index = mRNA_rRNA_structure["MinStructureID"]
                #mRNA_rRNA_structure.export_PDF(index, name = self.name + ": Before standby site", filename =  self.name + "_Before_Standby_rRNA.pdf", program = "subopt")

                #index = corrected_structure["MinStructureID"]
                #corrected_structure.export_PDF(index, name = self.name + ": After standby site", filename =  self.name + "_After_Standby_rRNA.pdf", program = "subopt")

                #print "3 state indicator = ", dG_mRNA_rRNA_nospacing - dG_mRNA + ddG_SD_open

                #print "C = dG_mRNA_rRNA (no spacing) - dG_mRNA + ddG_open = ", dG_mRNA_rRNA_nospacing - dG_mRNA + ddG_SD_open
                #print "mRNA structure = ", mRNA_structure  

            except CalcError, msg:
                print msg
                
                self.codon_list.append(codon)
                self.mRNA_structure_list.append([])
                self.mRNA_rRNA_uncorrected_structure_list.append([])
                self.mRNA_rRNA_corrected_structure_list.append([])

                self.most_5p_mRNA_list.append(self.infinity)
                self.dG_start_energy_list.append(self.infinity)
                self.dG_mRNA_list.append(self.infinity)
                self.dG_mRNA_rRNA_list.append(self.infinity)
                self.dG_SD_aSD_list.append(self.infinity)
                self.dG_spacing_list.append(self.infinity)
                self.dG_standby_site_list.append(self.infinity)
                self.dG_total_list.append(self.infinity)

                self.min_bp_prob_list.append(self.infinity)
                self.kinetic_score_list.append(self.infinity)
                self.three_state_indicator_list.append(self.infinity)
                
                self.start_pos_list.append(start_pos)
                self.Expression_list.append(0)
                
                # Added for UTR Designer
                self.dG_direct_list.append(self.infinity)
                self.dG_indirect_list.append(self.infinity)
                
        self.has_run = 1

        end = self.cpu_time()
        self.run_time = end - start

    def calc_expression_level(self, dG):

        import math
        return RBS_Calculator.K * math.exp(-dG / RBS_Calculator.RT_eff)

    def output(self):
        """Syntax: mRNA = RBS_Calculator.output(). In the future, apply RBS_Calculator on input mRNA class instance."""
        mRNA = Class_RNA()
        mRNA.sequence = self.mRNA_input
        
        mRNA.RBS_list = []
        
        for i in range(len(self.start_pos_list)):
        
            RBS = Class_Ribosome_Binding_Site()
            
            RBS.sequence = self.mRNA_input[max(0,self.start_pos_list[i] - self.cutoff):min(len(self.mRNA_input), self.start_pos_list[i] + self.cutoff)]
            RBS.pre = self.mRNA_input[max(0,self.start_pos_list[i] - self.cutoff):self.start_pos_list[i]]
            RBS.post = self.mRNA_input[self.start_pos_list[i]+3:min(len(self.mRNA_input), self.start_pos_list[i] + 3 + self.cutoff)]
            RBS.rRNA = self.rRNA
            
            RBS.start_codon = self.codon_list[i]
            RBS.start_position = self.start_pos_list[i]
            RBS.dG_total = self.dG_total_list[i]
            RBS.dG_mRNA = self.dG_mRNA_list[i]
            RBS.dG_mRNA_rRNA = self.dG_mRNA_rRNA_list[i]
            RBS.dG_SD_aSD = self.dG_SD_aSD_list[i]
            RBS.dG_start = self.dG_start_energy_list[i]
            RBS.dG_spacing = self.dG_spacing_list[i]
            RBS.dG_standby = self.dG_standby_site_list[i]
            RBS.tir = self.Expression_list[i]
            
            # Added for UTR Designer
            RBS.dG_direct = self.dG_direct_list[i]
            RBS.dG_indirect = self.dG_indirect_list[i]
            
            try:
                RBS.initial_structure = {'mRNA':self.mRNA_structure_list[i]["mRNA"],
                                   'bp_x':self.mRNA_structure_list[i]["bp_x"],
                                   'bp_y':self.mRNA_structure_list[i]["bp_y"],
                                   }
            except:
                #print RBS.sequence
                #print "No initial mRNA structure found. Saving absence of structure."
                RBS.initial_structure = {'mRNA':RBS.sequence,
                                   'bp_x':[],
                                   'bp_y':[],
                                   }
                       
            try:
                RBS.final_structure = {'mRNA':self.mRNA_rRNA_corrected_structure_list[i]["mRNA"],
                                   'rRNA':self.rRNA,
                                   'bp_x':self.mRNA_rRNA_corrected_structure_list[i]["bp_x"],
                                   'bp_y':self.mRNA_rRNA_corrected_structure_list[i]["bp_y"],
                                  }
            except:
                #print RBS.sequence
                #print "No final mRNA structure found. Saving absence of structure."
                RBS.final_structure = {'mRNA':RBS.sequence,
                                   'rRNA':self.rRNA,
                                   'bp_x':[],
                                   'bp_y':[],
                                  }
                
            RBS.kinetic_score = self.kinetic_score_list[i]
            RBS.three_state_indicator_list = self.three_state_indicator_list[i]
            
            RBS.warnings = []
            if len(RBS.post) < self.cutoff:
                RBS.warnings.append(self.Warning_List['SHORT_CDS'])
            
            if RBS.kinetic_score > 0.50:
                RBS.warnings.append(self.Warning_List['NON_EQUILIBRIUM'])
            
            if (i > 1) and (RBS.start_position  - self.start_pos_list[i-1] < 6):
                RBS.warnings.append(self.Warning_List['CLOSE_START_CODONS'])
            elif (i < len(self.start_pos_list)-1) and (self.start_pos_list[i+1]  - RBS.start_position < 6):
                RBS.warnings.append(self.Warning_List['CLOSE_START_CODONS'])                
            
            if len(RBS.warnings) == 0: RBS.warnings.append(self.Warning_List['NONE'])
            
            mRNA.RBS_list.append(RBS)
            
        return mRNA
        
    def print_dG(self,max_dG = 1e12,brief = 0, return_string = False, print_expression = False):
        '''Print out useful information about the mRNA sequence'''

        import math

        print_string = ""
        if self.has_run == 1:

            print "mRNA Sequence: ", self.name
            if return_string: print_string = print_string + "mRNA Sequence: " + self.name + "\n"

            if len(self.start_position_list) == 0:
                print "No start codons found in input mRNA sequence"
                print "----------------------------------------------------------------------------------------"
                if return_string: print_string = print_string + "No start codons found in input mRNA sequence" + "\n" + "----------------------------------------------------------------------------------------" + "\n"
            elif len( [dG_total for dG_total in self.dG_total_list if dG_total < max_dG] ) == 0:
                print "No RBSs found with dG <", str(max_dG), "."
                print "----------------------------------------------------------------------------------------" + "\n"
                if return_string: print_string = print_string + "No RBSs found with dG <" + str(max_dG) + "." + "\n" + "----------------------------------------------------------------------------------------" + "\n"

            else:

                if print_expression:
                    Headers = ("Start", "(pos)", "Expression level","Kinetic score")
                    format = "%4s %4s %15s %15s"
                else:
                    Headers = ("Start","(pos)", "dG total", "dG (rRNA:mRNA)", "dG (mRNA)", "dG (spacing)", "dG (standby)", "Kinetic Score")
                    format = "%4s %4s %12s %12s %12s %12s %12s %12s"

                print format % Headers

                if return_string: print_string = print_string + format % Headers + "\n"

                for (start,counter) in zip(self.start_position_list,range(len(self.start_position_list))):

                    dG_total = self.dG_total_list[counter]
                    Expression = RBS_Calculator.K * math.exp(-dG_total / RBS_Calculator.RT_eff)

                    if (dG_total < max_dG):

                        if (dG_total > self.infinity):
                            dG_total = "Inf"
                            Expression = 0.0

                        dG_mRNA = self.dG_mRNA_list[counter]
                        if (dG_mRNA > self.infinity): dG_mRNA = "Inf"

                        dG_mRNA_rRNA = self.dG_mRNA_rRNA_list[counter]
                        if (dG_mRNA_rRNA > self.infinity): dG_mRNA_rRNA = "Inf"

                        dG_spacing = self.dG_spacing_list[counter]
                        if (dG_spacing > self.infinity): dG_spacing = "Inf"

                        dG_standby_site = self.dG_standby_site_list[counter]
                        if (dG_standby_site > self.infinity): dG_standby_site = "Inf"

                        start_codon = self.start_codon_list[counter]

                        dG_rRNA = self.dG_rRNA

                        kinetic_score = self.kinetic_score_list[counter]

                        if print_expression:
                            print format % (start_codon, str(start), str(round(Expression,2)), str(round(kinetic_score,2)))

                        else:
                            print format % (start_codon, str(start), str(dG_total), str(dG_mRNA_rRNA), str(dG_mRNA), str(dG_spacing), str(dG_standby_site), str(round(kinetic_score,2)))

                        if return_string: print_string = print_string + format % (start_codon, str(start), str(dG_total), str(dG_mRNA_rRNA), str(dG_mRNA), str(dG_spacing), str(dG_standby_site), str(round(kinetic_score,2))) + "\n"

                print "----------------------------------------------------------------------------------------"
                if return_string: print_string = print_string + "----------------------------------------------------------------------------------------" + "\n"

                #print "Computation Time: ", str(self.run_time), " seconds."

                if return_string: return print_string
        else:
            raise RuntimeError("The RBS Calculator has not been run yet. Call the 'run' method.")

    def save_data(self, handle, header = False):

        infinity = 1e6

        if self.has_run == 1:

            if (header):
                #Print header information --> all the parameters
                parameters = "RNA model: \t%s \nDangles: \t%s \nTemperature: \t%s \nOptimal Spacing \t%s \nSpacing constant (push) \t%s\nSpacing constant (pull) \t%s\nCutoff \t%s\nStandby Site Length \t%s\nEnergy Cutoff \t%s\nrRNA sequence \t%s\n" % (self.RNA_model, str(self.dangles), str(self.temp), str(self.optimal_spacing), str(self.dG_spacing_constant_push), str(self.dG_spacing_constant_pull), str(self.cutoff), str(self.standby_site_length), str(self.energy_cutoff), self.rRNA)

                handle.writelines(parameters)

            #Name, dG total, dG rRNA:mRNA, dG mRNA, dG spacing, dG standby, dG start codon, kinetic score, longest helix, longest loop, start position
            #format = "%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\n"
            format = "%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\n"

            for counter in range(len(self.start_position_list)):

                dG_total = self.dG_total_list[counter]
                if (dG_total >= infinity): dG_total = "Inf"

                dG_mRNA = self.dG_mRNA_list[counter]
                if (dG_mRNA >= infinity): dG_mRNA = "Inf"

                dG_mRNA_rRNA = self.dG_mRNA_rRNA_list[counter]
                if (dG_mRNA_rRNA >= infinity): dG_mRNA_rRNA = "Inf"

                dG_spacing = self.dG_spacing_list[counter]
                if (dG_spacing >= infinity): dG_spacing = "Inf"

                dG_standby_site = self.dG_standby_site_list[counter]
                if (dG_standby_site >= infinity): dG_standby_site = "Inf"

                kinetic_score = self.kinetic_score_list[counter]

                most_5p_mRNA = self.most_5p_mRNA_list[counter]

                three_state = self.three_state_indicator_list[counter]

                handle.writelines(format % (self.name.split(" ")[0],dG_total, dG_mRNA_rRNA,dG_mRNA,dG_spacing,dG_standby_site, self.start_codon_energies[self.start_codon_list[counter]],round(kinetic_score,2), str(self.start_position_list[counter]), str(most_5p_mRNA), str(three_state) ))
        else:
            raise RuntimeError("The RBS Calculator has not been run yet. Call the 'run' method.")

#----------------------------------------------------------------------------------------------------------
#End RBS_Calculator class
#----------------------------------------------------------------------------------------------------------
def calc_dG_from_file(handle, output, verbose = True, parameters = {}):
    from Bio import SeqIO
    from RBS_Calculator import RBS_Calculator

    records = SeqIO.parse(handle,"fasta")
    First = True
    export_PDF = True

    #for i in range(30):
    #    records.next()

    for record in records:

        mRNA = record.seq.tostring().upper()

        #Set any defaults
        start_range = [0, len(mRNA)]
        name = record.description.split(" ")[0]

        #Create instance of RBS Calculator
        test = RBS_Calculator(mRNA, start_range, name)

        #Examine kvars dictionary and pull out any options. Assign them to instanced class.

        for (key,value) in parameters.items():

            if key == "cutoff":
                test.cutoff = value
            elif key == "start_range":
                test.start_range = value
            elif key == "rRNA":
                test.rRNA = value
            elif key == "energy_cutoff":
                test.energy_cutoff = value
            elif key == "standby_site_length":
                test.standby_site_length = value
            elif key == "dangles":
                test.dangles = value
            elif key == "pre_window_length":
                test.pre_window_length = value
            elif key == "post_window_length":
                test.post_window_length = value
            elif key == "export_PDF":
                export_PDF = value

        test.run()
        test.print_dG(test.infinity,print_expression=verbose)

        test.save_data(output, First)
        if First:
            First = False

        if export_PDF:
            num_structs = len(test.mRNA_rRNA_uncorrected_structure_list)
            for (structure,counter) in zip(test.mRNA_rRNA_uncorrected_structure_list,range(num_structs)):
                index = structure["MinStructureID"]
                structure.export_PDF(index, name, filename =  name + "_rRNA" + "_" + str(counter) + ".pdf", program = "subopt")

            num_structs = len(test.mRNA_structure_list)
            for (structure,counter) in zip(test.mRNA_structure_list,range(num_structs)):
                structure.export_PDF(0, name, filename = name + "_mRNA" + "_" + str(counter) + ".pdf")


    output.close()

def calc_dG_pre_post_RBS(pre_list,post_list,RBS_list,name_list,output,verbose = True, parameters = {}):

    from RBS_Calculator import RBS_Calculator

    First = True

    for (pre,post,RBS,name) in zip(pre_list,post_list,RBS_list,name_list):

        mRNA = pre + RBS + post

        start_range = [0, len(mRNA)]

        #Create instance of RBS Calculator
        test = RBS_Calculator(mRNA, start_range, name)

        #Examine kvars dictionary and pull out any options. Assign them to instanced class.

        for (key,value) in parameters.items():

            if key == "cutoff":
                test.cutoff = value
            elif key == "start_range":
                test.start_range = value
            elif key == "rRNA":
                test.rRNA = value
            elif key == "energy_cutoff":
                test.energy_cutoff = value
            elif key == "standby_site_length":
                test.standby_sitRBSe_length = value
            elif key == "dangles":
                test.dangles = value
            elif key == "pre_window_length":
                test.pre_window_length = value
            elif key == "post_window_length":
                test.post_window_length = value
            elif key == "export_PDF":
                export_PDF = value

        test.run()
        if verbose: test.print_dG(test.infinity)

        test.save_data(output, First)
        if First:
            First = False

        #index = test.mRNA_rRNA_corrected_structure_list[0]["MinStructureID"]

        #test.mRNA_rRNA_corrected_structure_list[0].export_PDF(index, name, filename =  name + "_rRNA" + ".pdf")

        #test.mRNA_structure_list[0].export_PDF(0, name, filename = name + "_mRNA" + ".pdf")


    output.close()

if __name__ == "__main__":
   
    name = "Test"
    
    mRNA = "GCCCACCGATTGATTCGTCTTCACCAGATGATTTAATTCACTGTTTCTCACGAACTCAGTCAAAGCGCTA"
    start_range = [0, len(mRNA)]
    test = RBS_Calculator(mRNA, start_range, name)
    test.run()
    test.print_dG(test.infinity)
            
    export_PDF = False

    if export_PDF:
        num_structs = len(test.mRNA_rRNA_corrected_structure_list)
        for (structure,counter) in zip(test.mRNA_rRNA_corrected_structure_list,range(num_structs)):

            index = structure["MinStructureID"]
            structure.export_PDF(index, name, filename =  name + "_rRNA" + "_" + str(counter) + ".pdf", program = "energy")

        for (structure,counter) in zip(test.mRNA_rRNA_uncorrected_structure_list,range(num_structs)):

            index = structure["MinStructureID"]
            structure.export_PDF(index, name, filename =  name + "_rRNA_uncorrected" + "_" + str(counter) + ".pdf", program = "energy")

        num_structs = len(test.mRNA_structure_list)
        for (structure,counter) in zip(test.mRNA_structure_list,range(num_structs)):
            structure.export_PDF(0, name, filename = name + "_mRNA" + "_" + str(counter) + ".pdf")
