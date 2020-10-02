#Python wrapper for the Vienna RNA Package by Andreas R. Gruber, Ronny Lorenz, Stephan H. Bernhart, Richard Neub?ck, and Ivo L. Hofacker (NAR, 2008).

#This file is part of the Ribosome Binding Site Calculator.

#The Ribosome Binding Site Calculator is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#The Ribosome Binding Site Calculator is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with Ribosome Binding Site Calculator.  If not, see <http://www.gnu.org/licenses/>.

#This Python wrapper is written by Howard Salis. Copyright 2008-2009. All rights reserved. :)
#Use at your own risk.

import os.path
import os, time
from subprocess import Popen, PIPE, STDOUT

import os, popen2, random, string

import cStringIO


current_dir = os.path.dirname(os.path.abspath(__file__)) + "/tmp"
if not os.path.exists(current_dir): os.mkdir(current_dir)

debug=0

#Class that encapsulates all of the functions from NuPACK 2.0
class ViennaRNA(dict):

    debug_mode = 0
    RT = 0.61597 #gas constant times 310 Kelvin (in units of kcal/mol)

    def __init__(self,Sequence_List,material = "rna37"):

        self.ran = 0

        import re
        import random
        import string

        exp = re.compile('[ATGCU]',re.IGNORECASE)

        for seq in Sequence_List:
            if exp.match(seq) == None:
                error_string = "Invalid letters found in inputted sequences. Only ATGCU allowed. \n Sequence is \"" + str(seq) + "\"."
                raise ValueError(error_string)

        self["sequences"] = Sequence_List
        self["material"] = material

        random.seed(time.time())
        long_id = "".join([random.choice(string.letters + string.digits) for x in range(10)])
        self.prefix = current_dir + "/temp_" + long_id

    def mfe(self, strands,Temp = 37.0, dangles = "all",outputPS = True):
        
        self["mfe_composition"] = strands
        
        if Temp <= 0: raise ValueError("The specified temperature must be greater than zero.")

        seq_string = "&".join(self["sequences"])
        input_string = seq_string + "\n"
              
        #handle = open(self.prefix,"w")
        #handle.write(input_string)
        #handle.close()

        #Set arguments
        material = self["material"]
        if dangles is "none":
            dangles = "-d0"
        elif dangles is "some":
            dangles = "-d1"
        elif dangles is "all":
            dangles = "-d2"

        if outputPS:
            outputPS_str = ""
        else:
            outputPS_str = "-noPS"

        #Call ViennaRNA C programs
        
        #========================================================================================================================================        
            
        if len(strands) == 1:
            CommandArgs=['RNAfold']+[dangles]+[outputPS_str]
            
        elif len(strands) >= 1:
            CommandArgs=['RNAcofold']+[dangles]+[outputPS_str]

        #========================================================================================================================================        
        
        p = Popen(CommandArgs, stdout=PIPE, stdin=PIPE, stderr=STDOUT)
        grep_stdout = p.communicate(input=input_string)
        VRNAResult=grep_stdout[0]
        output = cStringIO.StringIO(VRNAResult)
        line = output.readline()
        line = output.readline()
        
        words = line.split(" ")
        bracket_string = words[0]
        (strands,bp_x, bp_y) = self.convert_bracket_to_numbered_pairs(bracket_string)

        energy = float(words[len(words)-1].replace(")","").replace("(","").replace("\n",""))
          
        # ========================================================================================================
        # ========================================================================================================
        if len(strands) == 2:
            energy += -2.4809999999999999
        
        # ========================================================================================================
        # ========================================================================================================
     
        
        #self._cleanup()
        self["program"] = "mfe"
        self["mfe_basepairing_x"] = [bp_x]
        self["mfe_basepairing_y"] = [bp_y]
        self["mfe_energy"] = [energy]
        self["totalnt"]=strands

        #print "Minimum free energy secondary structure has been calculated."

    def subopt(self, strands,energy_gap,Temp = 37.0, dangles = "all", outputPS = False):

        self["subopt_composition"] = strands

        if Temp <= 0: raise ValueError("The specified temperature must be greater than zero.")

        seq_string = "&".join(self["sequences"])
        input_string = seq_string + "\n"

        #handle = open(self.prefix,"w")
        #handle.write(input_string)
        #handle.close()

        #Set arguments
        material = self["material"]
        if dangles is "none":
            dangles = "-d0"
        elif dangles is "some":
            dangles = "-d1"
        elif dangles is "all":
            dangles = "-d2"

        # FIXED FOR COMPATIBILITY WITH v2.2 - REIS 
        #if outputPS:
            #outputPS_str = ""
        #else:
            #outputPS_str = "-noPS"
        outputPS_str = ""

        #Call ViennaRNA C programs
        CommandArgs=['RNAsubopt']+["-e" , str(energy_gap) ,"-s", outputPS_str,dangles]
        p = Popen(CommandArgs, stdout=PIPE, stdin=PIPE, stderr=STDOUT)
        grep_stdout = p.communicate(input=input_string)
        VRNAResult=grep_stdout[0]
        output = cStringIO.StringIO(VRNAResult)
        line = output.readline()



        self["subopt_basepairing_x"] = []
        self["subopt_basepairing_y"] = []
        self["subopt_energy"] = []
        self["totalnt"]=[]
        counter=0

        while len(line)>0:            
            line = output.readline()
            if len(line) > 0:
                counter+=1
                words = line.split(" ")
                bracket_string = words[0]
                energy = float(words[len(words)-1].replace("\n",""))

                (strands,bp_x, bp_y) = self.convert_bracket_to_numbered_pairs(bracket_string)
                
                               
                # ========================================================================================================
                # ========================================================================================================
                if len(strands) == 2:
                    energy += -2.4809999999999999
                
                # ========================================================================================================
                # ========================================================================================================
 
                self["subopt_energy"].append(energy)
                self["subopt_basepairing_x"].append(bp_x)
                self["subopt_basepairing_y"].append(bp_y)
        
        self["subopt_NumStructs"] = counter

        #self._cleanup()
        self["program"] = "subopt"

        #print "Minimum free energy and suboptimal secondary structures have been calculated."

    def energy(self, strands, base_pairing_x, base_pairing_y, Temp = 37.0, dangles = "all"):

        self["energy_composition"] = strands

        if Temp <= 0: raise ValueError("The specified temperature must be greater than zero.")

        seq_string = "&".join(self["sequences"])
        strands = [len(seq) for seq in self["sequences"]]
        bracket_string = self.convert_numbered_pairs_to_bracket(strands,base_pairing_x,base_pairing_y)
        input_string = seq_string + "\n" + bracket_string + "\n"

        #handle = open(self.prefix,"w")
        #handle.write(input_string)
        #handle.close()

        #Set arguments
        material = self["material"]
        if dangles is "none":
            dangles = "-d0"
        elif dangles is "some":
            dangles = "-d1"
        elif dangles is "all":
            dangles = "-d2"


        self["energy_energy"] = []

        #Call ViennaRNA C programs

        CommandArgs=['RNAeval']+[dangles]
        p = Popen(CommandArgs, stdout=PIPE, stdin=PIPE, stderr=STDOUT)
        grep_stdout = p.communicate(input=input_string)
        VRNAResult=grep_stdout[0]
        
        output = cStringIO.StringIO(VRNAResult)
        line = output.readline()
        #print line
        
        line = output.readline()
        #print line
        
        words = line.split(" ")
        #print words

        line = output.readline()
        #print line
        
        energy = float(words[len(words)-1].replace("(","").replace(")","").replace("\n",""))
               
        # ========================================================================================================
        # ========================================================================================================
        if len(strands) == 2:
            energy += -2.4809999999999999
        
        # ========================================================================================================
        # ========================================================================================================
     

        self["program"] = "energy"
        self["energy_energy"].append(energy)
        self["energy_basepairing_x"] = [base_pairing_x]
        self["energy_basepairing_y"] = [base_pairing_y]
        #self._cleanup()

        return energy

    def export_PDF(self, complex_ID, name = "", filename = "temp.pdf", program = None):
        """Uses Zuker's sir_graph_ng and ps2pdf to convert a secondary structure described in .ct format
        to a PDF of the RNA"""

        if program is None:
            program = self["program"]

        inputfile = "temp.ct"

        self.Convert_to_ct(complex_ID,name,inputfile,program)

        cmd = "sir_graph" #Assumes it's on the path
        args = "-p" #to PostScript file
        output = popen2.Popen3(cmd + " " + args + " " + inputfile,"r")

        #print "Ran sir_graph"
        
        output.wait()
        if debug == 1: print output.fromchild.read()

        inputfile = inputfile[0:len(inputfile)-2] + "ps"

        #cmd = "ps2pdf" #Assumes it's on the path
        cmd = "gs -q -dBATCH  -dAutoFilterColorImages=false -sColorImageFilter=FlateEncode -dNOPAUSE -sPAPERSIZE=a3 -sDEVICE=pdfwrite -sOutputFile=temp.pdf" #Assumes it's on the path
        output = popen2.Popen3(cmd + " " + inputfile,"r")
        output.wait()

        if debug == 1: print output.fromchild.read()

        outputfile = inputfile[0:len(inputfile)-2] + "pdf"

        #Remove the temporary file "temp.ct" if it exists
        if os.path.exists("temp.ct"): os.remove("temp.ct")

        #Remove the temporary Postscript file if it exists
        if os.path.exists(inputfile): os.remove(inputfile)

        #Rename the output file to the desired filename.
        if os.path.exists(outputfile): os.rename(outputfile,filename)
        #Done!

    def Convert_to_ct(self,complex_ID,name,filename = "temp.ct",program = "ordered"):
        """Converts the secondary structure of a single complex into the .ct file format, which is used
        with sir_graph_ng (or other programs) to create an image of the secondary structure."""

        #hacksy way of reading from data produced by 'complex', by 'mfe', or by 'subopt'
        data_x = program + "_basepairing_x"
        data_y = program + "_basepairing_y"
        mfe_name = program + "_energy"
        composition_name = program + "_composition"

        #print data_x
        #print data_y
        #print mfe_name
        #print composition_name
        #print self[data_x][0]
        #print self[data_y][0]
        #print self[mfe_name][0]
        #Format of .ct file

        #Header: <Total # nt> \t dG = <# mfe> kcal/mol \t <name of sequence>
        #The Rest:
        #<nt num> \t <bp letter> \t <3' neighbor> \t <5' neighbor> \t <# of bp'ing, 0 if none> \t ...
        #<strand-specific nt num> \t <3' neighbor if connected by helix> \t <5' neighbor if connected by helix>

        #Extract the data for the desired complex using complex_ID
        bp_x = self[data_x][complex_ID]
        bp_y = self[data_y][complex_ID]
        mfe = self[mfe_name][complex_ID]

        if program == "mfe" or program == "subopt" or program == "energy":
            composition = self[composition_name]
        elif program == "ordered" or program == "unordered":
            composition = self[composition_name][complex_ID]


        #Determine concatenated sequence of all strands, their beginnings, and ends
        allseq = ""
        strand_begins = []
        strand_ends = []

        #Seemingly, the format of the composition is different for the program complex vs. mfe/subopt
        #for mfe/subopt, the composition is the list of strand ids
        #for complex, it is the number of each strand (in strand id order) in the complex
        #for mfe/subopt, '1 2 2 3' refers to 1 strand of 1, 2 strands of 2, and 1 strand of 3.
        #for complex, '1 2 2 3' refers to 1 strand of 1, 2 strands of 2, 2 strands of 3, and 3 strands of 4'.
        #what a mess.

        if program == "mfe" or program == "subopt" or program == "energy":
            for strand_id in composition:
                strand_begins.append(len(allseq) + 1)
                allseq = allseq + self["sequences"][strand_id-1]
                strand_ends.append(len(allseq))

        else:
            for (num_strands,strand_id) in zip(composition,range(len(composition))):
                for j in range(num_strands):
                    strand_begins.append(len(allseq) + 1)
                    allseq = allseq + self["sequences"][strand_id]
                    strand_ends.append(len(allseq))

        seq_len = len(allseq)
##============================================================
        #print "Seq Len = ", seq_len, "  Composition = ", composition
        #print "Sequence = ", allseq
        #print "Base pairing (x) = ", bp_x
        #print "Base pairing (y) = ", bp_y
        #print "Strand begins =", strand_begins
        #print "Strand ends =", strand_ends
#=============================================================

        #Create the header
        header = str(seq_len) + "\t" + "dG = " + str(mfe) + " kcal/mol" + "\t" + name + "\n"

        #Open the file
        handle = open(filename,"w")

        #Write the header
        handle.write(header)

        #Write a line for each nt in the secondary structure
        for i in range(1,seq_len+1):


            for (nt,pos) in zip(strand_begins,range(len(strand_begins))):
                if i >= nt:
                    strand_id = pos


            #Determine 3' and 5' neighbor
            #If this is the beginning of a strand, then the 3' neighbor is 0
            #If this is the end of a strand, then the 5' neighbor is 0

            if i in strand_begins:
                nb_5p = 0
            else:
                nb_5p = i - 1

            if i in strand_ends:
                nb_3p = 0
            else:
                nb_3p = i + 1


            if i in bp_x or i in bp_y:
                if i in bp_x: nt_bp = bp_y[bp_x.index(i)]
                if i in bp_y: nt_bp = bp_x[bp_y.index(i)]
            else:
                nt_bp = 0

            #Determine strand-specific counter
            strand_counter = i - strand_begins[strand_id] + 1

            #Determine the 3' and 5' neighbor helical connectivity
            #If the ith nt is connected to its 3', 5' neighbor by a helix, then include it
            #Otherwise, 0
            #Helix connectivity conditions:
            #The 5' or 3' neighbor is connected via a helix iff:
            #a) helix start: i not bp'd, i+1 bp'd, bp_id(i+1) - 1 is bp'd, bp_id(i+1) + 1 is not bp'd
            #b) helix end: i not bp'd, i-1 bp'd, bp_id(i-1) - 1 is not bp'd, bp_id(i-1) + 1 is bp'd
            #c) helix continued: i and bp_id(i)+1 is bp'd, 5' helix connection is bp_id(bp_id(i)+1)
            #d) helix continued: i and bp_id(i)-1 is bp'd, 3' helix connection is bp_id(bp_id(i)-1)
            #Otherwise, zero.

            #Init
            hc_5p = 0
            hc_3p = 0

            if i in bp_x or i in bp_y:  #helix continued condition (c,d)
                if i in bp_x: bp_i = bp_y[bp_x.index(i)]
                if i in bp_y: bp_i = bp_x[bp_y.index(i)]

                if bp_i+1 in bp_x or bp_i+1 in bp_y: #helix condition c
                    if bp_i+1 in bp_x: hc_3p = bp_y[bp_x.index(bp_i+1)]
                    if bp_i+1 in bp_y: hc_3p = bp_x[bp_y.index(bp_i+1)]

                if bp_i-1 in bp_x or bp_i-1 in bp_y: #helix condition d
                    if bp_i-1 in bp_x: hc_5p = bp_y[bp_x.index(bp_i-1)]
                    if bp_i-1 in bp_y: hc_5p = bp_x[bp_y.index(bp_i-1)]

            else: #helix start or end (a,b)

                if i+1 in bp_x or i+1 in bp_y: #Start, condition a
                    if i+1 in bp_x: bp_3p = bp_y[bp_x.index(i+1)]
                    if i+1 in bp_y: bp_3p = bp_x[bp_y.index(i+1)]

                    if bp_3p + 1 not in bp_x and bp_3p + 1 not in bp_y:
                        hc_3p = i + 1

                if i-1 in bp_x or i-1 in bp_y: #End, condition b
                    if i-1 in bp_x: bp_5p = bp_y[bp_x.index(i-1)]

                    if i-1 in bp_y: bp_5p = bp_x[bp_y.index(i-1)]

                    if bp_5p - 1 not in bp_x and bp_5p - 1 not in bp_y:
                        hc_5p = i - 1


            line = str(i) + "\t" + allseq[i-1] + "\t" + str(nb_5p) + "\t" + str(nb_3p) + "\t" + str(nt_bp) + "\t" + str(strand_counter) + "\t" + str(hc_5p) + "\t" + str(hc_3p) + "\n"

            handle.write(line)

        #Close the file. Done.
        handle.close()

    def convert_bracket_to_numbered_pairs(self,bracket_string):

        all_seq_len = len(bracket_string)
        bp_x = []
        bp_y = []
        strands = []

        for y in range(bracket_string.count(")")):
            bp_y.append([])

        last_nt_x_list = []
        counter=0
        num_strands=0
        for (pos,letter) in enumerate(bracket_string[:]):
            if letter is ".":
                counter += 1

            elif letter is "(":
                bp_x.append(pos-num_strands)
                last_nt_x_list.append(pos-num_strands)
                counter += 1

            elif letter is ")":
                nt_x = last_nt_x_list.pop()
                nt_x_pos = bp_x.index(nt_x)
                bp_y[nt_x_pos] = pos-num_strands
                counter += 1

            elif letter is "&":
                strands.append(counter)
                counter=0
                num_strands+=1

            else:
                print "Error! Invalid character in bracket notation."

        if len(last_nt_x_list) > 0:
            print "Error! Leftover unpaired nucleotides when converting from bracket notation to numbered base pairs."

        strands.append(counter)
        bp_x = [pos+1 for pos in bp_x[:]] #Shift so that 1st position is 1
        bp_y = [pos+1 for pos in bp_y[:]] #Shift so that 1st position is 1

        return (strands,bp_x, bp_y)

    def convert_numbered_pairs_to_bracket(self,strands,bp_x,bp_y):

        bp_x = [pos-1 for pos in bp_x[:]] #Shift so that 1st position is 0
        bp_y = [pos-1 for pos in bp_y[:]] #Shift so that 1st position is 0

        bracket_notation = []
        counter=0
        for (strand_number,seq_len) in enumerate(strands):
            if strand_number > 0: bracket_notation.append("&")
            for pos in range(counter,seq_len+counter):
                if pos in bp_x:
                    bracket_notation.append("(")
                elif pos in bp_y:
                    bracket_notation.append(")")
                else:
                    bracket_notation.append(".")
            counter+=seq_len

        return "".join(bracket_notation)

    def _cleanup(self):

        if os.path.exists(self.prefix): os.remove(self.prefix)
        return

if __name__ == "__main__":

    from datetime import datetime
    import time

    #sequences = ["AGGGGGGATCTCCCCCCAAAAAATAAGAGGTACACATGACTAAAACTTTCAAAGGCTCAGTATTCCCACT"] #,"acctcctta"]   
    sequences = ["AGATGGAGAAGCCATTAATATTAAGGAAGGTAACAATGGCGAGCTCTGAAGACGTTATCAAAGAGTTCAT"]
    sum1=0;

    t1= time.time()
    test = ViennaRNA(sequences,material = "rna37")
    for i in range(0,1):            
      
      test.mfe([1],Temp = 37.0, dangles = "none")
      bp_x = test["mfe_basepairing_x"][0]
      bp_y = test["mfe_basepairing_y"][0]
      strands = test["totalnt"]
      bracket_string = test.convert_numbered_pairs_to_bracket(strands,bp_x,bp_y)      

      (strands,bp_x, bp_y) = test.convert_bracket_to_numbered_pairs(bracket_string)
      test.subopt(strands,3.5,dangles = "all")
      energy_energy=test.energy(strands, bp_x, bp_y, dangles = "all")
      #print bracket_string
      #print "Strands = ", strands
      #print "bp_x = ", bp_x
      #print "bp_y = ", bp_y      
      
      
      #print i,q
      #sum1+=q
      
    t2= time.time()
    ConsumedTime=t2-t1
    print test["mfe_energy"][0], test["subopt_NumStructs"],energy_energy
    print "ConsumedTime=",ConsumedTime
