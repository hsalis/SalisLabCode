"""
Initialize database and add datasets

Copyright 2017 Alexander C. Reis, Howard M. Salis, all rights reserved.

"""

import re
import numpy as np
import pandas as pd
import cPickle as pickle
from openpyxl import load_workbook
import xlrd

def add_dataset(db,datasets):

    '''Amin Espah Borujeni, Anirudh S. Channarasappa, & Howard M. Salis
    Translation rate is controlled by coupled trade-offs between site accessibility, select RNA unfolding and sliding at upstream standby sites
    Nucleic Acids Research, 2014, Vol. 24, No. 4; doi: 10.1093/nar/gkt1139''' 
    paper = 'EspahBorujeni_NAR_2013'

    if paper in datasets:
        path = 'datasets/{}.xls'.format(paper)
        wb = xlrd.open_workbook(path,'r')
        sheet = wb.sheet_by_index(0)

        ds = {
            "5pUTR"     : sheet.col_values(colx=5, start_rowx=5, end_rowx=141),
            "CDS"       : sheet.col_values(colx=6, start_rowx=5, end_rowx=141),
            "PROT.MEAN" : sheet.col_values(colx=8, start_rowx=5, end_rowx=141),
            "PROT.STD"  : sheet.col_values(colx=9, start_rowx=5, end_rowx=141),
            "PROTEIN"   : "RFP",
            "ORGANISM"  : "Escherichia coli str. K-12 substr. DH10B",
            "METHOD"    : "Individually Characterized",
            "TEMP"      : 37.0,
            "DATASET"   : paper
        }
        
        # Add extended dataset
        paperext = 'EspahBorujeni_NAR_2013_extended'
        path = 'datasets/{}.xls'.format(paperext)
        wb = xlrd.open_workbook(path,'r')
        sheet = wb.sheet_by_index(0)

        ds["5pUTR"]     += sheet.col_values(colx=3, start_rowx=3, end_rowx=42)
        ds["CDS"]       += sheet.col_values(colx=5, start_rowx=3, end_rowx=42)
        ds["PROT.MEAN"] += sheet.col_values(colx=8, start_rowx=3, end_rowx=42)
        ds["PROT.STD"]  += sheet.col_values(colx=9, start_rowx=3, end_rowx=42)

        ds["SEQUENCE"] = [UTR+CDS for UTR,CDS in zip(ds["5pUTR"],ds["CDS"])]
        ds["STARTPOS"] = [len(UTR) for UTR in ds["5pUTR"]]

        df = pd.DataFrame(ds)
        db = db.append(df, ignore_index=True)


    '''Amin Espah Borujeni, Bennis M. Mishler, Jingzhi Wang, Walker Huso, & Howard M. Salis
    Automated Physics-Based Design of Synthetic Riboswitches from Diverse RNA Aptamers
    Nucleic Acids Research, 2015, doi: 10.1093/nar/gkv1289'''
    paper = 'EspahBorujeni_NAR_2015'
    
    if paper in datasets:
        path = 'datasets/{}.xls'.format(paper)
        wb = xlrd.open_workbook(path,'r')
        sheet = wb.sheet_by_index(0)

        ds = {
            "PRE.APTAMER"   : sheet.col_values(colx=3, start_rowx=3, end_rowx=75),
            "APTAMER"       : sheet.col_values(colx=4, start_rowx=3, end_rowx=75),
            "POST.APTAMER"  : sheet.col_values(colx=5, start_rowx=3, end_rowx=75),
            "CDS"           : sheet.col_values(colx=6, start_rowx=3, end_rowx=75),
            "PROT.MEAN"     : sheet.col_values(colx=16, start_rowx=3, end_rowx=75),
            "PROT.STD"      : sheet.col_values(colx=17, start_rowx=3, end_rowx=75),
            "PROTEIN"       : sheet.col_values(colx=1, start_rowx=3, end_rowx=75),
            "ORGANISM"      : "Escherichia coli str. K-12 substr. DH10B",
            "METHOD"        : "Individually Characterized",
            "TEMP"          : 37.0,
            "DATASET"       : paper
        }

        ds["5pUTR"] = ["{}{}{}".format(pre,aptamer,post) for pre,aptamer,post \
                        in zip(ds["PRE.APTAMER"],ds["APTAMER"],ds["POST.APTAMER"])]
        ds["SEQUENCE"] = [UTR+CDS for UTR,CDS in zip(ds["5pUTR"],ds["CDS"])]
        ds["STARTPOS"] = [len(UTR) for UTR in ds["5pUTR"]]

        df = pd.DataFrame(ds)
        db = db.append(df, ignore_index=True)


    '''Amin Espah Borujeni, Howard M. Salis
    Translation Initiation is Controlled by RNA Folding Kinetics via a Ribosome Drafting Mechanism
    Journal of the American Chemical Society (JACS), 2016'''
    paper = 'EspahBorujeni_JACS_2016'

    if paper in datasets:
        path = 'datasets/{}.xls'.format(paper)
        wb = xlrd.open_workbook(path,'r')
        sheet = wb.sheet_by_index(0)

        ds = {
            "5pUTR"     : sheet.col_values(colx=3, start_rowx=6, end_rowx=42),
            "CDS"       : sheet.col_values(colx=4, start_rowx=6, end_rowx=42),
            "PROT.MEAN" : sheet.col_values(colx=9, start_rowx=6, end_rowx=42),
            "PROT.STD"  : sheet.col_values(colx=10, start_rowx=6, end_rowx=42),
            "PROTEIN"   : "RFP",
            "ORGANISM"  : "Escherichia coli str. K-12 substr. DH10B",
            "METHOD"    : "Individually Characterized",
            "TEMP"      : 37.0,
            "DATASET"   : paper
        }

        ds["SEQUENCE"] = [UTR+CDS for UTR,CDS in zip(ds["5pUTR"],ds["CDS"])]
        ds["STARTPOS"] = [len(UTR) for UTR in ds["5pUTR"]]

        df = pd.DataFrame(ds)
        db = db.append(df, ignore_index=True)        


    '''Amin Espah Borujeni, Daniel P. Cetnar, Howard M. Salis
    Precise Quantification of Translation Inhibition by RNA structures that Overlap with the Ribosome Footprint at N-terminal Coding Sections
    Nucleic Acids Research, 2017'''
    paper = 'EspahBorujeni_Footprint'

    if paper in datasets:
        path = 'datasets/{}.xls'.format(paper)
        wb = xlrd.open_workbook(path,'r')
        sheet = wb.sheet_by_index(0)

        ds = {
            "5pUTR"     : sheet.col_values(colx=3, start_rowx=5, end_rowx=32),
            "CDS"       : sheet.col_values(colx=4, start_rowx=5, end_rowx=32),
            "PROT.MEAN" : sheet.col_values(colx=10, start_rowx=5, end_rowx=32),
            "PROT.STD"  : sheet.col_values(colx=11, start_rowx=5, end_rowx=32),
            "PROTEIN"   : "RFP",
            "ORGANISM"  : "Escherichia coli str. K-12 substr. DH10B",
            "METHOD"    : "Individually Characterized",
            "TEMP"      : 37.0,
            "DATASET"   : paper
        }

        ds["SEQUENCE"] = [UTR+CDS for UTR,CDS in zip(ds["5pUTR"],ds["CDS"])]
        ds["STARTPOS"] = [len(UTR) for UTR in ds["5pUTR"]]

        df = pd.DataFrame(ds)
        db = db.append(df, ignore_index=True)


    '''
    Amin Espah Borujeni, Manish Kushwaha, and Howard M Salis
    Reparameterizing the Biophysical Model of Translation Initiation inside Bacillus subtilis
    Unpublished. Reported in Thesis of Amin Espah Borujeni Penn State University 2016'''
    paper = 'EspahBorujeni_Bsubtilis_2016'

    if paper in datasets:
        path = 'datasets/{}.xls'.format(paper)
        wb = xlrd.open_workbook(path,'r')
        sheet = wb.sheet_by_index(0)

        ds = {
            "5pUTR"       : sheet.col_values(colx=1, start_rowx=1, end_rowx=48),
            "CDS"         : sheet.col_values(colx=2, start_rowx=1, end_rowx=48),
            "PROT.MEAN"   : sheet.col_values(colx=3, start_rowx=1, end_rowx=48),
            "PROT.STD"    : sheet.col_values(colx=4, start_rowx=1, end_rowx=48),
            "PROTEIN"     : "RFP",
            "PLASMID"     : "pDG1661",
            "PROMOTER.ID" : "pVeg",
            "ORGANISM"    : "Bacillus subtilis subsp. subtilis str. 168",
            "METHOD"      : "Individually Characterized",
            "TEMP"        : 30.0,
            "DATASET"     : paper
        }

        ds["SEQUENCE"] = [UTR+CDS for UTR,CDS in zip(ds["5pUTR"],ds["CDS"])]
        ds["STARTPOS"] = [len(UTR) for UTR in ds["5pUTR"]]

        df = pd.DataFrame(ds)
        db = db.append(df, ignore_index=True)


    '''Howard M. Salis, Ethan A. Mirsky, & Christopher A. Voigt
    Automated design of synthetic ribosome binding sites to control protein expression
    Nature Biotechnology, 2009, Vol. 27, No. 10; doi: 10.1038/nbt.1568'''
    paper = 'Salis_Nat_Biotech_2009'

    if paper in datasets:
        path = 'datasets/{}.xls'.format(paper)
        wb = xlrd.open_workbook(path,'r')
        sheet = wb.sheet_by_index(0)

        # Need to parse this dataset for the start positions
        RFP1_CDS = sheet.cell_value(1,3)
        excel_seqs = sheet.col_values(colx=3, start_rowx=3, end_rowx=135)        

        SacI = 'GAGCTC'
        SacI_positions = [seq.find(SacI) for seq in excel_seqs]
        start_positions = [p - 5 for p in SacI_positions]
        
        for i in range(len(excel_seqs)):
            seq = excel_seqs[i]
            start_pos = start_positions[i]
            start = seq[start_pos:start_pos+3]
            while start not in ['ATG','GTG','CTG','TTG']:
                start_pos -= 3
                start = seq[start_pos:start_pos+3]
            start_positions[i] = start_pos

        ds = {
            "5pUTR"     : [seq[0:start_pos] for seq,start_pos in zip(excel_seqs,start_positions)],
            "CDS"       : [seq[start_pos:SacI_pos] + RFP1_CDS for seq,start_pos,SacI_pos in \
                            zip(excel_seqs,start_positions,SacI_positions)],
            "PROT.MEAN" : sheet.col_values(colx=4, start_rowx=3, end_rowx=135),
            "PROT.STD"  : sheet.col_values(colx=5, start_rowx=3, end_rowx=135),
            "PROTEIN"   : "RFP",
            "ORGANISM"  : "Escherichia coli str. K-12 substr. DH10B",
            "METHOD"    : "Individually Characterized",
            "TEMP"      : 37.0,
            "DATASET"   : paper
        }

        ds["SEQUENCE"] = [UTR+CDS for UTR,CDS in zip(ds["5pUTR"],ds["CDS"])]
        ds["STARTPOS"] = [len(UTR) for UTR in ds["5pUTR"]]

        df = pd.DataFrame(ds)
        db = db.append(df, ignore_index=True)


    '''Iman Farasat, Manish Kushwaha, Jason Collens, Michael Easterbrook, Matthew Guido, & Howard M. Salis
    Efficient search, mapping, and optimization of multi-protein genetic systems in diverse bacteria
    Molecular Systems Biology. 2014 Jun 21; 10:731. doi: 10.15252/msb.20134955'''
    paper = 'Farasat_MSB_2014'

    if paper in datasets:
        path = 'datasets/{}.xls'.format(paper)
        wb = xlrd.open_workbook(path,'r')
        sheet = wb.sheet_by_index(0)

        ds = {
            "PRESEQ"    : sheet.col_values(colx=5, start_rowx=1, end_rowx=146),
            "RBS"       : sheet.col_values(colx=6, start_rowx=1, end_rowx=146),
            "CDS"       : sheet.col_values(colx=7, start_rowx=1, end_rowx=146),
            "PROT.MEAN" : sheet.col_values(colx=9, start_rowx=1, end_rowx=146),
            "PROT.STD"  : sheet.col_values(colx=10, start_rowx=1, end_rowx=146),
            "PROTEIN"   : sheet.col_values(colx=3, start_rowx=1, end_rowx=146),
            "ORGANISM"  : sheet.col_values(colx=2, start_rowx=1, end_rowx=146),
            "METHOD"    : "Individually Characterized",
            "TEMP"      : 37.0,
            "DATASET"   : paper
        }

        ds["5pUTR"] = ["{}{}".format(preseq,RBS) for preseq,RBS in zip(ds['PRESEQ'],ds['RBS'])]
        ds["SEQUENCE"] = [UTR+CDS for UTR,CDS in zip(ds["5pUTR"],ds["CDS"])]
        ds["STARTPOS"] = [len(UTR) for UTR in ds["5pUTR"]]

        df = pd.DataFrame(ds)
        db = db.append(df, ignore_index=True)        


    '''Tian Tian, & Howard M. Salis
    A Predictive Biophysical Model of Translational Coupling to Coordinate and Control Protein expression in Bacterial Operons
    Nucleic Acids Research, 2015 doi: 10.1093/nar/gkv635'''
    paper = 'Tian_NAR_2015'

    if paper in datasets:
        path = 'datasets/{}.xls'.format(paper)
        wb = xlrd.open_workbook(path,'r')
        sheet = wb.sheet_by_index(0)

        ds = {
            "5pUTR"     : sheet.col_values(colx=4, start_rowx=2, end_rowx=26),
            "CDS"       : sheet.col_values(colx=5, start_rowx=2, end_rowx=26),
            "PROT.MEAN" : sheet.col_values(colx=12, start_rowx=2, end_rowx=26),
            "PROT.STD"  : sheet.col_values(colx=13, start_rowx=2, end_rowx=26),
            "PROTEIN"   : sheet.col_values(colx=3, start_rowx=2, end_rowx=26),
            "ORGANISM"  : sheet.col_values(colx=2, start_rowx=2, end_rowx=26),
            "METHOD"    : "Individually Characterized",
            "TEMP"      : 37.0,
            "DATASET"   : paper
        }

        ds["SEQUENCE"] = [UTR+CDS for UTR,CDS in zip(ds["5pUTR"],ds["CDS"])]
        ds["STARTPOS"] = [len(UTR) for UTR in ds["5pUTR"]]

        df = pd.DataFrame(ds)
        db = db.append(df, ignore_index=True)


    '''Mark Mimee, Alex C. Tucker, Christopher A. Voigt, and Timothy K. Lu
    Programming a Human Commensal Bacterium, Bacteroides thetaiotaomicron, to Sense and Respond to Stimuli in the Murine Gut Microbiota
    Cell Systems 1, 62-71, July 29, 2015'''
    paper = 'Mimee_Cell_Sys_2015'

    if paper in datasets:
        path = 'datasets/{}.xls'.format(paper)
        wb = xlrd.open_workbook(path,'r')
        sheet = wb.sheet_by_index(0)

        ds = {
            "5pUTR"     : sheet.col_values(colx=6, start_rowx=3, end_rowx=146),
            "CDS"       : sheet.col_values(colx=7, start_rowx=3, end_rowx=146),
            "PROT.MEAN" : sheet.col_values(colx=8, start_rowx=3, end_rowx=146),
            "PROT.STD"  : sheet.col_values(colx=9, start_rowx=3, end_rowx=146),
            "PROTEIN"   : "NanoLuc",
            "ORGANISM"  : "Bacteroides thetaiotaomicron VPI-5482",
            "METHOD"    : "Individually Characterized",
            "TEMP"      : 37.0,
            "DATASET"   : paper
        }

        ds["SEQUENCE"] = [UTR+CDS for UTR,CDS in zip(ds["5pUTR"],ds["CDS"])]
        ds["STARTPOS"] = [len(UTR) for UTR in ds["5pUTR"]]

        df = pd.DataFrame(ds)
        db = db.append(df, ignore_index=True)


    '''Mads T Bonde, Margit Pederse, Michael S Klausen, Sheila I Jensen, Tune Wulff, Scott Harrison, Alex T Nielsen, Markus J Herrgard, Morten O A Sommer
    Predictable tuning of protein expression in bacteria
    Nature Methods, 13: 3, 233-236, 2016'''
    # Only including the IC data for now from Bonde et al., Nat Methods
    paper = 'Bonde_NatMethods_IC_2016'

    if paper in datasets:
        path = 'datasets/{}.xls'.format(paper)
        wb = xlrd.open_workbook(path,'r')
        sheet = wb.sheet_by_index(0)

        ds = {
            "5pUTR"     : sheet.col_values(colx=5, start_rowx=1, end_rowx=107),
            "CDS"       : sheet.col_values(colx=6, start_rowx=1, end_rowx=107),
            "PROT.MEAN" : sheet.col_values(colx=11, start_rowx=1, end_rowx=107),
            "PROT.STD"  : None,
            "PROTEIN"   : "sfGFP",
            "ORGANISM"  : "Escherichia coli str. K-12 substr. MG1655",
            "METHOD"    : "Individually Characterized",
            "TEMP"      : 37.0,
            "DATASET"   : paper
        }

        ds["SEQUENCE"] = [UTR+CDS for UTR,CDS in zip(ds["5pUTR"],ds["CDS"])]
        ds["STARTPOS"] = [len(UTR) for UTR in ds["5pUTR"]]

        df = pd.DataFrame(ds)
        db = db.append(df, ignore_index=True)


    '''Robert G Egbert and Eric Klavins
    Fine-tuning gene networks using simple sequence repeats
    PNAS, 109: 42, 16817-16822, 2012'''
    paper = 'Egbert_PNAS_2012'

    if paper in datasets:
        path = 'datasets/{}.xls'.format(paper)
        wb = xlrd.open_workbook(path,'r')
        sheet = wb.sheet_by_index(0)

        ds = {
            "5pUTR"     : sheet.col_values(colx=9, start_rowx=3, end_rowx=45),
            "CDS"       : sheet.cell_value(47,1),
            "PROT.MEAN" : sheet.col_values(colx=6, start_rowx=3, end_rowx=45),
            "PROT.STD"  : sheet.col_values(colx=7, start_rowx=3, end_rowx=45),
            "PROTEIN"   : "GFP",
            "ORGANISM"  : "Escherichia coli str. K-12 substr. MG1655",
            "METHOD"    : "Individually Characterized",
            "TEMP"      : 37.0,
            "DATASET"   : paper
        }

        ds["SEQUENCE"] = [UTR+ds["CDS"] for UTR in ds["5pUTR"]]
        ds["STARTPOS"] = [len(UTR) for UTR in ds["5pUTR"]]

        df = pd.DataFrame(ds)
        db = db.append(df, ignore_index=True)

    '''Ariel Hecht, Jeff Glasgow, Paul R. Jaschke, Lukmaan A. Bawazer, Drew Endy, Marc Salit
    Measurements of translation initiation from all 64 codons in E. coli
    Nucleic Acids Research, 45: 7, 3615-3626, 2017'''
    paper = 'Hecht_NAR_2017'

    if paper in datasets:
        path = 'datasets/{}.xls'.format(paper)
        wb = xlrd.open_workbook(path,'r')
        sheet = wb.sheet_by_name('Figure 3B')
        ds = {
            "PROMOTER"    : sheet.cell_value(17,1),
            "5pUTR"       : sheet.cell_value(18,1),
            "CDS"         : sheet.cell_value(19,1),
            "START.CODON" : sheet.col_values(colx=0,  start_rowx=2, end_rowx=14),
            "PROT.MEAN"   : sheet.col_values(colx=10, start_rowx=2, end_rowx=14),
            "PROT.STD"    : sheet.col_values(colx=11, start_rowx=2, end_rowx=14),
            "PROTEIN"     : "NanoLuc",
            "PLASMID"     : "p15A",
            "ORGANISM"    : "Escherichia coli BL21(DE3)",
            "METHOD"      : "Individually Characterized",
            "TEMP"        : 37.0,
            "DATASET"     : paper
        }

        ds["SEQUENCE"] = [ds["5pUTR"]+START+ds["CDS"] for START in ds["START.CODON"]]
        ds["STARTPOS"] = len(ds["5pUTR"])

        df = pd.DataFrame(ds)
        db = db.append(df, ignore_index=True)

        sheet = wb.sheet_by_name('Figure 3C')
        ds = {
            "PROMOTER"    : sheet.cell_value(17,1),
            "5pUTR"       : sheet.cell_value(18,1),
            "CDS"         : sheet.cell_value(19,1),
            "START.CODON" : sheet.col_values(colx=0,  start_rowx=2, end_rowx=14),
            "PROT.MEAN"   : sheet.col_values(colx=10, start_rowx=2, end_rowx=14),
            "PROT.STD"    : sheet.col_values(colx=11, start_rowx=2, end_rowx=14),
            "PROTEIN"     : "NanoLuc",
            "PLASMID"     : "oriS/BAC",
            "ORGANISM"    : "Escherichia coli BL21(DE3)",
            "METHOD"      : "Individually Characterized",
            "TEMP"        : 37.0,
            "DATASET"     : paper
        }

        ds["SEQUENCE"] = [ds["5pUTR"]+START+ds["CDS"] for START in ds["START.CODON"]]
        ds["STARTPOS"] = len(ds["5pUTR"])        

        df = pd.DataFrame(ds)
        db = db.append(df, ignore_index=True)

        sheet = wb.sheet_by_name('Figure 2')
        ds = {
            "PROMOTER"    : sheet.cell_value(111,1),
            "5pUTR"       : sheet.cell_value(112,1),
            "CDS"         : sheet.cell_value(113,1),
            "START.CODON" : sheet.col_values(colx=0, start_rowx=96, end_rowx=107),
            "PROT.MEAN"   : sheet.col_values(colx=5, start_rowx=96, end_rowx=107),
            "PROT.STD"    : sheet.col_values(colx=6, start_rowx=96, end_rowx=107),
            "PROTEIN"     : "sfGFP",
            "PLASMID"     : "pET20b",
            "ORGANISM"    : "Escherichia coli BL21(DE3)",
            "METHOD"      : "Individually Characterized",
            "TEMP"        : 37.0,
            "DATASET"     : paper
        }

        ds["SEQUENCE"] = [ds["5pUTR"]+START+ds["CDS"] for START in ds["START.CODON"]]
        ds["STARTPOS"] = len(ds["5pUTR"])        

        df = pd.DataFrame(ds)
        db = db.append(df, ignore_index=True)

    '''Heather J. Beck, Ian M. C. Fleming, Gary R. Janssen
    5'-Terminal AUGs in Escherichia coli mRNAs with Shine-Dalgarno Sequences:
    Identification and Analysis of Their Roles in Non-Canonical Translation Initiation
    PLoS ONE 11(7), 2016'''
    paper = 'Beck_PLoS_2016'

    if paper in datasets:
        path = 'datasets/{}.xls'.format(paper)
        wb = xlrd.open_workbook(path,'r')
        sheet = wb.sheet_by_name('Sheet2')

        SalI = sheet.cell_value(43,2)
        end_row = 34 # removed the last 2 sequences b/c they are too long!
        ds = {
            "CDS"         : sheet.cell_value(44,2),
            "GENE"        : sheet.col_values(colx=1, start_rowx=1, end_rowx=end_row),
            "ORF"         : sheet.col_values(colx=3, start_rowx=1, end_rowx=end_row),
            "PROT.MEAN"   : sheet.col_values(colx=6, start_rowx=1, end_rowx=end_row),
            "PROT.STD"    : sheet.col_values(colx=7, start_rowx=1, end_rowx=end_row),
            "PROTEIN"     : "LacZ-fusion",
            "ORGANISM"    : "Escherichia coli str. K-12 substr. MG1655",
            "METHOD"      : "Individually Characterized",
            "TEMP"        : 30.0, # educated guess
            "DATASET"     : paper
        }

        ds["SEQUENCE"] = [ORF+SalI+ds["CDS"] for ORF in ds["ORF"]]
        ds["5pUTR"] = [ORF+SalI for ORF in ds["ORF"]]
        ds["STARTPOS"] = [len(UTR) for UTR in ds["ORF"]]

        df = pd.DataFrame(ds)
        db = db.append(df, ignore_index=True)

    '''Sriram Kosuri, Daniel B. Goodman, George M. Church
    Composability of regulatory sequences controlling transcription and translation in Escherichia coli
    Proc Natl Acad Sci USA, 2013, Vol. 110 no. 34'''
    paper = 'Kosuri_PNAS_2013'

    if paper in datasets:
        path = 'datasets/{}.xls'.format(paper)
        wb = xlrd.open_workbook(path,'r')

        # import promoter and RBS information
        sheet = wb.sheet_by_name("Promoters")
        promoter_dict = {}
        for i in range(1,113):
            promoter = sheet.cell_value(i,0)
            TSS = int(sheet.cell_value(i,7))
            promoter_dict[promoter] = {}
            promoter_dict[promoter]["SEQ"] = re.sub('[ \"]','',sheet.cell_value(i,9))
            promoter_dict[promoter]["TSS"] = TSS

        bad_promoters = ['"BBa_J23113"','"BBa_J23109"','"em7"']
        for p in bad_promoters:
            promoter_dict[p] = {'SEQ':"",'TSS':0}

        sheet = wb.sheet_by_name("RBSs")
        RBS_dict = {}
        for i in range(1,112):
            RBS = sheet.cell_value(i,0)
            RBS_dict[RBS] = (sheet.cell_value(i,9).replace(' ','').replace('"',''))[:-3]        

        # add template RBS
        RBS_dict['"em7"'] = ""

        # import Flow-seq required data to replicate calculations
        sheet = wb.sheet_by_name("NGS counts")

        data = {}
        data['total.rna.a'] = sheet.cell_value(0,1)
        data['total.rna.b'] = sheet.cell_value(1,1)
        data['total.dna.a'] = sheet.cell_value(2,1)
        data['total.dna.b'] = sheet.cell_value(3,1)
        
        bin_list = ['BIN'+str(num) for num in range(1,13)]
        data['bin_list'] = bin_list    
        
        sheet = wb.sheet_by_name("bin_vals")
        data['bins'] = {b: {} for b in bin_list}
        for b,col in zip(bin_list,range(1,13)):
            data['bins'][b]['cell_fraction'] = sheet.cell_value(2,col)
            data['bins'][b]['fluo'] = sheet.cell_value(4,col)

        # import dataset
        sheet = wb.sheet_by_index(0)
        rowf = 12656
        ds = {
            "PROMOTER.ID"   : sheet.col_values(colx= 0, start_rowx=1, end_rowx=rowf),
            "RBS.ID"        : sheet.col_values(colx= 1, start_rowx=1, end_rowx=rowf),
            "CDS"           : sheet.cell_value(1,85),
            "COUNT.PROTEIN" : sheet.col_values(colx= 3, start_rowx=1, end_rowx=rowf),
            "COUNT.A.RNA"   : sheet.col_values(colx=58, start_rowx=1, end_rowx=rowf),
            "COUNT.B.RNA"   : sheet.col_values(colx=59, start_rowx=1, end_rowx=rowf),
            "COUNT.RNA"     : sheet.col_values(colx=60, start_rowx=1, end_rowx=rowf),
            "COUNT.A.DNA"   : sheet.col_values(colx=34, start_rowx=1, end_rowx=rowf),
            "COUNT.B.DNA"   : sheet.col_values(colx=35, start_rowx=1, end_rowx=rowf),
            "COUNT.DNA"     : sheet.col_values(colx=36, start_rowx=1, end_rowx=rowf),
            "PROTEIN"       : "sfGFP",
            "PLASMID"       : "pEGRC",
            "ORGANISM"      : "Escherichia coli str. K-12 substr. MG1655",
            "METHOD"        : "Flow-seq",
            "TEMP"          : 30.0,
            "DATASET"       : paper
        }

        # get contig counts from Flow-seq
        for b,col in zip(bin_list,range(4,16)):
            ds[b] = sheet.col_values(colx=col, start_rowx=1, end_rowx=rowf)

        ds["PROMOTER"] = [promoter_dict[ID]['SEQ'] for ID in ds["PROMOTER.ID"]]
        ds["TSS"]      = [promoter_dict[ID]['TSS'] for ID in ds["PROMOTER.ID"]]
        ds["RBS"]      = [RBS_dict[ID] for ID in ds['RBS.ID']]
        ds["5pUTR"]    = [p[TSS:]+RBS for p,TSS,RBS in zip(ds['PROMOTER'],ds['TSS'],ds['RBS'])]
        ds["SEQUENCE"] = [UTR+ds["CDS"] for UTR in ds["5pUTR"]]
        ds["STARTPOS"] = [len(UTR) for UTR in ds["5pUTR"]]

        # convert to pandas dataframe
        df = pd.DataFrame(ds)

        # remove entries corresponding to J23113, J23109, and EM7
        df = df[~df['PROMOTER.ID'].isin(bad_promoters)]

        # convert to numeric values, and remove entries that have NaNs
        numbers = ["COUNT.PROTEIN","COUNT.A.RNA","COUNT.B.RNA","COUNT.RNA",
                   "COUNT.A.DNA", "COUNT.B.DNA", "COUNT.DNA"]
        numbers += bin_list
        for label in numbers:
            df[label] = pd.to_numeric(df[label],errors='coerce')
        df = df.dropna(axis=0, how='any')

        df = calc_Flowseq(df,data)
        df = filter_Flowseq(df)
        db = db.append(df, ignore_index=True)


    '''Daniel B. Goodman, George M. Church, Sriram Kosuri
    Causes and Effects of N-Terminal Codon Bias in Bacterial Genes
    Science, 2013, Vol. 342'''
    paper = 'Goodman_Science_2013'

    if paper in datasets:
        path = 'datasets/{}.xls'.format(paper)
        wb = xlrd.open_workbook(path,'r')

        # import Flow-seq required data to replicate calculations
        sheet = wb.sheet_by_name("NGS counts")
        data = {}
        data['total.rna.a'] = sheet.cell_value(0,1)
        data['total.rna.b'] = sheet.cell_value(1,1)
        data['total.dna.a'] = sheet.cell_value(2,1)
        data['total.dna.b'] = sheet.cell_value(3,1)

        bin_list = ['BIN'+str(num) for num in range(1,13)]
        data['bin_list'] = bin_list    

        sheet = wb.sheet_by_name("bin_vals")
        data['bins'] = {b: {} for b in bin_list}
        for b,col in zip(bin_list,range(1,13)):
            data['bins'][b]['cell_fraction'] = sheet.cell_value(2,col)
            data['bins'][b]['fluo'] = sheet.cell_value(4,col)        

        # import dataset
        sheet = wb.sheet_by_name("Goodman_data")
        rowf = 14235
        ds = {
            "PROMOTER"      : sheet.col_values(colx=27, start_rowx=1, end_rowx=rowf),
            "TSS"           : [int(round(x)) if isinstance(x,float) else None for x in sheet.col_values(colx=24, start_rowx=1, end_rowx=rowf)],
            "RBS"           : sheet.col_values(colx=28, start_rowx=1, end_rowx=rowf),
            "CDS"           : sheet.col_values(colx=76, start_rowx=1, end_rowx=rowf),
            "N.TERMINAL.CDS": sheet.col_values(colx= 1, start_rowx=1, end_rowx=rowf),
            "COUNT.PROTEIN" : sheet.col_values(colx=23, start_rowx=1, end_rowx=rowf),
            "COUNT.DNA"     : sheet.col_values(colx= 3, start_rowx=1, end_rowx=rowf),
            "COUNT.RNA"     : sheet.col_values(colx= 4, start_rowx=1, end_rowx=rowf),
            "COUNT.A.DNA"   : sheet.col_values(colx= 5, start_rowx=1, end_rowx=rowf),
            "COUNT.A.RNA"   : sheet.col_values(colx= 6, start_rowx=1, end_rowx=rowf),
            "COUNT.B.DNA"   : sheet.col_values(colx= 7, start_rowx=1, end_rowx=rowf),
            "COUNT.B.RNA"   : sheet.col_values(colx= 8, start_rowx=1, end_rowx=rowf),
            "PROTEIN"       : "sfGFP",
            "PLASMID"       : "pEGRC",
            "ORGANISM"      : "Escherichia coli str. K-12 substr. MG1655",
            "METHOD"        : "Flow-seq",
            "TEMP"          : 30.0,
            "DATASET"       : paper
        }

        # get contig counts from Flow-seq
        for b,col in zip(bin_list,range(4,16)):
            ds[b] = np.array(sheet.col_values(colx=col, start_rowx=1, end_rowx=rowf))

        ds["5pUTR"] = [p[TSS:]+RBS+N+"CAT" for p,TSS,RBS,N in \
                                 zip(ds['PROMOTER'],ds['TSS'],ds['RBS'],ds['N.TERMINAL.CDS'])]
        ds["SEQUENCE"] = [UTR+CDS for UTR,CDS in zip(ds["5pUTR"],ds["CDS"])]
        ds["STARTPOS"] = [len(UTR) for UTR in ds["5pUTR"]]

        # convert to pandas dataframe
        df = pd.DataFrame(ds)

        # remove NA TSS entries
        df = df[~df["TSS"].isnull()]

        # convert to numeric values, and remove entries that have NaNs
        numbers = ["COUNT.PROTEIN","COUNT.A.RNA","COUNT.B.RNA","COUNT.RNA",
                   "COUNT.A.DNA", "COUNT.B.DNA", "COUNT.DNA"]
        numbers += bin_list
        for label in numbers:
            df[label] = pd.to_numeric(df[label],errors='coerce')
        df = df.dropna(axis=0, how='any')
        
        df = calc_Flowseq(df,data)
        df = filter_Flowseq(df)
        db = db.append(df, ignore_index=True)


    '''Mads T Bonde, Margit Pederse, Michael S Klausen, Sheila I Jensen, Tune Wulff, Scott Harrison, Alex T Nielsen, Markus J Herrgard, Morten O A Sommer
    Predictable tuning of protein expression in bacteria
    Nature Methods, 13: 3, 233-236, 2016'''
    # FS data from Bonde et al., Nat Methods
    paper = 'Bonde_NatMethods_FS_2016'

    if paper in datasets:
        path = 'datasets/{}.xls'.format(paper)
        wb = xlrd.open_workbook(path,'r')
        sheet = wb.sheet_by_index(0)

        ds = {
            "SD sequence": sheet.col_values(colx=0, start_rowx=1, end_rowx=3864),
            "5pUTR"      : sheet.col_values(colx=4, start_rowx=1, end_rowx=3864),
            "CDS"        : sheet.col_values(colx=5, start_rowx=1, end_rowx=3864),
            "PROT.MEAN"  : sheet.col_values(colx=6, start_rowx=1, end_rowx=3864),
            "Total count": sheet.col_values(colx=3, start_rowx=1, end_rowx=3864),
            "PROT.STD"   : None,
            "PROTEIN"    : "sfGFP",
            "ORGANISM"   : "Escherichia coli str. K-12 substr. MG1655",
            "METHOD"     : "Flow-seq",
            "TEMP"       : 37.0,
            "PLASMID"    : "pMA1",
            "DATASET"    : paper
        }

        ds["SEQUENCE"] = [UTR+CDS for UTR,CDS in zip(ds["5pUTR"],ds["CDS"])]
        ds["STARTPOS"] = [len(UTR) for UTR in ds["5pUTR"]]

        df = pd.DataFrame(ds)

        # only take variants with >50 reads
        df = df[df['Total count']>50]

        # compute "normalized" expression values
        Exp_TTGGGC = float(df.loc[df['SD sequence']=='TTGGGC']['PROT.MEAN'])
        Exp_AGGAGA = float(df.loc[df['SD sequence']=='AGGAGA']['PROT.MEAN'])

        print Exp_TTGGGC
        print Exp_AGGAGA

        df['PROT.PERCENT'] = (df['PROT.MEAN']-Exp_TTGGGC)/(Exp_AGGAGA-Exp_TTGGGC)

        db = db.append(df, ignore_index=True)

        print db.head()

    '''Data collected for this work'''
    for paper in ['Reis1_2018','Reis2_2018']:

        if paper in datasets:
            path = 'datasets/{}.csv'.format(paper)
            tempdf = pd.read_csv(path)

            with open('datasets/Reis_mRFP1.fasta') as handle:
                for line in handle:
                    CDS = str(line.rstrip()).upper()

            df = pd.DataFrame()
            df['5pUTR']     = tempdf['UTR']
            df['CDS']       = tempdf['CDS'].str.upper() + CDS
            df['PROT.MEAN'] = tempdf['RFP.mean']
            df['PROT.STD']  = tempdf['RFP.std']
            df['PROTEIN']   = 'RFP-modified'
            df['ORGANISM']  = 'Escherichia coli str. K-12 substr. DH10B'
            df['METHOD']    = 'Individually Characterized'
            df['TEMP']      = 37.0
            df['DATASET']   = paper

            df['SEQUENCE'] = df['5pUTR'] + df['CDS']
            df['STARTPOS'] = df['5pUTR'].map(len)

            db = db.append(df, ignore_index=True)

    # Clean up, define categories based on organism/host/dataset
    db = _make_categories(db)

    # Remove unicode strings
    db = _remove_unicode(db)
    
    return db

def _make_categories(db):
    # categories save some memory

    db["PROTEIN"]  = db["PROTEIN"].astype('category')
    db["ORGANISM"] = db["ORGANISM"].astype('category')
    db["METHOD"]   = db["METHOD"].astype('category')
    db["DATASET"]  = db["DATASET"].astype('category')

    # And let's define sub-groups of sequences categorized:
    # At the same time, in the same organism, with the same promoter, and same experimental conditions
    info = ["{}+{}+{}+{}".format(d,o,g,p) for d,o,g,p in zip(db["DATASET"],db["ORGANISM"],db["PROTEIN"],db["PLASMID"])]
    db["SUBGROUP"] = pd.Series(info, dtype="category")
    return db

def _remove_unicode(db):
    for label in db.keys():
        if isinstance(db[label][0],unicode):
            db[label] = db[label].str.encode('utf-8')
    return db

def filter_Flowseq(df):

    # Remove poor quality data from the Flow-seq datasets
    size0 = len(df)

    lowDNA = (df['COUNT.A.DNA']<10) | (df['COUNT.B.DNA']<10)
    lowRNA = (df['COUNT.A.RNA']<10) | (df['COUNT.B.RNA']<10)
    lowSEQ = (df['COUNT.DNA']<50) & (df['COUNT.RNA']<20)
    lowPRO = df['COUNT.PROTEIN']<1000
    loPROT = df['PROT.MEAN'] < 1584
    hiPROT = df['BIN12']/df['COUNT.PROTEIN'] > 0.95
    # loPROT = df['BIN1']/df['COUNT.PROTEIN']>0.025
    # hiPROT = df['BIN12']/df['COUNT.PROTEIN']>0.025
    df = df[~(lowDNA | lowRNA | lowSEQ | lowPRO | loPROT | hiPROT)]

    sizef = len(df)
    print "Filtered out {} many poor quality sequences. Final dataset size: {}.".format(size0-sizef,sizef)

    return df

def calc_Flowseq(ds,data):
    
    # Calculate RNA levels from NGS data
    # Two replicates, A & B
    ds['RNA.A'] = (ds['COUNT.A.RNA']/data['total.rna.a']) \
                 /(ds['COUNT.A.DNA']/data['total.dna.a'])
    
    ds['RNA.B'] = (ds['COUNT.B.RNA']/data['total.rna.a']) \
                 /(ds['COUNT.B.DNA']/data['total.dna.b'])
    
    ds['RNA'] = (ds['RNA.A'] + ds['RNA.B'])/2.0
    
    ds['RNA.VAR'] = ((ds['RNA.A']-ds['RNA'])**2.0 + \
                     (ds['RNA.B']-ds['RNA'])**2.0) / 2.0
    
    ds['RNA.STD'] = np.sqrt(ds['RNA.VAR'])

    # Calculate the reconstructed protein levels
    # calculate the normalized fractional contribution of each bin j per sequence i (a_ij)
    a = {}
    num_seqs = len(ds["PROMOTER"])
    denominator = np.zeros(num_seqs)
    for b in data['bin_list']:
        denominator += data['bins'][b]['cell_fraction']*ds[b]/ds['COUNT.PROTEIN']
    
    # Geometric mean to calculate protein level
    ds["PROT.MEAN"] = np.ones(num_seqs)
    for b in data['bin_list']:
        a[b] = (data['bins'][b]['cell_fraction']*ds[b]/ds['COUNT.PROTEIN'])/denominator
        ds["PROT.MEAN"] *= np.exp( a[b]*np.log(data['bins'][b]['fluo']) )
    
    '''
    # Arithmetic mean of protein level
    ds["PROT.MEAN"] = np.zeros(num_seqs)
    for b in data['bin_list']:
        ds["PROT.MEAN"] += a[b]*data['bins'][b]['fluo']
    '''

    # Variance & standard deviation of linear-scaled protein data
    var = 0.0
    for b in data['bin_list']:
        var += a[b]*( data['bins'][b]['fluo'] - ds["PROT.MEAN"] )**2.0
    ds["PROT.VAR"] = var
    ds["PROT.STD"] = np.sqrt(ds["PROT.VAR"])
    
    '''
    # Calculate variance on a log-scale
    var2 = 0.0
    for b in data['bin_list']:
        var2 += a[b]*( np.log(data['bins'][b]['fluo']) - np.log(ds["PROT.MEAN"]) )**2.0

    print len(ds["PROT.MEAN"]),np.mean(ds["PROT.MEAN"]),np.sqrt(np.mean(ds["PROT.VAR"])),np.sqrt(np.mean(var2))
    '''

    # Calculate apparent translation rate
    ds["TRANSL.RATE"] = ds["PROT.MEAN"]/ds["RNA"]
    ds["TRANSL.RATE.VAR"] = approx_var_ratio(ds["PROT.MEAN"],ds["PROT.VAR"],ds["RNA"],ds["RNA.VAR"])
    ds["TRANSL.RATE.STD"] = np.sqrt(ds["TRANSL.RATE.VAR"])
    
    # Final steps
    for key,val in data.iteritems():
        if isinstance(val,np.ndarray):
            data[key] = list(val)

    return ds

def approx_var_ratio(mu1,var1,mu2,var2):
    # assumed cov(x,y) = 0
    var = (mu2**2.0 * var1**2.0 + mu1**2.0 * var2**2.0) / mu2**2.0
    return var


if __name__ == "__main__":

    datasets = ['EspahBorujeni_NAR_2013',   # StandbySite
                'EspahBorujeni_NAR_2015',   # Riboswitch
                'EspahBorujeni_JACS_2016',  # RibosomeDrafting
                'EspahBorujeni_Footprint',  # Footprint
                'EspahBorujeni_Bsubtilis_2016', # B. subtilis
                'Salis_Nat_Biotech_2009',   # RBSCalculator
                'Farasat_MSB_2014',         # dRBS
                'Tian_NAR_2015',            # TranslationCoupling
                'Mimee_Cell_Sys_2015',      # Bthetaiotaomicron
                'Bonde_NatMethods_IC_2016', # EMOPEC
                'Egbert_PNAS_2012',         # SpacerComposition
                'Hecht_NAR_2017',           # HechtStartCodons
                'Beck_PLoS_2016',           # Beck Leaderless mRNA
                'Kosuri_PNAS_2013',         # Kosuri
                'Goodman_Science_2013',      # Goodman
                'Reis1_2018',
                'Reis2_2018'
                ]

    # datasets = ['Bonde_NatMethods_FS_2016']

    db = pd.DataFrame()
    db = add_dataset(db,datasets)
    
    handle = open('geneticsystems.db','w')
    pickle.dump(db,handle,protocol=2)
    handle.close()