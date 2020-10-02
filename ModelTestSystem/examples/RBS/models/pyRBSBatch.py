
import sys,os,re
import cStringIO
import multiprocessing as mp
import pickle
import subprocess
from subprocess import Popen, PIPE, STDOUT
import uuid

import numpy as np

def RBSBatch(data):
    # Windows RBSBatch (RBS Designer) wrapper
    
    mRNA,fivepUTR,CDS,rRNA,paper,i,j = data
    
    ppath=os.path.normpath("C:/Program Files (x86)/UNAFold/bin")
    fpath = os.getcwd()
    
    filein = "{}.txt".format(uuid.uuid4())
    fileout = "{}.txt".format(uuid.uuid4())
    
    os.chdir(fpath)
    f = open(filein,"w")
    f.write("name=temp\r\n")
    f.write("antiSD={}\r\n".format(rRNA))
    f.write("utr={}\r\n".format(fivepUTR))
    f.write("cds={}".format(CDS))
    f.close()
    
    os.chdir(ppath)
    program = "C:/Program Files (x86)/UNAFold/bin/RBSBatch.exe"

    stdin  = "C:/Users/areis_000/Dropbox/penn_state/hsalis_lab/proj_mRNA_database/RBS Designer/{}".format(filein)
    stdout = "C:/Users/areis_000/Dropbox/penn_state/hsalis_lab/proj_mRNA_database/RBS Designer/{}".format(fileout)
    
    # Call RBSBatch @ command line
    args = [os.path.normpath(program), os.path.normpath(stdin), os.path.normpath(stdout)]
    _ = subprocess.call(args);
    
    os.chdir(fpath)
    f = open(fileout,"r")
    _ = f.readline()
    strs = f.readline()
    output = re.split(r'\t',strs)
    f.close()
    
    data = {}
    data['translation efficiency'] = float(output[1])
    data['Expression'] = float(output[1])*1e5
    data['exposure probability'] = float(output[2])
    data['mRNA-ribosome probability'] = float(output[3])
    data['SD'] = output[4]
    data['dG_SD:aSD'] = float(output[5])
    data['spacer length'] = float(output[6].rstrip())
    
    print paper,i,j
    
    os.remove(filein)
    os.remove(fileout)
    
    return data
    
def main():
    
    # open mRNA_database.p
    handle = open('mRNA_database.p','r')
    db = pickle.load(handle)
    handle.close()
    
    # open RBSDesigner dictionary with all mRNA structures
    # handle = open('RBSDesigner_ALL.p','r')
    # RBSDesigner = pickle.load(handle)
    # handle.close()
    RBSDesigner = {}
    
    # for-loop; run RBSBatch for each mRNA sequence in database
    # save information to dictionary with each key as the 5'UTR+CDS
    
    poolin = []
    for paper in db.keys():
        for i in range(len(db[paper]['xls_ID'])):
            for j in range(len(db[paper]['xls_ID'][i])):
                
                mRNA = db[paper]['mRNA'][i][j]
                fivepUTR = db[paper]['fivepUTR'][i][j]
                CDS = db[paper]['CDS'][i][j]
                rRNA = db[paper]['rRNA'][i][j]
                
                if CDS[0:3] == 'ATG':
                    poolin.append((mRNA,fivepUTR,CDS,rRNA,paper,i,j))
    
    # Run with multiprocessing 
    pool = mp.Pool(processes=mp.cpu_count())
    output = pool.map(RBSBatch,poolin)
    pool.close()
    pool.join()
    
    for i in range(len(poolin)):
        mRNA = poolin[i][0]
        RBSDesigner[mRNA] = output[i]
    
    handle = open('RBSDesigner_ALL.p','wb')
    pickle.dump(RBSDesigner,handle)
    handle.close()


if __name__ == "__main__":
    
    # NEED TO RUN CMD AS ADMIN
    # Make sure admin account is setup with password
    # run Command Prompt (as Admin)
    # runas /user:Administrator cmd
    # net user Administrator [password] (make sure to set an Admin pw)
    
    # NEED TO MAKE SURE NUMPY IS INSTALLED ON WINDOWS
    # Make sure to get binary corresponding to your system & Python install
    # http://www.lfd.uci.edu/~gohlke/pythonlibs/#numpy
    # python -m pip install --upgrade pip
    # python -m pip install numpy-1.10.4+mkl-cp27-cp27m-win32.whl
    
    main()
