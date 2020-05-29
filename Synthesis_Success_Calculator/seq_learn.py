import sys, os, random, traceback, math, itertools, matplotlib

import numpy as np
import pandas as pd
import sklearn as sk
import scipy as sp
from scipy import stats
from ast import literal_eval
import pickle as pkl

from sklearn.linear_model import LinearRegression

from sklearn.ensemble              import RandomForestClassifier     as RFC 
from sklearn.ensemble              import AdaBoostClassifier         as Ada
from sklearn.model_selection       import train_test_split, RandomizedSearchCV,StratifiedShuffleSplit,GridSearchCV
from sklearn.metrics               import f1_score, accuracy_score, precision_score, recall_score, make_scorer, matthews_corrcoef, cohen_kappa_score, roc_auc_score
from sklearn.metrics import confusion_matrix

from sklearn.ensemble.forest import _generate_unsampled_indices
from sklearn.base import clone
from multiprocessing import Pool 

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from IPython.display import display, HTML
from collections import OrderedDict

from seq_eval import *

matplotlib.rc('font', **{'family': 'Arial', 'size': 12})
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['axes.facecolor'] = 'white'

sns.set(font_scale=1.5)
sns.set_style("whitegrid")

def extract_to_csv_train(directory,source):
    frame = pd.read_csv(directory+source+".csv")
    frame = frame[frame.Status!="Other"]
    results_dict = OrderedDict()
    repeat_dict = OrderedDict()
    hairpin_dict = OrderedDict()
    nucleotide_dict = OrderedDict()
    sequence_dict = OrderedDict()
    for j in range(len(frame["Name"])):
        print "Running sequence ", str(j)
        i=frame['Name'][j]
        sequence_dict[i] = [frame['Sequence'][j]]
        results_dict[i] = seq_evaluate(frame['Sequence'][j])
        repeat_dict[i] = results_dict[i]["r_metrics"]
        hairpin_dict[i] = results_dict[i]["h_metrics"]
        nucleotide_dict[i] = results_dict[i]["n_metrics"]
    for i,r in frame.iterrows():
        if r["Status"] == "Delayed":
            results_dict[r["Name"]]["Status"]= "Canceled"
            repeat_dict[r["Name"]]["Status"]= "Canceled"
            hairpin_dict[r["Name"]]["Status"]= "Canceled"
            nucleotide_dict[r["Name"]]["Status"]= "Canceled"
        elif r["Status"] == "Synthesized":
            results_dict[r["Name"]]["Status"]= "Synthesized"
            repeat_dict[r["Name"]]["Status"]= "Synthesized"
            hairpin_dict[r["Name"]]["Status"]= "Synthesized"
            nucleotide_dict[r["Name"]]["Status"]= "Synthesized"
        elif r["Status"] == "Canceled":
            results_dict[r["Name"]]["Status"]= "Canceled"
            repeat_dict[r["Name"]]["Status"]= "Canceled"
            hairpin_dict[r["Name"]]["Status"]= "Canceled"
            nucleotide_dict[r["Name"]]["Status"]= "Canceled"
        else:
            continue

    sequence_frame = pd.DataFrame(sequence_dict).T
    repeat_frame = pd.DataFrame(repeat_dict).T
    hairpin_frame = pd.DataFrame(hairpin_dict).T
    nucleotide_frame = pd.DataFrame(nucleotide_dict).T

    sequence_frame.to_csv(source+"_Sequences.csv")
    repeat_frame.to_csv(source+"_Analysis_Repeat.csv")
    hairpin_frame.to_csv(source+"_Analysis_Hairpin.csv")
    nucleotide_frame.to_csv(source+"_Analysis_Nucleotide.csv")

def extract_to_csv_test(directory,source_list):
    try:
        
        #sets up multiprocessing via MPI
        from mpi4py import MPI
        from MPI_pool import Pool as MPIPool
        import copy_reg
        import cPickle as pickle

        # Using copy_reg to allow the pickle module to instance methods
        # Replicates multiprocessing module's ForkingPickler

        #copy_reg.pickle(types.MethodType, _reduce_method)
        print "test"
        pool = MPIPool(MPI.COMM_WORLD)
        print "MPI Pool Started"

    except:
        pool = None
        print traceback.format_exc()
        print "Could not start MPI Pool. Using serial map instead."

    if pool is not None:
        pool.start()
        if pool.is_master():
            map_function = pool.map
            print 'INFO: Using MPI Pool'
    else:
        map_function = map
        print 'INFO: Using Serial Mode'

    if pool is None or pool.is_master():
        for source in source_list:
            frame = pd.read_csv(directory+source+".csv")
            if "Status" in frame.keys():
                frame = frame[frame.Status!="Other"]
            results_dict = OrderedDict()
            repeat_dict = OrderedDict()
            hairpin_dict = OrderedDict()
            nucleotide_dict = OrderedDict()
            sequence_dict = OrderedDict()
            input_tups = [(frame["Sequence"][j],j) for j in range(len(frame))]

            output_list=map_function(extraction_subroutine,input_tups)

            for o in output_list:
                i =frame["Name"][o[0]]
                results_dict[i] = o[2]
                sequence_dict[i]=[o[1]]
                repeat_dict[i] = results_dict[i]["r_metrics"]
                hairpin_dict[i] = results_dict[i]["h_metrics"]
                nucleotide_dict[i] = results_dict[i]["n_metrics"]
            if "Status" in frame.keys():
                for i,r in frame.iterrows():
                    if r["Status"] == "Delayed":
                        results_dict[r["Name"]]["Status"]= "Canceled"
                        repeat_dict[r["Name"]]["Status"]= "Canceled"
                        hairpin_dict[r["Name"]]["Status"]= "Canceled"
                        nucleotide_dict[r["Name"]]["Status"]= "Canceled"
                    elif r["Status"] == "Synthesized":
                        results_dict[r["Name"]]["Status"]= "Synthesized"
                        repeat_dict[r["Name"]]["Status"]= "Synthesized"
                        hairpin_dict[r["Name"]]["Status"]= "Synthesized"
                        nucleotide_dict[r["Name"]]["Status"]= "Synthesized"
                    elif r["Status"] == "Canceled":
                        results_dict[r["Name"]]["Status"]= "Canceled"
                        repeat_dict[r["Name"]]["Status"]= "Canceled"
                        hairpin_dict[r["Name"]]["Status"]= "Canceled"
                        nucleotide_dict[r["Name"]]["Status"]= "Canceled"
                    else:
                        continue
            sequence_frame = pd.DataFrame(sequence_dict).T
            repeat_frame = pd.DataFrame(repeat_dict).T
            hairpin_frame = pd.DataFrame(hairpin_dict).T
            nucleotide_frame = pd.DataFrame(nucleotide_dict).T

            sequence_frame.to_csv(source+"_f_Sequences.csv")
            repeat_frame.to_csv(source+"_f_Analysis_Repeat.csv")
            hairpin_frame.to_csv(source+"_f_Analysis_Hairpin.csv")
            nucleotide_frame.to_csv(source+"_f_Analysis_Nucleotide.csv")
        pool.stop()
    else:
        pass

def extraction_subroutine((i,j)):
    print "Running sequence ", str(j)
    result = seq_evaluate(i)
    return (j,i,result)

def load_to_data_frames(source):

    #product_frame= pd.read_csv(source+"_Analysis.csv")
    repeat_frame=pd.read_csv(source+"_Analysis_Repeat.csv",index_col=[0])
    hairpin_frame=pd.read_csv(source+"_Analysis_Hairpin.csv",index_col=[0])
    nucleotide_frame=pd.read_csv(source+"_Analysis_Nucleotide.csv",index_col=[0])

    df_y = repeat_frame.loc[:,'Status']
    df_r_x = repeat_frame[[i for i in list(repeat_frame.columns) if i != 'Status']]
    #df_r_x = df_r_x.loc[:, (df_r_x != 0).any(axis=0)]
    df_h_x = hairpin_frame[[i for i in list(hairpin_frame.columns) if i != 'Status']]
    #df_h_x = df_h_x.loc[:, (df_h_x != 0).any(axis=0)]
    df_n_x = nucleotide_frame[[ i for i in list(nucleotide_frame.columns) if i != 'Status']]
    #df_n_x = df_n_x.loc[:, (df_n_x != 0).any(axis=0)]
    df_all_x = pd.concat([df_r_x,df_h_x,df_n_x],axis=1)
    #df_all_x = df_all_x.loc[:, (df_all_x != 0).any(axis=0)]

    #return df_y, df_r_x, df_h_x, df_n_x, df_all_x
    return df_y, df_all_x


def concat_sets(source_list):

    for j in xrange(len(source_list)):
        if j ==0:
            #product_frame= pd.read_csv(source_list[j]+"_Analysis.csv")
            repeat_frame=pd.read_csv(source_list[j]+"_Analysis_Repeat.csv")
            hairpin_frame=pd.read_csv(source_list[j]+"_Analysis_Hairpin.csv")
            nucleotide_frame=pd.read_csv(source_list[j]+"_Analysis_Nucleotide.csv")

            df_y = repeat_frame.loc[:,'Status']
            df_r_x = repeat_frame[[i for i in list(repeat_frame.columns) if i != 'Status']].iloc[:,1:]
            #df_r_x = df_r_x.loc[:, (df_r_x != 0).any(axis=0)]
            df_h_x = hairpin_frame[[i for i in list(hairpin_frame.columns) if i != 'Status']].iloc[:,1:]
            #df_h_x = df_h_x.loc[:, (df_h_x != 0).any(axis=0)]
            df_n_x = nucleotide_frame[[ i for i in list(nucleotide_frame.columns) if i != 'Status']].iloc[:,1:]
            #df_n_x = df_n_x.loc[:, (df_n_x != 0).any(axis=0)]
            df_all_x = pd.concat([df_r_x,df_h_x,df_n_x],axis=1)
            #df_all_x = df_all_x.loc[:, (df_all_x != 0).any(axis=0)]
            #print source_list[j]
            #print len(df_r_x)
            #print len(df_all_x)
            #print len(df_y)
        else:
            #product_frame= pd.read_csv(source_list[j]+"_Analysis.csv")
            repeat_frame=pd.read_csv(source_list[j]+"_Analysis_Repeat.csv")
            hairpin_frame=pd.read_csv(source_list[j]+"_Analysis_Hairpin.csv")
            nucleotide_frame=pd.read_csv(source_list[j]+"_Analysis_Nucleotide.csv")

            df_y_0 = repeat_frame.loc[:,'Status']
            df_r_x_0 = repeat_frame[[i for i in list(repeat_frame.columns) if i != 'Status']].iloc[:,1:]
            #df_r_x = df_r_x.loc[:, (df_r_x != 0).any(axis=0)]
            df_h_x_0 = hairpin_frame[[i for i in list(hairpin_frame.columns) if i != 'Status']].iloc[:,1:]
            #df_h_x = df_h_x.loc[:, (df_h_x != 0).any(axis=0)]
            df_n_x_0 = nucleotide_frame[[ i for i in list(nucleotide_frame.columns) if i != 'Status']].iloc[:,1:]
            #df_n_x = df_n_x.loc[:, (df_n_x != 0).any(axis=0)]
            df_all_x_0 = pd.concat([df_r_x_0,df_h_x_0,df_n_x_0],axis=1)
            #df_all_x = df_all_x.loc[:, (df_all_x != 0).any(axis=0)]
            #print source_list[j]
            #print len(df_r_x_0)
            #print len(df_all_x_0)
            #print len(df_y_0)

            df_y = pd.concat([df_y,df_y_0],axis=0)
            df_r_x = pd.concat([df_r_x,df_r_x_0],axis=0)
            #df_r_x = df_r_x.loc[:, (df_r_x != 0).any(axis=0)]
            df_h_x = pd.concat([df_h_x,df_h_x_0],axis=0)
            #df_h_x = df_h_x.loc[:, (df_h_x != 0).any(axis=0)]
            df_n_x = pd.concat([df_n_x,df_n_x_0],axis=0)
            #df_n_x = df_n_x.loc[:, (df_n_x != 0).any(axis=0)]
            df_all_x = pd.concat([df_all_x,df_all_x_0],axis=0)
            #df_all_x = df_all_x.loc[:, (df_all_x != 0).any(axis=0)]

    #df_r_x = df_r_x.loc[:, (df_r_x != 0).any(axis=0)]
    #df_h_x = df_h_x.loc[:, (df_h_x != 0).any(axis=0)]
    #df_n_x = df_n_x.loc[:, (df_n_x != 0).any(axis=0)]
    df_all_x = df_all_x.loc[:, (df_all_x != 0).any(axis=0)]

    return df_y, df_all_x

def split_shuffle_sets(source_list,test_frac=0.443, validation_size= 76,sublist=None):
    rand_state=stream_random_states()
    rand_state2=stream_random_states()
    #get final columns labels
    df_y_0, df_x_0 = concat_sets(source_list)
    total_set_size = len(df_y_0)
    print total_set_size
    if sublist:
        features = sublist
    else:
        features = df_x_0.keys()
    #make frame lists
    fy_list =[]
    fx_list = []
    val_size = []
    for j in xrange(len(source_list)):
        
        repeat_frame=pd.read_csv(source_list[j]+"_Analysis_Repeat.csv")
        hairpin_frame=pd.read_csv(source_list[j]+"_Analysis_Hairpin.csv")
        nucleotide_frame=pd.read_csv(source_list[j]+"_Analysis_Nucleotide.csv")
        df_y = repeat_frame.loc[:,'Status']
        df_r_x = repeat_frame[[i for i in list(repeat_frame.columns) if i != 'Status']].iloc[:,1:]
        #df_r_x = df_r_x.loc[:, (df_r_x != 0).any(axis=0)]
        df_h_x = hairpin_frame[[i for i in list(hairpin_frame.columns) if i != 'Status']].iloc[:,1:]
        #df_h_x = df_h_x.loc[:, (df_h_x != 0).any(axis=0)]
        df_n_x = nucleotide_frame[[ i for i in list(nucleotide_frame.columns) if i != 'Status']].iloc[:,1:]
        #df_n_x = df_n_x.loc[:, (df_n_x != 0).any(axis=0)]
        df_all_x = pd.concat([df_r_x,df_h_x,df_n_x],axis=1)
        df_all_x = df_all_x[[i for i in features]]
        #df_all_x['Random']=np.random.random(size=len(df_all_x))
        fy_list.append(df_y)
        fx_list.append(df_all_x)
        val_size.append(int(validation_size*len(df_y)/total_set_size))
    print val_size

    X_train_list = []
    X_test_list = []
    X_val_list = []
    y_val_list = []
    y_train_list = []
    y_test_list = []
    train_sizes = []
    for j in xrange(len(source_list)):
        if len(set(list(fy_list[j])))<2 :
            X_train, X_test, y_train, y_test = train_test_split(fx_list[j], fy_list[j], test_size=test_frac, random_state=rand_state.next(),shuffle=True)
        else:
            X_train, X_test, y_train, y_test = train_test_split(fx_list[j], fy_list[j], test_size=test_frac, random_state=rand_state.next(), stratify = fy_list[j],shuffle=True)
        print len(y_train)
        print len(y_test)
        if len(set(list(y_test)))<2 :
            X_val, X_test_f, y_val, y_test_f = train_test_split(X_test, y_test, test_size=.5, random_state=rand_state2.next(),shuffle=True)
        else:
            X_val, X_test_f, y_val, y_test_f = train_test_split(X_test, y_test, test_size=.5, random_state=rand_state2.next(), stratify =y_test,shuffle=True)
        train_sizes.append(len(y_train))
        X_train_list.append(X_train) 
        X_test_list.append(X_test_f) 
        X_val_list.append(X_val)
        y_val_list.append(y_val)
        y_train_list.append(y_train) 
        y_test_list.append(y_test_f)
    train_validate_sets=[]
    for j in xrange(len(source_list)):
        validate = [list(xrange(3*i, 3*i+ val_size[j])) for i in xrange(int(train_sizes[j]/3)-val_size[j])]
        print len(validate)
        train_final = []
        for v in validate:
            train_final.append([i for i in xrange(train_sizes[j]) if i not in v])
        train_validate_sets.append((train_final,validate))
        X_train_list[j].to_csv("{}_X_train.csv".format(source_list[j]))
        X_test_list[j].to_csv("{}_X_test_f.csv".format(source_list[j]))
        X_val_list[j].to_csv("{}_X_val.csv".format(source_list[j]))
        y_val_list[j].to_csv("{}_y_val.csv".format(source_list[j]))
        y_train_list[j].to_csv("{}_y_train.csv".format(source_list[j]))
        y_test_list[j].to_csv("{}_y_test_f.csv".format(source_list[j]))
    return X_train_list,  X_val_list, y_train_list,  y_val_list, train_validate_sets

def split_shuffle_sets2(source_list,test_frac=0.464,sublist=None):
    rand_state=stream_random_states()
    rand_state2=stream_random_states()
    #get final columns labels
    df_y_0, df_x_0 = concat_sets(source_list)
    total_set_size = len(df_y_0)
    print total_set_size
    if sublist:
        features = sublist
    else:
        features = df_x_0.keys()
    #make frame lists
    fy_list =[]
    fx_list = []
    val_size = []
    for j in xrange(len(source_list)):
        
        repeat_frame=pd.read_csv(source_list[j]+"_Analysis_Repeat.csv")
        hairpin_frame=pd.read_csv(source_list[j]+"_Analysis_Hairpin.csv")
        nucleotide_frame=pd.read_csv(source_list[j]+"_Analysis_Nucleotide.csv")
        df_y = repeat_frame.loc[:,'Status']
        df_r_x = repeat_frame[[i for i in list(repeat_frame.columns) if i != 'Status']].iloc[:,1:]
        #df_r_x = df_r_x.loc[:, (df_r_x != 0).any(axis=0)]
        df_h_x = hairpin_frame[[i for i in list(hairpin_frame.columns) if i != 'Status']].iloc[:,1:]
        #df_h_x = df_h_x.loc[:, (df_h_x != 0).any(axis=0)]
        df_n_x = nucleotide_frame[[ i for i in list(nucleotide_frame.columns) if i != 'Status']].iloc[:,1:]
        #df_n_x = df_n_x.loc[:, (df_n_x != 0).any(axis=0)]
        df_all_x = pd.concat([df_r_x,df_h_x,df_n_x],axis=1)
        df_all_x = df_all_x[[i for i in features]]
        #df_all_x['Random']=np.random.random(size=len(df_all_x))
        fy_list.append(df_y)
        fx_list.append(df_all_x)
        
    X_train_list = []
    X_test_list = []
    X_val_list = []
    y_val_list = []
    y_train_list = []
    y_test_list = []
    train_sizes = []
    for j in xrange(len(source_list)):
        if len(set(list(fy_list[j])))<2 :
            X_train, X_test, y_train, y_test = train_test_split(fx_list[j], fy_list[j], test_size=test_frac, random_state=rand_state.next(),shuffle=True)
        else:
            X_train, X_test, y_train, y_test = train_test_split(fx_list[j], fy_list[j], test_size=test_frac, random_state=rand_state.next(), stratify = fy_list[j],shuffle=True)
        print len(y_train)
        print len(y_test)
        if len(set(list(y_test)))<2 :
            X_val, X_test_f, y_val, y_test_f = train_test_split(X_test, y_test, test_size=.5, random_state=rand_state2.next(),shuffle=True)
        else:
            X_val, X_test_f, y_val, y_test_f = train_test_split(X_test, y_test, test_size=.5, random_state=rand_state2.next(), stratify =y_test,shuffle=True)
        train_sizes.append(len(y_train))
        X_train_list.append(X_train) 
        X_test_list.append(X_test_f) 
        X_val_list.append(X_val)
        y_val_list.append(y_val)
        y_train_list.append(y_train) 
        y_test_list.append(y_test_f)
        X_test.to_csv("{}_X_test.csv".format(source_list[j]))
        y_test.to_csv("{}_y_test.csv".format(source_list[j]))
        X_train_list[j].to_csv("{}_X_train.csv".format(source_list[j]))
        X_test_list[j].to_csv("{}_X_test_f.csv".format(source_list[j]))
        X_val_list[j].to_csv("{}_X_val.csv".format(source_list[j]))
        y_val_list[j].to_csv("{}_y_val.csv".format(source_list[j]))
        y_train_list[j].to_csv("{}_y_train.csv".format(source_list[j]))
        y_test_list[j].to_csv("{}_y_test_f.csv".format(source_list[j]))
    return X_train_list,  X_val_list, y_train_list,  y_val_list

def restore_train_test(source_list):

    #regenerating internal train val tests for local training 
    X_train_list = []
    X_test_list = []
    y_train_list = []
    y_test_list = []
    train_sizes = []
    val_size=[21, 26, 28]
    for j in xrange(len(source_list)):
        X_train =pd.read_csv(source_list[j]+"_X_train.csv",index_col=[0])
        X_test =pd.read_csv(source_list[j]+"_X_test.csv",index_col=[0])
        y_train =pd.read_csv(source_list[j]+"_y_train.csv",header=None,squeeze=True,index_col=[0])
        y_test =pd.read_csv(source_list[j]+"_y_test.csv",header=None,squeeze=True,index_col=[0])
        print y_train
        print len(y_train)
        print len(y_test)
        train_sizes.append(len(y_train))
        X_train_list.append(X_train) 
        X_test_list.append(X_test) 
        y_train_list.append(y_train) 
        y_test_list.append(y_test)
    train_validate_sets=[]
    for j in xrange(len(source_list)):
        validate = [list(xrange(3*i, 3*i+ val_size[j])) for i in xrange(int(train_sizes[j]/3)-val_size[j])]
        print len(validate)
        train_final = []
        for v in validate:
            train_final.append([i for i in xrange(train_sizes[j]) if i not in v])
        train_validate_sets.append((train_final,validate))
    return X_train_list,  X_test_list, y_train_list,  y_test_list, train_validate_sets

def quantileNormalize(df_input):
    df = df_input.copy()
    #compute rank
    dic = {}
    for col in df:
        dic.update({col : sorted(df[col])})
    sorted_df = pd.DataFrame(dic)
    rank = sorted_df.mean(axis = 1).tolist()
    #sort
    for col in df:
        t = np.searchsorted(np.sort(df[col]), df[col])
        df[col] = [rank[i] for i in t]
    return df

def load_to_data_frames_alt(source,model):

    rand_state = stream_random_states()

    #product_frame= pd.read_csv(source+"_Analysis.csv")
    repeat_frame=pd.read_csv(source+"_Analysis_Repeat.csv")
    hairpin_frame=pd.read_csv(source+"_Analysis_Hairpin.csv")
    nucleotide_frame=pd.read_csv(source+"_Analysis_Nucleotide.csv")

    df_y = repeat_frame.loc[:,'Status']
    df_r_x = repeat_frame[[i for i in list(repeat_frame.columns) if i != 'Status']].iloc[:,1:]
    df_r_x = df_r_x.loc[:, (df_r_x != 0).any(axis=0)]
    df_h_x = hairpin_frame[[i for i in list(hairpin_frame.columns) if i != 'Status']].iloc[:,1:]
    df_h_x = df_h_x.loc[:, (df_h_x != 0).any(axis=0)]
    df_n_x = nucleotide_frame[[ i for i in list(nucleotide_frame.columns) if i != 'Status']].iloc[:,1:]
    df_n_x = df_n_x.loc[:, (df_n_x != 0).any(axis=0)]
    df_all_x = pd.concat([df_r_x,df_h_x,df_n_x],axis=1)
    df_all_x = df_all_x.loc[:, (df_all_x != 0).any(axis=0)]

    class_index = {}
    class_vals  = np.array(df_y)

    X_train_full = None
    y_train_full = None
    X_test_full = None
    y_test_full = None
    for i in set(class_vals):
        class_index[i] = [x for x in np.where(class_vals == i)[0]]
        #print i, class_index[i]
        y = df_y[[j for j in class_index[i]]]
        #print y
        if model == "repeat":
            x = df_r_x.iloc[class_index[i]]
        elif model == "hairpin":
            x = df_h_x.iloc[class_index[i]]
        elif model == "nucleotide":
            x = df_n_x.iloc[class_index[i]]
        else:
            x = df_all_x.iloc[class_index[i]]
        #print x
        X_train, X_test, y_train, y_test = train_test_split(x, y, test_size=0.3, random_state=rand_state.next())
        
        if X_train_full is not None:
            X_train_full = pd.concat([X_train_full,X_train],axis=0)
            y_train_full = pd.concat([y_train_full,y_train],axis=0)
            X_test_full = pd.concat([X_test_full,X_test],axis=0)
            y_test_full = pd.concat([y_test_full,y_test],axis=0)
        else:
            X_train_full = X_train
            y_train_full = y_train
            X_test_full = X_test
            y_test_full = y_test        

    return y_train_full, X_train_full, y_test_full, X_test_full

def convert_to_SVM_inputs(source):
    #product_frame= pd.read_csv(source+"_Analysis.csv")
    repeat_frame=pd.read_csv(source+"_Analysis_Repeat.csv")
    hairpin_frame=pd.read_csv(source+"_Analysis_Hairpin.csv")
    nucleotide_frame=pd.read_csv(source+"_Analysis_Nucleotide.csv")

    df_y = repeat_frame.loc[:,'Status']
    df_r_x = repeat_frame[[i for i in list(repeat_frame.columns) if i != 'Status']].iloc[:,1:]
    df_r_x = df_r_x.loc[:, (df_r_x != 0).any(axis=0)]
    df_h_x = hairpin_frame[[i for i in list(hairpin_frame.columns) if i != 'Status']].iloc[:,1:]
    df_h_x = df_h_x.loc[:, (df_h_x != 0).any(axis=0)]
    df_n_x = nucleotide_frame[[ i for i in list(nucleotide_frame.columns) if i != 'Status']].iloc[:,1:]
    df_n_x = df_n_x.loc[:, (df_n_x != 0).any(axis=0)]
    df_all_x = pd.concat([df_r_x,df_h_x,df_n_x],axis=1)
    df_all_x = df_all_x.loc[:, (df_all_x != 0).any(axis=0)]

    df_r_x_n = df_r_x/df_r_x.max()
    df_h_x_n = df_r_x/df_r_x.max()
    df_n_x_n = df_r_x/df_r_x.max()
    df_all_x_n = df_r_x/df_r_x.max()

    return df_y, df_r_x_n, df_h_x_n, df_n_x_n, df_all_x_n

def stream_random_states():
    random.seed(1000)
    while True:
        yield random.randint(1, 1000000000)



def permutation_importances(rf, X_train, y_train, metric):

    baseline = metric(rf, X_train, y_train)

    imp = []

    for col in X_train.columns:

        save = X_train[col].copy()

        X_train[col] = np.random.permutation(X_train[col])

        m = metric(rf, X_train, y_train)

        X_train[col] = save


        imp.append(baseline - m)

    return np.array(imp)

def oob_classifier_accuracy(rf, X_train, y_train):
    """
    Compute out-of-bag (OOB) accuracy for a scikit-learn random forest
    classifier. We learned the guts of scikit's RF from the BSD licensed
    code:
    https://github.com/scikit-learn/scikit-learn/blob/a24c8b46/sklearn/ensemble/forest.py#L425
    """
    X = X_train.values
    y = y_train.values   
    n_samples = len(X)
    n_classes = len(np.unique(y))
    predictions = np.zeros((n_samples, n_classes))
    for tree in rf.estimators_:
        unsampled_indices = _generate_unsampled_indices(tree.random_state, n_samples)
        tree_preds = tree.predict_proba(X[unsampled_indices, :])
        predictions[unsampled_indices] += tree_preds    
    predicted_class_indexes = np.argmax(predictions, axis=1)
    predicted_classes = [rf.classes_[i] for i in predicted_class_indexes]    
    oob_score = np.mean(y == predicted_classes)
    return oob_score


def train_best_hyper_params(source,destination,seed=None, iterations=100,sublist=None):
    # import dataset
    if isinstance(source,list):
        for s in range(len(source)):
            if s==0:
                X_train = pd.read_csv(source[s]+"_X_train.csv",index_col=[0])
                X_test = pd.read_csv(source[s]+"_X_test.csv",index_col=[0])
                y_train = pd.read_csv(source[s]+"_y_train.csv",header=None,index_col=[0],squeeze=True)
                y_test = pd.read_csv(source[s]+"_y_test.csv",header=None,index_col=[0],squeeze=True)
            else:
                X_train = pd.concat([X_train, pd.read_csv(source[s]+"_X_train.csv",index_col=[0])],axis=0)
                X_test = pd.concat([X_test, pd.read_csv(source[s]+"_X_test.csv",index_col=[0])],axis=0)
                y_train = pd.concat([y_train, pd.read_csv(source[s]+"_y_train.csv",header=None,index_col=[0],squeeze=True)],axis=0)
                y_test = pd.concat([y_test, pd.read_csv(source[s]+"_y_test.csv",header=None,index_col=[0],squeeze=True)],axis=0)
    else:
        X_train =pd.read_csv(source+"_X_train.csv",index_col=[0])
        X_test =pd.read_csv(source+"_X_test.csv",index_col=[0])
        y_train =pd.read_csv(source+"_y_train.csv",header=None,squeeze=True,index_col=[0])
        y_test =pd.read_csv(source+"_y_test.csv",header=None,squeeze=True,index_col=[0])
    df_all_x = pd.concat([X_train,X_test],axis=0)
    df_y = pd.concat([y_train,y_test],axis=0)

    # sublist= ['longest_repeat',
    #             'repeat_20',
    #             'size',
    #             'longest_hairpin',
    #             'richest_hairpin',
    #             'wide_hairpins',
    #             'GC_long_h',
    #             'GC_long_l',
    #             'GC_short_h',
    #             'GC_short_l',
    #             'GC_term_l',
    #             'Tm_high',
    #             'Tm_low',
    #             'dGC',
    #             'dTm',
    #             'i_motifs',
    #             'total_GC']
    if sublist:
        features=sublist
    else:
        features= [f for f in df_all_x.columns]          

    df_all_x = df_all_x[[i for i in features]]
    # df_all_x['Random']=np.random.random(size=len(df_all_x))
    #df_all_x=quantileNormalize(df_all_x)
    name = destination
    features = list(df_all_x.columns)
    # initialize grid/bounds for random CV
    if not seed:
        n_estimators = [100,200,400,600,800,1000,2000,4000,6000]
        #n_estimators.append(5000)
        min_samples_leaf =[2]
        min_samples_leaf.extend([int(x) for x in np.linspace(start=4,stop=40,num=10)])
        min_samples_split = [4]
        min_samples_split.extend([int(x) for x in np.linspace(start=8,stop=80,num=10)])
        criterion = ['gini']
        max_features = ['auto']
        weight = [None,"balanced","balanced_subsample"]

        #round 0, filtering features

        # 100 rounds of 500 iteration random search 10 fold CV
        grid_c={'n_estimators':n_estimators,
                'min_samples_split':min_samples_leaf,
                'min_samples_leaf':min_samples_split,
                'criterion':criterion,
                'max_features': max_features,
                "class_weight":weight}

        # initialize best param dict and best oob variable
        best_oob_score = 0
        all_oobs=[]
        all_params=[]
        best_param_dicts = []

        RF_i = RFC(bootstrap =True, oob_score=True)

        rand_state=stream_random_states()

        export_dict={}

        n_features= len(list(df_all_x))
        #grid_c['max_features'] = ['auto',int(math.sqrt(n_features)/2.),int(math.sqrt(n_features)*2.)]
        for j in range(iterations):
            print "Round", j
            RCV = RandomizedSearchCV(estimator=RF_i, 
                            param_distributions=grid_c, 
                            n_iter=200, 
                            scoring=make_scorer(fmk_score),  
                            n_jobs=-1, iid=False, 
                            refit=True, 
                            cv=10, 
                            random_state=rand_state.next(),verbose=1)
            results=RCV.fit(X=df_all_x,y=df_y)
            RF_i.random_state=rand_state.next()
            all_oobs.append(results.best_estimator_.oob_score_)
            all_params.append(results.best_params_)
            if results.best_estimator_.oob_score_ > best_oob_score:
                print "replaced"
                best_oob_score = results.best_estimator_.oob_score_
                best_param_dicts = [results.best_params_]

            elif results.best_estimator_.oob_score == best_oob_score:
                print "appending"
                best_param_dicts.append(results.best_params_)

            print best_param_dicts
            print best_oob_score

            export_dict['best_oob_score']=best_oob_score
            export_dict['best_params']=best_param_dicts
            export_dict['params']=all_params
            export_dict['oob_score']=all_oobs
            pkl.dump(export_dict,open("{}_random_{}_hyperparams.pkl".format(name,iterations),"w"))

        return best_param_dicts, best_oob_score

    else:
        ne_seed = seed[0]['n_estimators']
        mss_seed = seed[0]['min_samples_split']
        msl_seed = seed[0]['min_samples_leaf']
        w_seed = seed[0]["class_weight"]
        c_seed = seed[0]["criterion"]
        mf_seed = seed[0]["max_features"]

        n_estimators=[int(math.ceil(x)) for x in np.logspace(start=math.log10(0.5*ne_seed),stop=math.log10(2*ne_seed),num=7)]
        if mss_seed>8:
            min_samples_split = [int(math.ceil(x)) for x in np.logspace(start=math.log10(0.5*mss_seed),stop=math.log10(2*mss_seed),num=7)]
        else:
            min_samples_split = [int(math.ceil(x)) for x in np.linspace(start=4,stop=10,num=7)]
        if msl_seed>4:
            min_samples_leaf = [int(math.ceil(x)) for x in np.logspace(start=math.log10(0.5*msl_seed),stop=math.log10(2*msl_seed),num=7)]
        else:
            min_samples_leaf = [int(math.ceil(x)) for x in np.linspace(start=2,stop=8,num=7)]
        class_weight = [w_seed]
        criterion = [c_seed]
        max_features = [mf_seed]

        grid_f = {'n_estimators':n_estimators,
                'min_samples_split':min_samples_split,
                'min_samples_leaf':min_samples_leaf,
                'criterion':criterion,
                'class_weight':class_weight,
                'max_features':max_features}


        RF_i = RFC(bootstrap =True, oob_score=True)

        # initialize best param dict and best oob variable
        best_oob_score = 0
        all_oobs=[]
        all_params=[]
        best_param_dicts = []

        # initialize RF Class
        RF_i = RFC(bootstrap =True, oob_score=True)

        rand_state=stream_random_states()

        export_dict={}

        n_features= len(list(df_all_x))
        #grid_c['max_features'] = ['auto',int(math.sqrt(n_features)/2.),int(math.sqrt(n_features)*2.)]
        for j in range(iterations):
            print "Round", j

            RF_i.random_state=rand_state.next()
            GSCV = GridSearchCV(estimator=RF_i, 
                                param_grid=grid_f, 
                                n_jobs=-1, 
                                cv=10, 
                                scoring=make_scorer(fmk_score),
                                refit=True,
                                verbose=1,iid=False)
            results = GSCV.fit(X=df_all_x,y=df_y)
            all_oobs.append(results.best_estimator_.oob_score_)
            all_params.append(results.best_params_)
            if results.best_estimator_.oob_score_ > best_oob_score:
                print "replaced"
                best_oob_score = results.best_estimator_.oob_score_
                best_param_dicts = [results.best_params_]
            elif results.best_estimator_.oob_score == best_oob_score:
                print "appending"
                best_param_dicts.append(results.best_params_)
        
            print best_param_dicts
            print best_oob_score
            
            export_dict['best_oob_score']=best_oob_score
            export_dict['best_params']=best_param_dicts
            export_dict['params']=all_params
            export_dict['oob_score']=all_oobs
            pkl.dump(export_dict,open("{}_grid_{}_hyperparams.pkl".format(name,iterations),"w"))
        
        return best_param_dicts,best_oob_score

def dropcol_importances(name,rf, X_train, y_train,X_test,y_test,scorer,seed):

    rf_ = clone(rf)
    baseline=[]
    dropped=[]
    for i in range(len(seed)):
        print i
        rf_.random_state = seed[i]

        rf_.fit(X_train, y_train)

        y_pred=rf_.predict(X_test)
    
        temp = scorer(y_test,y_pred)
        baseline.append(temp)
        #baseline=rf_.oob_score_
    
        dropped.append([])
    
        for col in X_train.columns:

            X = X_train.drop(col, axis=1)

            X_t =X_test.drop(col, axis=1)

            rf_ = clone(rf)

            rf_.random_state = seed[i]

            rf_.fit(X, y_train)

            y_pred=rf_.predict(X_t)

            offset = scorer(y_test,y_pred)
            
            dropped[i].append(offset)
            #o =rf_.oob_score_
            #print o
            #imp.append(baseline - o)
    baseline = np.array(baseline)
    dropped = np.array(dropped)
    
    super_structure=[]
    for row in baseline:
        b=row*np.ones((len(seed),1))
        o=np.copy(dropped)
        imp=np.subtract(b,o)
        print len(imp.T)
        super_structure.append(imp.T)
    pkl.dump(super_structure,open(name+"feature_reduction_data.pkl","w"))
    return super_structure
    #return I

def feature_reduction(source,destination,params):
    if isinstance(source,list):
        for s in range(len(source)):
            if s==0:
                X_train = pd.read_csv(source[s]+"_X_train.csv",index_col=[0])
                X_test = pd.read_csv(source[s]+"_X_test.csv",index_col=[0])
                Y_train = pd.read_csv(source[s]+"_y_train.csv",header=None,index_col=[0],squeeze=True)
                Y_test = pd.read_csv(source[s]+"_y_test.csv",header=None,index_col=[0],squeeze=True)
            else:
                X_train = pd.concat([X_train, pd.read_csv(source[s]+"_X_train.csv",index_col=[0])],axis=0)
                X_test = pd.concat([X_test, pd.read_csv(source[s]+"_X_test.csv",index_col=[0])],axis=0)
                Y_train = pd.concat([Y_train, pd.read_csv(source[s]+"_y_train.csv",header=None,index_col=[0],squeeze=True)],axis=0)
                Y_test = pd.concat([Y_test, pd.read_csv(source[s]+"_y_test.csv",header=None,index_col=[0],squeeze=True)],axis=0)
    else:
        X_train = pd.read_csv(source+"_X_train.csv",index_col=[0])
        X_test = pd.read_csv(source+"_X_test.csv",index_col=[0])
        Y_train = pd.read_csv(source+"_y_train.csv",header=None,index_col=[0],squeeze=True)
        Y_test = pd.read_csv(source+"_y_test.csv",header=None,index_col=[0],squeeze=True)
    print len(X_train)
    print len(X_test)
    RF=RFC(n_jobs=-1,bootstrap =True,oob_score=True,criterion=params[0]['criterion'],max_features=params[0]['max_features'],n_estimators=params[0]['n_estimators'],min_samples_split= params[0]['min_samples_split'],  min_samples_leaf=params[0]['min_samples_leaf'])
    rand_state=stream_random_states()
    seed_list = [rand_state.next() for i in range(100)]
    importances=dropcol_importances(destination, RF, X_train, Y_train,X_test,Y_test,fmk_score,seed_list)
    imps=np.array([np.array([importances[j][i]for j in range(len(importances))]).flatten() for i in range(len(importances[0]))])
    frame=pd.DataFrame(imps, index=X_train.columns).T
    f_dict={"means":frame.mean(axis=0),
            "std":frame.std(axis=0)}
    df = pd.DataFrame(f_dict, columns=["means","std"], index=X_train.columns ).sort_values(by=['means'],ascending=True)
    print df
    fig1 = plt.figure(figsize=(10,10))
    #cmap = sns.color_palette("Reds",n_colors=1)+sns.color_palette("Reds",n_colors=8,desat=0.5)+sns.color_palette("Blues",n_colors=7)+sns.color_palette("Greens",n_colors=8)+sns.color_palette("Greens",n_colors=7,desat=0.5)

    plt.barh(y=df.index,width=df['means'],xerr=df['std'],capsize=3)
    #plt.xticks(rotation=90)
    plt.savefig(destination+"feature_reduction.eps",format='eps',dpi=600)
    plt.show()
    print df[df.means>0]
    return df[df.means>0].T.columns

def feature_reduction_data_prep(source_list,subset,test_frac=0.20):
    rand_state=stream_random_states()
    #get final columns labels
    df_y_0, df_x_0 = concat_sets(source_list)
    total_set_size = len(df_y_0)
    #print total_set_size
    features = df_x_0.keys()
    #make frame lists
    fy_list =[]
    fx_list = []
    for j in xrange(len(source_list)):
        repeat_frame=pd.read_csv(source_list[j]+"_Analysis_Repeat.csv")
        hairpin_frame=pd.read_csv(source_list[j]+"_Analysis_Hairpin.csv")
        nucleotide_frame=pd.read_csv(source_list[j]+"_Analysis_Nucleotide.csv")
        df_y = repeat_frame.loc[:,'Status']
        df_r_x = repeat_frame[[i for i in list(repeat_frame.columns) if i != 'Status']].iloc[:,1:]
        #df_r_x = df_r_x.loc[:, (df_r_x != 0).any(axis=0)]
        df_h_x = hairpin_frame[[i for i in list(hairpin_frame.columns) if i != 'Status']].iloc[:,1:]
        #df_h_x = df_h_x.loc[:, (df_h_x != 0).any(axis=0)]
        df_n_x = nucleotide_frame[[ i for i in list(nucleotide_frame.columns) if i != 'Status']].iloc[:,1:]
        #df_n_x = df_n_x.loc[:, (df_n_x != 0).any(axis=0)]
        df_all_x = pd.concat([df_r_x,df_h_x,df_n_x],axis=1)
        df_all_x = df_all_x[[i for i in features]]
        if subset[j]>0:
            if len(set(list(df_y)))<2 :
                X_train, df_all_x, y_train, df_y = train_test_split(df_all_x, df_y, test_size=subset[j], random_state=rand_state.next())
            else:
                X_train, df_all_x, y_train, df_y = train_test_split(df_all_x, df_y, test_size=subset[j], random_state=rand_state.next(), stratify = df_y)
            #fy_list.append(y_test)
            #fx_list.append(X_test)
        # else:
        #     fy_list.append(df_y)
        #     fx_list.append(df_all_x)
        # df_X_full=pd.concat(fx_list,axis=0)
        # df_Y_full=pd.concat(fy_list,axis=0)
       
        if len(set(list(df_y)))<2 :
            X_train, X_test, y_train, y_test = train_test_split(df_all_x, df_y, test_size=test_frac, random_state=rand_state.next())
        else:
            X_train, X_test, y_train, y_test = train_test_split(df_all_x, df_y, test_size=test_frac, random_state=rand_state.next(), stratify = df_y)
            print len(y_train)
            print len(y_test)
        X_train.to_csv(source_list[j]+"_fr_X_train.csv")
        X_test.to_csv(source_list[j]+"_fr_X_test.csv")
        y_train.to_csv(source_list[j]+"_fr_y_train.csv")
        y_test.to_csv(source_list[j]+"_fr_y_test.csv")

def feature_dataset_preparation(source,split=.20):
    if isinstance(source,list):
        try:
            for s in range(len(source)):
                if s==0:
                    X_train = pd.read_csv(source[s]+"_fr_X_train.csv",index_col=[0])
                    X_test = pd.read_csv(source[s]+"_fr_X_test.csv",index_col=[0])
                    Y_train = pd.read_csv(source[s]+"_fr_y_train.csv",header=None,index_col=[0],squeeze=True)
                    Y_test = pd.read_csv(source[s]+"_fr_y_test.csv",header=None,index_col=[0],squeeze=True)
                else:
                    X_train = pd.concat([X_train, pd.read_csv(source[s]+"_fr_X_train.csv",index_col=[0])],axis=0)
                    X_test = pd.concat([X_test, pd.read_csv(source[s]+"_fr_X_test.csv",index_col=[0])],axis=0)
                    Y_train = pd.concat([Y_train, pd.read_csv(source[s]+"_fr_y_train.csv",header=None,index_col=[0],squeeze=True)],axis=0)
                    Y_test = pd.concat([Y_test, pd.read_csv(source[s]+"_fr_y_test.csv",header=None,index_col=[0],squeeze=True)],axis=0)
            return True
        except:
            print "failed to find", source[s]
            try:
                for s in range(len(source)):
                    feature_reduction_data_prep(source,[0 for i in range(len(source))],test_frac=split)
                return True
            except:
                return False

    else:
        try:
            train_X =  pd.read_csv(source+"_fr_X_train.csv",index_col=[0])
            test_X = pd.read_csv(source+"_fr_X_test.csv",index_col=[0])
            train_y = pd.read_csv(source+"_fr_y_train.csv",header=None,index_col=[0],squeeze=True)
            test_y = pd.read_csv(source+"_fr_y_test.csv",header=None,index_col=[0],squeeze=True)
            return True
        except:
            try:
                feature_reduction_data_prep([source],[0],test_frac=split)
                return True
            except:
                return False



def feature_reduction_pipeline(source,destination,split=.20,param_dict=None):
    flag=feature_dataset_preparation(source,split)
    
    if flag:
        if isinstance(source,list):
            source2 = [s+"_fr" for s in source]
        else:
            source2 = source+"_fr"
    else:
        print "One or more datasets in this pipeline is missing feature data, run extract_to_csv_train to enable its use in this pipeline"
    if param_dict==None:
        paramsf,best_oobf=train_best_hyper_params(source2,destination,iterations=50)
        #params2,best_oob2=train_best_hyper_params(source2,destination,seed=params1,iterations=50)
        #paramsf,best_oobf=train_best_hyper_params(source2,destination,seed=params2,iterations=100)
    else:
        paramsf=param_dict
    feature_results=feature_reduction(source2,destination,paramsf)
    return feature_results
   

def forest_training_pipeline(source,destination,split=.464,sublist=None,param_dict=None):
    if isinstance(source,list):
        X_train_list,  X_val_list, y_train_list,  y_val_list = split_shuffle_sets2(source,test_frac=split,sublist=sublist)
    else:
        X_train_list,  X_val_list, y_train_list,  y_val_list = split_shuffle_sets2([source],test_frac=split,sublist=sublist)
    #except:
    #    print "One or more datasets in this pipeline is missing feature data, run extract_to_csv_train to enable its use in this pipeline"
    if param_dict:
        paramsf=param_dict
    else:
        params1,best_oob1=train_best_hyper_params(source,destination,iterations=100)
        params2,best_oob2=train_best_hyper_params(source,destination,seed=params1,iterations=50)
        paramsf,best_oobf=train_best_hyper_params(source,destination,seed=params2,iterations=100)
    
    if isinstance(paramsf,list):
        results_dict=train_best_forest_true(destination,X_train_list,y_train_list,X_val_list,y_val_list,paramsf[0])
    else:
        results_dict=train_best_forest_true(destination,X_train_list,y_train_list,X_val_list,y_val_list,paramsf)
    return results_dict
   

def forest_analysis_pipeline(source_list,path,destination,paramsf):
    #sorting trained forests
    print "loading_results", destination
    results_dict_2 = pkl.load(open(path+destination+"_parallel.pkl",'r'))
    frame_alt = pd.DataFrame(results_dict_2['validation_scores'],columns=results_dict_2['validation_scores'].keys())
    if "seeds" not in frame_alt.keys():
        print "missing seeds, repopulating"
        random_state=stream_random_states()
        frame_alt['seeds'] = [random_state.next() for i in range(len(frame_alt))]
    frame_alt = frame_alt.sort_values(by=['v_f1'])
    seeds = list(frame_alt.seeds[-100:])
    seed= seeds[-1]
    print "calculating learning iterations"
    best_f1=[]
    best_mcc=[]
    best_kappa=[]
    best_seed=[]
    best_indices=[]
    best_oob=[]
    current_oob=0
    current_f1=0
    current_mcc=0
    current_kappa=0
    current_seed=None
    for i in range(len(frame_alt)):
        if frame_alt.v_f1[i]>=current_f1 and frame_alt.v_mcc[i]>=current_mcc and frame_alt.v_kappa[i]>=current_kappa:
            #print frame_alt.v_f1[i]
            current_f1=frame_alt.v_f1[i]
            current_mcc=frame_alt.v_mcc[i]
            current_kappa=frame_alt.v_kappa[i]
            current_seed=frame_alt.seeds[i]
            #current_index=frame_alt.indices[i]
            current_oob=frame_alt.v_oob[i]
        best_f1.append(current_f1)
        best_mcc.append(current_mcc)
        best_kappa.append(current_kappa)
        #best_seed.append(current_seed)
        # best_indices.append(current_index)
        best_oob.append(current_oob)
    
    #examining forests training results
    fig=plt.figure(figsize=(9,6))
    plot=plt.plot(#xrange(len(frame_alt)),frame_alt['t_accuracy'],
                  #xrange(len(frame_alt)),frame_alt[['t_fmk']],
                  xrange(len(frame_alt)),best_f1,
                  xrange(len(frame_alt)),best_mcc,
                  xrange(len(frame_alt)),best_kappa,
                  xrange(len(frame_alt)),best_oob)
    plt.xlabel("Iterations")
    plt.ylabel("Score")
    plt.legend(["F1","MCC","CK"],bbox_to_anchor=(1, 1))
    #plt.savefig("Training_Results_gini_reduced.eps",format='eps',dpi=600)
    plt.savefig(path+"Training_Results_{}.eps".format(destination),format='eps',dpi=600)

    print "loading training data"
    for j in range(len(source_list)):
        if j ==0:
            X_train = pd.read_csv(path+source_list[j]+"_X_train.csv",index_col=[0])
            y_train = pd.read_csv(path+source_list[j]+"_y_train.csv",index_col=[0], header= None,squeeze=True)
        else:
            X_train = pd.concat([X_train,pd.read_csv(path+source_list[j]+"_X_train.csv",index_col=[0])],axis=0)
            y_train = pd.concat([y_train,pd.read_csv(path+source_list[j]+"_y_train.csv",index_col=[0], header= None,squeeze=True)],axis=0)
    print X_train

    forest= RFC(n_jobs=1,
                bootstrap =True, 
                oob_score=True,
                criterion=paramsf['criterion'],
                max_features=paramsf['max_features'],
                n_estimators=paramsf['n_estimators'],
                min_samples_split= paramsf['min_samples_split'],  
                min_samples_leaf= paramsf['min_samples_leaf'],
                random_state=seed)
    print "training_best_forest"
    best_forest=forest.fit(X_train,y_train)

    features=X_train.columns

    trees=best_forest.estimators_
    forest_imports=best_forest.feature_importances_
    print forest_imports
    tree_weights=[[i if i>0 else None for i in t.feature_importances_ ] for t in trees]
    tree_frame=pd.DataFrame(tree_weights, columns=features)

    tree_imports=[[i for i in t.feature_importances_ ] for t in trees]
    timports_frame=pd.DataFrame(tree_imports, columns=features)
    cmap=sns.color_palette("Reds",n_colors=5)+sns.color_palette("Blues",n_colors=2)+sns.color_palette("Greens",n_colors=8)+sns.color_palette("Reds",n_colors=3,desat=0.5)

    stdev=np.std(tree_imports,axis=0)
    print tree_frame.count(axis=0)
    print "plotting forest importances"


    fig,ax4=plt.subplots(figsize=(10,3))
    plot2=ax4.bar(x=features,height=tree_frame.count(axis=0),color=cmap,align='center')
    plt.xticks(rotation=90)
    plt.ylim([0,best_forest.n_estimators])
    rect_labels = []
    # Lastly, write in the ranking inside each bar to aid in interpretation
    for rect in plot2:
        # Rectangle widths are already integer-valued but are floating
        # type, so it helps to remove the trailing decimal point and 0 by
        # converting width to int type
        width = round(float(rect.get_height()),2)
        rankStr=str(width)
        # The bars aren't wide enough to print the ranking inside
        # Shift the text to the left side of the right edge
        yloc = width+.2
        # White on magenta
        clr = 'black'
        align = 'right'

        # Center the text vertically in the bar
        xloc = rect.get_x() + rect.get_width()/2.0 +0.1
        label = ax4.text(xloc, yloc, rankStr, horizontalalignment='center',
                         verticalalignment='top', color=clr,clip_on=True,rotation=90)
        rect_labels.append(label)
    plt.savefig(path+"{}_participation_bar_plot.eps".format(destination),format='eps',dpi=600)

    fig=plt.figure(figsize=(10,3))

    plt.bar(x=features,height=forest_imports,color=cmap,align='center',yerr=stdev,capsize=3)
    plt.xticks(rotation=90)
    plt.ylim([-.025,0.5])
    plt.savefig(path+"{}_inital_forest_imps_error_bars.eps".format(destination),format='eps',dpi=600)
    #plt.show()
    print "plotting permuation importances"
    imp= permutation_importances(best_forest, X_train, y_train, make_scorer(fmk_score))

    fig=plt.figure(figsize=(10,3))

    plt.bar(x=features,height=imp,color=cmap,align='center')
    plt.xticks(rotation=90)
    plt.ylim([-.025,0.50])
    plt.savefig(path+"{}_initial_forest_imps_permutation.eps",format='eps',dpi=600)
    
    print "generating inital metrics for valiation set"
    for j in range(len(source_list)):
        if j ==0:
            X_retrain = pd.read_csv(path+source_list[j]+"_X_val.csv",index_col=[0])
            y_retrain = pd.read_csv(path+source_list[j]+"_y_val.csv",index_col=[0], header= None,squeeze=True)
        else:
            X_retrain = pd.concat([X_retrain,pd.read_csv(path+source_list[j]+"_X_val.csv",index_col=[0])],axis=0)
            y_retrain = pd.concat([y_retrain,pd.read_csv(path+source_list[j]+"_y_val.csv",index_col=[0], header= None,squeeze=True)],axis=0)

    print len(X_retrain)
    print len(y_retrain)
    y_pred=best_forest.predict(X_retrain)

    f1 = f1_score(y_retrain,y_pred,average="micro")
    print "f1",f1
    ac = accuracy_score(y_retrain,y_pred)
    print "ac",ac
    pr = precision_score(y_retrain,y_pred,average="micro")
    print "pr",pr
    re = recall_score(y_retrain,y_pred,average="micro")
    print "re",re
    mcc = matthews_corrcoef(y_retrain,y_pred)
    print "mcc",mcc
    ck = cohen_kappa_score(y_retrain,y_pred)
    print "ck",ck
    #print best_forest.predict_proba(X_test)[:,1]
    #print list(Y_test[1])
    roc = roc_auc_score([1 if i =="Synthesized" else 0 for i in list(y_retrain)],best_forest.predict_proba(X_retrain)[:,1],average="micro")
    print "auc_roc",roc
    cf=confusion_matrix(y_retrain,y_pred)
    print cf
    print 

    export_dict={"forest":best_forest,
                 "features":features,
                 "seed":seed,
                 "params":paramsf,
                 "importances":{features[i]:forest_imports[i] for i in range(len(features))},
                 "metrics":{"f1":f1,
                            "ac":ac,
                            "pr":pr,
                            "re":re,
                            "mcc":mcc,
                            "ck":ck,
                            "roc":roc},
                 "confusion":cf}

    pkl.dump(export_dict,open(path+"{}_inital_forest_results.pkl".format(destination),"w"))
    #regenerating internal train val tests for local training 



def retrain_functions(source_list,path,destination,oversample=[0.6,0.25]):
    parent_dict=pkl.load(open(path+"{}_inital_forest_results.pkl".format(destination),"r"))
    seed=parent_dict["seed"]
    forest=parent_dict["forest"]
    features=parent_dict["features"]
    paramsf=parent_dict["params"]

    for j in xrange(len(source_list)):
        if j == 0:
            X_train = pd.read_csv("{}_X_train.csv".format(path+source_list[j]),index_col=[0])
            X_val = pd.read_csv("{}_X_val.csv".format(path+source_list[j]),index_col=[0])
            Y_train = pd.read_csv("{}_y_train.csv".format(path+source_list[j]),header=None,index_col=[0],squeeze=True)
            Y_val = pd.read_csv("{}_y_val.csv".format(path+source_list[j]),header=None,index_col=[0],squeeze=True)
        else:
            X_train = pd.concat([X_train,pd.read_csv("{}_X_train.csv".format(path+source_list[j]),index_col=[0])],axis=0)
            X_val = pd.concat([X_val,pd.read_csv("{}_X_val.csv".format(path+source_list[j]),index_col=[0])],axis=0)
            Y_train = pd.concat([Y_train,pd.read_csv("{}_y_train.csv".format(path+source_list[j]),index_col=[0],header=None,squeeze=True)],axis=0)
            Y_val = pd.concat([Y_val,pd.read_csv("{}_y_val.csv".format(path+source_list[j]),index_col=[0],header=None,squeeze=True)],axis=0)

    forest= RFC(n_jobs=1,
                bootstrap =True, 
                oob_score=True,
                criterion=paramsf['criterion'],
                max_features=paramsf['max_features'],
                n_estimators=paramsf['n_estimators'],
                min_samples_split= paramsf['min_samples_split'],  
                min_samples_leaf= paramsf['min_samples_leaf'],
                random_state=seed)
    forest.fit(X_train,Y_train)

    print "Train"
    y_pred = forest.predict(X_train)
    f1 = f1_score(Y_train,y_pred,average="micro")
    print "f1",f1
    ac = accuracy_score(Y_train,y_pred)
    print "ac",ac
    pr = precision_score(Y_train,y_pred,average="micro")
    print "pr",pr
    re = recall_score(Y_train,y_pred,average="micro")
    print "re",re
    mcc = matthews_corrcoef(Y_train,y_pred)
    print "mcc",mcc
    ck = cohen_kappa_score(Y_train,y_pred)
    print "ck",ck
    roc = roc_auc_score([1 if i =="Synthesized" else 0 for i in list(Y_train)],forest.predict_proba(X_train)[:,1],average="micro")
    print "auc_roc",roc

    print "Validate"

    y_pred = forest.predict(X_val)
    f1 = f1_score(Y_val,y_pred,average="micro")
    print "f1",f1
    ac = accuracy_score(Y_val,y_pred)
    print "ac",ac
    pr = precision_score(Y_val,y_pred,average="micro")
    print "pr",pr
    re = recall_score(Y_val,y_pred,average="micro")
    print "re",re
    mcc = matthews_corrcoef(Y_val,y_pred)
    print "mcc",mcc
    ck = cohen_kappa_score(Y_val,y_pred)
    print "ck",ck
    roc = roc_auc_score([1 if i =="Synthesized" else 0 for i in list(Y_val)],forest.predict_proba(X_val)[:,1],average="micro")
    print "auc_roc",roc

    cf=confusion_matrix(Y_val,y_pred)
    print cf

    confusion={}
    FN={}
    FP={}
    for j in range(len(source_list)):
        confusion[path+source_list[j]]={}
        FN[path+source_list[j]]={}
        FP[path+source_list[j]]={}
        X_train = pd.read_csv("{}_X_train.csv".format(path+source_list[j]),index_col=[0])
        X_val = pd.read_csv("{}_X_val.csv".format(path+source_list[j]),index_col=[0])
        Y_train = pd.read_csv("{}_y_train.csv".format(path+source_list[j]),header=None,index_col=[0],squeeze=True)
        Y_val = pd.read_csv("{}_y_val.csv".format(path+source_list[j]),index_col=[0],header=None,squeeze=True)

        yprtr=forest.predict(X_train)
        yprts=forest.predict(X_val)

        ypbtr=forest.predict_proba(X_train)[:,1]
        #ypbv=forest.predict_proba(X_val)[:,1]
        ypbts=forest.predict_proba(X_val)[:,1]
        c_matrix=[[0,0],[0,0]]
        #print X_val
        #print Y_val.iloc[0]
        #print len(Y_val)
        FN[path+source_list[j]]["val"]=[]
        FP[path+source_list[j]]["val"]=[]
        for i in range(len(Y_val)):
            #print Y_val.index[i]
            if Y_val.iloc[i] == "Synthesized":
                if yprts[i]==Y_val.iloc[i]:
                    c_matrix[0][0]+=1
                else:
                    #print "FN",Y_val.index[i], ypbts[i]
                    #print Y_val.index[i]
                    #print
                    FN[path+source_list[j]]["val"].append((Y_val.index[i],ypbts[i]))
                    c_matrix[0][1]+=1
            if Y_val.iloc[i] == "Canceled":
                if yprts[i]==Y_val.iloc[i]:
                    c_matrix[1][1]+=1
                else:
                    #print "FP", Y_val.index[i], ypbts[i]
                    #print Y_val.index[i]
                    #print
                    FP[path+source_list[j]]["val"].append((Y_val.index[i],ypbts[i]))
                    c_matrix[1][0]+=1
        #print confusion_matrix
        confusion[path+source_list[j]]["val"]=c_matrix
        c_matrix=[[0,0],[0,0]]
        #print "Y_val"
        FN[path+source_list[j]]["train"]=[]
        FP[path+source_list[j]]["train"]=[]
        for i in range(len(Y_train)):
            if Y_train.iloc[i] == "Synthesized":
                if yprtr[i]==Y_train.iloc[i]:
                    c_matrix[0][0]+=1
                else:
                    #print "FN",Y_val.index[i], ypbv[i]
                    #print Y_train.index[i]
                    FN[path+source_list[j]]["train"].append((Y_train.index[i],ypbtr[i]))
                    c_matrix[0][1]+=1
            if Y_train.iloc[i] == "Canceled":
                if yprtr[i]==Y_train.iloc[i]:
                    c_matrix[1][1]+=1
                else:
                    #print "FP", Y_val.index[i], ypbv[i]
                    #print Y_train.index[i]
                    FP[path+source_list[j]]["train"].append((Y_train.index[i],ypbtr[i]))
                    c_matrix[1][0]+=1
        #print confusion_matrix
        confusion[path+source_list[j]]["train"]=c_matrix
        
    
    # optimization of oversampling

    print confusion
    s_list=[path+s for s in source_list]
    for s in s_list:
        X_train = pd.read_csv("{}_X_train.csv".format(s),index_col=[0])
        #X_val = pd.read_csv("{}_X_val_0.csv".format(s),index_col=[0])
        X_val = pd.read_csv("{}_X_val.csv".format(s),index_col=[0])
        Y_train = pd.read_csv("{}_y_train.csv".format(s),header=None,index_col=[0],squeeze=True)
        #Y_val = pd.read_csv("{}_y_val_0.csv".format(s),header=None,index_col=[0],squeeze=True)
        Y_val = pd.read_csv("{}_y_val.csv".format(s),header=None,index_col=[0],squeeze=True)
        
        X_problems=[]
        Y_problems=[]
        X_retest=[]
        Y_retest=[]
        print len(X_train)
        print len(X_val)
        for tup in FP[s]["val"]:
            if tup[1]>oversample[0]:
                #print X_val.loc[tup[0],:]
                #print X_val
                rowx = pd.DataFrame(X_val.loc[tup[0],:]).T
                #rowx.set_index(pd.Index([tup[0]+2000]))
                rowy = Y_val.loc[tup[0]]
                print "val",tup

                X_problems.append(rowx)
                Y_problems.append(rowy)

        for tup in FP[s]["train"]:
            if tup[1]>oversample[0]:
                #print X_val.loc[tup[0],:]
                rowx = pd.DataFrame(X_train.loc[tup[0],:]).T
                rowy = Y_train.loc[tup[0]]
                X_problems.append(rowx)
                Y_problems.append(rowy)
                X_retest.append(rowx)
                Y_retest.append(rowy)
                print "train",tup

        
        for tup in FN[s]["val"]:
            if tup[1]<oversample[1]:
                #print X_val.loc[tup[0],:]
                rowx = pd.DataFrame(X_val.loc[tup[0],:]).T
                rowy = Y_val.loc[tup[0]]
                #X_train=pd.concat([X_retrain,rowx],axis=0)
                #Y_retrain[tup[0]]=rowy
                X_problems.append(rowx)
                Y_problems.append(rowy)
                print "val",tup
                #X_val=X_val.drop(tup[0])
                #Y_val=Y_val.drop(tup[0])
        for tup in FN[s]["train"]:
            if tup[1]<oversample[1]:
                #print X_val.loc[tup[0],:]
                rowx = pd.DataFrame(X_train.loc[tup[0],:]).T
                rowy = Y_train.loc[tup[0]]
                #X_retrain=pd.concat([X_retrain,rowx],axis=0)
                #Y_retrain[tup[0]]=rowy
                X_problems.append(rowx)
                Y_problems.append(rowy)
                X_retest.append(rowx)
                Y_retest.append(rowy)
                print "train",tup
                #X_val=X_val.drop(tup[0])
                #Y_val=Y_val.drop(tup[0])

        if len(X_problems)>0:
            X_retrain=pd.concat(X_problems,axis=0)
            Y_retrain=pd.Series(Y_problems,index=X_retrain.index)
            X_retrain.to_csv("{}_X_retrain_test.csv".format(s))
            Y_retrain.to_csv("{}_y_retrain_test.csv".format(s))

    # intial forest on test set
    for j in xrange(len(source_list)):
        if j == 0:
            X_train_2 = pd.read_csv("{}_X_train.csv".format(path+source_list[j]),index_col=[0])
            X_test_2 = pd.read_csv("{}_X_test_f.csv".format(path+source_list[j]),index_col=[0])
            Y_train_2 = pd.read_csv("{}_y_train.csv".format(path+source_list[j]),header=None,index_col=[0],squeeze=True)
            Y_test_2 = pd.read_csv("{}_y_test_f.csv".format(path+source_list[j]),index_col=[0],header=None,squeeze=True)
            #print X_test_2
            #print len(Y_test_2)
            
        else:
            X_train_2 = pd.concat([X_train_2,pd.read_csv("{}_X_train.csv".format(path+source_list[j]),index_col=[0])],axis=0)
            X_test_2 = pd.concat([X_test_2,pd.read_csv("{}_X_test_f.csv".format(path+source_list[j]),index_col=[0])],axis=0)
            Y_train_2 = pd.concat([Y_train_2,pd.read_csv("{}_y_train.csv".format(path+source_list[j]),index_col=[0],header=None,squeeze=True)],axis=0)
            Y_test_2 = pd.concat([Y_test_2,pd.read_csv("{}_y_test_f.csv".format(path+source_list[j]),index_col=[0],header=None,squeeze=True)],axis=0)
    print "inital"

    y_pred = forest.predict(X_test_2)

    f11 = f1_score(Y_test_2,y_pred,average="micro")
    print "f1",f11
    ac1 = accuracy_score(Y_test_2,y_pred)
    print "ac",ac1
    pr1 = precision_score(Y_test_2,y_pred,average="micro")
    print "pr",pr1
    re1 = recall_score(Y_test_2,y_pred,average="micro")
    print "re",re1
    mcc1 = matthews_corrcoef(Y_test_2,y_pred)
    print "mcc",mcc1
    ck1 = cohen_kappa_score(Y_test_2,y_pred)
    print "ck",ck1
    roc1 = roc_auc_score([1 if i =="Synthesized" else 0 for i in list(Y_test_2)],forest.predict_proba(X_test_2)[:,1],average="micro")
    print "auc_roc",roc1
    cf1=confusion_matrix(Y_test_2,y_pred)
    print cf1
    print len(Y_test_2)


    print "oversample"

    
    for j in xrange(len(source_list)):
        if j == 0:
            X_train_2 = pd.read_csv("{}_X_train.csv".format(path+source_list[j]),index_col=[0])
            X_test_2 = pd.read_csv("{}_X_test_f.csv".format(path+source_list[j]),index_col=[0])
            Y_train_2 = pd.read_csv("{}_y_train.csv".format(path+source_list[j]),header=None,index_col=[0],squeeze=True)
            Y_test_2 = pd.read_csv("{}_y_test_f.csv".format(path+source_list[j]),header=None,index_col=[0],squeeze=True)
            #print Y_train
        else:
            X_train_2 = pd.concat([X_train_2,pd.read_csv("{}_X_train.csv".format(path+source_list[j]),index_col=[0])],axis=0)
            X_test_2 = pd.concat([X_test_2,pd.read_csv("{}_X_test_f.csv".format(path+source_list[j]),index_col=[0])],axis=0)
            Y_train_2 = pd.concat([Y_train_2,pd.read_csv("{}_y_train.csv".format(path+source_list[j]),index_col=[0],header=None,squeeze=True)],axis=0)
            Y_test_2 = pd.concat([Y_test_2,pd.read_csv("{}_y_test_f.csv".format(path+source_list[j]),header=None,index_col=[0],squeeze=True)],axis=0)
            #print Y_train
    #for j in ["./training_datasets/reduced_gini/bulk_training_set","./training_datasets/reduced_gini/designed_failures"]:
    for j in s_list:
        if os.path.exists(j+"_X_retrain_test.csv"):
            X_train_2 = pd.concat([X_train_2,pd.read_csv("{}_X_retrain_test.csv".format(j),index_col=[0])],axis=0)
            Y_train_2 = pd.concat([Y_train_2,pd.read_csv("{}_y_retrain_test.csv".format(j),header=None,index_col=[0],squeeze=True)],axis=0)
        #X_train_2 = pd.concat([X_train_2,pd.read_csv("{}_X_retrain.csv".format(j),index_col=[0])],axis=0)
        #Y_train_2 = pd.concat([Y_train_2,pd.read_csv("{}_y_retrain.csv".format(j),index_col=[0],squeeze=True)],axis=0)
        #X_test_2 = pd.concat([X_test_2,pd.read_csv("{}_X_retrain_test.csv".format(j),index_col=[0])],axis=0)
        #Y_test_2 = pd.concat([Y_test_2,pd.read_csv("{}_y_retrain_test.csv".format(j),index_col=[0],header=None,squeeze=True)],axis=0)

    print len(X_train_2)
    print len(X_test_2)
    print len(X_train_2.columns)
    #new_forest = RFC(n_jobs=-1,bootstrap =True, oob_score=True,criterion='gini',max_features='auto',n_estimators=1000,min_samples_split= 32,  min_samples_leaf= 2,random_state=seed)
    new_forest_3 = RFC(n_jobs=1,
                bootstrap =True, 
                oob_score=True,
                criterion=paramsf['criterion'],
                max_features=paramsf['max_features'],
                n_estimators=paramsf['n_estimators'],
                min_samples_split= paramsf['min_samples_split'],  
                min_samples_leaf= paramsf['min_samples_leaf'],
                random_state=seed)
    new_forest_3.fit(X_train_2,Y_train_2)

    y_pred = new_forest_3.predict(X_test_2)

    f12 = f1_score(Y_test_2,y_pred,average="micro")
    print "f1",f12
    ac2 = accuracy_score(Y_test_2,y_pred)
    print "ac",ac2
    pr2 = precision_score(Y_test_2,y_pred,average="micro")
    print "pr",pr2
    re2 = recall_score(Y_test_2,y_pred,average="micro")
    print "re",re2
    mcc2 = matthews_corrcoef(Y_test_2,y_pred)
    print "mcc",mcc2
    ck2 = cohen_kappa_score(Y_test_2,y_pred)
    print "ck",ck2
    roc2 = roc_auc_score([1 if i =="Synthesized" else 0 for i in list(Y_test_2)],new_forest_3.predict_proba(X_test_2)[:,1],average="micro")
    print "auc_roc",roc2
    cf2=confusion_matrix(Y_test_2,y_pred)
    print cf2
    print len(Y_test_2)

    if f12>=f11:
        best_forest=forest
    else:
        best_forest=new_forest_3
    # stats for best forest
    trees=best_forest.estimators_
    forest_imports=best_forest.feature_importances_
    tree_weights=[[i if i>0 else None for i in t.feature_importances_ ] for t in trees]
    tree_frame=pd.DataFrame(tree_weights, columns=features)

    tree_imports=[[i for i in t.feature_importances_ ] for t in trees]
    timports_frame=pd.DataFrame(tree_imports, columns=features)

    stdev=np.std(tree_imports,axis=0)
    print tree_frame.count(axis=0)

    fig=plt.figure(figsize=(10,10))
    #cmap=sns.color_palette("Reds",n_colors=2)+sns.color_palette("Blues",n_colors=4)+sns.color_palette("Greens",n_colors=5)+sns.color_palette("Greens",n_colors=6,desat=0.5)
    #cmap=sns.color_palette("Reds",n_colors=8)+sns.color_palette("Reds",n_colors=8,desat=0.5)+sns.color_palette("Blues",n_colors=7)+sns.color_palette("Greens",n_colors=8)+sns.color_palette("Greens",n_colors=7,desat=0.5)
    cmap=sns.color_palette("Reds",n_colors=5)+sns.color_palette("Blues",n_colors=2)+sns.color_palette("Greens",n_colors=8)+sns.color_palette("Reds",n_colors=3,desat=0.5)

    plot=sns.boxplot(data=tree_frame,palette=cmap,orient="h")
    #plt.xlim([0,38])
    plt.xlim([0,0.75])
    #plt.xscale('log')
    plt.ylabel('Features')
    plt.xlabel('Feature Importance')
    #plt.legend(features)
    plt.yticks(rotation=0)
    #plt.savefig("trees_box_plot_LINEAR.eps",format='eps',dpi=600)
    #plt.show()
    fig,ax4=plt.subplots(figsize=(10,3))
    plot2=ax4.bar(x=features,height=tree_frame.count(axis=0),color=cmap,align='center')
    plt.xticks(rotation=90)
    plt.ylim([0,best_forest.n_estimators])
    rect_labels = []
    # Lastly, write in the ranking inside each bar to aid in interpretation
    for rect in plot2:
        # Rectangle widths are already integer-valued but are floating
        # type, so it helps to remove the trailing decimal point and 0 by
        # converting width to int type
        width = round(float(rect.get_height()),2)
        rankStr=str(width)
        # The bars aren't wide enough to print the ranking inside
        # Shift the text to the left side of the right edge
        yloc = width+.2
        # White on magenta
        clr = 'black'
        align = 'right'

        # Center the text vertically in the bar
        xloc = rect.get_x() + rect.get_width()/2.0 +0.1
        label = ax4.text(xloc, yloc, rankStr, horizontalalignment='center',
                         verticalalignment='top', color=clr,clip_on=True,rotation=90)
        rect_labels.append(label)
    plt.savefig(path+"{}_final_forest_bar_plot.eps".format(destination),format='eps',dpi=600)

    

    fig=plt.figure(figsize=(10,3))

    plt.bar(x=features,height=forest_imports,color=cmap,align='center', yerr=stdev,capsize=3)
    plt.xticks(rotation=90)
    plt.ylim([-.025,0.5])
    plt.savefig(path+"{}_final_forest_gini_importances.eps".format(destination),format='eps',dpi=600)

    imp = permutation_importances(best_forest, X_train_2, Y_train_2, make_scorer(fmk_score))

    fig=plt.figure(figsize=(10,3))

    plt.bar(x=features,height=imp,color=cmap,align='center')
    plt.xticks(rotation=90)
    plt.ylim([-.025,0.50])
    plt.savefig(path+"{}_final_forest_permutation_imps.eps".format(destination),format='eps',dpi=600)


    export_dict={"forest":best_forest,
                 "features":features,
                 "seed":seed,
                 "params":paramsf,
                 "importances":{features[i]:forest_imports[i] for i in range(len(features))},
                 "metrics":{"f1":f1,
                            "ac":ac,
                            "pr":pr,
                            "re":re,
                            "mcc":mcc,
                            "ck":ck,
                            "roc":roc},
                 "confusion":cf}

    pkl.dump(export_dict,open(path+"{}_final_forest_results.pkl".format(destination),"w"))

def fmk_score(y_t,y_p):
    #print y_t
    #print y_p
    f1=f1_score(y_t, y_p, average='micro')
    #print "f1",f1,type(f1)
    if len(set(list(y_p)))<2 :
        matt = 0.0
    elif len(set(list(y_t)))<2 :
        matt = 0.0
    else:
        matt = matthews_corrcoef(np.array(y_t),np.array(y_p))
    kappa = cohen_kappa_score(y_t,y_p)
    #print "kappa",kappa,type(kappa)
    return (f1*matt*kappa)**(1./3.)
    '''
    except:
        f1=f1_score(y_test, y_pred, average='micro')
        print f1
        kappa = cohen_kappa_score(y_test,y_pred)
        print kappa
        return (f1*kappa)**(1./2.)
    '''

def train_best_forest_true(name,x_train,y_train,x_val,y_val,param_dict):
    for i in xrange(len(y_train)):
        if i == 0:
            X_train = x_train[i]
            Y_train = y_train[i]
            X_val = x_val[i]
            Y_val = y_val[i]
        else:
            X_train = pd.concat([X_train,x_train[i]],axis=0)
            Y_train = pd.concat([Y_train,y_train[i]],axis=0)
            X_val = pd.concat([X_val,x_val[i]],axis=0)
            Y_val = pd.concat([Y_val,y_val[i]],axis=0)

    export_dict = OrderedDict([('validation_sets',[]),
                               ('seeds',[]),
                               ('feature_importances',[]),
                               ('perm_importances',[]),
                               ('validation_scores',OrderedDict([('v_accuracy',[]),
                                                                 ('v_fmk',[]),
                                                                 ('v_f1',[]),
                                                                 ('v_mcc',[]),
                                                                 ('v_kappa',[]),
                                                                 ('v_oob',[])])),
                               ("best_scores",OrderedDict([('v_accuracy',[]),
                                                                 ('v_fmk',[]),
                                                                 ('v_f1',[]),
                                                                 ('v_mcc',[]),
                                                                 ('v_kappa',[]),
                                                                 ('v_oob',[])]))])
    rand_state=stream_random_states()
    best_fmk_score=0
    best_oob_score=0
    best_accuracy_score=0
    best_estimator = []
    p=Pool(21)
    sub_inputs =[(k,X_train,Y_train,X_val,Y_val, rand_state.next(),param_dict) for k in xrange(40000)]
    #sub_inputs =[(i,j,k,x_train,y_train,X_val,Y_val,,rand_state.next(),param_dict) for k in xrange(len(train_validate_sets[2][0]))]
    output_dict = p.map(training_subroutine,sub_inputs)
    for o in output_dict:
        for key in ['feature_importances','perm_importances','seeds']:
            export_dict[key].extend(o[key])
        for key1 in ["validation_scores"]:
            for key2 in o[key1].keys():
                export_dict[key1][key2].extend(o[key1][key2]) 
    best_f1=[]
    best_mcc=[]
    best_kappa=[]
    best_seed=[]
    current_f1=0
    current_mcc=0
    current_kappa=0
    current_seed=None
    for i in xrange(len(export_dict['seeds'])):
        if export_dict['validation_scores']['v_f1'][i]>=current_f1 and export_dict['validation_scores']['v_mcc'][i]>=current_mcc and export_dict['validation_scores']['v_kappa'][i]>=current_kappa:
            #print frame_alt.v_f1[i]
            current_f1=export_dict['validation_scores']['v_f1'][i]
            current_mcc=export_dict['validation_scores']['v_mcc'][i]
            current_kappa=export_dict['validation_scores']['v_kappa'][i]
            current_seed=export_dict['seeds'][i]
    best_f1.append(current_f1)
    best_mcc.append(current_mcc)
    best_kappa.append(current_kappa)
    best_seed.append(current_seed)

                
    pkl.dump(export_dict,open("{}_parallel.pkl".format(name),"w"))
    #return export_dict

def training_subroutine((k,X_train,Y_train, X_val,Y_val, seed,params)):
    local_dict=OrderedDict([("seeds",[]),
                            ('feature_importances',[]),
                               ('perm_importances',[]),
                               ('validation_scores',OrderedDict([('v_accuracy',[]),
                                                                 ('v_fmk',[]),
                                                                 ('v_f1',[]),
                                                                 ('v_mcc',[]),
                                                                 ('v_kappa',[]),
                                                                 ('v_oob',[])]))])
                               


    print k
    forest_cand = None
    oob_cand = 0
    fmk_cand = 0
    #for r in range(4):
    RF_i = RFC(n_jobs=1,
                bootstrap =True, 
                oob_score=True,
                criterion=params['criterion'],
                max_features=params['max_features'],
                n_estimators=params['n_estimators'],
                min_samples_split= params['min_samples_split'],  
                min_samples_leaf= params['min_samples_leaf'],
                random_state=seed)
    #RF_i.random_state = rand_state.next()
    forest_cand = RF_i.fit(X_train,Y_train)
    #y_0 =forest.predict(X_val)
    #fmk  = fmk_score(Y_val,y_0)
    #if forest.oob_score_>oob_cand:# and fmk>fmk_cand:
    #oob_cand = forest.oob_score
    #fmk_cand = fmk
    #forest_cand=forest
    y_pred =forest_cand.predict(X_val)
    local_dict['seeds'].append(seed)
    local_dict['feature_importances'].append(forest_cand.feature_importances_)
    local_dict['perm_importances'].append(permutation_importances(forest_cand,X_train,Y_train,oob_classifier_accuracy))
    model_oob_score = forest_cand.oob_score_
    local_dict['validation_scores']['v_oob'].append(model_oob_score)
    model_accuracy_score = accuracy_score(Y_val,y_pred)
    local_dict['validation_scores']['v_accuracy'].append(model_accuracy_score)
    model_fmk_score = fmk_score(Y_val,y_pred)
    local_dict['validation_scores']['v_fmk'].append(model_fmk_score)
    model_f1_score = f1_score(Y_val,y_pred,average='micro')
    local_dict['validation_scores']['v_f1'].append(model_f1_score)
    model_mcc_score = matthews_corrcoef(Y_val,y_pred)
    local_dict['validation_scores']['v_mcc'].append(model_mcc_score)
    model_kappa_score = cohen_kappa_score(Y_val,y_pred)
    local_dict['validation_scores']['v_kappa'].append(model_kappa_score)
    
    
    return local_dict



if __name__=="__main__":

    #source_list = ["./bulk_training_set_f","./designed_failures_randomized_f","./derived_successes_randomized_f"]
    #feature_reduction_pipeline(source=source_list,destination="feature_reduction_pipeline_8_15_19")
    # sublist= ['longest_repeat',              
    #             'high_local_density_2',
    #             'high_specific_density', 
    #             'repeat_10',              
    #             'tandem',    
    #             'strong_hairpins',          
    #             'palindromes',  
    #             'GC_short_l',                  
    #             'GC_long_l',
    #             'GC_short_h',        
    #             'GC_long_h',
    #             'GC_term_h',
    #             'dGC',  
    #             'Tm_low',   
    #             'dTm',              
    #             'pattern_runs',
    #             'i_motifs',
    #             'size']


    sublist=['longest_repeat', 
                'total_GC',
                'GC_term_l',
                "repeat_9_metric",
                "Tm_low",
                "size",
                "Tm_high",
                "GC_short_l",
                "richest_hairpin"]

    #sublist=['longest_repeat', 
    #            "repeat_9_metric",
    #            "richest_hairpin",
    #            "GC_short_l",
    #            'GC_term_l',
    #            "Tm_low",
    #            "Tm_high",
    #            'total_GC',
    #            "size"]
    #         'repeat_9_metric', 
    #         'high_local_density_2',
    #         'high_total_density',
    #         'repeat_10', 
    #         'repeat_20', 
    #         'repeat_25',
    #         'repeat_40', 
    #         'repeat_large', 
    #         'tandem', 
    #         'wide_hairpins', 
    #         'strong_hairpins', 
    #         'palindromes', 
    #         'terminal_hairpins',
    #         'Tm_high', 
    #         'Tm_low',
    #         #'GC_short_l', 
    #         'GC_long_h',
    #         'GC_long_l',
    #         'GC_term_l',
    #         'total_GC',
    #         'pattern_runs',
    #         'poly_runs',
    #         'g_quad_motifs',
    #         'size']


    #param_dict={'min_samples_leaf': 2, 'n_estimators': 319, 'min_samples_split': 5, 'criterion': 'gini', 'max_features': 'auto', 'class_weight': None}

    #param_dict={'min_samples_leaf': 2, 'n_estimators': 700, 'min_samples_split': 4, 'criterion': 'gini', 'max_features': 'auto', 'class_weight': None}

    #param_dict ={'min_samples_leaf': 3, 'n_estimators': 80, 'min_samples_split': 4, 'criterion': 'gini', 'max_features': 'auto', 'class_weight': None}

    source_list = ["SupplementaryData1"]
#    source_list = ["designed_failures_randomized","derived_successes_randomized","Compiled_IDT_Data_2019","designed_failures","derived_successes"]
    param_dict={'min_samples_leaf': 2, 'n_estimators': 320, 'min_samples_split': 5, 'criterion': 'gini', 'max_features': 'auto', 'class_weight': None}

    forest_training_pipeline(source = source_list, destination="pipeline_full_8_16_19", sublist=sublist)
    forest_analysis_pipeline(source_list = source_list, path='./', destination = "pipeline_full_7_14_19", paramsf = param_dict)
    retrain_functions(source_list=source_list,path='./',destination="pipeline_full_7_14_19",oversample=[0.6,0.25])

    
    # X_train_list,  X_test_list, y_train_list,  y_test_list, train_validate_sets = split_shuffle_sets(source_list,sublist=sublist)
    # X_train_list,  X_test_list, y_train_list,  y_test_list, train_validate_sets = restore_train_test(source_list)
    