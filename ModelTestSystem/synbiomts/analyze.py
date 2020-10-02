"""
Main classes for SynBioMTS module

Copyright 2017 Alexander C. Reis, Howard M. Salis, all rights reserved.

"""

import os
import copy_reg
import types
import cPickle as pickle
from itertools import product
import multiprocessing as mp
import scipy
import numpy as np
import pandas
import shelve
import dbms
import stats

# Using copy_reg to allow the pickle module to instance methods
# Replicates multiprocessing module's ForkingPickler
def _reduce_method(m):
    if m.im_self is None:
        return getattr, (m.im_class, m.im_func.func_name)
    else:
        return getattr, (m.im_self, m.im_func.func_name)
copy_reg.pickle(types.MethodType, _reduce_method)


class ModelTest(object):

    # identifiers is used to uniquely identify sequence entries
    identifiers = ["SEQUENCE","SUBGROUP"]

    def __init__(self,models,dbfilename,filters={},recalc=False,add_data=True,nprocesses=(mp.cpu_count()-1),verbose=False):
        '''Inputs:
        models (interface.Container)  = see interface.Container
        dbfilename (string)           = filename of the geneticsystems database
        filters (dictionary)          = dictionary to filter dataset using dbms
        nprocesses (int)              = number of processes to use with
                                        multiprocessing if 1, ModelTest
                                        does not use multiprocessing
        recalc (bool)                 = boolean to tell the testsystem
                                        to recalcualte model predictions
                                        on existing datasets
        verbose (bool)                = talk to me'''

        if models.__name__ != "Models":
            raise Exception("Not an interface.Container object: {}.".format(models))
        
        assert nprocesses > 0,            "nprocesses should be an int > 0"
        assert isinstance(recalc,bool),   "recalc should be boolean"
        assert isinstance(add_data,bool), "add_data should be boolean"
        assert isinstance(filters,dict),  "filters should be a dictionary"

        self.models      = models
        self.dbfilename  = dbfilename
        self.recalc      = recalc
        self.add_data    = add_data
        self.nprocesses  = nprocesses
        self.verbose     = verbose
        self.filters     = filters
        self.predictions = {}
        self.statistics  = {}

        # import sequences from genetic systems database
        # based on specified dbfilename and filters
        self._update_database()        

    def run(self,calcsFilename=None,statsFilename=None):
        '''Use run to (1) run model calculations and (2) model statistics
        Inputs:
        calcsFilename (str) = If not None, saves model predictions to a shelve persistance object
        statsFilename (str) = If not None, saves model statistics to a shelve persistance object'''

        self.predict(calcsFilename)
        self.calc_stats(statsFilename)

    def predict(self,filename=None):
        '''Use predict to (1) run model calculations without model statistics
        Inputs:
        filename (str) = If not None, saves model predictions to a shelve persistance object'''

        db = self.database

        # Remove sequences that have already been calculated if self.recalc is False,
        # convert pandas dataframe into entries, a list of records (dictionaries)
        # then bundle these entries with the models that are available
        n_entries = []
        if (not self.recalc) and (not filename is None):
            bundles = []
            d = shelve.open(filename)
            for model in self.models.available:
                if model in d.keys():
                    kargs = {i: d[model][i] for i in self.identifiers}
                    db.remove_data(kargs,ordered=True)
                    dict_list = db.get_entries()
                    n_entries.append((model,len(dict_list)))
                    bundles += [(model,entry) for entry in dict_list]
                else:
                    dict_list = db.get_entries()
                    bundles += product([model],dict_list)
                    n_entries.append((model,len(dict_list)))
            d.close()
        else:
            dict_list = db.get_entries()
            bundles = product(self.models.available,dict_list)
            n = len(dict_list)
            n_entries = [(m,n) for m in self.models.available]

        # Call multiprocessing (or MPI) to run model predictions
        if self.nprocesses > 1:
            pool = mp.Pool(processes=self.nprocesses)
            output = pool.map(self._wrap,bundles)
            pool.close()
            pool.join()
        else:
            output = [self._wrap(bundle) for bundle in bundles]

        # Convert model predictions (list of dictionaries) to pandas dataframes
        db.data.reset_index(drop=True, inplace=True)
        total = 0
        for model,n in n_entries:
            modelcalcs = pandas.DataFrame(output[total:total+n])
            if self.add_data:
                dfsave = pandas.concat([db.data, modelcalcs], axis=1)
            else:
                dfsave = pandas.concat([db.data[self.identifiers], modelcalcs], axis=1)
            self.predictions[model] = dfsave
            total += n



    def _wrap(self,bundle):
        ''' _wrap interprets the inputs of the interface.Model and pulls those values
        from the database; this method requires a tuple input for Python's map() function.'''

        (name,entry) = bundle
        # entry = {'ORGANISM':"Escherichia coli",'SEQUENCE':"ACTCGATCTT",...}

        # dev notes
        #('ACTGTAC',) # args
        #{'organism': 'E. coli', 'temp': 37.0} # keywords
        #['sequence', 'organism', 'temp', 'startpos'] # variables

        # Remove args and keywords from variables list
        vrs = self.models[name].variables[len(self.models[name].args):]
        vrs = [k for k in vrs[:] if k not in self.models[name].keywords.keys()]

        # Exception handling when data is not available
        if any(k.upper() not in entry.keys() for k in vrs):
            err = "One of {}'s arguments is not in the database.".format(name)
            print "Model requested arguments: " + str(vrs)
            print "Database available values: " + str(entry.keys())
            raise KeyError(err)
        else:
            kargs = {k: entry[k.upper()] for k in vrs}

        if self.verbose:
            fmt = "model={m:s}, SUBGROUP={sbgrp:s}, seq={seq:s}..."
            print fmt.format(m=name,sbgrp=entry['SUBGROUP'],seq=entry['SEQUENCE'][:50])

        # Run model
        return self.models[name](**kargs)

    def _update_database(self):
        ''' _update_database is run on __init__ and anytime the database, datasets,
        or filters are updated with any of the following methods:
        add_datasets(), remove_datasets(), remove_filters().'''

        try:
            handle = open(self.dbfilename,'r')
            database = pickle.load(handle)
            handle.close()
            assert isinstance(database,pandas.DataFrame)
        except:
            raise Exception("Database filename: {} is not valid.".format(self.dbfilename))

        for i in self.identifiers:
            assert i in database.keys(), "{} isn't a database label.".format(i)

        # dbms.get_indexes checks that filters are labels in database
        # This code block code be removed if dbms.get_indexes gets no values
        if "DATASET" in self.filters:
            listed = database["DATASET"].cat.categories
            unlisted = [x for x in self.filters["DATASET"] if x not in listed]
            if unlisted:
                error = "These datasets are unlisted: " + ", ".join(unlisted)
                raise ValueError(error)

        # Use filters specified by user on database
        # if self.filters:
        #     database.filter(self.filters,False)

        if self.filters:
            for key,valuelist in self.filters.iteritems():
                database = database[database[key].isin(valuelist)].reset_index()

        self.database = dbms.DataBase(database)

    def calc_stats(self,filename=None):
        ''' calc_stats runs stats.linear_complete for models with defined 
        functional forms; calc_stats runs the statistics on each of the subgroups
        as well as the full dataset defined by the database any any provided filters.
        Inputs:
        filename = If not None, saves model statistics to a shelve persistance object'''

        for m in self.models.available:
            df = self.predictions[m]
            allError = np.nan*np.ones(len(df))
            allPredicted = np.nan*np.ones(len(df))
            entries = []
            yScale = self.models[m].yScale
            xScale = self.models[m].xScale

            for subgroup in df["SUBGROUP"].unique():

                # Extract subgroup data & predictions
                indx = np.array(df["SUBGROUP"] == subgroup)
                x = np.array(df[self.models[m].x][indx])
                y = np.array(df[self.models[m].y][indx])
                std = np.array(df[self.models[m].std][indx])

                # filter out no-prediction sequences (nan and inf values)
                setsize = len(x)
                invalid = np.isnan(x) + np.isinf(x) + (x == 0.0)
                count_invalid = np.sum(invalid)
                x_valid = x[~invalid]
                y_valid = y[~invalid]
                std_valid = std[~invalid]

                # Insufficient number of sequences in dataset or
                # model was unable to predict more than two sequences
                if (setsize - count_invalid < 2):
                    data = {}
                    data["valid dataset"] = False
                    yErrorAll = np.nan*np.ones(setsize)
                    yPredicted = np.nan*np.ones(setsize)
                    
                else:
                    # Run statistics and information theory calcs
                    # At least two data points required to run stats
                    data,yError = stats.linear_complete(x_valid,y_valid,std_valid,xScale,yScale,self.models[m].a1)
                    data["valid dataset"] = True

                    # "place" yError into correct size array
                    yErrorAll = np.nan*np.ones(setsize)
                    np.place(yErrorAll,~invalid,yError)

                    # Calculate yPredicted based on linear model fit
                    if   xScale == 'ln':     x = np.log(x)
                    elif xScale == 'log10':  x = np.log10(x)
                    elif xScale == 'linear': pass
                    else: raise Exception('Bad xScale set for {}'.format(m))

                    yPredicted = data['slope']*x + data['intercept']
                    if   yScale == 'ln':     yPredicted = np.exp(yPredicted)
                    elif yScale == 'log10':  yPredicted = np.power(10,yPredicted)
                    elif yScale == 'linear': pass
                    else: raise Exception('Bad yScale set for {}'.format(m))

                    # Add yError and yPredicted to model predictions
                    allError[indx] = yErrorAll
                    allPredicted[indx] = yPredicted
                    
                # Define subgroup, number of non-predicted sequence, and sequence entropy
                data["SUBGROUP"] = subgroup
                data["Count.Invalid"] = count_invalid
                data["Sequence entropy"],_ = stats.sequence_entropy(df["SEQUENCE"][indx],\
                                                          positions=df["STARTPOS"][indx])
                
                if data["valid dataset"]:
                    data["MC"] = data["Sequence entropy"]*(data["N-states"]-1)*data["RIG"]
                else:
                    data["MC"] = 0.0

                # Append data to entries
                entries.append(data)

            # Calculate statistics for model {m} on ALL data
            x = allPredicted
            y = np.array(df[self.models[m].y])
            std = np.array(df[self.models[m].std])

            setsize = len(x)
            invalid = np.isnan(x) + np.isinf(x) + (x == 0.0)
            count_invalid = np.sum(invalid)
            x_valid = x[~invalid]
            y_valid = y[~invalid]
            std_valid = std[~invalid]
                
            if (setsize - count_invalid > 1):
                data,_ = stats.linear_simple(x_valid,y_valid,std_valid,xScale=yScale,yScale=yScale)
                data["valid dataset"] = True
            else:
                data = {}
                data["valid dataset"] = False

            data["SUBGROUP"] = "ALL"
            data["Count.Nan"] = count_invalid
            data["Sequence entropy"],_ = stats.sequence_entropy(df["SEQUENCE"],\
                                                      positions=df["STARTPOS"])
            if data["valid dataset"]:
                data["MC"] = data["Sequence entropy"]*(data["N-states"]-1)*data["RIG"]
                print m, data["MC"]
            else:
                data["MC"] = 0.0
            entries.append(data)

            self.statistics[m] = pandas.DataFrame(entries)
            self.predictions[m]['yPredicted'] = pandas.Series(allPredicted, index=df.index)
            self.predictions[m]['yError'] = pandas.Series(allError, index=df.index)

        # write statistics to shelve if filename given
        if not filename is None:
            pass
    
    def compare2models(self,modelNames=[],modelCalcs=[]):
        '''Compares error distributions of two models by computing 
        F-test and two-sample t-test for equal variance and means respectively.
        This method calls functions from the stats module.
        Input:
        modelNames (list) = one or two modelNames of model(s) that were predicted
        modelNames (list) = one or two model prediction dataframes
        *input should be one of either, or two of one only!
        Output:
        IN BETA... currently prints a few values to bash

        NOTE: SECOND MODEL IS CONSIDERED THE NULL HYPOTHESIS'''

        error = "Make sure to provide a list for modelNames for compare2models!"
        assert isinstance(modelNames,list),error
        error = "Make sure to provide a list for modelCalcs for compare2models!"
        assert isinstance(modelCalcs,list),error
        
        # Compare to internal models just run with the test system
        if len(modelNames) == 2:
            error = "To compare2models, you should run calc_stats first!"
            assert all([name in self.statistics.keys() for name in modelNames]),error

        # Compare a model in the test system to an external one
        elif len(modelNames) == 1 and len(modelCalcs) == 1:
            error = "modelName, {}, provided was not analyzed with calc_stats".format(modelNames[0])
            assert modelNames[0] in self.statistics.keys(),error
            
            error = "Looking for a pandas dataframe as the item in modelCalcs, not {}.".format(type(modelCalcs[0]))
            assert isinstance(modelCalcs[0],pandas.DataFrame),error
            

            self.predictions['Model2'] = modelCalcs[0]
            modelNames.append('Model2')

        # Compare to external models from outside of the test system
        elif len(modelCalcs) == 2:
            error = "You should provide pandas dataframes in a list for compare2models."
            assert all([isinstance(df,pandas.DataFrame) for df in modelCalcs]),error
            error = "yError, or the model error, has not been computed for this dataframe."
            assert all(['yError' in df.keys() for df in modelCalcs]),error
            
            self.predictions['Model1'] = modelCalcs[0]
            self.predictions['Model2'] = modelCalcs[1]
            modelNames = ['Model1','Model2']

        else:
            raise Exception("The number of items for the argument lists in compare2models \
                is {} which exceeds 2.".format(len(modelNames)+len(modelCalcs)))

        # Let's calculate the F-test and the t-tests for these two model error distributions
        x = self.predictions[modelNames[0]]['yError']
        y = self.predictions[modelNames[1]]['yError']
        x = x[~np.isnan(x)]
        y = y[~np.isnan(y)]

        (_,F,F_pval) = stats.vartest2(x,y,logNormal=True,test="F")
        (_,t,t_pval) = stats.ttest2(x,y)
        
        R2_1 = list(self.statistics[modelNames[0]]['Pearson R-squared'])[-1]
        R2_2= list(self.statistics[modelNames[1]]['Pearson R-squared'])[-1]
        mean1 = np.mean(x)
        mean2 = np.mean(y)
        var1 = np.exp(np.var(np.log(x)))
        var2 = np.exp(np.var(np.log(y)))

        output = {
                'R2_1': R2_1,
                'R2_2': R2_2,
                'mean1': mean1,
                'mean2': mean2,
                'var1': var1,
                'var2': var2,
                'F-statistic': F,
                'pval(F)': F_pval,
                't-statistic': t,
                'pval(t)': t_pval}

        return output

    def to_shelve(self,filename):
        '''Export model predictions to shelve.
        Input:
        filename (str) = name of shelve object to create
        Output:
        A python shelve with model names as keys, and pandas dataframes as values'''

        d = shelve.open(filename)
        if not d: # if empty
            d.update(self.predictions)
        else:
            for model in self.models.available:
                pass
                # d[model] = dbms.udpate(d[model],self.predictions[model],self.identifiers)
                
                # ERROR, CODE NEEDS TO BE UPDATED HERE:
                '''d[model] = dbms.udpate(d[model],self.predictions[model],self.identifiers)
                AttributeError: 'module' object has no attribute 'udpate' '''
        d.close()

    def to_excel(self,filename,predictColumns=[],statsColumns=[],models=[]):
        '''Export model predictions and statistics to an Excel workbook with one
        worksheet for each model. Preferred method of creating readable output.
        Input:
        filename (str)        = name of Excel workbook to create
        predictColumns (list) = labels from the pandas dataframes to override
                                automatic alphabetization of all dataframe labels (default behavior)
        statsColumns (list)   = labels from the stats pandas dataframes to override write
                                (same behavior as predictColumns)
        models (list)         = models to export, if [], to_excel writes all predicted models by default
        Output:
        A Excel workbook with model predictions and statistics.'''

        assert isinstance(filename,str), "Filename provided, {}, for export() needs to be a string.".format(filename)
        if not models:
            models = self.predictions.keys()
        else:
            for m in models:
                assert m in self.predictions.keys(), "Model {} was not tested.".format(m)
        if filename[-4:] == ".xlsx":
            fn = filename
        else:
            fn = filename + ".xlsx"
        writer = pandas.ExcelWriter(fn)
        if predictColumns:
            for model in models:
                self.predictions[model].to_excel(writer,sheet_name=model,columns=predictColumns)
        else:
            for model in models:
                self.predictions[model].to_excel(writer,sheet_name=model)
            
        if statsColumns:        
            for model in models:
                if model in self.statistics.keys(): 
                    self.statistics[model].to_excel(writer,sheet_name="{}-stats".format(model),columns=statsColumns)
                else:
                    print "Functional form for {} was not specified. Not writing stats.".format(model)
        else:
            for model in models:
                if model in self.statistics.keys():
                    print "WRITING STATISTICS"
                    self.statistics[model].to_excel(writer,sheet_name="{}-stats".format(model))
                else:
                    print "Functional form for {} was not specified. Not writing stats.".format(model)

        writer.save()

    def add_datasets(self,datasets):
        self.datasets = list(set(self.datasets+datasets))
        self._update_database()

    def remove_datasets(self,datasets):
        self.datasets = [ds for ds in self.datasets[:] if ds not in datasets]
        self._update_database()

    def add_filters(self,filters={}):
        self.filters.update(filters)

    def remove_filters(sef,filters=[]):
        for key in filters:
            self.filters.pop(key,None)
        self._update_database()

    def add_model(self,alias,model,args,kargs):
        self.models.add(alias, model, args, kargs)

    def remove_model(self,alias):
        self.models.remove(alias)


if __name__ == "__main__":
    pass