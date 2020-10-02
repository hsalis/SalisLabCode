"""
Custom database management system built by subclassing pandas classes

Copyright 2017 Alexander C. Reis, Howard M. Salis, all rights reserved.

"""

import cPickle as pickle
import pandas as pd

'''Load a database from file, creates and returns a DataBase instance
Inputs: filename (string) :: the filename with the extension
Output: DB (DataBase)     :: a DataBase instance with data from filename'''
def load(filename):
    
    if filename.endswith('.p'):
        with open(filename,'rb') as handle:
            data = pickle.load(handle)

    # load from 
    elif filename.endswith('.csv'):
        data = pd.read_csv(filename)

    else:
        raise Exception('Filename, {}, should be a pickled file (.p) or a csv (.csv)'.format(filename))

    DB = DataBase(data)
    return DB


class DataBase(object):

    def __init__(self,data=None):

        if isinstance(data,dict):
            data = remove_unicode(data)
            self.data = pd.DataFrame(data)

        elif isinstance(data,pd.DataFrame):
            data = remove_unicode(data)
            self.data = data

        else:
            self.data = pd.DataFrame()

    def __repr__(self):
        return str(self.data)

    def __len__(self):
        return len(self.data)

    ''' Alternative method to add_data(), other is either a dict or a DataFrame'''
    def __add__(self,other):
        self.add_data(other)
        return self

    ''' Alternative method to remove_data(), other is either a dict or a DataFrame'''
    def __sub__(self,other):
        self.remove_data(other)
        return self

    ''' Alternative method to filter_data(), other is either a dict or a DataFrame'''
    def __mul__(self,other):
        self.filter_data(other)
        return self

    '''Use add_data() to append data in form of DataFrame or dict
    Input: data :: Either a dictionary or a DataFrame'''
    def add_data(self,data):
        assert isinstance(data,(dict,pd.DataFrame))
        if isinstance(data,dict):
            data = pd.DataFrame(data)
        data = remove_unicode(data)
        self.data = self.data.append(data,ignore_index=True)

    '''Use remove() to remove entries that match kargs
    Input:  database (pandas DataFrame)
            kargs (dictionary)  :: keys=database label, values=a list of filter values
            ordered (bool)      :: specifies if kargs contains lists with 
                                   of record attributes that are ordered            
    Output: Database filtered to remove kargs'''
    def remove_data(self,kargs,ordered=False):
        self.data = self.data[~self.get_indexes(kargs,ordered)]

    '''Use filter() to filter to remove all but a subset defined by kargs
    Input:  kargs (dictionary)  :: keys=database label, values=a list of filter values
            ordered (bool)      :: specifies if kargs contains lists with 
                                   of record attributes that are ordered
    Output: Database filtered to only include kargs'''
    def filter_data(self,kargs,ordered=False):
        self.data = self.data[get_indexes(kargs,ordered)]

    # Gets indexes for values in database that match kargs
    def get_indexes(self,kargs,ordered=False,allQueries=False):
        assert isinstance(ordered,bool), "Argument, ordered, should be a boolean. Type given = {}".format(type(ordered))
        kargs = {k.upper(): v for k,v in kargs.iteritems()}
        keylist = []
        for key in kargs.keys():
            assert key in self.data.keys(), "{} is not a label in the database".format(key)
            keylist.append(key)
        if ordered:
            indexes = []
            for valuetup in zip(kargs.itervalues()):
                boolarray = np.all([self.data[keylist[i]]==valuetup[i] for i in range(len(valuetup))],axis=1)
                indx = [i for i,x in enumerate(boolarray) if x]
                if not indx:
                    indexes.append(None)
                else:
                    indexes += indx
        else:
            indexes = self.data[kargs.keys()].isin(kargs).all(1)
        if not allQueries:
            indexes = [indx for indx in indexes[:] if not indx is None]
        return indexes

    '''Return Database as a list of dictionaries, where each entry (row)
    is a dictionary with the DataFrame labels as keys (sorted by index).'''
    def get_entries(self):
        asDict = self.data.T.to_dict()
        indxs = sorted(asDict.keys())
        return [asDict[i] for i in indxs]

    '''Save the DataBase.data to ``filename`` for persistance;
    Input:  filename (string) :: Name of the file without extension to save to
            type (string) :: Can be ``pickle`` (default) or ``csv`` '''
    def save(self,filename,type='pickle'):

        if filename.endswith('.csv'): fn = filename[:-4]
        elif filename.endswith('.p'): fn = filename[:-2]
        else:                         fn = filename

        # pickle data
        if type == 'pickle':
            with open(fn+'.p','wb') as handle:
                pickle.dump(self.data,handle,protocol=2)

        # write to human readable csv
        elif type == 'csv':
            self.data.to_csv(fn+'.csv')

        else:
            raise Exception('Bad argument value for ``type`` in method DataBase.save().')


'''Use to remove unicode from a dict or a DataFrame'''
def remove_unicode(d):
    for label in d.keys():
        if isinstance(d[label],list):
            if isinstance(d[label][0],unicode):
                d[label] = [s.encode('utf-8') for s in d[label]]
        elif isinstance(d[label],unicode):
            d[label] = d[label].str.encode('utf-8')
    return d

if __name__ == "__main__":
    DB = Database()