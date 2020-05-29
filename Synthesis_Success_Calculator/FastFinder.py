import os
import subprocess
import uuid
import atexit
import RepeatFinder
import cPickle as pickle

from itertools import chain

cleanup_files = []

class FastFinder(object):

    def __init__(self, str_type='DNA'):
        self.stype   = [1, 2][str_type == 'RNA']

    def bash_command(self, cmd):
        try:
            subprocess.call(cmd)
        except Exception as e:
            print 'Unexpected error when running RepeatFinder.py: {}'.format(e)
            return False
        return True

    def execute_RepeatFinder(self, seq_list, k_low, dictype, verbose, k_high):
        # Discover RepeatFinder.py locally
        # cwd = os.getcwd()
        RepeatFinder_directory = '/'.join(RepeatFinder.__file__.split('/')[:-1])
        
        # print "cwd: ", cwd
        # print "cng: ", RepeatFinder_directory
        # print RepeatFinder.__file__
        
        #os.chdir(RepeatFinder_directory)

        # Setup IO files
        infile  = RepeatFinder_directory + '/seq_list_{}.txt'.format(str(uuid.uuid4()))
        with open(infile, 'w') as out_file:
            for seq in seq_list:
                out_file.write('{}\n'.format(seq))
        outfile = RepeatFinder_directory + '/repeat_dict_{}.dict'.format(str(uuid.uuid4()))
        global cleanup_files
        cleanup_files.extend([infile, outfile])

        # Call RepeatFinder via system's commandline
        cmd_1 = map(str, ['-u {}'.format(k_high), ''][k_high == None])
        cmd_2 = [['', '-v'][verbose == True]]
        cmd_3 = map(str, ['pypy', RepeatFinder_directory + '/RepeatFinder.py', '-i', infile, '-o', outfile, '-k', k_low, '-s', self.stype, '-d', dictype])
        cmd   = []
        for cmd_el in chain(cmd_3, cmd_2, cmd_1):
            if cmd_el != '':
                cmd.append(cmd_el)

        # Setup Execution and Parsing
        execution   = self.bash_command(cmd)
        repeat_dict = {}

        if execution:
            # Parse back and return repeat_dict
            with open(outfile, 'r') as in_file:
                repeat_dict = pickle.load(in_file)

        # Cleanups and resets
        delete_dicts()

        return repeat_dict

    def get_repeat_dict(self, seq_list, k_low, verbose=False, k_high=None):
        return self.execute_RepeatFinder(seq_list, k_low, dictype=1, verbose=verbose, k_high=k_high)

    def get_alt_repeat_dict(self, seq_list, k_low, verbose=False, k_high=None):
        return self.execute_RepeatFinder(seq_list, k_low, dictype=2, verbose=verbose, k_high=k_high)

@atexit.register
def delete_dicts():
    global cleanup_files
    for file in cleanup_files:
        try:
            os.remove(file)
        except:
            pass
