#import mpi4py
#mpi4py.rc.threaded = False
import sys
from mpi4py import MPI

class Pool(object):

    def __init__(self, comm, master=0):
        assert comm.size > 1
        assert 0 <= master < comm.size
        self.comm = comm
        self.master = master
        self.workers = set(range(comm.size))
        self.workers.discard(self.master)

    def is_master(self):
        return self.master == self.comm.rank

    def is_worker(self):
        return self.comm.rank in self.workers

    def map(self, function, iterable):
        assert self.is_master()

        comm = self.comm
        workerset = self.workers.copy()
        tasklist = [(tid, (function, arg)) for tid, arg in enumerate(iterable)]
        resultlist = [None] * len(tasklist)
        pending = len(tasklist)

        while pending:

            if workerset and tasklist:
                worker = workerset.pop()
                taskid, task = tasklist.pop()
                comm.send(task, dest=worker, tag=taskid)

            if tasklist:
                flag = comm.Iprobe(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG)
                if not flag: continue
            else:
                comm.Probe(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG)

            status = MPI.Status()
            result = comm.recv(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=status)
            worker = status.source
            workerset.add(worker)
            taskid = status.tag
            resultlist[taskid] = result
            pending -= 1

        return resultlist

    def start(self):
        if not self.is_worker(): return
        comm = self.comm
        master = self.master
        status = MPI.Status()
        while True:
            task = comm.recv(source=master, tag=MPI.ANY_TAG, status=status)
            if task is None:
                print "Worker %s received MPI_Pool.stop(). Exiting." % self.comm.rank
                break
            function, arg = task
            result = function(arg)
            comm.ssend(result, master, status.tag)
        sys.exit()

    def stop(self):
        if not self.is_master(): return
        for worker in self.workers:
            self.comm.send(None, worker, 0)


if __name__ == '__main__':

    pool = Pool(MPI.COMM_WORLD)

    def sq(x): return x*x

    pool.start()

    if pool.is_master():

        tic = MPI.Wtime()
        res = pool.map(sq, range(100))
        toc = MPI.Wtime()

        for y, x in zip(res, range(100)):
            assert y ==  sq(x)

    pool.stop()
