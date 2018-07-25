cimport client_lib
from libc.stdlib cimport free

import time, multiprocessing
try:
    import queue
except ImportError:
    import Queue as queue
from multiprocessing.managers import SyncManager

cdef client_lib.node_t * cfl_list

PORTNUM = 5000
AUTHKEY = b'abc'
IP = ''
NUM_BLOCKS = multiprocessing.cpu_count()/2 if multiprocessing.cpu_count() > 1 else 1 #3-4 ottimo su una macchina(4 core)
MAX_EMPTY_RECEIVED = 10
RESET_EMPTY_AFTER = 500
C = 10

#Settaggio parametri
print("Benvenuto nel modulo client di Spire\n")
IP = input("Prego fornisca indirizzo IP del server (Tra apici \". Se è questa macchina inserisca solo due apici)\n")
C = input("Fornisca la dimensione massima di ciascun fattore\n")

print("Inizio ricezione dei blocchi. Client in attesa\n")

#Funzioni
def call_c_CFL(str):
    cdef char *factorization_c
    cfl_list = client_lib.CFL(str)
    #client_lib.print_list_reverse(cfl_list)
    factorization_c = client_lib.list_to_string(cfl_list, 0)
    #free fact created by malloc in c function (other free need import module level from cpython.mem cimport PyMem_Free)
    try:
        factorization = <bytes> factorization_c
    finally:
        free(factorization_c)
    return factorization


def call_c_CFL_icfl(str, c):
    cdef char *factorization_c
    cfl_list = client_lib.CFL_icfl(str, c)
    #client_lib.print_list_reverse(cfl_list)
    factorization_c = client_lib.list_to_string(cfl_list, 0)
    #free fact created by malloc in c function (other free need import module level from cpython.mem cimport PyMem_Free)
    try:
        factorization = <bytes> factorization_c
    finally:
        free(factorization_c)
    return factorization


def call_c_ICFL(str):
    cdef char *factorization_c
    cfl_list = client_lib.ICFL_recursive(str)
    #client_lib.print_list_reverse(cfl_list)
    factorization_c = client_lib.list_to_string(cfl_list, 1)
    #free fact created by malloc in c function (other free need import module level from cpython.mem cimport PyMem_Free)
    try:
        factorization = <bytes> factorization_c
    finally:
        free(factorization_c)
    return factorization


def call_c_ICFL_cfl(str, c):
    cdef char *factorization_c
    cfl_list = client_lib.ICFL_cfl(str, c)
    #client_lib.print_list_reverse(cfl_list)
    factorization_c = client_lib.list_to_string(cfl_list, 1)
    #free fact created by malloc in c function (other free need import module level from cpython.mem cimport PyMem_Free)
    try:
        factorization = <bytes> factorization_c
    finally:
        free(factorization_c)
    return factorization


def factorizer_worker(job_q, result_q):
    """ A worker function to be launched in a separate process. Takes jobs from
        job_q - each job a list of numbers to factorize. When the job is done,
        the result is placed into result_q. Runs until job_q is empty.
    """
    empty_received = 0
    no_empty_count = 0 # reset empty count
    #to do: provare break e non return eccezioni perche sembra che non esce dalla join
    while True:
        try:
            block = job_q.get_nowait()
            #if block is not None:
            if block == 'fine':
                #try:
                job_q.put_nowait(block)
                return
                #except queue.Full:
                    #return
            read_id = 0
            for i in range(len(block)):
                if i % 2 == 0:
                    read_id = block[i]
                else:
                    block[i] = call_c_CFL(block[i])
            result_q.put(block)
            if no_empty_count == RESET_EMPTY_AFTER:
                no_empty_count = 0
                empty_received = 0
            else:
                no_empty_count += 1
        except queue.Empty:
            time.sleep(1)
            #try:
            if job_q.empty():
                empty_received += 1
                time.sleep(2)
                no_empty_count = 0
            if empty_received > MAX_EMPTY_RECEIVED:
                print('client received too many empty, sleep')
                time.sleep(5)
                #return
            #except (EOFError, IOError) as exception:
             #   print('queue closed, shutting down...')
              #  return
        #gestisce perdita di read in caso di kill da user
        except KeyboardInterrupt: #fattorizza blocco gia preso, poi termina
            print('pressed CTRL+C, waiting for last block factorizations..')
            is_busy = True
            while is_busy: #necessario se utente continua a cercare di chiudere il processo, altrimenti eseguito una volta
                try:
                    while i < (len(block)):
                        if i % 2 == 0:
                            read_id = block[i]
                        else:
                            block[i] = call_c_CFL(block[i])
                        i += 1
                    result_q.put(block)
                    is_busy = False
                except KeyboardInterrupt:
                    print('another CTRL+C, waiting for last block factorizations..')
            return
        #except (EOFError, IOError) as exception:
         #   print('queue closed, shutting down..')
          #  return


def mp_factorizer(shared_job_q, shared_result_q, nprocs):
    """ Split the work with jobs in shared_job_q and results in
        shared_result_q into several processes. Launch each process with
        factorizer_worker as the worker function, and wait until all are
        finished.
    """
    procs = []
    for i in range(nprocs):
        p = multiprocessing.Process(
                target=factorizer_worker,
                args=(shared_job_q, shared_result_q))
        procs.append(p)
        p.start()

    for p in procs:
        p.join()


def make_client_manager(ip, port, authkey):
    """ Create a manager for a client. This manager connects to a server on the
        given address and exposes the get_job_q and get_result_q methods for
        accessing the shared queues from the server.
        Return a manager object.
    """
    class ServerQueueManager(SyncManager):
        pass

    ServerQueueManager.register('get_job_q')
    ServerQueueManager.register('get_result_q')

    manager = ServerQueueManager(address=(ip, port), authkey=authkey)
    manager.connect()

    print('Client connected to %s:%s' % (ip, port))
    return manager


def runclient():
    manager = make_client_manager(IP, PORTNUM, AUTHKEY)
    job_q = manager.get_job_q()
    result_q = manager.get_result_q()
    mp_factorizer(job_q, result_q, NUM_BLOCKS) # num of blocks presi dalla coda


runclient()
