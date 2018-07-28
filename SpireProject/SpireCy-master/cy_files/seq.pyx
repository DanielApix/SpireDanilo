#!/usr/bin/python

import datetime, os

cimport client_lib
from libc.stdlib cimport free

cdef client_lib.node_t * cfl_list

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

def compute_fingerprint_by_fasta(name_file, second_format):
  print("Compute fingerprint...")
  max_length_fact = 0
  factorization_file = open(name_file)
  if (second_format):
    fingerprint_second_file = open(name_file.replace("results", "fingerprint_second_format"), "w")
  fingerprint_file = open(name_file.replace("results", "fingerprint"), "w")
  txt_fingerprint = ""
  id_read = ''
  s = ' '

  while True:
    if len(s) == 0:
      break

    while True:
      s = factorization_file.readline().rstrip()
      if s == "":
        break
      if s[0] == '>':
        id_read = s
        break
      else:
        s = s.replace('"', '')
        s = s.replace("[", "")
        s = s.replace("]", "")
        list_fact = s.split()
      print ("step1")
      second_txt_fingerprint = ""
      if (second_format):
        second_txt_fingerprint = second_txt_fingerprint + str(id_read) + " $ "
      txt_fingerprint = txt_fingerprint + str(id_read) + " "

      for i in range(len(list_fact)):

        factor = list_fact[i]

        if factor == "<<" or factor == ">>":
          txt_fingerprint = txt_fingerprint + factor
        else:
          txt_fingerprint = txt_fingerprint + str(len(factor))
          if max_length_fact < len(factor):
            max_length_fact = len(factor)

        if i == len(list_fact) - 1:
          txt_fingerprint = txt_fingerprint + "\n"
        else:
          txt_fingerprint = txt_fingerprint + ","

      txt_fingerprint = txt_fingerprint.replace("<<", "0")
      txt_fingerprint = txt_fingerprint.replace(">>", "-1")
      print ("step2")
      if (second_format):
        fingerprint_second_file.write(second_txt_fingerprint + " $ " + s)
      fingerprint_file.write(txt_fingerprint)
      txt_fingerprint = ""

  fingerprint_file.close()
  factorization_file.close()
  return max_length_fact

def check_and_execute_factorization(choice, read, c):
  if (choice == 1):
    return call_c_CFL(read)
  elif (choice == 2):
    return call_c_ICFL(str)
  elif (choice == 3):
    return call_c_CFL_icfl(read, c)
  elif (choice == 4):
    return call_c_ICFL_cfl(read, c)


#Settaggio dei parametri

print('Benvenuto nel programma Spire sequenziale\n')

C = input('Prego, fornisca la dimensione massima di ciascun fattore\n')

dir_path_experiment = input("Fornisca il percorso della cartella dei file Fasta (tra apici \")\n")

fact_choice = input("Indichi quale funzione di fattorizzazione adottare:\n1.CFL\n2.ICFL\n3.CFL con sottofattorizzazioni ICFL\n4.ICFL con sottofattorizzazzioni CFL\n")

#Inizio del processamento dei file

#import main_fingerprint
list_runs = os.listdir(dir_path_experiment)

for run in list_runs:
    run_path = dir_path_experiment + "/" + run
    list_fasta = os.listdir(run_path)
    list_fasta = [file for file in list_fasta if file.endswith('.fasta')]

    if len(list_fasta) == 0:#or len(list_fasta) > 1:
        print("La directory deve contenere un file .fasta")
        continue

    fasta = open(run_path + "/" + list_fasta[0], 'r')
    pos = fasta.tell()  # check file
    row = fasta.readline()
    if row == "" or row[0] != '>':
        print("Errore file fasta")
    fasta.seek(pos)

    filename = '/results_' + run + '.txt'
    mode = 'w' # make a new file if not
    results = open(run_path + filename, mode)

    first = True
    last_block_size = -1
    part = ' '  #
    while True:
        if part == ' ':  #first time#
            results.write(fasta.readline().rstrip() + '\n')  #primo id su file
	#check read and id
        part = ' '
        while part[0] != '>':  # o last_block_size cambiato
            if first:
                read = fasta.readline().rstrip()
                first = False
            else:
                part = fasta.readline().rstrip()
                if part == "":
                    #fact = call_c_CFL(read)
                    fact = check_and_execute_factorization(fact_choice, read, C)
                    #fact = call_c_CFL_icfl(read, C)
                    #fact = apply_factorization(read, C, call_c_CFL_icfl)
                    results.write(str(fact))
                    last_block_size = 10
                    break
                elif part[0] == '>':
                    #part(id) su file e fact su file
                    #fact = call_c_CFL(read)
                    fact = check_and_execute_factorization(fact_choice, read, C)
                    #fact = call_c_CFL_icfl(read, C)
                    #fact = apply_factorization(read, C, call_c_CFL_icfl)
                    results.write(str(fact) + '\n' + str(part) + '\n')
                    first = True
                else:
                    read += part.rstrip()

        if last_block_size != -1:  # end for
            break
    results.close()
    fasta.close()

    list_fasta = os.listdir(run_path)
    list_fasta = [file for file in list_fasta if file.startswith('results')]

    if (fact_choice == 3 or fact_choice == 4):
      #compute_fingerprint_by_fasta(dir_path_experiment + "/" + run + "/" + list_fasta[0], False)
      compute_fingerprint_by_fasta(dir_path_experiment + "/" + run + "/" + list_fasta[0], True)
    else:
      compute_fingerprint_by_fasta(dir_path_experiment + "/" + run + "/" + list_fasta[0], False)

    #list_fasta = os.listdir(run_path)
    #list_fasta = [file for file in list_fasta if file.startswith('fingerprint')]

    #main_fingerprint.kfingerprint(dir_path_experiment + "/" + list_fasta[0], 4)


    print(run + ' processato')
print('Tutti i file sono stati processati')
input() #Per non causare la chiusura immediata della console
