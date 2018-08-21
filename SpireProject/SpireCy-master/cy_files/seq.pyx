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
    fingerprint_second_file = open(name_file.replace("results", "fingerprintsecond_format"), "w")
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
    return call_c_ICFL(read)
  elif (choice == 3):
    return call_c_CFL_icfl(read, c)
  elif (choice == 4):
    return call_c_ICFL_cfl(read, c)

def kfingerprint(name_file, k):
  print("Compute k-fingerprint...")
  fingerprint_file = open(name_file)
  kfinger_file = open(name_file.replace("fingerprint", "k-finger"), 'w')

  while True:
    s = fingerprint_file.readline()
    if len(s) == 0 or s == "\n":
      break

    position = 0
    inner_fact = 0 # 0=fatt.principale, 1=fatt.secondaria
    fact_list = [] # lista lunghezze fattori effettivi
    binary_list = [] # lista fattorizzazioni applicate

    splitted = s.split()
    s = splitted[1].split(',')

    # creazione liste fact e binary
    for i in range(len(s)):
      if s[i] == '0':
        inner_fact = 1
        continue
      if s[i] == '-1':
        inner_fact = 0
        continue

      fact_list.append(int(s[i]))
      binary_list.append(inner_fact)


    if len(fact_list) < k: # se numero dei fattori piu' piccolo della finestra
      diff = k - len(fact_list)
      row = ' '.join(map(str, fact_list)) + " " + ' '.join(['*' for h in range(diff)])
      row += " $ "
      row += ' '.join(map(str, binary_list)) + " " + ' '.join(['*' for h in range(diff)])
      row += " $ "
      row += str(splitted[0]) + " " + str(position)

      kfinger_file.write(row + "\n")

    else: # scorrimento finestre
      for j in range(len(fact_list) - k + 1):
        k_finger = fact_list[j:j+k]
        b_finger = binary_list[j:j+k]

        row = ' '.join(map(str, k_finger))
        row += " $ "
        row += ' '.join(map(str, b_finger))
        row += " $ "
        row += str(splitted[0]) + " " + str(position)  # contiene id

        position += k_finger[0]

        kfinger_file.write(row + "\n")
  kfinger_file.close()



#Settaggio dei parametri

print('Benvenuto nel programma Spire sequenziale\n')

C = input('Prego, fornisca la dimensione massima di ciascun fattore\n')

client_lib.communicate_max_fact_length(C)

dir_path_experiment = input("Fornisca il percorso della cartella dei file Fasta (tra apici \")\n")

fact_choice = input("Indichi quale funzione di fattorizzazione adottare:\n1.CFL\n2.ICFL\n3.CFL con sottofattorizzazioni ICFL\n4.ICFL con sottofattorizzazzioni CFL\n")

k = input("Fornisca la dimensione k per la creazione del file k-fingerprint\n");

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
    for fasta_name in list_fasta:
      fasta = open(run_path + "/" + fasta_name, 'r')
      pos = fasta.tell()  # check file
      row = fasta.readline()
      if row == "" or row[0] != '>':
          print("Errore file fasta")
      fasta.seek(pos)

      filename = '/results_' + fasta_name + '.txt'
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

    for fasta_name in list_fasta:
      if (fact_choice == 3 or fact_choice == 4):
        compute_fingerprint_by_fasta(dir_path_experiment + "/" + run + "/" + fasta_name, False)
        compute_fingerprint_by_fasta(dir_path_experiment + "/" + run + "/" + fasta_name, True)
      else:
        compute_fingerprint_by_fasta(dir_path_experiment + "/" + run + "/" + fasta_name, False)

    list_fasta = os.listdir(run_path)
    list_fasta = [file for file in list_fasta if file.startswith('fingerprint_')]

    for fasta_name in list_fasta:
      kfingerprint(dir_path_experiment + "/" + run + "/" + fasta_name, k)


    print(run + ' processato')
print('Tutti i file sono stati processati')
client_lib.print_statistics()
input() #Per non causare la chiusura immediata della console
