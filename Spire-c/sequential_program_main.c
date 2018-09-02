/*
  Main program of Sequential Spire-c. It allows the user to factorize fasta files containing reads.
  Author: Danilo Apicella
  Contact: danielapix@hotmail.it
*/

#include <stdio.h>
#include <string.h>
#include <dirent.h>
#include <errno.h>
#include <stdlib.h>
#include "factorizations.h"
#include "utils.h"

int fact_choice;

char root_path[100];     //...of the directory to process
int max_fact_length = 4; //arbitrary chosen

FILE *factorization_file, *fingerprint_file, *kfingerprint_file;


void for_each_element_in(char* directory_path,  void (*apply_function) (struct dirent *, char *));

char* append_filename_to_path(char* path, char *name);
char* process_and_write_in_file(char* to_process, char* (*process_function) (), FILE* file_to_write, char* path);
char* create_fingerprint(char* factorized_genom);
char* create_k_fingerprint(char* fingerprint);

char *apply_factorization(char *genom);

void open_towrite_file(char *name, char *fasta_name, char *directory_path);
void process_fasta(struct dirent *file_description, char *directory_path);
void process_all_fasta_files(struct dirent *subdirectory_description, char* current_path);

/*
Struttura generale del progamma

int file_aperti

apri_in_scrittura(char* filename)


Chiedi input dall'utente
Per ogni sottodirectory nella root
  per ogni file nella sottodirectory
    se è un file fasta
      apri in lettura il file
      se l'apertura è fallita
        termina esecuzione funzione
        stampa errore
      setta nome results
      apri in scrittura("results")
      apri in scrittura("fingerprint")
      apri in scrittura("k-fingerprint")
      per ogni read
        se è una read valida
          char* result = processa_e_salva_nel_file(read, processa_read, file_fattorizzazione)
          processa_e_salva_nel_file(result, processa_fingerprint, file_fingerprint)
          read = processa_e_salva_nel_file(read, processa_read, file_kfingerprint)
        altrimenti
          esci dal ciclo
stampa le statistiche
*/

int main() {

  char path[100];
  DIR *root_directory;
  DIR *subdirectory;
  struct dirent *dir;
  struct dirent *inner_file; // cointained in one of the subdirectories

  char str[100];

  /*Catching input from the user*/

  printf("Benvenuto nel programma sequenziale Spire\n\n");

  printf("Fornisca la directory dei file fasta\n");
  scanf("%s", root_path); //warning: not possible to enter more than 100 characters and with spaces

  do {
    printf("Fornisca una delle seguenti opzioni di fattorizzazione\n");
    printf("1.CFL\n");
    printf("2.ICFL\n");
    printf("3.CFL seguito da ICFL\n");
    printf("4.ICFL seguito da CFL\n");
    printf("Risposta:\n");
    scanf("%d", &fact_choice);
    if (fact_choice >= 1 && fact_choice <= 4)
      break;
    else
      printf("Spiacente, la risposta data non è nelle opzioni possibili\n\n");
  } while (1);

  if (fact_choice > 2) {
    printf("Fornisca dimensione massima di ciascun fattore\n");
    scanf("%d", &max_fact_length);
  }

  for_each_element_in(root_path, process_all_fasta_files);

  return 0;
}

/*
  param apply_function: the first param of this function is referred to the description of the current examined file and
      the second one for the path of the current examined file
  pre-condition: path refers to an existing directory location
*/
void for_each_element_in(char* directory_path,  void (*apply_function) (struct dirent *, char *)) {

  DIR *file = opendir(directory_path);
  struct dirent *inner_file;

  if (file) {
     while ((inner_file = readdir(file)) != NULL) {
       printf("arguments passed: %s, %s\n", directory_path, inner_file->d_name);
       (*apply_function)(inner_file, append_filename_to_path(directory_path, inner_file->d_name));
     }
  }
  else {
    printf("Errore nella scansione della cartella\nNon è stato possibile aprire il file %s\n%s\n", directory_path, strerror(errno));
  }
  closedir(file);
}

/*adding '/' to create a correct path*/
char* append_filename_to_path(char* path, char *name) {
  char *new_path = malloc(strlen(path) + strlen(name) + 3);

  strcpy(new_path, path);
  strcat(new_path, "/");
  strcat(new_path, name);
  return new_path;
}


void process_all_fasta_files(struct dirent *subdirectory_description, char* current_path) {
  DIR *subdirectory;

  if ((strcmp(subdirectory_description->d_name, ".") != 0) && (strcmp(subdirectory_description->d_name, "..") != 0)) {
    printf("processing %s\n", current_path);
    for_each_element_in(current_path, process_fasta);
  }
}


/*processes fasta file specified by the file description in the directory specified by the path*/
void process_fasta(struct dirent *file_description, char *path) {

  char header_read[300];
  char genom_read[300];   //also body of read
  char directory_path[100];
  char str[300];
  FILE *fasta_file;
  int error;
  char *result;
  char *result2;

  if (strlen(file_description->d_name) > strlen(".fasta")) {
    if (strstr(file_description->d_name, ".fasta") != NULL) {   //if it has .fasta extention

      strncpy(directory_path, path, strlen(path) - strlen(file_description->d_name) - 1);
      open_towrite_file("factorization", file_description->d_name, directory_path);
      open_towrite_file("fingerprint", file_description->d_name, directory_path);
      open_towrite_file("kfingerprint", file_description->d_name, directory_path);

      if ((fasta_file = fopen(path, "r")) == NULL) {
        printf("error: %s\n\n", strerror(errno));
      }
      /*Processing of fasta file*/
      if (fasta_file != NULL) {
        while (fscanf(fasta_file, "%s %s", header_read, genom_read) != EOF) {
          if (factorization_file != NULL) {
            result = apply_factorization(genom_read);
            error = fprintf(factorization_file, "%s\n%s\n", header_read, result);
            result = create_fingerprint(result);
            fprintf(fingerprint_file, "%s %s\n", header_read, result);
 //           printf("factorization: %s\n", result);
           // fprintf(fingerprint_file, "%s %s\n", header_read, result);
          }
        }
        fclose(factorization_file);
      }
      else {
        printf("Non è stato possibile aprire il file fasta %s\nErrore: %s\n", file_description->d_name, strerror(errno));
      }
      /**/
    }
  }
}


void open_towrite_file(char *name, char *fasta_name, char *directory_path) {

  FILE *file_to_open;
  char *filename;

  filename = malloc(strlen(name) + strlen(fasta_name));
  strncpy(filename, fasta_name, strlen(fasta_name) - strlen(".fasta"));
  strcat(filename, "-");
  strcat(filename, name);

  if (strcmp(name, "factorization") == 0) {
    file_to_open = factorization_file = fopen(append_filename_to_path(directory_path, filename), "w");
  }
  else if (strcmp(name, "fingerprint") == 0) {
    file_to_open = fingerprint_file = fopen(append_filename_to_path(directory_path, filename), "w");
  }
  else if (strcmp(name, "kfingerprint") == 0) {
    file_to_open = kfingerprint_file = fopen(append_filename_to_path(directory_path, filename), "w");
  }

  if (file_to_open == NULL) {
    printf("Non è stato possibile aprire in scrittura il file %s: %s\n\n", filename, strerror(errno));
  }
}

char *apply_factorization(char *genom) {

  node_t *factorized_genom;
  int second_parameter_value;
  char *factorized_genom_string;

  switch (fact_choice) {
    case 1:
      factorized_genom = CFL(genom);
      second_parameter_value = 0;
      break;
    case 2:
      factorized_genom = ICFL_recursive(genom);
      second_parameter_value = 1;
      break;
    case 3:
      factorized_genom = CFL_icfl(genom, max_fact_length);
      second_parameter_value = 0;
      break;
    case 4:
      factorized_genom = ICFL_cfl(genom, max_fact_length);
      second_parameter_value = 1;
      break;
  }
  return list_to_string(factorized_genom, second_parameter_value);
}

/*path is referred to the location of the file to write*/
char* process_and_write_in_file(char* to_process, char* (*process_function) (), FILE* file_to_write, char* path) {

  char* p = NULL;
  return p;
}

/*pre-condition: param must be the result of factorize_read or format equivalent*/
char* create_fingerprint(char* factorized_genom) {

  int i = 2;
  int j = 0;
  int cont = 0;
  int dim = strlen(factorized_genom);
  char *fingerprint = malloc(dim);
  char converted_number[5];

  while (1) {

    switch (factorized_genom[i]) {
      case '<':
        fingerprint[j++] = '-';
        fingerprint[j++] = '1';
        fingerprint[j++] = ',';
        i += 4;
        //printf("next character: %c\n", factorized_genom[i]);
        break;
      case '>':
        fingerprint[j++] = '0';
        fingerprint[j++] = ',';
        i += 4;
        break;
      case ']':
        j--;
        fingerprint[j] = '\0';
        return fingerprint;
//        printf("found '['\n");
        break;
      case '\"':
        if (cont > 0) {   //if it defines the end of a factor
          sprintf(converted_number, "%d", cont);
          strcat(fingerprint, converted_number);
          j += strlen(converted_number);
          fingerprint[j++] = ',';
          cont = 0;
        }
        else {
          i++;
        }
        break;
      default:
        i++;
        cont++;   //cause it is one of the characters
    }
    fingerprint[j] = '\0';
  }
}

/*pre-condition: param must be the result of create_fingerprint or format equivalent*/
char* create_k_fingerprint(char* fingerprint) {
}

