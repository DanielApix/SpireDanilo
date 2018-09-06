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
#include <time.h>
#include "factorizations.h"
#include "utils.h"

int fact_choice;

char root_path[100];     //...of the directory to process
int max_fact_length = 0; //arbitrary chosen and requested to the user

FILE *factorization_file, *fingerprint_file, *kfingerprint_file, *oneformat_file;

char header_read[300];  //refers to the current read

/*these refer to the k-fingerprint calculus. To better understand refer to the documentation.*/
int *zero_one_tail, *finger_tail;
int end_window_limit = -1;
int window_dimension;
int is_one = 0;   //boolean equivalent to determine if current finger refers to a second time factorized factor
int start_window_limit = 0;
int cont_shift = 0;  //number of shitfts from the start of the read dealing with the current window
int has_not_been_filled_once = 1;  //if it hasn't been filled at least once for the current read
int is_last_fingerprint = 0;
int recived_exactly_k_fingers = 0;
/**/


void for_each_element_in(char* directory_path,  void (*apply_function) (struct dirent *, char *));

char* append_filename_to_path(char* path, char *name);
char* process_and_write_in_file(char* to_process, char* (*process_function) (), FILE* file_to_write, char* path);
char* create_fingerprint(char* factorized_genom);

/*all these are exclusively referred to the k-fingerprint calculus*/

void pop_tail(int fingerprint_number);
void flush();
void initialize_tail();
void fill_k_fingerprint(int fingerprint_number);
void initialize_k_finger();
/**/

char *apply_factorization(char *genom);

void open_towrite_file(char *name, char *fasta_name, char *directory_path); void process_fasta(struct dirent *file_description, char 
*directory_path); void process_all_fasta_files(struct dirent *subdirectory_description, char* current_path);

/* Struttura generale del progamma

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

  communicate_max_fact_length(max_fact_length);

  printf("fornisca il numero di elementi per ciascuna finestra per le k-fingerprint\n");
  scanf("%d", &window_dimension);

  time_t m;
  time_t now = time(NULL);

  initialize_k_finger();
  for_each_element_in(root_path, process_all_fasta_files);

  print_statistics();

  m = difftime(time(NULL), now);
  printf("tempo totale in secondi: %ld\n",m);
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
    //   printf("arguments passed: %s, %s\n", directory_path, inner_file->d_name);
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
    //printf("processing %s\n", current_path);
    for_each_element_in(current_path, process_fasta);
  }
}


/*processes fasta file specified by the file description in the directory specified by the path*/
void process_fasta(struct dirent *file_description, char *path) {

  char genom_read[300];   //also body of read
  char directory_path[100];
  char str[300];
  FILE *fasta_file;
  int error;
  char *result1;
  char *result2;

  int i;  //can be deleted

  if (strlen(file_description->d_name) > strlen(".fasta")) {
    if (strstr(file_description->d_name, ".fasta") != NULL) {   //if it has .fasta extention

      strncpy(directory_path, path, strlen(path) - strlen(file_description->d_name) - 1);
      printf("d_name: %s\n", file_description->d_name);
      open_towrite_file("factorization", file_description->d_name, directory_path);
      open_towrite_file("fingerprint", file_description->d_name, directory_path);
      open_towrite_file("kfingerprint", file_description->d_name, directory_path);
      open_towrite_file("oneformat", file_description->d_name, directory_path);

      if ((fasta_file = fopen(path, "r")) == NULL) {
        printf("error: %s\n\n", strerror(errno));
      }
      /*Processing of fasta file*/
      if (fasta_file != NULL) {
        while (fscanf(fasta_file, "%s %s", header_read, genom_read) != EOF) {
          if (factorization_file != NULL) {
            result1 = apply_factorization(genom_read);
            //printf("worked\n");
            error = fprintf(factorization_file, "%s\n%s\n", header_read, result1);
            //printf("not seen\n");
            result2 = create_fingerprint(result1);
            fprintf(fingerprint_file, "%s %s\n", header_read, result2);
            fprintf(oneformat_file, "%s %c %s %c %s\n", header_read, '$', result2, '$', result1);
          }
        }
        fclose(factorization_file);
        fclose(fingerprint_file);
        fclose(kfingerprint_file);
        fclose(oneformat_file);
      }
    }
    else {
      printf("Non è stato possibile aprire il file fasta %s\nErrore: %s\n", file_description->d_name, strerror(errno));
    }
  }
}


void open_towrite_file(char *name, char *fasta_name, char *directory_path) {

  FILE *file_to_open;
  char *filename;

  filename = malloc(strlen(name) + strlen(fasta_name));
  strncpy(filename, fasta_name, strlen(fasta_name) - strlen(".fasta"));
  filename[strlen(fasta_name) - strlen(".fasta")] = '\0';
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
  else if (strcmp(name, "oneformat") == 0) {
    file_to_open = oneformat_file = fopen(append_filename_to_path(directory_path, filename), "w");
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
//    printf("found %c\n", factorized_genom[i]);
//    printf("blocked with %s at position %d\n", factorized_genom, i);
    switch (factorized_genom[i]) {
      case '<':
        fingerprint[j++] = '-';
        fingerprint[j++] = '1';
        fingerprint[j++] = ',';
        i += 4;
        fill_k_fingerprint(-1);
        //printf("next character: %c\n", factorized_genom[i]);
        break;
      case '>':
        fingerprint[j++] = '0';
        fingerprint[j++] = ',';
        fill_k_fingerprint(0);
        i += 4;
        break;
      case ']':
        j--;
        fingerprint[j] = '\0';
        fill_k_fingerprint(-2);
        return fingerprint;
//        printf("found '['\n");
        break;
      case '\"':
        if (cont > 0) {   //if it defines the end of a factor
          sprintf(converted_number, "%d", cont);
          strcat(fingerprint, converted_number);
          i += 2;
          j += strlen(converted_number);
          fingerprint[j++] = ',';
//          printf("probable end: %c\n", factorized_genom[i + 3]);
          fill_k_fingerprint(cont);
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

/*Adds a new element to the fingerprint tail*/
void fill_k_fingerprint(int fingerprint_number) {

  //printf("got %d\n", fingerprint_number);
  if (fingerprint_number == -2) {
    if (!recived_exactly_k_fingers)
      flush();
    initialize_tail();
    return;
  }

  //printf("inizio processamento k-fingerprint\n");
  if ( (fingerprint_number == -1) || (fingerprint_number == 0) ) {
    is_one = (fingerprint_number == -1);
//    printf("revealed is_one: %d\n", is_one);
  }
  else {
    if (has_not_been_filled_once) {
      end_window_limit++;
//      printf("has_not_been_filled_once\n");
      pop_tail(fingerprint_number);
      if (end_window_limit == window_dimension - 1) {
//        printf("exiting the has_not_been_filled_once\n");
        flush();
        recived_exactly_k_fingers = 1;
        end_window_limit = 0;
        start_window_limit = (start_window_limit + 1) % window_dimension;
        has_not_been_filled_once = 0;
      }
    }
    else {
//      printf("from tail as empty\n");
      cont_shift += finger_tail[start_window_limit];
      pop_tail(fingerprint_number);
      flush();
      start_window_limit = (start_window_limit + 1) % window_dimension;
      end_window_limit = (end_window_limit + 1) % window_dimension;
    }
  }
}

void pop_tail(int fingerprint_number) {

  finger_tail[end_window_limit] = fingerprint_number;
  zero_one_tail[end_window_limit] = is_one;
}

void flush() {

  const int FINGERPRINT_MAX_DIMENSION = 3;

  int i, number_of_iterations = 0;
  /*n window element * FINGERPRINT_MAX_DIMENSION + number of spaces(n - 1).
    result = FINGERPRINT_MAX_DIMENSION * n + n  - 1 = (FINGERPRINT_MAX_DIMENSION + 1)n - 1. Considering
    than the end string character ('\0'), we have result + 1 =
    = (FINGERPRINT_MAX_DIMENSION + 1) * n - 1 + 1 = (FINGERPRINT_MAX_DIMENSION + 1) * n as resulting dimension*/
  int string_dimension = ((FINGERPRINT_MAX_DIMENSION + 1) * window_dimension) ;
  char zero_one_result[string_dimension], fingerprint_result[string_dimension];

  char converted_number[FINGERPRINT_MAX_DIMENSION + 1]; char converted_bit[3 + 1];  //considering the end character
  int current_length = 0;
  int current_length2 = 0;

//  printf("start of flush\n");
  strcpy(zero_one_result, "");
  strcpy(fingerprint_result, "");

/*  printf("dealing with the following values\n");
  printf("about fingerprints\n");
  printf("%d %d %d %d\n", finger_tail[0], finger_tail[1], finger_tail[2], finger_tail[3]);
  printf("start_window_limit: %d, end_window_limit: %d, window_dimension = %d\n", start_window_limit, end_window_limit, window_dimension); */
  for (int i = start_window_limit; number_of_iterations < window_dimension; i = (i + 1) % window_dimension) {
//    printf("i: %d\n", i);
    sprintf(converted_number, "%d ", finger_tail[i]);
    sprintf(converted_bit, "%d ", zero_one_tail[i]);
    //printf("first section\n");
    current_length += strlen(converted_number);
    strcat(fingerprint_result, converted_number);
   // printf("second section\n");
    fingerprint_result[current_length] = '\0';   //userful to prevent dirty characters forward program failure
    current_length2 = strlen(zero_one_result);
   // printf("third section\n");
    strcat(zero_one_result, converted_bit);
    zero_one_result[current_length2 + 3] = '\0'; //userful to prevent dirty characters forward program failure
    number_of_iterations++;
  }

  fprintf(kfingerprint_file, "%s %c %s %c %s %d\n", fingerprint_result, '$', zero_one_result, '$', header_read, cont_shift);
}

void initialize_tail() {

  int i = 0;

  cont_shift = 0;
  has_not_been_filled_once = 1;
  end_window_limit = -1;
  start_window_limit = 0;
  is_one = 0;
  is_last_fingerprint = 0;
  recived_exactly_k_fingers = 0;

  for (i = 0; i < window_dimension; i++) {
    finger_tail[i] = zero_one_tail[i] = 0;
  }
}

/*Sets the number of elements for each window and sets the tails*/
void initialize_k_finger() {
  zero_one_tail = calloc(window_dimension, sizeof(int));
  finger_tail = calloc(window_dimension, sizeof(int));
}

