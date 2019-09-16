/*
  Main program of Sequential Spire-c. It allows the user to process fasta files containing reads.
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
#include <assert.h>   //the last 4 are used for testing
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <pthread.h>
#include <sys/wait.h>

#define NUM_PROCESSES 4

typedef struct {
  int *finger_tail;
  int *zero_one_tail;
  int end_window_limit;
  int window_dimension;
  int is_one;   //boolean equivalent to determine if current finger refers to a second time factorized factor
  int start_window_limit;
  int cont_shift;  //number of shitfts from the start of the read dealing with the current window
  int has_not_been_filled_once;  //if it hasn't been filled at least once for the current read
  int recived_exactly_k_fingers;
} tail;


struct stat st = {0};

int fact_choice;

int max_fact_length = 0; //arbitrary chosen and requested to the user


char *header_read;  //refers to the current read
char *genom_read;
int current_header_size = 0;
int current_genom_size = 0;

int window_dimension;

time_t time_spent_to_read_file = 0;
time_t time_spent_to_write_in_file = 0;
FILE *log_file = NULL;

/*necessary to set the dimension of the string returned by list_to_string (efficiency reasons)*/
int get_number_of_factors();
int get_number_of_delimeters();
void set_number_of_elements(int num);
void set_read_dimension(int value);
void communicate_max_fact_length(int c);
void print_statistics();

void for_each_element_in(char* directory_path,  void (*apply_function) (struct dirent *, char *));

char* append_filename_to_path(char* path, char *name);
char* process_and_write_in_file(char* to_process, char* (*process_function) (), FILE* file_to_write, char* path);
char* create_fingerprint(char* factorized_genom, FILE *file_to_write, char *header);

/*all these are exclusively referred to the k-fingerprint calculus*/

void pop_tail(int fingerprint_number, tail *t);
void flush(FILE *file_to_write, tail *t, char* header);
tail *initialize_tail(tail *t);
void fill_k_fingerprint(int fingerprint_number, FILE *file_to_write, tail *t, char* header_to_pass);
void test_initialize_tail();
/**/

char *apply_factorization(char *genom);

FILE *open_towrite_file(char *name, char *fasta_name, char *directory_path);
void process_fasta(struct dirent *file_description, char *directory_path);
void process_all_fasta_files(struct dirent *subdirectory_description, char* current_path);

/* Struttura generale del progamma
avvia_unit_test()  //può essere evitato
avvia_integration_test()  //può essere evitato

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
      apri in scrittura("oneformat");
      per ogni read
          char* result = processa_e_salva_nel_file(read, fattorizza_read, file_fattorizzazione)
          processa_e_salva_nel_file(result, crea_fingerprint, file_fingerprint)  //causa la creazione di k-fingerprint
          crea_one_format_con_dati_precedenti();
stampa le statistiche
*/

/*Reads an entire stream line until it reaches the \n character and records it in a pointer. If there is no
  line to be read or error occurs, NULL is returned.
  if current_size is zero, that is a buffer is considered as it hadn't been previousely allocated, another one will be.
  NOTE: buffer doesn't store the result at the end of the processing.
  pre-condition buffer: buffer must point to a previousely dinamically allocated memory or to NULL
  pre-condition current_size: must point to a variable whose value is the dimension of buffer (== 0 if buffer is NULL)
  pre-condition stream: must point to an existing file
  return: next string read up to the next '\n' reached or NULL if reading fails
*/
char* safe_fgets(char* buffer, int *current_size, FILE *stream) {


  char* tmp;
  char* aux_tmp;
  const int STARTING_DIMENSION = 5;
  int actual_dimension;

  tmp = NULL;

  if (*current_size == 0) {
    buffer = malloc(STARTING_DIMENSION * sizeof(char));
    *current_size = STARTING_DIMENSION;
  }
  if (fgets(buffer, *current_size, stream) == NULL) {
    *current_size = 0;
    return NULL;
  }
  while (buffer[strlen(buffer) - 1] != '\n') {
    *current_size = (*current_size * 2) + 1;
    actual_dimension = *current_size * sizeof(char);
    tmp = (tmp == NULL) ? malloc(actual_dimension) : realloc(tmp, actual_dimension);
    strcpy(tmp, buffer);
    tmp[strlen(buffer)] = '\0';
    fgets(buffer, strlen(buffer), stream);
    strcat(tmp, buffer);

    aux_tmp = buffer;
    buffer = tmp;
    tmp = aux_tmp;
  }
  buffer[strlen(buffer) - 1] = '\0';
  buffer[strlen(buffer)] = '\0';
  if (tmp != NULL)
    free(tmp);
  return buffer;
}

/*
  It scans the directory of the given path and applies the function defined to each contained file or subdirectory
  passing as argument relevant information as description of the file and it's path.
  param directory_path: must contain the name of the directory in addition to the location description
  param apply_function: the first param of this function is referred to the description of the current examined file and
      the second one to the path of the current examined file
  pre-condition directory_path: refers to an existing directory location
  pre-condition apply_function: != NULL
  post-condition: apply_function is applied to each file or subdirectory in the directory specified by the path
*/
void for_each_element_in(char* directory_path,  void (*apply_function) (struct dirent *, char *)) {

  DIR *file = opendir(directory_path);
  struct dirent *inner_file;

  if (file) {
     while ((inner_file = readdir(file)) != NULL) {
       (*apply_function)(inner_file, append_filename_to_path(directory_path, inner_file->d_name));
     }
  }
  else {
    printf("Errore nella scansione della cartella\nNon è stato possibile aprire il file %s\n%s\n", directory_path, strerror(errno));
  }
  closedir(file);
}

/*...adding '/' to create a correct path.
  pre-condition path: != NULL
  pre-condition name: != NULL
  return: *path + '/' + *name
*/
char* append_filename_to_path(char* path, char *name) {
  char *new_path = malloc(strlen(path) + strlen(name) + 3);

  strcpy(new_path, path);
  strcat(new_path, "/");
  strcat(new_path, name);
  return new_path;
}

/*
  pre-condition subdirectory_description: != NULL
  pre-condition current_path: must refer to the file descripted by subdirectory_description
  post-condition: all fasta files in the specified directory will be processed
*/
void process_all_fasta_files(struct dirent *subdirectory_description, char* current_path) {

  if ((strcmp(subdirectory_description->d_name, ".") != 0) && (strcmp(subdirectory_description->d_name, "..") != 0)) {
      for_each_element_in(current_path, process_fasta);
  }
}

void print_error(const char *message, const char *value) {
  
  if (log_file == NULL) {
    log_file = fopen("log_file.txt", "w");
  }

  if (log_file != NULL)
    fprintf(log_file, "%s %s\n", message, value);
}

/*If the file descripted is a fasta file, factorization, fingerprint, kfingerprint and the one format of the first three
  will be created and saved in the same directory containing the file to be processed
  pre-condition file_description: must descript an existing file with reads that respect the correct format
  pre-condition path: must refer to the descripted file
  pre-condition: fact_choice >= 1 || fact_choice <= 4

  post-condition: given the name nam of the fasta file without ".fasta" at the end, nam-factorization, nam-fingerprint,
      nam-kfingerprint, nam-oneformat will be created in the same directory of the fasta file with the respective output inside
*/
void process_fasta(struct dirent *file_description, char *path) {

  char directory_path[100];
  FILE *fasta_file;
  char *result1;
  char *result2;
  char filename[255];
  char *header;
  char *genom;
  int current_header_size = 0;
  int current_genom_size = 0;
  FILE *factorization_file;
  FILE *fingerprint_file;
  FILE *kfingerprint_file;
  FILE *oneformat_file;

  if (strlen(file_description->d_name) > strlen(".fasta")) {
    if (strstr(file_description->d_name, ".fasta") != NULL) {   //if it has .fasta extention
      directory_path[strlen(path) - strlen(file_description->d_name) - 1] = '\0';
      strncpy(directory_path, path, strlen(path) - strlen(file_description->d_name) - 1);
      strcpy(filename, file_description->d_name);  //useful cause sometimes file_description->d_name strangely changes in opening files
      factorization_file = open_towrite_file("factorization", filename, directory_path);
      fingerprint_file = open_towrite_file("fingerprint", filename, directory_path);
      kfingerprint_file = open_towrite_file("kfingerprint", filename, directory_path);
      oneformat_file = open_towrite_file("oneformat", filename, directory_path);

      if ((fasta_file = fopen(path, "r")) == NULL) {
        printf("error: %s\n\n", strerror(errno));
      }
      /*Processing of fasta file*/
      if (fasta_file != NULL) {
next:   while ( (header = safe_fgets(header, &current_header_size, fasta_file)) != NULL ) {
          genom = safe_fgets(genom, &current_genom_size, fasta_file);

          if (genom == NULL) {
            printf("got null\n");
            print_error("Incorrect read genom for", header);
            goto next;
          } else {
            int error = 0, i;
            for(i = 0; i < current_genom_size; i++) {
              switch(genom[i]) {
                printf("%c\n", genom[i]);
                case 'A': break;
                case 'G': break;
                case 'T': break;
                case 'C': break;
                case '\0': goto after_cycle;
                default: error = 1;
              }
              if(error == 1) {
                print_error("Incorrect read genom for", header);
                goto next;
              }
            }
          }
after_cycle: result1 = apply_factorization(genom);
          fprintf(factorization_file, "%s\n%s\n", header, result1);
          result2 = create_fingerprint(result1, kfingerprint_file, header);
          fprintf(fingerprint_file, "%s %s\n", header, result2);
          fprintf(oneformat_file, "%s %c %s %c %s\n", header, '$', result2, '$', result1);
        }
      }
      fclose(kfingerprint_file);
      fclose(factorization_file);
      fclose(fingerprint_file);
      fclose(oneformat_file);
      printf("%s processed\n", path);
    }
    else {
      print_error("Impossible to open file ", file_description->d_name);
    }
  }
}

/*
  Thread that perform all the factorization and fingerprint calculus operations on a signle file
*/
void *process_file(void* arg) {
  return NULL;
}

/*
  It returns a pointer to a file in write mode with the name made up of two parts: fasta-name and name.
  param name: name to concatenate with the fasta file one that will be the name of the file to be opened
  param directory_path: directory path in which the file to be created will be
  pre-condition name: != NULL
  pre-condition fasta_name: must contain ".fasta" at the end of its name
  pre-condition directory_path: must be the path of an existing directory
  post-condition: an opened file with the descripted name will be returned
*/
FILE* open_towrite_file(char *name, char *fasta_name, char *directory_path) {

  FILE *file_to_open;
  char *filename;

  if (name == NULL || fasta_name == NULL) {
    printf("one of the names are NULL\n");
    return NULL;
  }

  if (strstr(fasta_name, ".fasta") == NULL) {
    printf("fasta name must contain .fasta suffix\n");
    return NULL;
  }

  filename = malloc(strlen(name) + strlen(fasta_name));
  strncpy(filename, fasta_name, strlen(fasta_name) - strlen(".fasta"));
  filename[strlen(fasta_name) - strlen(".fasta")] = '\0';
  strcat(filename, "-");
  strcat(filename, name);

  file_to_open = fopen(append_filename_to_path(directory_path, filename), "w");

  if (file_to_open == NULL) {
    printf("Non è stato possibile aprire in scrittura il file %s: %s\n\n", filename, strerror(errno));
  }

  return file_to_open;
}

/*
  pre-condition: fact_choice >= 1 && fact_choice <= 4
  return: factorized genom
*/
char *apply_factorization(char *genom) {

  node_t *factorized_genom;
  int second_parameter_value;
  char* result;
  int number_of_factors;

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

  /*needed to set the size of the string returned by list_to_string*/
  get_number_of_delimeters();
  number_of_factors = get_number_of_factors();
  set_number_of_elements(number_of_factors);
  set_number_of_elements(strlen(genom));
  /**/

  result = list_to_string(factorized_genom, second_parameter_value);

  number_of_factors = get_number_of_factors();
  set_number_of_elements(number_of_factors);

  return result;
}


/* It creates the fingerprint and k-fingerprint of the factorized genom. The k-fingerprints will be stored in
       the file pointed by kfingerprint_file variable.
  pre-condition factorized_genom: must be the result of factorize_read or format equivalent
  pre-condition file_to_write: must point a file opened in write mode
  pre-condition header: != NULL
  pre-condition file_to_write must be closed to save the results in the file permanently
  post-condition: kfingerprints will be written in the correct format
  return: fingerprint of the factorized genom
*/
char* create_fingerprint(char* factorized_genom, FILE *file_to_write, char *header) {

  int i = 2;
  int j = 0;
  int cont = 0;
  int dim = strlen(factorized_genom);
  char *fingerprint = malloc(dim);
  char converted_number[5];
  char fact_gen[strlen(factorized_genom)];
  tail *t = NULL;

  strcpy(fact_gen, factorized_genom);  //just to avoid external changes on the location pointed by factorized_genom. Necessary cause this happens
  factorized_genom = fact_gen;
  t = initialize_tail(t);
  while (1) {
    switch (factorized_genom[i]) {
      case '<':
        fingerprint[j++] = '-';
        fingerprint[j++] = '1';
        fingerprint[j++] = ',';
        i += 4;
        fill_k_fingerprint(-1, file_to_write, t, header);
        break;
      case '>':
        fingerprint[j++] = '0';
        fingerprint[j++] = ',';
        fill_k_fingerprint(0, file_to_write, t, header);
        i += 4;
        break;
      case ']':
        j--;
        fingerprint[j] = '\0';
        fill_k_fingerprint(-2, file_to_write, t, header);
        return fingerprint;
        break;
      case '\"':
        if (cont > 0) {   //if it defines the end of a factor
          sprintf(converted_number, "%d", cont);
          strcat(fingerprint, converted_number);
          i += 2;
          j += strlen(converted_number);
          fingerprint[j++] = ',';
          fill_k_fingerprint(cont, file_to_write, t, header);
          cont = 0;
        }
        else {
          i++;
        }
        break;
      default:
        cont++;   //cause it is one of the characters
        i++;
    }
    fingerprint[j] = '\0';
  }
}

/*It adds a new element to the fingerprint tail and in case it is full or it isn't and the fingerprints are ended,
      writes the next kfingerprint line.
  param fingerprint_number: if it is -2, all the fingerprints have been created and the tail must deal with that,
      and if it is -1 or 0, it means start and end of the second time factorized factors are found respectively.
      Other values are just normal fingerprints
  param *t: is referred to the rispective tail
  pre-condition fingerprint_number: >= -2
  pre-condition file_to_write: must be opened in write mode
  pre-condition *t: must be previousely initialized by initialize_tail function or simply set in a proper way
  pre-condition *header_to_pass: must not be null
  pre-condition: after calling this function, when result file is needed to be saved, fclose(file_to_write) must be called
  post_condition: if the got informations are enough, new kfingerprints are written in the file_to_write file
*/
void fill_k_fingerprint(int fingerprint_number, FILE *file_to_write, tail *t, char* header_to_pass) {

  if (fingerprint_number == -2) {
    if (!t->recived_exactly_k_fingers) {
      flush(file_to_write, t, header_to_pass);
    }
    if (t != NULL)
      t = initialize_tail(t);
    return;
  }

  if ( (fingerprint_number == -1) || (fingerprint_number == 0) ) {
    t->is_one = (fingerprint_number == -1);
  }
  else {
    if (t->has_not_been_filled_once) {
      t->end_window_limit++;
      pop_tail(fingerprint_number, t);
      if (t->end_window_limit == window_dimension - 1) {
        flush(file_to_write, t, header_to_pass);
        t->cont_shift += t->finger_tail[t->start_window_limit];
        t->recived_exactly_k_fingers = 1;
        t->end_window_limit = 0;
        t->start_window_limit = (t->start_window_limit + 1) % window_dimension;
        t->has_not_been_filled_once = 0;
      }
    }
    else {
      pop_tail(fingerprint_number, t);
      flush(file_to_write, t, header_to_pass);
      t->cont_shift += t->finger_tail[t->start_window_limit];
      t->start_window_limit = (t->start_window_limit + 1) % window_dimension;
      t->end_window_limit = (t->end_window_limit + 1) % window_dimension;
    }
  }
}

/* It adds a new fingerprint and the bool value about the fact it is the result of a second factorization
   at the head of the tails respectively in the position of end_window_limit variable (this one is not updated).
   The end_window_limit zero_one_tail element is set to the bool is_one that defines if the current element is
   the result of sub-factorization or not.
   pre-condiion: finger_tail and zero_one_tail must be synchronized (zero_one_tail[i] must regard finger_tail[i])
*/
void pop_tail(int fingerprint_number, tail *t) {

  t->finger_tail[t->end_window_limit] = fingerprint_number;
  t->zero_one_tail[t->end_window_limit] = t->is_one;
}

/* It records the tails content in the k-fingerprint file as the right format if you close the file_to_write after called this function.
   pre-condition file_to_write: must be opened in write mode;
   pre-condition t: must be previousely initialized and only manipulated by fill_k_fingerprint function or modified in proper way
   pre-condition header != NULL
   pre-condition: after calling this function fclose(file_to_write) must be called
*/
void flush(FILE *file_to_write, tail *t, char* header) {

  const int FINGERPRINT_MAX_DIMENSION = 3;

  int number_of_iterations = 0;
  int string_dimension = ((FINGERPRINT_MAX_DIMENSION + 1) * window_dimension) ;
  char zero_one_result[string_dimension], fingerprint_result[string_dimension];

  char converted_number[FINGERPRINT_MAX_DIMENSION + 1]; char converted_bit[3 + 1];  //considering the end character
  int current_length = 0;
  int current_length2 = 0;

  strcpy(zero_one_result, "");
  strcpy(fingerprint_result, "");


  for (int i = t->start_window_limit; number_of_iterations < window_dimension; i = (i + 1) % window_dimension) {
    number_of_iterations++;
    sprintf(converted_number, "%d ", t->finger_tail[i]);
    sprintf(converted_bit, "%d ", t->zero_one_tail[i]);
    current_length += strlen(converted_number);
    strcat(fingerprint_result, converted_number);
    fingerprint_result[current_length] = '\0';   //userful to prevent dirty characters forward program failure
    current_length2 = strlen(zero_one_result);
    strcat(zero_one_result, converted_bit);
    zero_one_result[current_length2 + 3] = '\0'; //userful to prevent dirty characters forward program failure
  }
  fingerprint_result[current_length - 1] = '\0';
  zero_one_result[strlen(zero_one_result) - 1] = '\0';

  int err;
  err = fprintf(file_to_write, "%s %c %s %c %s %d\n", fingerprint_result, '$', zero_one_result, '$', header, t->cont_shift);
  if (err < 0) {
    printf("impossible to write next k-fingerprints\n");
    printf("error: %s\n", strerror(errno));
    exit(1);
  }
}

/*
  if parameter t is NULL, t will point to a new tail. If t is not NULL, it will be initialized.
  pre-condition: window_dimension > 0
  pre-condition: if t is not NULL, zero_one_tail && finger_tail must have window_dimension as dimension
  post-condition: t points to a tail ready to be used
*/
tail *initialize_tail(tail *t) {

  int i = 0;

  if (window_dimension <= 0) {
    printf("window_dimension is not set correctly: %d\n", window_dimension);
    exit(1);
  }

  if (t == NULL) {
    t = (tail *) malloc(sizeof(tail));
    t->finger_tail = calloc(window_dimension, sizeof(int));
    t->zero_one_tail = calloc(window_dimension, sizeof(int));
  }

  t->cont_shift = 0;
  t->has_not_been_filled_once = 1;
  t->end_window_limit = -1;
  t->start_window_limit = 0;
  t->is_one = 0;
  t->recived_exactly_k_fingers = 0;

  if ( (t->finger_tail == NULL) || (t->zero_one_tail == NULL) ) {
    printf("error creating the tails. See initialite_tail function\n");
    exit(1);
  }
  for (i = 0; i < window_dimension; i++) {
    t->finger_tail[i] = t->zero_one_tail[i] = 0;
  }
  return t;
}

