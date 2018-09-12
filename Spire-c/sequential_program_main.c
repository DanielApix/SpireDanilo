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

struct stat st = {0};

int fact_choice;

char *root_path;     //...of the directory to process
int max_fact_length = 0; //arbitrary chosen and requested to the user

FILE *factorization_file, *fingerprint_file, *kfingerprint_file, *oneformat_file;

char *header_read;  //refers to the current read
char *genom_read;
int current_header_size = 0;
int current_genom_size = 0;

/*these refer to the k-fingerprint calculus. To better understand refer to the documentation.*/
int *zero_one_tail, *finger_tail;
int end_window_limit = -1;
int window_dimension;
int is_one = 0;   //boolean equivalent to determine if current finger refers to a second time factorized factor
int start_window_limit = 0;
int cont_shift = 0;  //number of shitfts from the start of the read dealing with the current window
int has_not_been_filled_once = 1;  //if it hasn't been filled at least once for the current read
int recived_exactly_k_fingers = 0;
/**/

char *inputString(FILE* fp, size_t size, char ending_character);
char* safe_fgets(char* buffer, int *current_size, FILE *stream);

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

  do_unit_testing();
  do_integration_testing();
  /*Catching input from the user*/

  printf("Benvenuto nel programma sequenziale Spire\n\n");

  printf("Fornisca la directory dei file fasta\n");
  root_path = inputString(stdin, 100, '\n');

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
/*
  fact_choice = 3;
  strcpy(root_path, "/home/danilo/Scrivania/example-c");
  window_dimension = 4;
  max_fact_length = 30;
  communicate_max_fact_length(max_fact_length);
*/
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

/*Used as scanf handling in the case in which the defined dimension is not enough.
  pre-condition fp: == stdio
  pre-condition size: > 0
  pre-condition ending_character: != '\0'
  return: input read from the start to the ending_character
*/
char *inputString(FILE* fp, size_t size, char ending_character){
//The size is extended by the input with the value of the provisional
    char *str;
    int ch;
    size_t len = 0;
    str = realloc(NULL, sizeof(char)*size);//size is start size
    if(!str)return str;
    while(EOF!=(ch=fgetc(fp)) && ch != '\n'){
        str[len++]=ch;
        if(len==size){
            str = realloc(str, sizeof(char)*(size+=16));
            if(!str)return str;
        }
    }
    str[len++]='\0';

    return realloc(str, sizeof(char)*len);
}

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
  int i;

  tmp = NULL;

  if (*current_size == 0) {
    buffer = malloc(STARTING_DIMENSION * sizeof(char));
    *current_size = STARTING_DIMENSION;
  }
  if (fgets(buffer, *current_size, stream) == NULL) {
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


/*If the file descripted is a fasta file, factorization, fingerprint, kfingerprint and the one format of the first three
  will be created and saved in the same directory containing the file to be processed
  pre-condition file_description: must descript an existing file
  pre-condition path: must refer to the descripted file
  pre-condition: current_header_size and current_genom_size must be 0
  post-condition: given the name nam of the fasta file without ".fasta" at the end, nam-factorization.txt, nam-fingerprint.txt,
      nam-kfingerprint, nam-oneformat.txt will be created in the same directory of the fasta file with the respective output inside*/
void process_fasta(struct dirent *file_description, char *path) {

  char directory_path[100];
  FILE *fasta_file;
  char *result1;
  char *result2;

  if (strlen(file_description->d_name) > strlen(".fasta")) {
    if (strstr(file_description->d_name, ".fasta") != NULL) {   //if it has .fasta extention
      directory_path[strlen(path) - strlen(file_description->d_name) - 1] = '\0';
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
        while ( (header_read = safe_fgets(header_read, &current_header_size, fasta_file)) != NULL ) {
          genom_read = safe_fgets(genom_read, &current_genom_size, fasta_file);

          result1 = apply_factorization(genom_read);
          fprintf(factorization_file, "%s\n%s\n", header_read, result1);
          result2 = create_fingerprint(result1);
          fprintf(fingerprint_file, "%s %s\n", header_read, result2);
          fprintf(oneformat_file, "%s %c %s %c %s\n", header_read, '$', result2, '$', result1);
        }
      }
  //    if (factorization_file == NULL)
  //      printf("factorization_file lost\n");
      fclose(factorization_file);
      fclose(fingerprint_file);
      fclose(oneformat_file);
      current_header_size = 0;
      current_genom_size = 0;
    }
    else {
      printf("Non è stato possibile aprire il file fasta %s\nErrore: %s\n", file_description->d_name, strerror(errno));
    }
  }
}

/*
  It creates the file in the specified directory and opens it into the respective variable in write mode.
  param name: name to concatenate with the fasta file one that will be the name of the file to be opened
  param directory_path: directory path in which the file to be created will be
  pre-condition name: name == "factorization" || name == "fingerprint" || name == "kfingerprint" || name == "oneformat"
  pre-condition fasta_name: must contain ".fasta" at the end of its name
  pre-condition directory_path: must be the path of an existing directory
  post-condition one of the variables factorization_file, fingerprint_file, kfingerprint_file, oneformat_file will contain
      a pointer to an opened file in write mode with the name descripted below depending on the name parameter value respectively
*/
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

/*
  pre-condition: fact_choice >= 1 && fact_choice <= 4
  return: factorized genom
*/
char *apply_factorization(char *genom) {

  node_t *factorized_genom;
  int second_parameter_value;

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

/* It creates the fingerprint and k-fingerprint of the factorized genom. The k-fingerprints will be stored in
       the file pointed by kfingerprint_file variable.
  pre-condition factorized_genom: must be the result of factorize_read or format equivalent
  pre-condition: kfingerprint_file must point to an existing opened file in writing mode
  post-condition: kfingerprints will be written in the correct format in the file pointed by the kfingerprint_file variable
  return: fingerprint of the factorized genom
*/
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
        fill_k_fingerprint(-1);
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
        break;
      case '\"':
        if (cont > 0) {   //if it defines the end of a factor
          sprintf(converted_number, "%d", cont);
          strcat(fingerprint, converted_number);
          i += 2;
          j += strlen(converted_number);
          fingerprint[j++] = ',';
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

/*It adds a new element to the fingerprint tail and in case it is full or it isn't and the fingerprints are ended,
      writes the next kfingerprint line.
  param fingerprint_number: if it is -2, all the fingerprints have been created and the tail must deal with that,
      and if it is -1 or 0, it means start and end of the second time factorized factors are found respectively.
      Other values are just normal fingerprints
  pre-condition fingerprint_number: >= -2
  pre-condition: all the values must be initialized before the start of each series of fingerprint to be scanned
  post_condition: if the got informations are enough, new kfingerprints are written in the k-fingerprint file
*/
void fill_k_fingerprint(int fingerprint_number) {

  if (fingerprint_number == -2) {
    if (!recived_exactly_k_fingers)
      flush();
    initialize_tail();
    return;
  }

  if ( (fingerprint_number == -1) || (fingerprint_number == 0) ) {
    is_one = (fingerprint_number == -1);
  }
  else {
    if (has_not_been_filled_once) {
      end_window_limit++;
      pop_tail(fingerprint_number);
      if (end_window_limit == window_dimension - 1) {
        flush();
        recived_exactly_k_fingers = 1;
        end_window_limit = 0;
        start_window_limit = (start_window_limit + 1) % window_dimension;
        has_not_been_filled_once = 0;
      }
    }
    else {
      cont_shift += finger_tail[start_window_limit];
      pop_tail(fingerprint_number);
      flush();
      start_window_limit = (start_window_limit + 1) % window_dimension;
      end_window_limit = (end_window_limit + 1) % window_dimension;
    }
  }
}

/* It adds a new fingerprint and the bool value about the fact it is the result of a second factorization
   at the head of the tails respectively in the position of end_window_limit variable (this one is not updated).
   The end_window_limit zero_one_tail element is set to the bool is_one that defines if the current element is
   the result of sub-factorization or not.
   pre-condiion: finger_tail and zero_one_tail must be synchronized (zero_one_tail[i] must regard finger_tail[i])
*/
void pop_tail(int fingerprint_number) {

  finger_tail[end_window_limit] = fingerprint_number;
  zero_one_tail[end_window_limit] = is_one;
}

/* It records the tails content in the k-fingerprint file as the right format if you close the kfingerprint_file after called this function.
   pre-condition: window_dimension must match the distance between start_window_limit and end_window_limit in modulo arithmetic
   pre-condtion: zero_one_tail and finger_tail must be integer arrays
   pre-condition: zero_one_tail and finger_tail must have window_dimension as dimension
   pre-condition: kfingerprint_file must point to an opened file in write mode
   pre-condition: header_read != NULL and must be a string
   pre-condition: after calling this function fclose(kfingerprint_file) must be called
*/
void flush() {

  const int FINGERPRINT_MAX_DIMENSION = 3;

  int number_of_iterations = 0;
  /*n window element * FINGERPRINT_MAX_DIMENSION + number of spaces(n - 1).
    result = FINGERPRINT_MAX_DIMENSION * n + n  - 1 = (FINGERPRINT_MAX_DIMENSION + 1)n - 1. Considering
    than the end string character ('\0'), we have result + 1 =
    = (FINGERPRINT_MAX_DIMENSION + 1) * n - 1 + 1 = (FINGERPRINT_MAX_DIMENSION + 1) * n as resulting dimension*/
  int string_dimension = ((FINGERPRINT_MAX_DIMENSION + 1) * window_dimension) ;
  char zero_one_result[string_dimension], fingerprint_result[string_dimension];

  char converted_number[FINGERPRINT_MAX_DIMENSION + 1]; char converted_bit[3 + 1];  //considering the end character
  int current_length = 0;
  int current_length2 = 0;

  strcpy(zero_one_result, "");
  strcpy(fingerprint_result, "");

  for (int i = start_window_limit; number_of_iterations < window_dimension; i = (i + 1) % window_dimension) {
    number_of_iterations++;
    sprintf(converted_number, "%d ", finger_tail[i]);
    sprintf(converted_bit, "%d ", zero_one_tail[i]);
    current_length += strlen(converted_number);
    strcat(fingerprint_result, converted_number);
    fingerprint_result[current_length] = '\0';   //userful to prevent dirty characters forward program failure
    current_length2 = strlen(zero_one_result);
    strcat(zero_one_result, converted_bit);
    zero_one_result[current_length2 + 3] = '\0'; //userful to prevent dirty characters forward program failure
  }
  fingerprint_result[current_length - 1] = '\0';
  zero_one_result[strlen(zero_one_result) - 1] = '\0';
  fprintf(kfingerprint_file, "%s %c %s %c %s %d\n", fingerprint_result, '$', zero_one_result, '$', header_read, cont_shift);
}

/*
  pre-condition: variables to be initialized must mach the compatible types
  post-condition: tails are set to a corret state configuration
*/
void initialize_tail() {

  int i = 0;

  cont_shift = 0;
  has_not_been_filled_once = 1;
  end_window_limit = -1;
  start_window_limit = 0;
  is_one = 0;
  recived_exactly_k_fingers = 0;

  for (i = 0; i < window_dimension; i++) {
    finger_tail[i] = zero_one_tail[i] = 0;
  }
}

/*Sets the number of elements for each window and sets the tails
  pre-condition: window_dimension > 0
*/
void initialize_k_finger() {
  zero_one_tail = calloc(window_dimension, sizeof(int));
  finger_tail = calloc(window_dimension, sizeof(int));
}

void do_unit_testing() {
  printf("start of unit test\n\n");
  test_safe_fgets();
  test_foreach_element_in();
  test_initialize_k_finger();
  test_flush();
  test_open_towrite_file();
  printf("unit test successfully passed\n\n");
}

void do_integration_testing() {
  printf("start of integration test\n");
  test_create_fingerprint();
  printf("integration test successfully completed\n");
  printf("program works as it should\n");
}

/*Reads an entire stream line until it reaches the \n character and records it in a pointer. If there is no
  line to be read or error occurs, NULL is returned.
  if current_size is zero, that is a buffer is considered as it hadn't been previousely allocated, another one will be.
  NOTE: buffer doesn't store the result at the end of the processing.
  pre-condition buffer: buffer must point to a previousely dinamically allocated memory or to NULL
  pre-condition current_size: must point to a variable whose value is the dimension of buffer (== 0 if buffer is NULL)
  pre-condition stream: must point to an existing file
  return: next string read up to the next '\n' reached or NULL if reading fails

char* safe_fgets(char* buffer, int *current_size, FILE *stream) */

void test_safe_fgets() {

  FILE *test_file;
  char *buffer;
  int current_size;
  char *t;

  printf("\n\nstart of test_safe_fgets test\n");

  test_file = fopen("test_file.txt", "w");
  if (test_file == NULL) {
    printf("test_safe_fgets: test couldn't be executed cause test file cannot be opened in writing mode\n");
    printf("error: %s\n", strerror(errno));
    return;
  }

  fprintf(test_file, "%s\n%s\n%s", "this is the first sentence", "this is the second one", "this is the third");
  fclose(test_file);

  test_file = fopen("test_file.txt", "r");

  if (test_file == NULL) {
    printf("test_safe_fgets: test couldn't be executed cause test file cannot be opened in reading mode\n");
    printf("error: %s\n", strerror(errno));
    return;
  }

  //printf("test conditions: buffer dimension and variable current_size are unsynchronized\n");

  /*buffer = NULL;
  current_size = 3;

  printf("test case: buffer = NULL and current_dimension = 3\n");
  buffer = safe_fgets(buffer, &current_size, test_file);
  assert(strcmp(buffer, "this is the first sentence") == 0);
  printf("test case passed\n");
*/
/*
  printf("buffer dimension is not current_size\n");
  printf("buffer dimension is 3 and current_size is 5\n");

  buffer = (char *) malloc(3 * sizeof(char));
  current_size = 5;

  strcpy(buffer, "ads");
  buffer = safe_fgets(buffer, &current_size, test_file);
  assert(strcmp(buffer, "this is the first sentence") == 0);

*/
  printf("test conditions: buffer dimension and variable current_size are synchronized\n");
  printf("test case: buffer is NULL and current_size is 0\n");

  buffer = NULL;
  current_size = 0;

  buffer = safe_fgets(buffer, &current_size, test_file);
  assert(strcmp(buffer, "this is the first sentence") == 0);

  printf("test case passed\n");

  printf("test case: buffer is of dimension current_size = 26 and second line is going to be read\n");
  buffer = safe_fgets(buffer, &current_size, test_file);
  assert(strcmp(buffer, "this is the second one") == 0);

  printf("test case passed\n");

  free(buffer);
  fclose(test_file);

  test_file = fopen("test_file.txt", "r");
  if (test_file == NULL) {
    printf("test_safe_fgets: test couldn't be executed cause test file cannot be opened in writing mode\n");
    printf("error: %s\n", strerror(errno));
    return;
  }

  printf("test conditions: dimension of buffer < dimension of the string to be read\n");
  printf("test case: buffer is dimension current_size = 5 and the sentence to be read si of dimension 26\n");

  buffer = (char *) malloc(5 * sizeof(char));
  current_size = 5;
  buffer = safe_fgets(buffer, &current_size, test_file);
  assert(strcmp(buffer, "this is the first sentence") == 0);

  printf("test case passed\n");

  free(buffer);

  printf("safe_fgets test passed\n");
  fclose(test_file);
  remove("test_file.txt");
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

void for_each_element_in(char* directory_path,  void (*apply_function) (struct dirent *, char *)) {
*/


char path_test[255];  //watch out: path can't be more the 100 long or test fails and not because of program to be tested
int directories_found_test[4] = {0};

void print_names (struct dirent * element, char *path) {

  char* name = element->d_name;
  char path_to_compare[255];
  int b;

  b = ((strcmp(name, ".") == 0) || (strcmp(name, "..") == 0) || (strcmp(name, "directory1") == 0) || (strcmp(name, "directory2") == 0));

  strcpy(path_to_compare, path_test);
  strcat(path_to_compare, "/");
  strcat(path_to_compare, name);
  assert(strcmp(path, path_to_compare) == 0);
  assert(b);

  if (!directories_found_test[0])
    directories_found_test[0] =  (strcmp(name, ".") == 0);
  if (!directories_found_test[1])
    directories_found_test[1] =  (strcmp(name, "..") == 0);
  if (!directories_found_test[2])
    directories_found_test[2] =  (strcmp(name, "directory1") == 0);
  if (!directories_found_test[3])
    directories_found_test[3] =  (strcmp(name, "directory2") == 0);
}

void test_foreach_element_in() {

  DIR *file;
  struct dirent *inner_file;

  printf("\n\nStart of for_each_element_in test\n");
  if (stat("directory-test", &st) != -1) {
    system("rm -r directory-test");
  }

  mkdir("directory-test", 0700);
  mkdir("directory-test/directory1", 0700);
  mkdir("directory-test/directory2", 0700);

  if (getcwd(path_test, 255) == NULL) {
    printf("for_each_element_in test cannot be run cause it hasn't been possible to find the current path\n");
    printf("error: %s\n", strerror(errno));
    return;
  }

  printf("test conditions: scanning of 4 created subdirectories");

  strcat(path_test, "/");
  strcat(path_test, "directory-test");

  file = opendir("/home/danilo/Scrivania/Spire-c/directory-test/directory1");
  if (file == NULL) {
    printf("for_each_element test cannot be completed cause created directory1 cannot be open\n");
    return;
  }
  inner_file = readdir(file);
  if (inner_file == NULL) {
    printf("for_each_element test cannot be completed cause something gone wrong in catching a file\n");
    return;
  }
  //print_names(inner_file, path_test);

  for_each_element_in(path_test, print_names);

  assert(directories_found_test[0] != 0);
  assert(directories_found_test[1] != 0);
  assert(directories_found_test[2] != 0);
  assert(directories_found_test[3] != 0);


  system("rm -r directory-test");

  printf("test case passed\n");
  printf("for_each_element_in test passed\n");
}


/*Sets the number of elements for each window and sets the tails
  pre-condition: window_dimension > 0

void initialize_k_finger() */

void test_initialize_k_finger() {

  int i;

  printf("\n\nStart of initialize_k_finger test\n");

  printf("test conditions: window dimension = 4\n");

  window_dimension = 4;
  initialize_k_finger();

  for (i = 0; i < window_dimension; i++) {  //it stops the program if pre-condition are unsatisfied
    zero_one_tail[i] = 0;
    finger_tail[i] = 0;
  }

  free(zero_one_tail);
  free(finger_tail);

  printf("initialize_k_finger test passed\n");
}

/* It records the tails content in the k-fingerprint file as the right format if you close the kfingerprint_file after called this function.
   pre-condition: window_dimension must match the distance between start_window_limit and end_window_limit in modulo arithmetic
   pre-condtion: zero_one_tail and finger_tail must be integer arrays
   pre-condition: zero_one_tail and finger_tail must have window_dimension as dimension
   pre-condition: kfingerprint_file must point to an opened file in write mode
   pre-condition: header_read != NULL and must be a string
   pre-condition: after calling this function fclose(kfingerprint_file) must be called

void flush()
*/

void test_flush() {

  FILE *test_file;
  char s[300], c1, c2, header[100];
  int f1, f2, f3, f4, z1, z2, z3,z4, cont;

  printf("\n\nStart of flush test\n");
  printf("test conditions: window_dimension = 5, start_window_limit = 0, end_window_limit = 4\n");

  window_dimension = 4;
  zero_one_tail = calloc(window_dimension, sizeof(int));
  finger_tail = calloc(window_dimension, sizeof(int));

  start_window_limit = 0;
  end_window_limit = 3;

  kfingerprint_file = fopen("test_file.txt", "w");
  if (kfingerprint_file == NULL) {
    printf("test cannot be completed cause has not been possible to open kfingerprint_file\n");
    return;
  }
  header_read = (char *) malloc(100 * sizeof(char));
  strcpy(header_read, "header");

  zero_one_tail[0] = 1;
  zero_one_tail[1] = 1;
  zero_one_tail[2] = 0;
  zero_one_tail[3] = 0;

  finger_tail[0] = 2;
  finger_tail[1] = 3;
  finger_tail[2] = 4;
  finger_tail[3] = 5;

  flush();

  fclose(kfingerprint_file);
  test_file = fopen("test_file.txt", "r");
  if (test_file == NULL) {
    printf("test cannot be completed cause has not been possible to open test_file\n");
    return;
  }
  if(fgets(s, 300, test_file) == NULL)
    assert(0);

  sscanf(s, "%d %d %d %d %c %d %d %d %d %c %s %d", &f1, &f2, &f3, &f4, &c1, &z1, &z2, &z3, &z4, &c2, header, &cont);

  assert(f1 == 2);
  assert(f2 == 3);
  assert(f3 == 4);
  assert(f4 == 5);
  assert(c1 == '$');
  assert(z1 == 1);
  assert(z2 == 1);
  assert(z3 == 0);
  assert(z4 == 0);
  assert(c2 == '$');
  assert(strcmp(header, "header") == 0);

  fclose(test_file);
  system("rm -r test_file.txt");
  free(zero_one_tail);
  free(finger_tail);
  free(header_read);

  printf("test case passed\n");

  printf("test conditions: window_dimension = 4, start_window_limit = 0, end_window_limit = 4\n");

  window_dimension = 4;
  zero_one_tail = calloc(window_dimension, sizeof(int));
  finger_tail = calloc(window_dimension, sizeof(int));

  start_window_limit = 0;
  end_window_limit = 3;

  kfingerprint_file = fopen("test_file.txt", "w");
  if (kfingerprint_file == NULL) {
    printf("test cannot be completed cause has not been possible to open kfingerprint_file\n");
    return;
  }
  header_read = (char *) malloc(100 * sizeof(char));
  strcpy(header_read, "header");

  zero_one_tail[0] = 1;
  zero_one_tail[1] = 1;
  zero_one_tail[2] = 0;
  zero_one_tail[3] = 0;

  finger_tail[0] = 2;
  finger_tail[1] = 3;
  finger_tail[2] = 4;
  finger_tail[3] = 5;

  flush();

  fclose(kfingerprint_file);
  test_file = fopen("test_file.txt", "r");
  if (test_file == NULL) {
    printf("test cannot be completed cause has not been possible to open test_file\n");
    return;
  }
  if(fgets(s, 300, test_file) == NULL)
    assert(0);

  printf("%s\n", s);
  sscanf(s, "%d %d %d %d %c %d %d %d %d %c %s %d", &f1, &f2, &f3, &f4, &c1, &z1, &z2, &z3, &z4, &c2, header, &cont);

  assert(f1 == 2);
  assert(f2 == 3);
  assert(f3 == 4);
  assert(f4 == 5);
  assert(c1 == '$');
  assert(z1 == 1);
  assert(z2 == 1);
  assert(z3 == 0);
  assert(z4 == 0);
  assert(c2 == '$');
  assert(strcmp(header, "header") == 0);

  fclose(test_file);
  system("rm -r test_file.txt");
  free(zero_one_tail);
  free(finger_tail);
  free(header_read);

  printf("test case passed\n");

  printf("test conditions: window_dimension = 4, start_window_limit = 3, end_window_limit = 2\n");

  window_dimension = 4;
  zero_one_tail = calloc(window_dimension, sizeof(int));
  finger_tail = calloc(window_dimension, sizeof(int));

  kfingerprint_file = fopen("test_file.txt", "w");
  if (kfingerprint_file == NULL) {
    printf("test cannot be completed cause has not been possible to open kfingerprint_file\n");
    return;
  }
  header_read = (char *) malloc(100 * sizeof(char));
  strcpy(header_read, "header");

  zero_one_tail[3] = 1;
  zero_one_tail[0] = 1;
  zero_one_tail[1] = 0;
  zero_one_tail[2] = 0;

  finger_tail[3] = 2;
  finger_tail[0] = 3;
  finger_tail[1] = 4;
  finger_tail[2] = 5;

  start_window_limit = 3;
  end_window_limit = 2;

  flush();

  fclose(kfingerprint_file);
  test_file = fopen("test_file.txt", "r");
  if (test_file == NULL) {
    printf("test cannot be completed cause has not been possible to open test_file\n");
    return;
  }
  if(fgets(s, 300, test_file) == NULL)
    assert(0);

  printf("%s\n", s);
  sscanf(s, "%d %d %d %d %c %d %d %d %d %c %s %d", &f1, &f2, &f3, &f4, &c1, &z1, &z2, &z3, &z4, &c2, header, &cont);

  assert(f1 == 2);
  assert(f2 == 3);
  assert(f3 == 4);
  assert(f4 == 5);
  assert(c1 == '$');
  assert(z1 == 1);
  assert(z2 == 1);
  assert(z3 == 0);
  assert(z4 == 0);
  assert(c2 == '$');
  assert(strcmp(header, "header") == 0);

  fclose(test_file);
  system("rm -r test_file.txt");
  free(zero_one_tail);
  free(finger_tail);
  free(header_read);

  printf("test case passed\n");

  printf("flush test passed\n");
}

/*
  It creates the file in the specified directory and opens it into the respective variable in write mode.
  param name: name to concatenate with the fasta file one that will be the name of the file to be opened
  param directory_path: directory path in which the file to be created will be
  pre-condition name: name == "factorization" || name == "fingerprint" || name == "kfingerprint" || name == "oneformat"
  pre-condition fasta_name: must contain ".fasta" at the end of its name
  pre-condition directory_path: must be the path of an existing directory
  post-condition one of the variables factorization_file, fingerprint_file, kfingerprint_file, oneformat_file will contain
      a pointer to an opened file in write mode with the name descripted below depending on the name parameter value respectively

void open_towrite_file(char *name, char *fasta_name, char *directory_path)
*/

void test_open_towrite_file() {

  FILE *test_file;
  char path_test[300];

  printf("\n\nStart of open_towrite_file test\n");
  printf("test conditions: fasta_name = \".fasta\" for factorization, fingerprint, kfingerprint and oneformat\n");

  if (getcwd(path_test, 255) == NULL) {
    printf("for_each_element_in test cannot be run cause it hasn't been possible to find the current path\n");
    printf("error: %s\n", strerror(errno));
    return;
  }

  open_towrite_file("factorization", ".fasta", path_test);
  open_towrite_file("fingerprint", ".fasta", path_test);
  open_towrite_file("kfingerprint", ".fasta", path_test);
  open_towrite_file("oneformat", ".fasta", path_test);

  if (factorization_file == NULL)
    assert(0);
  if (fingerprint_file == NULL)
    assert(0);
  if (kfingerprint_file == NULL)
    assert(0);
  if (oneformat_file == NULL)
    assert(0);

  fclose(factorization_file);
  fclose(fingerprint_file);
  fclose(kfingerprint_file);
  fclose(oneformat_file);

  test_file = fopen("-factorization", "r");
  if (test_file == NULL) {
    printf("-factorization.txt not created\n");
    assert(0);
  }

  test_file = fopen("-fingerprint", "r");
  if (test_file == NULL) {
    printf("-fingerprint.txt not created\n");
    assert(0);
  }

  test_file = fopen("-kfingerprint", "r");
  if (test_file == NULL) {
    printf("-kfingerprint.txt not created\n");
    assert(0);
  }

  test_file = fopen("-oneformat", "r");
  if (test_file == NULL) {
    printf("-oneformat.txt not created\n");
    assert(0);
  }

  remove("-factorization");
  remove("-fingerprint");
  remove("-kfingerprint");
  remove("-oneformat");

  printf("test case passed\n");

  printf("test conditions: fasta_name = \"normal_name.fasta\" for factorization, fingerprint, kfingerprint and oneformat\n");

  if (getcwd(path_test, 255) == NULL) {
    printf("for_each_element_in test cannot be run cause it hasn't been possible to find the current path\n");
    printf("error: %s\n", strerror(errno));
    return;
  }

  open_towrite_file("factorization", "normal_name.fasta", path_test);
  open_towrite_file("fingerprint", "normal_name.fasta", path_test);
  open_towrite_file("kfingerprint", "normal_name.fasta", path_test);
  open_towrite_file("oneformat", "normal_name.fasta", path_test);

  if (factorization_file == NULL)
    assert(0);
  if (fingerprint_file == NULL)
    assert(0);
  if (kfingerprint_file == NULL)
    assert(0);
  if (oneformat_file == NULL)
    assert(0);

  fclose(factorization_file);
  fclose(fingerprint_file);
  fclose(kfingerprint_file);
  fclose(oneformat_file);

  test_file = fopen("normal_name-factorization", "r");

  test_file = fopen("normal_name-fingerprint", "r");
  if (test_file == NULL) {
    printf("-fingerprint.txt not created\n");
    assert(0);
  }

  test_file = fopen("normal_name-kfingerprint", "r");
  if (test_file == NULL) {
    printf("-kfingerprint.txt not created\n");
    assert(0);
  }

  test_file = fopen("normal_name-oneformat", "r");
  if (test_file == NULL) {
    printf("-oneformat.txt not created\n");
    assert(0);
  }

  remove("normal_name-factorization");
  remove("normal_name-fingerprint");
  remove("normal_name-kfingerprint");
  remove("normal_name-oneformat");

  printf("test case passed\n");

  printf("open_towrite_file test passed\n");


  /*file name of length next to 255 is not considered here cause of time constraints*/
}

/* It creates the fingerprint and k-fingerprint of the factorized genom. The k-fingerprints will be stored in
       the file pointed by kfingerprint_file variable.
  pre-condition factorized_genom: must be the result of factorize_read or format equivalent
  pre-condition: kfingerprint_file must point to an existing opened file in writing mode
  post-condition: kfingerprints will be written in the correct format in the file pointed by the kfingerprint_file variable
  return: fingerprint of the factorized genom

char* create_fingerprint(char* factorized_genom)
*/

void test_create_fingerprint() {

}
