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

#define PROCESSORS_NUMBER 3

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

typedef struct {
  char directory_path[100];
  char filename[255];
  int deleted;
} filesystem_node;

struct stat st = {0};

pthread_t thread_ids[PROCESSORS_NUMBER];
filesystem_node files_to_process[PROCESSORS_NUMBER];
int available_threads[PROCESSORS_NUMBER];
int end_of_processing = 0;  //used to comunicate all the threads to finish
int fact_choice;

char *root_path;     //...of the directory to process
int max_fact_length = 0; //arbitrary chosen and requested to the user


/*char *header_read;  //refers to the current read
char *genom_read;
int current_header_size = 0;
int current_genom_size = 0;
*/
int window_dimension;

time_t time_spent_to_read_file = 0;
time_t time_spent_to_write_in_file = 0;

int test_finish = PROCESSORS_NUMBER;
int test_activate_process_file_log = 0;

/*necessary to set the dimension of the string returned by list_to_string (efficiency reasons)*/
int get_number_of_factors();
int get_number_of_delimeters();
void set_number_of_elements(int num);
void set_read_dimension(int value);
void communicate_max_fact_length(int c);
void print_statistics();

char *inputString(FILE* fp, size_t size, char ending_character);
char* safe_fgets(char* buffer, int *current_size, FILE *stream);

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
void *process_file(void* arg);
void process_fasta(struct dirent *file_description, char *directory_path);
void process_all_fasta_files(struct dirent *subdirectory_description, char* current_path);

void test_safe_fgets();
void test_foreach_element_in();
void test_flush();
void test_open_towrite_file();
void do_unit_testing();

void  test_fill_k_finger();
void test_create_fingerprint();
void test_process_file();
void test_process_fasta();
void do_integration_testing();

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

  //do_unit_testing();
  //do_integration_testing();
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

  communicate_max_fact_length(max_fact_length);
  printf("fornisca il numero di elementi per ciascuna finestra per le k-fingerprint\n");
  scanf("%d", &window_dimension);

  time_t m;
  time_t now = time(NULL);

  /*thread initialization*/
  int i, *j;
  pthread_attr_t *attr;

  end_of_processing = 0;
  for(i = 0; i < PROCESSORS_NUMBER; i++) {
    available_threads[i] = 1;
    files_to_process[i].deleted = 1;
    j = malloc(sizeof(int));
    *j = i;
    attr = malloc(sizeof(pthread_attr_t));
    pthread_attr_init(attr);
    pthread_create(&thread_ids[i], attr, process_file, j);
  }
  /**/

  for_each_element_in(root_path, process_all_fasta_files);
  end_of_processing = 1;

  int threads_have_done;
  while (!threads_have_done) {
    threads_have_done = 1;
    for (i = 0; i < PROCESSORS_NUMBER; i++)
      threads_have_done = threads_have_done && available_threads[i];  //if at least one thread hasn't finished, continue waiting
  }

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
  PROCESSORS_NUMBER > 0
  post-condition: all fasta files in the specified directory will be processed
*/
void process_all_fasta_files(struct dirent *subdirectory_description, char* current_path) {

  if ((strcmp(subdirectory_description->d_name, ".") != 0) && (strcmp(subdirectory_description->d_name, "..") != 0)) {
      for_each_element_in(current_path, process_fasta);
  }
}


/*If the file descripted is a fasta file, factorization, fingerprint, kfingerprint and the one format of the first three
  will be created and saved in the same directory containing the file to be processed
  pre-condition file_description: must descript an existing file with reads that respect the correct format
  pre-condition path: must refer to the descripted file
  pre-condition: fact_choice >= 1 || fact_choice <= 4
  pre-condition: PROCESSORS_NUMBER >= 0;
  pre-condition: PROCESSORS_NUMBER threads must be created before calling this function initialializing available_threads, thread_ids and files_to_process (.deleted set to 1)

  post-condition: given the name nam of the fasta file without ".fasta" at the end, nam-factorization, nam-fingerprint,
      nam-kfingerprint, nam-oneformat will be created in the same directory of the fasta file with the respective output inside
*/
void process_fasta(struct dirent *file_description, char *path) {

  int i;
  char fasta_path[strlen(path)];

  if (strlen(file_description->d_name) > strlen(".fasta")) {
    if (strstr(file_description->d_name, ".fasta") != NULL) {   //if it has .fasta extention
      i = 0;
      while (1) {
        if(available_threads[i]) {
          available_threads[i] = 0;
          break;
        }
        i = (i + 1) % PROCESSORS_NUMBER;
      }
      strncpy(files_to_process[i].directory_path, path, strlen(path) - strlen(file_description->d_name) - 1);
      strcpy(files_to_process[i].filename, file_description->d_name);  //useful cause sometimes file_description->d_name strangely changes in opening files
      strcpy(fasta_path, path);
      files_to_process[i].deleted = 0;
    }
  }
  test_finish--;
}

int running = 0;

/*
  It is the thread function that processes the file given by the array files_to_process at index given by arg, setting
  available_threads at same index to true. If no file is given, that is if the filesystem_node of files_to_process at the same
  index has variable deleted to true, it stays in wait mode until variable end_of_processing is set to true.
  param: integer that is the index of the interested informations in thread_ids, files_to_process, available_threads
  pre-condition arg: must be an integer between 0 and NUMBER_PROCESSORS
*/
void *process_file(void* arg) {

  FILE *fasta_file;
  char *result1;
  char *result2;
  char *header;
  char *genom;
  int current_header_size = 0;
  int current_genom_size = 0;
  /*used to not loose the allocated memory in case NULL is returned by safe_fgets cause EOF has been reached*/
//  char *header_read_beckup;
//  int current_header_size_beckup;
  /**/
  FILE *factorization_file;
  FILE *fingerprint_file;
  FILE *kfingerprint_file;
  FILE *oneformat_file;
  char *directory_path;
  char *filename;
  const int *identity_in_arrays = arg;
  while(1) {
  if (!files_to_process[*identity_in_arrays].deleted) {   //if there is a file exclusively to process by you
    if (test_activate_process_file_log) {
      printf("running state\n");
    }
    running = 1;
        available_threads[*identity_in_arrays] = 0;
        filesystem_node *file = &files_to_process[*identity_in_arrays];
        directory_path = file->directory_path;
        filename = file->filename;
        char path[strlen(directory_path) + 1 + strlen(filename)];
        strcpy(path, directory_path);
        strcat(path, "/");
        strcat(path, filename);

        factorization_file = open_towrite_file("factorization", filename, directory_path);
        fingerprint_file = open_towrite_file("fingerprint", filename, directory_path);
        kfingerprint_file = open_towrite_file("kfingerprint", filename, directory_path);
        oneformat_file = open_towrite_file("oneformat", filename, directory_path);
        if ((fasta_file = fopen(path, "r")) == NULL) {
         printf("error: %s\n\n", strerror(errno));
        }
       /*Processing of fasta file*/
        if (fasta_file != NULL) {
          while ( (header = safe_fgets(header, &current_header_size, fasta_file)) != NULL ) {
            genom = safe_fgets(genom, &current_genom_size, fasta_file);
            result1 = apply_factorization(genom);
            fprintf(factorization_file, "%s\n%s\n", header, result1);
            result2 = create_fingerprint(result1, kfingerprint_file, header);
            fprintf(fingerprint_file, "%s %s\n", header, result2);
            fprintf(oneformat_file, "%s %c %s %c %s\n", header, '$', result2, '$', result1); 
            /*  header_read_beck = header;
            current_header_size_beckup = current_header_size; */
          }
          /*   header_read = header_read_beckup;
          current_header_size = current_header_size_beckup;  */
          fclose(kfingerprint_file);
          fclose(factorization_file);
          fclose(fingerprint_file);
          fclose(oneformat_file);
          printf("%s processed\n", path);
          test_finish--;  //used only for test
        }
        else {
          printf("Non è stato possibile aprire il file fasta %s\nErrore: %s\n", path, strerror(errno));
        }
        available_threads[*identity_in_arrays] = 1;
        files_to_process[*identity_in_arrays].deleted = 1;
      }
    else {   
      if(end_of_processing) {  //if you have not to wait and there is nothing else to process
         pthread_exit(0);
        if (test_activate_process_file_log) {
          printf("exit state\n");
        }
      }
      else {
        if (test_activate_process_file_log) {
          printf("wait state\n");
        }
      }
    }
  }
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
 /* get_number_of_delimeters();
  number_of_factors = get_number_of_factors();
  set_number_of_elements(number_of_factors);
  set_number_of_elements(strlen(genom));*/
  /**/

  result = list_to_string(factorized_genom, second_parameter_value);

/*  number_of_factors = get_number_of_factors();
  set_number_of_elements(number_of_factors);
*/  return result;
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

  tail *t2 = malloc(sizeof(tail));
//a strange bug happens here causing t to have dirty values after fprintf. allocating another tail fixes it

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

void do_unit_testing() {
  printf("start of unit test\n");
  test_safe_fgets();
  test_foreach_element_in();
  test_initialize_tail();
  test_open_towrite_file();
  printf("unit test successfully passed\n\n");
}

void do_integration_testing() {
  printf("start of integration test\n");
  test_flush();
  test_fill_k_finger();
  test_create_fingerprint();
  test_process_file();
  test_process_fasta();
  //test of process_all_fasta_files has not been accomplished cause of time constraints
  printf("\nintegration test successfully completed\n\n");
  printf("running the program...\n\n");
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

  printf("\n\nstart of test_safe_fgets test\n");

  test_file = fopen("test_file.txt", "w");
  if (test_file == NULL) {
    printf("test_safe_fgets: test couldn't be executed cause test file cannot be opened in writing mode\n");
    printf("error: %s\n", strerror(errno));
    exit(1);
  }

  fprintf(test_file, "%s\n%s\n%s", "this is the first sentence", "this is the second one", "this is the third");
  fclose(test_file);

  test_file = fopen("test_file.txt", "r");

  if (test_file == NULL) {
    printf("test_safe_fgets: test couldn't be executed cause test file cannot be opened in reading mode\n");
    printf("error: %s\n", strerror(errno));
    exit(1);
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
    exit(1);
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
    exit(1);
  }

  printf("test conditions: scanning of 4 created subdirectories");

  strcat(path_test, "/");
  strcat(path_test, "directory-test");

  char subdir_path[strlen(path_test) + strlen("directory1") + 1];
  strcpy(subdir_path, path_test);
  strcat(subdir_path, "/directory1");
  file = opendir(subdir_path);
  if (file == NULL) {
    printf("for_each_element test cannot be completed cause created directory1 cannot be open\n");
    exit(1);
  }
  inner_file = readdir(file);
  if (inner_file == NULL) {
    printf("for_each_element test cannot be completed cause something gone wrong in catching a file\n");
    exit(1);
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


/*
  if parameter t is NULL, t will point to a new tail. If t is not NULL, it will be initialized.
  pre-condition: window_dimension > 0
  pre-condition: if t is not NULL, zero_one_tail && finger_tail must have window_dimension as dimension
  post-condition: t points to a tail ready to be used

  tail *intialize_tail(tail *t)
*/


void test_initialize_tail() {

  printf("\n\nStart of initialize_tail test\n");

  printf("test conditions: window dimension = 4 and tail = NULL\n");

  window_dimension = 4;
  tail *t = NULL;

  initialize_tail(t);

  printf("before initialization\n");
  int i;
  for (i = 0; i < window_dimension; i++) {
    t->zero_one_tail[i];
    t->finger_tail[i];
  }
  printf("test case passed\n");

  printf("test conditions: window dimension = 4 and tail != NULL\n");

  initialize_tail(t);

  window_dimension = 4;
  for (i = 0; i < window_dimension; i++) {
    t->zero_one_tail[i];
    t->finger_tail[i];
  }
  printf("test case passed\n");

  printf("initialize tail test correctly passed\n");
}

/* It records the tails content in the k-fingerprint file as the right format if you close the file_to_write after called this function.
   pre-condition file_to_write: must be opened in write mode;
   pre-condition t: must be previousely initialized and only manipulated by fill_k_fingerprint function
   pre-condition header != NULL
   pre-condition: after calling this function fclose(file_to_write) must be called

void flush(FILE *file_to_write, tail *t, char* header)
*/

void test_flush() {

  FILE *test_file;
  char s[300], c1, c2, header[100];
  int f1, f2, f3, f4, z1, z2, z3,z4, cont;
  tail *t = NULL;
  char *header_to_pass;

  printf("\n\nStart of flush test\n");
  printf("test conditions: window_dimension = 5, start_window_limit = 0, end_window_limit = 4\n");

  window_dimension = 4;

  t = initialize_tail(t);

  if (t == NULL)
    printf("problem is t");
  printf("correctly initialized\n");
  t->start_window_limit = 0;
  t->end_window_limit = 3;

  printf("before opening\n");
  test_file = fopen("test_file.txt", "w");
  if (test_file == NULL) {
    printf("test cannot be completed cause has not been possible to open kfingerprint_file\n");
    exit(1);
  }
  printf("building header\n");
  header_to_pass = (char *) malloc(100 * sizeof(char));
  strcpy(header_to_pass, "header");

  t->zero_one_tail[0] = 1;
  t->zero_one_tail[1] = 1;
  t->zero_one_tail[2] = 0;
  t->zero_one_tail[3] = 0;

  t->finger_tail[0] = 2;
  t->finger_tail[1] = 3;
  t->finger_tail[2] = 4;
  t->finger_tail[3] = 5;

  printf("before calling flush\n");
  flush(test_file, t, header_to_pass);

  printf("after\n");
  fclose(test_file);

  test_file = fopen("test_file.txt", "r");
  if (test_file == NULL) {
    printf("test cannot be completed cause has not been possible to open test_file\n");
    exit(1);
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
  free(t->zero_one_tail);
  free(t->finger_tail);
  free(t);
  free(header_to_pass);

  printf("test case passed\n");

  printf("test conditions: window_dimension = 4, start_window_limit = 0, end_window_limit = 4\n");

  window_dimension = 4;
  t = NULL;
  t = initialize_tail(t);

  t->start_window_limit = 0;
  t->end_window_limit = 3;

  test_file = fopen("test_file.txt", "w");
  if (test_file == NULL) {
    printf("test cannot be completed cause has not been possible to open kfingerprint_file\n");
    exit(1);
  }
  header_to_pass = (char *) malloc(100 * sizeof(char));
  strcpy(header_to_pass, "header");

  t->zero_one_tail[0] = 1;
  t->zero_one_tail[1] = 1;
  t->zero_one_tail[2] = 0;
  t->zero_one_tail[3] = 0;

  t->finger_tail[0] = 2;
  t->finger_tail[1] = 3;
  t->finger_tail[2] = 4;
  t->finger_tail[3] = 5;

  flush(test_file, t, header_to_pass);

  fclose(test_file);
  test_file = fopen("test_file.txt", "r");
  if (test_file == NULL) {
    printf("test cannot be completed cause has not been possible to open test_file\n");
    exit(1);
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
  free(t->zero_one_tail);
  free(t->finger_tail);
  free(t);
  free(header_to_pass);

  printf("test case passed\n");

  printf("test conditions: window_dimension = 4, start_window_limit = 3, end_window_limit = 2\n");

  window_dimension = 4;
  t = NULL;
  t = initialize_tail(t);

  test_file = fopen("test_file.txt", "w");
  if (test_file == NULL) {
    printf("test cannot be completed cause has not been possible to open kfingerprint_file\n");
    exit(1);
  }
  header_to_pass = (char *) malloc(100 * sizeof(char));
  strcpy(header_to_pass, "header");

  t->zero_one_tail[3] = 1;
  t->zero_one_tail[0] = 1;
  t->zero_one_tail[1] = 0;
  t->zero_one_tail[2] = 0;

  t->finger_tail[3] = 2;
  t->finger_tail[0] = 3;
  t->finger_tail[1] = 4;
  t->finger_tail[2] = 5;

  t->start_window_limit = 3;
  t->end_window_limit = 2;

  flush(test_file, t, header_to_pass);

  fclose(test_file);
  test_file = fopen("test_file.txt", "r");
  if (test_file == NULL) {
    printf("test cannot be completed cause has not been possible to open test_file\n");
    exit(1);
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
  free(t->zero_one_tail);
  free(t->finger_tail);
  free(t);
  free(header_to_pass);

  printf("test case passed\n");

  printf("flush test passed\n");
}

/*
  It returns a pointer to a file in write mode with the name made up of two parts: fasta-name and name.
  param name: name to concatenate with the fasta file one that will be the name of the file to be opened
  param directory_path: directory path in which the file to be created will be
  pre-condition name: != NULL
  pre-condition fasta_name: must contain ".fasta" at the end of its name
  pre-condition directory_path: must be the path of an existing directory
  post-condition: an opened file with the descripted name will be returned

FILE* open_towrite_file(char *name, char *fasta_name, char *directory_path) {
*/

void test_open_towrite_file() {

  FILE *test_file;
  char path_test[300];

  printf("\n\nStart of open_towrite_file test\n");
  printf("test conditions: fasta_name = \".fasta\" for factorization, fingerprint, kfingerprint and oneformat\n");

  if (getcwd(path_test, 255) == NULL) {
    printf("for_each_element_in test cannot be run cause it hasn't been possible to find the current path\n");
    printf("error: %s\n", strerror(errno));
    exit(1);
  }

  test_file = open_towrite_file("factorization", ".fasta", path_test);

  if (test_file == NULL)
    assert(0);

  fclose(test_file);

  test_file = fopen("-factorization", "r");
  if (test_file == NULL) {
    printf("-factorization.txt not created\n");
    assert(0);
  }

  remove("-factorization");

  printf("test case passed\n");

  printf("test conditions: fasta_name = \"normal_name.fasta\" for factorization, fingerprint, kfingerprint and oneformat\n");

  if (getcwd(path_test, 255) == NULL) {
    printf("for_each_element_in test cannot be run cause it hasn't been possible to find the current path\n");
    printf("error: %s\n", strerror(errno));
    exit(1);
  }

  test_file = open_towrite_file("factorization", "normal_name.fasta", path_test);

  if (test_file == NULL)
    assert(0);

  fclose(test_file);

  test_file = fopen("normal_name-factorization", "r");

  if (test_file == NULL) {
    printf("normal_name_factorization not created\n");
    assert(0);
  }

  remove("normal_name-factorization");

  printf("test case passed\n");

  printf("open_towrite_file test passed\n");


  /*file name of length next to 255 is not considered here cause of time constraints*/
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

void fill_k_fingerprint(int fingerprint_number, FILE *file_to_write, tail *t, char* header_to_pass) {
*/

void test_fill_k_finger() {

  char s[300], c1, c2, header[100];
  int f1, f2, f3, f4, z1, z2, z3,z4, cont;
  FILE *test_file;
  char *header_to_pass;
  tail *t;

  printf("\n\nStart of fill_k_finger test\n");
  printf("test conditions: number of fingerprints all less then the actual size of the tails not considering sub-factorization values\n");

  window_dimension = 4;

  test_file = fopen("test_file", "w");
  if (test_file == NULL) {
    printf("test couldn't be completed cause kfingerprint file cannot be opened in writing mode\n");
    exit(1);
  }

  header_to_pass = (char *) malloc(100 * sizeof(char));
  if (header_to_pass == NULL) {
    printf("test couldn't be completed cause malloc returned NULL on header_read\n");
  }
  strcpy(header_to_pass, "header1");

  t = NULL;
  t = initialize_tail(t);

  fill_k_fingerprint(1, test_file, t, header_to_pass);
  fill_k_fingerprint(2, test_file, t, header_to_pass);
  fill_k_fingerprint(-2, test_file, t, header_to_pass);

  fclose(test_file);

  test_file = fopen("test_file", "r");
  if (test_file == NULL) {
    printf("test couldn't be completed cause test_file cannot be opened in reading mode\n");
    remove("test_file");
    exit(1);
  }

  if(fgets(s, 300, test_file) == NULL)
    assert(0);

  sscanf(s, "%d %d %d %d %c %d %d %d %d %c %s %d", &f1, &f2, &f3, &f4, &c1, &z1, &z2, &z3, &z4, &c2, header, &cont);

  assert(f1 == 1);
  assert(f2 == 2);
  assert(f3 == 0);
  assert(f4 == 0);
  assert(c1 == '$');
  assert(z1 == 0);
  assert(z2 == 0);
  assert(z3 == 0);
  assert(z4 == 0);
  assert(c2 == '$');
  assert(strcmp(header, "header1") == 0);
  assert(cont == 0);

  if(fgets(s, 300, test_file) != NULL)
    assert(0);

  system("remove -r test_file");

  free(t->zero_one_tail);
  free(t->finger_tail);
  free(t);

  printf("test case passed\n");

  printf("test conditions: number of fingerprints equals the actual size of the tails not considering sub-factorization values\n");

  window_dimension = 4;

  test_file = fopen("test_file", "w");
  if (test_file == NULL) {
    printf("test couldn't be completed cause kfingerprint file cannot be opened in writing mode\n");
    exit(1);
  }

  header_to_pass = (char *) malloc(100 * sizeof(char));
  if (header_to_pass == NULL) {
    printf("test couldn't be completed cause malloc returned NULL on header_read\n");
    exit(1);
  }
  strcpy(header_to_pass, "header1");

  t = NULL;
  t = initialize_tail(t);

  fill_k_fingerprint(1, test_file, t, header_to_pass);
  fill_k_fingerprint(2, test_file, t, header_to_pass);
  fill_k_fingerprint(3, test_file, t, header_to_pass);
  fill_k_fingerprint(4, test_file, t, header_to_pass);
  fill_k_fingerprint(-2, test_file, t, header_to_pass);

  fclose(test_file);

  test_file = fopen("test_file", "r");
  if (test_file == NULL) {
    printf("test couldn't be completed cause test_file cannot be opened in reading mode\n");
    remove("test_file");
    exit(1);
  }

  if(fgets(s, 300, test_file) == NULL)
    assert(0);

  sscanf(s, "%d %d %d %d %c %d %d %d %d %c %s %d", &f1, &f2, &f3, &f4, &c1, &z1, &z2, &z3, &z4, &c2, header, &cont);

  assert(f1 == 1);
  assert(f2 == 2);
  assert(f3 == 3);
  assert(f4 == 4);
  assert(c1 == '$');
  assert(z1 == 0);
  assert(z2 == 0);
  assert(z3 == 0);
  assert(z4 == 0);
  assert(c2 == '$');
  assert(strcmp(header, "header1") == 0);
  assert(cont == 0);

  if(fgets(s, 300, test_file) != NULL)
    assert(0);

  remove("test_file");

  free(t->zero_one_tail);
  free(t->finger_tail);
  free(t);

  printf("test case passed\n");

  printf("test conditions: number of fingerprints are more then the actual size of the tails not considering sub-factorization values\n");

  window_dimension = 4;

  test_file = fopen("test_file", "w");
  if (test_file == NULL) {
    printf("test couldn't be completed cause kfingerprint file cannot be opened in writing mode\n");
    exit(1);
  }

  header_to_pass = (char *) malloc(100 * sizeof(char));
  if (header_to_pass == NULL) {
    printf("test couldn't be completed cause malloc returned NULL on header_read\n");
    exit(1);
  }
  strcpy(header_to_pass, "header1");

  t = NULL;
  t = initialize_tail(t);

  fill_k_fingerprint(1, test_file, t, header_to_pass);
  fill_k_fingerprint(2, test_file, t, header_to_pass);
  fill_k_fingerprint(3, test_file, t, header_to_pass);
  fill_k_fingerprint(4, test_file, t, header_to_pass);
  fill_k_fingerprint(5, test_file, t, header_to_pass);
  fill_k_fingerprint(6, test_file, t, header_to_pass);

  fclose(test_file);

  test_file = fopen("test_file", "r");
  if (test_file == NULL) {
    printf("test couldn't be completed cause test_file cannot be opened in reading mode\n");
    remove("test_file");
    exit(1);
  }

  if(fgets(s, 300, test_file) == NULL)
    assert(0);

  sscanf(s, "%d %d %d %d %c %d %d %d %d %c %s %d", &f1, &f2, &f3, &f4, &c1, &z1, &z2, &z3, &z4, &c2, header, &cont);

  assert(f1 == 1);
  assert(f2 == 2);
  assert(f3 == 3);
  assert(f4 == 4);
  assert(c1 == '$');
  assert(z1 == 0);
  assert(z2 == 0);
  assert(z3 == 0);
  assert(z4 == 0);
  assert(c2 == '$');
  assert(strcmp(header, "header1") == 0);
  assert(cont == 0);

  if(fgets(s, 300, test_file) == NULL)
    assert(0);

  sscanf(s, "%d %d %d %d %c %d %d %d %d %c %s %d", &f1, &f2, &f3, &f4, &c1, &z1, &z2, &z3, &z4, &c2, header, &cont);

  assert(f1 == 2);
  assert(f2 == 3);
  assert(f3 == 4);
  assert(f4 == 5);
  assert(c1 == '$');
  assert(z1 == 0);
  assert(z2 == 0);
  assert(z3 == 0);
  assert(z4 == 0);
  assert(c2 == '$');
  assert(strcmp(header, "header1") == 0);
  assert(cont == 1);

  if(fgets(s, 300, test_file) == NULL)
    assert(0);

  sscanf(s, "%d %d %d %d %c %d %d %d %d %c %s %d", &f1, &f2, &f3, &f4, &c1, &z1, &z2, &z3, &z4, &c2, header, &cont);

  assert(f1 == 3);
  assert(f2 == 4);
  assert(f3 == 5);
  assert(f4 == 6);
  assert(c1 == '$');
  assert(z1 == 0);
  assert(z2 == 0);
  assert(z3 == 0);
  assert(z4 == 0);
  assert(c2 == '$');
  assert(strcmp(header, "header1") == 0);
  assert(cont == 3);

  if(fgets(s, 300, test_file) != NULL)
    assert(0);

  remove("test_file");

  free(t->zero_one_tail);
  free(t->finger_tail);
  free(t);

  printf("test case passed\n");

  printf("test conditions: number of fingerprints is less than the actual size of the tails considering sub-factorization values at extrems\n");

  window_dimension = 4;

  test_file = fopen("test_file", "w");
  if (test_file == NULL) {
    printf("test couldn't be completed cause kfingerprint file cannot be opened in writing mode\n");
    exit(1);
  }

  header_to_pass = (char *) malloc(100 * sizeof(char));
  if (header_to_pass == NULL) {
    printf("test couldn't be completed cause malloc returned NULL on header_read\n");
    exit(1);
  }
  strcpy(header_to_pass, "header1");

  t = NULL;
  t = initialize_tail(t);

  fill_k_fingerprint(-1, test_file, t, header_to_pass);
  fill_k_fingerprint(2, test_file, t, header_to_pass);
  fill_k_fingerprint(3, test_file, t, header_to_pass);
  fill_k_fingerprint(0, test_file, t, header_to_pass);
  fill_k_fingerprint(-2, test_file, t, header_to_pass);

  fclose(test_file);

  test_file = fopen("test_file", "r");
  if (test_file == NULL) {
    printf("test couldn't be completed cause test_file cannot be opened in reading mode\n");
    remove("test_file");
    exit(1);
  }

  if(fgets(s, 300, test_file) == NULL)
    assert(0);

  sscanf(s, "%d %d %d %d %c %d %d %d %d %c %s %d", &f1, &f2, &f3, &f4, &c1, &z1, &z2, &z3, &z4, &c2, header, &cont);

  assert(f1 == 2);
  assert(f2 == 3);
  assert(f3 == 0);
  assert(f4 == 0);
  assert(c1 == '$');
  assert(z1 == 1);
  assert(z2 == 1);
  assert(z3 == 0);
  assert(z4 == 0);
  assert(c2 == '$');
  assert(strcmp(header, "header1") == 0);
  assert(cont == 0);

  if(fgets(s, 300, test_file) != NULL)
    assert(0);

  remove("test_file");

  free(t->zero_one_tail);
  free(t->finger_tail);
  free(t);

  printf("test case passed\n");

  printf("test conditions: number of fingerprints equals the actual size of the tails considering sub-factorization values placed at extrems\n");

  window_dimension = 4;

  test_file = fopen("test_file", "w");
  if (test_file == NULL) {
    printf("test couldn't be completed cause kfingerprint file cannot be opened in writing mode\n");
    exit(1);
  }

  header_to_pass = (char *) malloc(100 * sizeof(char));
  if (header_to_pass == NULL) {
    printf("test couldn't be completed cause malloc returned NULL on header_read\n");
    exit(1);
  }
  strcpy(header_to_pass, "header1");

  t = NULL;
  t = initialize_tail(t);

  fill_k_fingerprint(-1, test_file, t, header_to_pass);
  fill_k_fingerprint(2, test_file, t, header_to_pass);
  fill_k_fingerprint(3, test_file, t, header_to_pass);
  fill_k_fingerprint(4, test_file, t, header_to_pass);
  fill_k_fingerprint(5, test_file, t, header_to_pass);
  fill_k_fingerprint(0, test_file, t, header_to_pass);
  fill_k_fingerprint(-2, test_file, t, header_to_pass);

  fclose(test_file);

  test_file = fopen("test_file", "r");
  if (test_file == NULL) {
    printf("test couldn't be completed cause test_file cannot be opened in reading mode\n");
    remove("test_file");
    exit(1);
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
  assert(z3 == 1);
  assert(z4 == 1);
  assert(c2 == '$');
  assert(strcmp(header, "header1") == 0);
  assert(cont == 0);

  if(fgets(s, 300, test_file) != NULL)
    assert(0);

  remove("test_file");

  free(t->zero_one_tail);
  free(t->finger_tail);
  free(t);

  printf("test case passed\n");
  printf("test conditions: number of fingerprints is more than the actual size of the tails considering sub-factorization values\n");

  window_dimension = 4;

  test_file = fopen("test_file", "w");
  if (test_file == NULL) {
    printf("test couldn't be completed cause kfingerprint file cannot be opened in writing mode\n");
    exit(1);
  }

  header_to_pass = (char *) malloc(100 * sizeof(char));
  if (header_to_pass == NULL) {
    printf("test couldn't be completed cause malloc returned NULL on header_read\n");
    exit(1);
  }
  strcpy(header_to_pass, "header1");

  t = NULL;
  t = initialize_tail(t);

  fill_k_fingerprint(-1, test_file, t, header_to_pass);
  fill_k_fingerprint(1, test_file, t, header_to_pass);
  fill_k_fingerprint(2, test_file, t, header_to_pass);
  fill_k_fingerprint(0, test_file, t, header_to_pass);
  fill_k_fingerprint(3, test_file, t, header_to_pass);
  fill_k_fingerprint(4, test_file, t, header_to_pass);
  fill_k_fingerprint(5, test_file, t, header_to_pass);
  fill_k_fingerprint(6, test_file, t, header_to_pass);

  fclose(test_file);

  test_file = fopen("test_file", "r");
  if (test_file == NULL) {
    printf("test couldn't be completed cause test_file cannot be opened in reading mode\n");
    remove("test_file");
    exit(1);
  }

  if(fgets(s, 300, test_file) == NULL)
    assert(0);

  sscanf(s, "%d %d %d %d %c %d %d %d %d %c %s %d", &f1, &f2, &f3, &f4, &c1, &z1, &z2, &z3, &z4, &c2, header, &cont);

  assert(f1 == 1);
  assert(f2 == 2);
  assert(f3 == 3);
  assert(f4 == 4);
  assert(c1 == '$');
  assert(z1 == 1);
  assert(z2 == 1);
  assert(z3 == 0);
  assert(z4 == 0);
  assert(c2 == '$');
  assert(strcmp(header, "header1") == 0);
  assert(cont == 0);

  if(fgets(s, 300, test_file) == NULL)
    assert(0);

  sscanf(s, "%d %d %d %d %c %d %d %d %d %c %s %d", &f1, &f2, &f3, &f4, &c1, &z1, &z2, &z3, &z4, &c2, header, &cont);

  assert(f1 == 2);
  assert(f2 == 3);
  assert(f3 == 4);
  assert(f4 == 5);
  assert(c1 == '$');
  assert(z1 == 1);
  assert(z2 == 0);
  assert(z3 == 0);
  assert(z4 == 0);
  assert(c2 == '$');
  assert(strcmp(header, "header1") == 0);
  assert(cont == 1);

  if(fgets(s, 300, test_file) == NULL)
    assert(0);

  sscanf(s, "%d %d %d %d %c %d %d %d %d %c %s %d", &f1, &f2, &f3, &f4, &c1, &z1, &z2, &z3, &z4, &c2, header, &cont);

  assert(f1 == 3);
  assert(f2 == 4);
  assert(f3 == 5);
  assert(f4 == 6);
  assert(c1 == '$');
  assert(z1 == 0);
  assert(z2 == 0);
  assert(z3 == 0);
  assert(z4 == 0);
  assert(c2 == '$');
  assert(strcmp(header, "header1") == 0);
  assert(cont == 3);

  if(fgets(s, 300, test_file) != NULL)
    assert(0);

  remove("test_file");

  free(t->zero_one_tail);
  free(t->finger_tail);
  free(t);

  printf("test case passed\n");

  printf("fill_k_finger test passed\n");
}

/* It creates the fingerprint and k-fingerprint of the factorized genom. The k-fingerprints will be stored in
       the file pointed by kfingerprint_file variable.
  pre-condition factorized_genom: must be the result of factorize_read or format equivalent
  pre-condition file_to_write: must point a file opened in write mode
  pre-condition header: != NULL
  pre-condition file_to_write must be closed to save the results in the file permanently
  pre-condition: window_dimension > 0
  post-condition: kfingerprints will be written in the correct format
  return: fingerprint of the factorized genom

char* create_fingerprint(char* factorized_genom, FILE *file_to_write, char *header) {
*/
void test_create_fingerprint() {

  char s[300], c1, c2, header[100], *finger_result;
  int f1, f2, f3, f4, z1, z2, z3,z4, cont;
  FILE *test_file;
  tail *t;
  char *header_to_pass;

  printf("\n\nStart of create_fingerprint test\n");
  printf("test conditions: number of fingerprints all less then the actual size of the tails with no sub_factorization\n");

  window_dimension = 4;

  test_file = fopen("test_file", "w");
  if (test_file == NULL) {
    printf("test couldn't be completed cause kfingerprint file cannot be opened in writing mode\n");
    exit(1);
  }

  header_to_pass = (char *) malloc(100 * sizeof(char));
  if (header_to_pass == NULL) {
    printf("test couldn't be completed cause malloc returned NULL on header_read\n");
    exit(1);
  }
  strcpy(header_to_pass, "header1");

  t = NULL;
  t = initialize_tail(t);
/*
  fill_k_fingerprint(1);
  fill_k_fingerprint(2);
  fill_k_fingerprint(-2);
*/
  finger_result = create_fingerprint("[ \"A\" \"GG\" ]", test_file, header_to_pass);
  printf("%s\n", finger_result);
  assert(strcmp("1,2", finger_result) == 0);
  fclose(test_file);
  test_file = fopen("test_file", "r");
  if (test_file == NULL) {
    printf("test couldn't be completed cause test_file cannot be opened in reading mode\n");
    remove("test_file");
    exit(1);
  }

  if(fgets(s, 300, test_file) == NULL)
    assert(0);

  sscanf(s, "%d %d %d %d %c %d %d %d %d %c %s %d", &f1, &f2, &f3, &f4, &c1, &z1, &z2, &z3, &z4, &c2, header, &cont);

  assert(f1 == 1);
  assert(f2 == 2);
  assert(f3 == 0);
  assert(f4 == 0);
  assert(c1 == '$');
  assert(z1 == 0);
  assert(z2 == 0);
  assert(z3 == 0);
  assert(z4 == 0);
  assert(c2 == '$');
  assert(strcmp(header, "header1") == 0);
  assert(cont == 0);

  if(fgets(s, 300, test_file) != NULL)
    assert(0);

  remove("test_file");

  printf("test case passed\n");

  printf("test conditions: number of fingerprints equals the actual size of the tails with no sub-factorization values\n");

  window_dimension = 4;

  test_file = fopen("test_file", "w");
  if (test_file == NULL) {
    printf("test couldn't be completed cause kfingerprint file cannot be opened in writing mode\n");
    exit(1);
  }

  header_to_pass = (char *) malloc(100 * sizeof(char));
  if (header_to_pass == NULL) {
    printf("test couldn't be completed cause malloc returned NULL on header_read\n");
    exit(1);
  }
  strcpy(header_to_pass, "header1");

  t = initialize_tail(t);

/*
  fill_k_fingerprint(1);
  fill_k_fingerprint(2);
  fill_k_fingerprint(3);
  fill_k_fingerprint(4);
  fill_k_fingerprint(-2);
*/
  finger_result = create_fingerprint("[ \"A\" \"GG\" \"TTT\" \"CCCC\" ]", test_file, header_to_pass);

  assert(strcmp("1,2,3,4", finger_result) == 0);

  fclose(test_file);
  test_file = fopen("test_file", "r");
  if (test_file == NULL) {
    printf("test couldn't be completed cause test_file cannot be opened in reading mode\n");
    remove("test_file");
    exit(1);
  }

  if(fgets(s, 300, test_file) == NULL)
    assert(0);

  sscanf(s, "%d %d %d %d %c %d %d %d %d %c %s %d", &f1, &f2, &f3, &f4, &c1, &z1, &z2, &z3, &z4, &c2, header, &cont);

  assert(f1 == 1);
  assert(f2 == 2);
  assert(f3 == 3);
  assert(f4 == 4);
  assert(c1 == '$');
  assert(z1 == 0);
  assert(z2 == 0);
  assert(z3 == 0);
  assert(z4 == 0);
  assert(c2 == '$');
  assert(strcmp(header, "header1") == 0);
  assert(cont == 0);

  if(fgets(s, 300, test_file) != NULL)
    assert(0);

  remove("test_file");

  printf("test case passed\n");

  printf("test conditions: number of fingerprints are more then the actual size of the tails not considering sub-factorization values\n");

  window_dimension = 4;

  test_file = fopen("test_file", "w");
  if (test_file == NULL) {
    printf("test couldn't be completed cause kfingerprint file cannot be opened in writing mode\n");
    exit(1);
  }

  header_to_pass = (char *) malloc(100 * sizeof(char));
  if (header_to_pass == NULL) {
    printf("test couldn't be completed cause malloc returned NULL on header_read\n");
    exit(1);
  }
  strcpy(header_to_pass, "header1");

  t = initialize_tail(t);
/*
  fill_k_fingerprint(1);
  fill_k_fingerprint(2);
  fill_k_fingerprint(3);
  fill_k_fingerprint(4);
  fill_k_fingerprint(5);
  fill_k_fingerprint(6);
*/
  finger_result = create_fingerprint("[ \"A\" \"GG\" \"CCC\" \"TTTT\" \"GGTTA\" \"AGATTC\" ]", test_file, header_to_pass);
  assert(strcmp("1,2,3,4,5,6", finger_result) == 0);

  fclose(test_file);
  test_file = fopen("test_file", "r");
  if (test_file == NULL) {
    printf("test couldn't be completed cause test_file cannot be opened in reading mode\n");
    remove("test_file");
    exit(1);
  }

  if(fgets(s, 300, test_file) == NULL)
    assert(0);

  sscanf(s, "%d %d %d %d %c %d %d %d %d %c %s %d", &f1, &f2, &f3, &f4, &c1, &z1, &z2, &z3, &z4, &c2, header, &cont);

  assert(f1 == 1);
  assert(f2 == 2);
  assert(f3 == 3);
  assert(f4 == 4);
  assert(c1 == '$');
  assert(z1 == 0);
  assert(z2 == 0);
  assert(z3 == 0);
  assert(z4 == 0);
  assert(c2 == '$');
  assert(strcmp(header, "header1") == 0);
  assert(cont == 0);

  if(fgets(s, 300, test_file) == NULL)
    assert(0);

  sscanf(s, "%d %d %d %d %c %d %d %d %d %c %s %d", &f1, &f2, &f3, &f4, &c1, &z1, &z2, &z3, &z4, &c2, header, &cont);

  assert(f1 == 2);
  assert(f2 == 3);
  assert(f3 == 4);
  assert(f4 == 5);
  assert(c1 == '$');
  assert(z1 == 0);
  assert(z2 == 0);
  assert(z3 == 0);
  assert(z4 == 0);
  assert(c2 == '$');
  assert(strcmp(header, "header1") == 0);
  assert(cont == 1);

  if(fgets(s, 300, test_file) == NULL)
    assert(0);

  sscanf(s, "%d %d %d %d %c %d %d %d %d %c %s %d", &f1, &f2, &f3, &f4, &c1, &z1, &z2, &z3, &z4, &c2, header, &cont);

  assert(f1 == 3);
  assert(f2 == 4);
  assert(f3 == 5);
  assert(f4 == 6);
  assert(c1 == '$');
  assert(z1 == 0);
  assert(z2 == 0);
  assert(z3 == 0);
  assert(z4 == 0);
  assert(c2 == '$');
  assert(strcmp(header, "header1") == 0);
  assert(cont == 3);

  if(fgets(s, 300, test_file) != NULL)
    assert(0);

  remove("test_file");

  printf("test case passed\n");

  printf("test conditions: number of fingerprints is less than the actual size of the tails considering sub-factorization values at extrems\n");

  window_dimension = 4;

  test_file = fopen("test_file", "w");
  if (test_file == NULL) {
    printf("test couldn't be completed cause kfingerprint file cannot be opened in writing mode\n");
    exit(1);
  }

  header_to_pass = (char *) malloc(100 * sizeof(char));
  if (header_to_pass == NULL) {
    printf("test couldn't be completed cause malloc returned NULL on header_read\n");
    exit(1);
  }
  strcpy(header_to_pass, "header1");

  t = initialize_tail(t);
/*
  fill_k_fingerprint(-1);
  fill_k_fingerprint(2);
  fill_k_fingerprint(3);
  fill_k_fingerprint(0);
  fill_k_fingerprint(-2);
*/

  finger_result = create_fingerprint("[ \"<<\" \"AA\" \"GGG\" \">>\" ]", test_file, header_to_pass);

  assert(strcmp("-1,2,3,0", finger_result) == 0);

  fclose(test_file);
  test_file = fopen("test_file", "r");
  if (test_file == NULL) {
    printf("test couldn't be completed cause test_file cannot be opened in reading mode\n");
    remove("test_file");
    exit(1);
  }

  if(fgets(s, 300, test_file) == NULL)
    assert(0);

  sscanf(s, "%d %d %d %d %c %d %d %d %d %c %s %d", &f1, &f2, &f3, &f4, &c1, &z1, &z2, &z3, &z4, &c2, header, &cont);

  assert(f1 == 2);
  assert(f2 == 3);
  assert(f3 == 0);
  assert(f4 == 0);
  assert(c1 == '$');
  assert(z1 == 1);
  assert(z2 == 1);
  assert(z3 == 0);
  assert(z4 == 0);
  assert(c2 == '$');
  assert(strcmp(header, "header1") == 0);
  assert(cont == 0);

  if(fgets(s, 300, test_file) != NULL)
    assert(0);

  remove("test_file");

  printf("test case passed\n");

  printf("test conditions: number of fingerprints equals the actual size of the tails considering sub-factorization values placed at extrems\n");

  window_dimension = 4;

  test_file = fopen("test_file", "w");
  if (test_file == NULL) {
    printf("test couldn't be completed cause kfingerprint file cannot be opened in writing mode\n");
    exit(1);
  }

  header_to_pass = (char *) malloc(100 * sizeof(char));
  if (header_to_pass == NULL) {
    printf("test couldn't be completed cause malloc returned NULL on header_read\n");
  }
  strcpy(header_to_pass, "header1");

  t = initialize_tail(t);
/*
  fill_k_fingerprint(-1);
  fill_k_fingerprint(2);
  fill_k_fingerprint(3);
  fill_k_fingerprint(4);
  fill_k_fingerprint(5);
  fill_k_fingerprint(0);
  fill_k_fingerprint(-2);
*/

  finger_result = create_fingerprint("[ \"<<\" \"AA\" \"GGG\" \"CCCC\" \"TTTTT\" \">>\" ]", test_file, header_to_pass);

  assert(strcmp("-1,2,3,4,5,0", finger_result) == 0);

  fclose(test_file);
  test_file = fopen("test_file", "r");
  if (test_file == NULL) {
    printf("test couldn't be completed cause test_file cannot be opened in reading mode\n");
    remove("test_file");
    exit(1);
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
  assert(z3 == 1);
  assert(z4 == 1);
  assert(c2 == '$');
  assert(strcmp(header, "header1") == 0);
  assert(cont == 0);

  if(fgets(s, 300, test_file) != NULL)
    assert(0);

  remove("test_file");

  printf("test case passed\n");
  printf("test conditions: number of fingerprints is more than the actual size of the tails considering sub-factorization values\n");

  window_dimension = 4;

  test_file = fopen("test_file", "w");
  if (test_file == NULL) {
    printf("test couldn't be completed cause kfingerprint file cannot be opened in writing mode\n");
    exit(1);
  }

  header_to_pass = (char *) malloc(100 * sizeof(char));
  if (header_to_pass == NULL) {
    printf("test couldn't be completed cause malloc returned NULL on header_read\n");
  }
  strcpy(header_to_pass, "header1");

  t = initialize_tail(t);

  fill_k_fingerprint(-1, test_file, t, header_to_pass);
  fill_k_fingerprint(1, test_file, t, header_to_pass);
  fill_k_fingerprint(2, test_file, t, header_to_pass);
  fill_k_fingerprint(0, test_file, t, header_to_pass);
  fill_k_fingerprint(3, test_file, t, header_to_pass);
  fill_k_fingerprint(4, test_file, t, header_to_pass);
  fill_k_fingerprint(5, test_file, t, header_to_pass);
  fill_k_fingerprint(6, test_file, t, header_to_pass);


  fclose(test_file);

  test_file = fopen("test_file", "r");
  if (test_file == NULL) {
    printf("test couldn't be completed cause test_file cannot be opened in reading mode\n");
    remove("test_file");
    exit(1);
  }

  if(fgets(s, 300, test_file) == NULL)
    assert(0);

  sscanf(s, "%d %d %d %d %c %d %d %d %d %c %s %d", &f1, &f2, &f3, &f4, &c1, &z1, &z2, &z3, &z4, &c2, header, &cont);

  assert(f1 == 1);
  assert(f2 == 2);
  assert(f3 == 3);
  assert(f4 == 4);
  assert(c1 == '$');
  assert(z1 == 1);
  assert(z2 == 1);
  assert(z3 == 0);
  assert(z4 == 0);
  assert(c2 == '$');
  assert(strcmp(header, "header1") == 0);
  assert(cont == 0);

  if(fgets(s, 300, test_file) == NULL)
    assert(0);

  sscanf(s, "%d %d %d %d %c %d %d %d %d %c %s %d", &f1, &f2, &f3, &f4, &c1, &z1, &z2, &z3, &z4, &c2, header, &cont);

  assert(f1 == 2);
  assert(f2 == 3);
  assert(f3 == 4);
  assert(f4 == 5);
  assert(c1 == '$');
  assert(z1 == 1);
  assert(z2 == 0);
  assert(z3 == 0);
  assert(z4 == 0);
  assert(c2 == '$');
  assert(strcmp(header, "header1") == 0);
  assert(cont == 1);

  if(fgets(s, 300, test_file) == NULL)
    assert(0);

  sscanf(s, "%d %d %d %d %c %d %d %d %d %c %s %d", &f1, &f2, &f3, &f4, &c1, &z1, &z2, &z3, &z4, &c2, header, &cont);

  assert(f1 == 3);
  assert(f2 == 4);
  assert(f3 == 5);
  assert(f4 == 6);
  assert(c1 == '$');
  assert(z1 == 0);
  assert(z2 == 0);
  assert(z3 == 0);
  assert(z4 == 0);
  assert(c2 == '$');
  assert(strcmp(header, "header1") == 0);
  assert(cont == 3);

  if(fgets(s, 300, test_file) != NULL)
    assert(0);

  remove("test_file");

  printf("test case passed\n");

  printf("create_fingerprint test passed\n");

}

/*
  It is the thread function that processes the file given by the array files_to_process at index given by arg, setting
  available_threads at same index to true. If no file is given, that is if the filesystem_node of files_to_process at the same
  index has variable deleted to true, it stays in wait mode until variable end_of_processing is set to true.
  param: integer that is the index of the interested informations in thread_ids, files_to_process, available_threads
  pre-condition: this function must be passed to pthread_create
  pre-condition arg: must be an integer between 0 and NUMBER_PROCESSORS
  pre-condition: window_dimension > 0;

void *process_file(void* arg)
*/
void test_process_file() {

  FILE *test_file;
  char path_test[300];
  DIR *file;
  struct dirent *inner_file;
  char header_test[300], genom[600];
  char s[300], c1, c2;
  int f1, f2, f3, f4, z1, z2, z3,z4, cont;
  pthread_t tid;
  pthread_attr_t att;
  int index = 0;

  printf("\n\nStart of process_file test\n");
  printf("test conditions: processing fasta with 3 reads, creting the filesystem_node after launching the thread\n");

   if (getcwd(path_test, 255) == NULL) {
    printf("process_fasta test cannot be run cause it hasn't been possible to find the current path\n");
    printf("error: %s\n", strerror(errno));
    exit(1);
  }

  test_file = fopen("test_file.fasta", "w");
  if (test_file == NULL) {
    printf("test couldn't be completed cause test_file cannot be opened in reading mode\n");
    remove("test_file");
    exit(1);
  }

  fprintf(test_file, "%s\n%s\n", ">header1", "CGTTGCGGAAAGGTC");
  fprintf(test_file, "%s\n%s\n", ">header2", "GTCCCCCAAAAGGGCTC");
  fprintf(test_file, "%s\n%s\n", ">header3", "GTCTCCCACCTCAG");

  fclose(test_file);

  file = opendir(path_test);
  if (path_test == NULL) {
    printf("process_file test could not be completed cause directory cannot be opened\n");
    exit(1);
  }

  inner_file = readdir(file);
  if (inner_file == NULL) {
    printf("process file test could not be completed cause directory cannot be read\n");
    exit(1);
  }

  while (strcmp(inner_file->d_name, "test_file.fasta") != 0){
    inner_file = readdir(file);
    if (inner_file == NULL) {
      printf("process file test could not be completed cause directory cannot be read\n");
      exit(1);
    }
  }
  closedir(file);

  end_of_processing = 0;
  files_to_process[0].deleted = 1;
  pthread_attr_init(&att);
  pthread_create(&tid, &att, process_file, &index);

  window_dimension = 4;
  strcpy(files_to_process[0].directory_path, path_test);
  strcpy(files_to_process[0].filename, "test_file.fasta");
  strcat(path_test, "/");
  strcat(path_test, "test_file.fasta");
  communicate_max_fact_length(0);
  fact_choice = 1;
  files_to_process[0].deleted = 0;
  
  end_of_processing = 1;
  pthread_join(tid, NULL);

  test_file = fopen("test_file-factorization", "r");
  if (test_file == NULL) {
    printf("test couldn't be completed cause test_file-factorization cannot be opened in reading mode\n");
    printf("error: %s\n", strerror(errno));
    remove("test_file-factorization");
    remove("test_file-fingerprint");
    remove("test_file-kfingerprint");
    remove("test_file-oneformat");
    remove("test_file.fasta");
    exit(1);
  }

  fgets(header_test, 300, test_file);
  fgets(genom, 300, test_file);
  header_test[strlen(header_test) - 1] = '\0';
  genom[strlen(genom) - 1] = '\0';
  assert(strcmp(header_test, ">header1") == 0);
  assert(strcmp(genom, "[ \"CGTTG\" \"CGG\" \"AAAGGTC\" ]") == 0);
  fgets(header_test, 300, test_file);
  fgets(genom, 300, test_file);
  header_test[strlen(header_test) - 1] = '\0';
  genom[strlen(genom) - 1] = '\0';
  assert(strcmp(header_test, ">header2") == 0);
  assert(strcmp(genom, "[ \"GT\" \"C\" \"C\" \"C\" \"C\" \"C\" \"AAAAGGGCTC\" ]") == 0);
  fgets(header_test, 300, test_file);
  fgets(genom, 300, test_file);
  header_test[strlen(header_test) - 1] = '\0';
  genom[strlen(genom) - 1] = '\0';
  assert(strcmp(header_test, ">header3") == 0);
  assert(strcmp(genom, "[ \"GT\" \"CT\" \"C\" \"C\" \"C\" \"ACCTCAG\" ]") == 0);

  if(fgets(s, 300, test_file) != NULL)
    assert(0);

  fclose(test_file);

  test_file = fopen("test_file-fingerprint", "r");
  if (test_file == NULL) {
    printf("test couldn't be completed cause test_file-fingerprint.txt cannot be opened in reading mode\n");
    printf("error: %s\n", strerror(errno));
    remove("test_file-factorization");
    remove("test_file-fingerprint");
    remove("test_file-kfingerprint");
    remove("test_file-oneformat");
    remove("test_file.fasta");
    exit(1);
  }

  fscanf(test_file, "%s %s\n", header_test, genom);
  assert(strcmp(header_test, ">header1") == 0);
  assert(strcmp(genom, "5,3,7") == 0);
  fscanf(test_file, "%s\n%s\n", header_test, genom);
  assert(strcmp(header_test, ">header2") == 0);
  assert(strcmp(genom, "2,1,1,1,1,1,10") == 0);
  fscanf(test_file, "%s\n%s\n", header_test, genom);
  assert(strcmp(header_test, ">header3") == 0);
  assert(strcmp(genom, "2,2,1,1,1,7") == 0);

  if(fgets(s, 300, test_file) != NULL)
    assert(0);
  fclose(test_file);

  test_file = fopen("test_file-kfingerprint", "r");
  if (test_file == NULL) {
    printf("test couldn't be completed cause test_file-kfingerprint.txt cannot be opened in reading mode\n");
    remove("test_file-factorization");
    remove("test_file-fingerprint");
    remove("test_file-kfingerprint");
    remove("test_file-oneformat");
    remove("test_file.fasta");
    exit(1);
  }

 if(fgets(s, 300, test_file) == NULL)
    assert(0);
  s[strlen(s) - 1] = '\0';

  sscanf(s, "%d %d %d %d %c %d %d %d %d %c %s %d", &f1, &f2, &f3, &f4, &c1, &z1, &z2, &z3, &z4, &c2, header_test, &cont);

  assert(f1 == 5);
  assert(f2 == 3);
  assert(f3 == 7);
  assert(f4 == 0);
  assert(c1 == '$');
  assert(z1 == 0);
  assert(z2 == 0);
  assert(z3 == 0);
  assert(z4 == 0);
  assert(c2 == '$');
  assert(strcmp(header_test, ">header1") == 0);
  assert(cont == 0);

 if(fgets(s, 300, test_file) == NULL)
    assert(0);
  s[strlen(s) - 1] = '\0';

  sscanf(s, "%d %d %d %d %c %d %d %d %d %c %s %d", &f1, &f2, &f3, &f4, &c1, &z1, &z2, &z3, &z4, &c2, header_test, &cont);

  assert(f1 == 2);
  assert(f2 == 1);
  assert(f3 == 1);
  assert(f4 == 1);
  assert(c1 == '$');
  assert(z1 == 0);
  assert(z2 == 0);
  assert(z3 == 0);
  assert(z4 == 0);
  assert(c2 == '$');
  assert(strcmp(header_test, ">header2") == 0);
  assert(cont == 0);

 if(fgets(s, 300, test_file) == NULL)
    assert(0);
  s[strlen(s) - 1] = '\0';

  sscanf(s, "%d %d %d %d %c %d %d %d %d %c %s %d", &f1, &f2, &f3, &f4, &c1, &z1, &z2, &z3, &z4, &c2, header_test, &cont);

  assert(f1 == 1);
  assert(f2 == 1);
  assert(f3 == 1);
  assert(f4 == 1);
  assert(c1 == '$');
  assert(z1 == 0);
  assert(z2 == 0);
  assert(z3 == 0);
  assert(z4 == 0);
  assert(c2 == '$');
  assert(strcmp(header_test, ">header2") == 0);
  assert(cont == 2);

 if(fgets(s, 300, test_file) == NULL)
    assert(0);
  s[strlen(s) - 1] = '\0';

  sscanf(s, "%d %d %d %d %c %d %d %d %d %c %s %d", &f1, &f2, &f3, &f4, &c1, &z1, &z2, &z3, &z4, &c2, header_test, &cont);

  assert(f1 == 1);
  assert(f2 == 1);
  assert(f3 == 1);
  assert(f4 == 1);
  assert(c1 == '$');
  assert(z1 == 0);
  assert(z2 == 0);
  assert(z3 == 0);
  assert(z4 == 0);
  assert(c2 == '$');
  assert(strcmp(header_test, ">header2") == 0);
  assert(cont == 3);

 if(fgets(s, 300, test_file) == NULL)
    assert(0);
  s[strlen(s) - 1] = '\0';

  sscanf(s, "%d %d %d %d %c %d %d %d %d %c %s %d", &f1, &f2, &f3, &f4, &c1, &z1, &z2, &z3, &z4, &c2, header_test, &cont);

  assert(f1 == 1);
  assert(f2 == 1);
  assert(f3 == 1);
  assert(f4 == 10);
  assert(c1 == '$');
  assert(z1 == 0);
  assert(z2 == 0);
  assert(z3 == 0);
  assert(z4 == 0);
  assert(c2 == '$');
  assert(strcmp(header_test, ">header2") == 0);
  assert(cont == 4);

 if(fgets(s, 300, test_file) == NULL)
    assert(0);
  s[strlen(s) - 1] = '\0';

  sscanf(s, "%d %d %d %d %c %d %d %d %d %c %s %d", &f1, &f2, &f3, &f4, &c1, &z1, &z2, &z3, &z4, &c2, header_test, &cont);

  assert(f1 == 2);
  assert(f2 == 2);
  assert(f3 == 1);
  assert(f4 == 1);
  assert(c1 == '$');
  assert(z1 == 0);
  assert(z2 == 0);
  assert(z3 == 0);
  assert(z4 == 0);
  assert(c2 == '$');
  assert(strcmp(header_test, ">header3") == 0);
  assert(cont == 0);

  if(fgets(s, 300, test_file) == NULL)
    assert(0);
  s[strlen(s) - 1] = '\0';

  sscanf(s, "%d %d %d %d %c %d %d %d %d %c %s %d", &f1, &f2, &f3, &f4, &c1, &z1, &z2, &z3, &z4, &c2, header_test, &cont);

  assert(f1 == 2);
  assert(f2 == 1);
  assert(f3 == 1);
  assert(f4 == 1);
  assert(c1 == '$');
  assert(z1 == 0);
  assert(z2 == 0);
  assert(z3 == 0);
  assert(z4 == 0);
  assert(c2 == '$');
  assert(strcmp(header_test, ">header3") == 0);
  assert(cont == 2);

  if(fgets(s, 300, test_file) == NULL)
    assert(0);
  s[strlen(s) - 1] = '\0';

  sscanf(s, "%d %d %d %d %c %d %d %d %d %c %s %d", &f1, &f2, &f3, &f4, &c1, &z1, &z2, &z3, &z4, &c2, header_test, &cont);

  assert(f1 == 1);
  assert(f2 == 1);
  assert(f3 == 1);
  assert(f4 == 7);
  assert(c1 == '$');
  assert(z1 == 0);
  assert(z2 == 0);
  assert(z3 == 0);
  assert(z4 == 0);
  assert(c2 == '$');
  assert(strcmp(header_test, ">header3") == 0);
  assert(cont == 4);

 if(fgets(s, 300, test_file) != NULL)
    assert(0);

  fclose(test_file);

  test_file = fopen("test_file-oneformat", "r");
  if (test_file == NULL) {
    printf("test couldn't be completed cause test_file-oneformat.txt cannot be opened in reading mode\n");
    printf("error: %s\n", strerror(errno));
    remove("test_file-factorization");
    remove("test_file-fingerprint");
    remove("test_file-kfingerprint");
    remove("test_file-oneformat");
    remove("test_file.fasta");
    exit(1);
  }

  if(fgets(s, 300, test_file) == NULL)
    assert(0);
  s[strlen(s) - 1] = '\0';
  assert(strcmp(s, ">header1 $ 5,3,7 $ [ \"CGTTG\" \"CGG\" \"AAAGGTC\" ]") == 0);

  if(fgets(s, 300, test_file) == NULL)
    assert(0);
  s[strlen(s) - 1] = '\0';
  assert(strcmp(s, ">header2 $ 2,1,1,1,1,1,10 $ [ \"GT\" \"C\" \"C\" \"C\" \"C\" \"C\" \"AAAAGGGCTC\" ]") == 0);

  if(fgets(s, 300, test_file) == NULL)
    assert(0);
  s[strlen(s) - 1] = '\0';
  assert(strcmp(s, ">header3 $ 2,2,1,1,1,7 $ [ \"GT\" \"CT\" \"C\" \"C\" \"C\" \"ACCTCAG\" ]") == 0);

  fclose(test_file);

  remove("test_file-factorization");
  remove("test_file-fingerprint");
  remove("test_file-kfingerprint");
  remove("test_file-oneformat");
  remove("test_file.fasta");

  printf("test_case_passed\n");

  printf("test conditions: processing fasta with 3 reads, creting the filesystem_node before launching the thread\n");

   if (getcwd(path_test, 255) == NULL) {
    printf("process_fasta test cannot be run cause it hasn't been possible to find the current path\n");
    printf("error: %s\n", strerror(errno));
    exit(1);
  }

  test_file = fopen("test_file.fasta", "w");
  if (test_file == NULL) {
    printf("test couldn't be completed cause test_file cannot be opened in reading mode\n");
    remove("test_file");
    exit(1);
  }

  fprintf(test_file, "%s\n%s\n", ">header1", "CGTTGCGGAAAGGTC");
  fprintf(test_file, "%s\n%s\n", ">header2", "GTCCCCCAAAAGGGCTC");
  fprintf(test_file, "%s\n%s\n", ">header3", "GTCTCCCACCTCAG");

  fclose(test_file);

  file = opendir(path_test);
  if (path_test == NULL) {
    printf("process_file test could not be completed cause directory cannot be opened\n");
    exit(1);
  }

  inner_file = readdir(file);
  if (inner_file == NULL) {
    printf("process file test could not be completed cause directory cannot be read\n");
    exit(1);
  }

  while (strcmp(inner_file->d_name, "test_file.fasta") != 0){
    inner_file = readdir(file);
    if (inner_file == NULL) {
      printf("process file test could not be completed cause directory cannot be read\n");
      exit(1);
    }
  }
  closedir(file);

  end_of_processing = 0;
  window_dimension = 4;
  strcpy(files_to_process[0].directory_path, path_test);
  strcpy(files_to_process[0].filename, "test_file.fasta");
  strcat(path_test, "/");
  strcat(path_test, "test_file.fasta");
  files_to_process[0].deleted = 0;

  printf("%s\n", files_to_process[0].directory_path);
  printf("%s\n", files_to_process[0].filename);

  communicate_max_fact_length(0);
  fact_choice = 1;

  pthread_attr_init(&att);
  pthread_create(&tid, &att, process_file, &index);
  end_of_processing = 1;
  pthread_join(tid, NULL);

  test_file = fopen("test_file-factorization", "r");
  if (test_file == NULL) {
    printf("test couldn't be completed cause test_file-factorization cannot be opened in reading mode\n");
    printf("error: %s\n", strerror(errno));
    remove("test_file-factorization");
    remove("test_file-fingerprint");
    remove("test_file-kfingerprint");
    remove("test_file-oneformat");
    remove("test_file.fasta");
    exit(1);
  }

  fgets(header_test, 300, test_file);
  fgets(genom, 300, test_file);
  header_test[strlen(header_test) - 1] = '\0';
  genom[strlen(genom) - 1] = '\0';
  assert(strcmp(header_test, ">header1") == 0);
  assert(strcmp(genom, "[ \"CGTTG\" \"CGG\" \"AAAGGTC\" ]") == 0);
  fgets(header_test, 300, test_file);
  fgets(genom, 300, test_file);
  header_test[strlen(header_test) - 1] = '\0';
  genom[strlen(genom) - 1] = '\0';
  assert(strcmp(header_test, ">header2") == 0);
  assert(strcmp(genom, "[ \"GT\" \"C\" \"C\" \"C\" \"C\" \"C\" \"AAAAGGGCTC\" ]") == 0);
  fgets(header_test, 300, test_file);
  fgets(genom, 300, test_file);
  header_test[strlen(header_test) - 1] = '\0';
  genom[strlen(genom) - 1] = '\0';
  assert(strcmp(header_test, ">header3") == 0);
  assert(strcmp(genom, "[ \"GT\" \"CT\" \"C\" \"C\" \"C\" \"ACCTCAG\" ]") == 0);

  if(fgets(s, 300, test_file) != NULL)
    assert(0);

  fclose(test_file);

  test_file = fopen("test_file-fingerprint", "r");
  if (test_file == NULL) {
    printf("test couldn't be completed cause test_file-fingerprint.txt cannot be opened in reading mode\n");
    printf("error: %s\n", strerror(errno));
    remove("test_file-factorization");
    remove("test_file-fingerprint");
    remove("test_file-kfingerprint");
    remove("test_file-oneformat");
    remove("test_file.fasta");
    exit(1);
  }

  fscanf(test_file, "%s %s\n", header_test, genom);
  assert(strcmp(header_test, ">header1") == 0);
  assert(strcmp(genom, "5,3,7") == 0);
  fscanf(test_file, "%s\n%s\n", header_test, genom);
  assert(strcmp(header_test, ">header2") == 0);
  assert(strcmp(genom, "2,1,1,1,1,1,10") == 0);
  fscanf(test_file, "%s\n%s\n", header_test, genom);
  assert(strcmp(header_test, ">header3") == 0);
  assert(strcmp(genom, "2,2,1,1,1,7") == 0);

  if(fgets(s, 300, test_file) != NULL)
    assert(0);
  fclose(test_file);

  test_file = fopen("test_file-kfingerprint", "r");
  if (test_file == NULL) {
    printf("test couldn't be completed cause test_file-kfingerprint.txt cannot be opened in reading mode\n");
    remove("test_file-factorization");
    remove("test_file-fingerprint");
    remove("test_file-kfingerprint");
    remove("test_file-oneformat");
    remove("test_file.fasta");
    exit(1);
  }

 if(fgets(s, 300, test_file) == NULL)
    assert(0);
  s[strlen(s) - 1] = '\0';

  sscanf(s, "%d %d %d %d %c %d %d %d %d %c %s %d", &f1, &f2, &f3, &f4, &c1, &z1, &z2, &z3, &z4, &c2, header_test, &cont);

  assert(f1 == 5);
  assert(f2 == 3);
  assert(f3 == 7);
  assert(f4 == 0);
  assert(c1 == '$');
  assert(z1 == 0);
  assert(z2 == 0);
  assert(z3 == 0);
  assert(z4 == 0);
  assert(c2 == '$');
  assert(strcmp(header_test, ">header1") == 0);
  assert(cont == 0);

 if(fgets(s, 300, test_file) == NULL)
    assert(0);
  s[strlen(s) - 1] = '\0';

  sscanf(s, "%d %d %d %d %c %d %d %d %d %c %s %d", &f1, &f2, &f3, &f4, &c1, &z1, &z2, &z3, &z4, &c2, header_test, &cont);

  assert(f1 == 2);
  assert(f2 == 1);
  assert(f3 == 1);
  assert(f4 == 1);
  assert(c1 == '$');
  assert(z1 == 0);
  assert(z2 == 0);
  assert(z3 == 0);
  assert(z4 == 0);
  assert(c2 == '$');
  assert(strcmp(header_test, ">header2") == 0);
  assert(cont == 0);

 if(fgets(s, 300, test_file) == NULL)
    assert(0);
  s[strlen(s) - 1] = '\0';

  sscanf(s, "%d %d %d %d %c %d %d %d %d %c %s %d", &f1, &f2, &f3, &f4, &c1, &z1, &z2, &z3, &z4, &c2, header_test, &cont);

  assert(f1 == 1);
  assert(f2 == 1);
  assert(f3 == 1);
  assert(f4 == 1);
  assert(c1 == '$');
  assert(z1 == 0);
  assert(z2 == 0);
  assert(z3 == 0);
  assert(z4 == 0);
  assert(c2 == '$');
  assert(strcmp(header_test, ">header2") == 0);
  assert(cont == 2);

 if(fgets(s, 300, test_file) == NULL)
    assert(0);
  s[strlen(s) - 1] = '\0';

  sscanf(s, "%d %d %d %d %c %d %d %d %d %c %s %d", &f1, &f2, &f3, &f4, &c1, &z1, &z2, &z3, &z4, &c2, header_test, &cont);

  assert(f1 == 1);
  assert(f2 == 1);
  assert(f3 == 1);
  assert(f4 == 1);
  assert(c1 == '$');
  assert(z1 == 0);
  assert(z2 == 0);
  assert(z3 == 0);
  assert(z4 == 0);
  assert(c2 == '$');
  assert(strcmp(header_test, ">header2") == 0);
  assert(cont == 3);

 if(fgets(s, 300, test_file) == NULL)
    assert(0);
  s[strlen(s) - 1] = '\0';

  sscanf(s, "%d %d %d %d %c %d %d %d %d %c %s %d", &f1, &f2, &f3, &f4, &c1, &z1, &z2, &z3, &z4, &c2, header_test, &cont);

  assert(f1 == 1);
  assert(f2 == 1);
  assert(f3 == 1);
  assert(f4 == 10);
  assert(c1 == '$');
  assert(z1 == 0);
  assert(z2 == 0);
  assert(z3 == 0);
  assert(z4 == 0);
  assert(c2 == '$');
  assert(strcmp(header_test, ">header2") == 0);
  assert(cont == 4);

 if(fgets(s, 300, test_file) == NULL)
    assert(0);
  s[strlen(s) - 1] = '\0';

  sscanf(s, "%d %d %d %d %c %d %d %d %d %c %s %d", &f1, &f2, &f3, &f4, &c1, &z1, &z2, &z3, &z4, &c2, header_test, &cont);

  assert(f1 == 2);
  assert(f2 == 2);
  assert(f3 == 1);
  assert(f4 == 1);
  assert(c1 == '$');
  assert(z1 == 0);
  assert(z2 == 0);
  assert(z3 == 0);
  assert(z4 == 0);
  assert(c2 == '$');
  assert(strcmp(header_test, ">header3") == 0);
  assert(cont == 0);

  if(fgets(s, 300, test_file) == NULL)
    assert(0);
  s[strlen(s) - 1] = '\0';

  sscanf(s, "%d %d %d %d %c %d %d %d %d %c %s %d", &f1, &f2, &f3, &f4, &c1, &z1, &z2, &z3, &z4, &c2, header_test, &cont);

  assert(f1 == 2);
  assert(f2 == 1);
  assert(f3 == 1);
  assert(f4 == 1);
  assert(c1 == '$');
  assert(z1 == 0);
  assert(z2 == 0);
  assert(z3 == 0);
  assert(z4 == 0);
  assert(c2 == '$');
  assert(strcmp(header_test, ">header3") == 0);
  assert(cont == 2);

  if(fgets(s, 300, test_file) == NULL)
    assert(0);
  s[strlen(s) - 1] = '\0';

  sscanf(s, "%d %d %d %d %c %d %d %d %d %c %s %d", &f1, &f2, &f3, &f4, &c1, &z1, &z2, &z3, &z4, &c2, header_test, &cont);

  assert(f1 == 1);
  assert(f2 == 1);
  assert(f3 == 1);
  assert(f4 == 7);
  assert(c1 == '$');
  assert(z1 == 0);
  assert(z2 == 0);
  assert(z3 == 0);
  assert(z4 == 0);
  assert(c2 == '$');
  assert(strcmp(header_test, ">header3") == 0);
  assert(cont == 4);

 if(fgets(s, 300, test_file) != NULL)
    assert(0);

  fclose(test_file);

  test_file = fopen("test_file-oneformat", "r");
  if (test_file == NULL) {
    printf("test couldn't be completed cause test_file-oneformat.txt cannot be opened in reading mode\n");
    printf("error: %s\n", strerror(errno));
    remove("test_file-factorization");
    remove("test_file-fingerprint");
    remove("test_file-kfingerprint");
    remove("test_file-oneformat");
    remove("test_file.fasta");
    exit(1);
  }

  if(fgets(s, 300, test_file) == NULL)
    assert(0);
  s[strlen(s) - 1] = '\0';
  assert(strcmp(s, ">header1 $ 5,3,7 $ [ \"CGTTG\" \"CGG\" \"AAAGGTC\" ]") == 0);

  if(fgets(s, 300, test_file) == NULL)
    assert(0);
  s[strlen(s) - 1] = '\0';
  assert(strcmp(s, ">header2 $ 2,1,1,1,1,1,10 $ [ \"GT\" \"C\" \"C\" \"C\" \"C\" \"C\" \"AAAAGGGCTC\" ]") == 0);

  if(fgets(s, 300, test_file) == NULL)
    assert(0);
  s[strlen(s) - 1] = '\0';
  assert(strcmp(s, ">header3 $ 2,2,1,1,1,7 $ [ \"GT\" \"CT\" \"C\" \"C\" \"C\" \"ACCTCAG\" ]") == 0);

  fclose(test_file);

  remove("test_file-factorization");
  remove("test_file-fingerprint");
  remove("test_file-kfingerprint");
  remove("test_file-oneformat");
  remove("test_file.fasta");

  printf("test_case_passed\n");

  printf("test conditions: processing fasta file with no reads\n");

  if (getcwd(path_test, 255) == NULL) {
    printf("process_fasta test cannot be run cause it hasn't been possible to find the current path\n");
    printf("error: %s\n", strerror(errno));
    exit(1);
  }

  test_file = fopen("test_file.fasta", "w");
  if (test_file == NULL) {
    printf("test couldn't be completed cause test_file cannot be opened in reading mode\n");
    remove("test_file");
    exit(1);
  }

  fclose(test_file);

  file = opendir(path_test);
  if (path_test == NULL) {
    printf("process_fasta test could not be completed cause directory cannot be opened\n");
    exit(1);
  }

  inner_file = readdir(file);
  if (inner_file == NULL) {
    printf("process fasta could not be completed cause directory cannot be read\n");
    exit(1);
  }

  while (strcmp(inner_file->d_name, "test_file.fasta") != 0){
    inner_file = readdir(file);
    if (inner_file == NULL) {
      printf("process fasta could not be completed cause directory cannot be read\n");
      exit(1);
    }
  }
  closedir(file);

  end_of_processing = 0;
  window_dimension = 4;
  strcpy(files_to_process[0].directory_path, path_test);
  strcpy(files_to_process[0].filename, "test_file.fasta");
  strcat(path_test, "/");
  strcat(path_test, "test_file.fasta");
  files_to_process[0].deleted = 0;

  communicate_max_fact_length(0);
  fact_choice = 1;

  pthread_attr_init(&att);
  pthread_create(&tid, &att, process_file, &index);
  end_of_processing = 1;
  pthread_join(tid, NULL);

  test_file = fopen("test_file-factorization", "r");
  if (test_file == NULL) {
    printf("test couldn't be completed cause test_file-factorization.txt cannot be opened in reading mode\n");
    printf("error: %s\n", strerror(errno));
    remove("test_file-factorization");
    remove("test_file-fingerprint");
    remove("test_file-kfingerprint");
    remove("test_file-oneformat");
    remove("test_file.fasta");
    exit(1);
  }

  if(fgets(s, 300, test_file) != NULL)
    assert(0);

  fclose(test_file);

  test_file = fopen("test_file-fingerprint", "r");
  if (test_file == NULL) {
    printf("test couldn't be completed cause test_file-fingerprint.txt cannot be opened in reading mode\n");
    printf("error: %s\n", strerror(errno));
    remove("test_file-factorization");
    remove("test_file-fingerprint");
    remove("test_file-kfingerprint");
    remove("test_file-oneformat");
    remove("test_file.fasta");
    exit(1);
  }

  if(fgets(s, 300, test_file) != NULL)
    assert(0);

  fclose(test_file);

  test_file = fopen("test_file-kfingerprint", "r");
  if (test_file == NULL) {
    printf("test couldn't be completed cause test_file-kfingerprint.txt cannot be opened in reading mode\n");
    remove("test_file-factorization");
    remove("test_file-fingerprint");
    remove("test_file-kfingerprint");
    remove("test_file-oneformat");
    remove("test_file.fasta");

    exit(1);
  }

  if(fgets(s, 300, test_file) != NULL)
    assert(0);

  fclose(test_file);

  test_file = fopen("test_file-oneformat", "r");
  if (test_file == NULL) {
    printf("test couldn't be completed cause test_file-oneformat.txt cannot be opened in reading mode\n");
    printf("error: %s\n", strerror(errno));
    remove("test_file-factorization");
    remove("test_file-fingerprint");
    remove("test_file-kfingerprint");
    remove("test_file-oneformat");
    remove("test_file.fasta");
    exit(1);

  }

  if(fgets(s, 300, test_file) != NULL)
    assert(0);

  fclose(test_file);

  remove("test_file-factorization");
  remove("test_file-fingerprint");
  remove("test_file-kfingerprint");
  remove("test_file-oneformat");

  remove("test_file.fasta");

  printf("test_case_passed\n");

  printf("test conditions: processing fasta for the second time\n");

  if (getcwd(path_test, 255) == NULL) {
    printf("process_fasta test cannot be run cause it hasn't been possible to find the current path\n");
    printf("error: %s\n", strerror(errno));
    exit(1);
  }

  test_file = fopen("test_file.fasta", "w");
  if (test_file == NULL) {
    printf("test couldn't be completed cause test_file cannot be opened in reading mode\n");
    remove("test_file");
    exit(1);
  }

  fprintf(test_file, "%s\n%s\n", ">header1", "CGTTGCGGAAAGGTC");
  fprintf(test_file, "%s\n%s\n", ">header2", "GTCCCCCAAAAGGGCTC");
  fprintf(test_file, "%s\n%s\n", ">header3", "GTCTCCCACCTCAG");

  fclose(test_file);

  file = opendir(path_test);
  if (path_test == NULL) {
    printf("process_fasta test could not be completed cause directory cannot be opened\n");
    exit(1);
  }

  inner_file = readdir(file);
  if (inner_file == NULL) {
    printf("process fasta could not be completed cause directory cannot be read\n");
    exit(1);
  }

  while (strcmp(inner_file->d_name, "test_file.fasta") != 0){
    inner_file = readdir(file);
    if (inner_file == NULL) {
      printf("process fasta could not be completed cause directory cannot be read\n");
      exit(1);
    }
  }
  closedir(file);

  fact_choice = 1;
  end_of_processing = 0;
  window_dimension = 4;
  strcpy(files_to_process[0].directory_path, path_test);
  strcpy(files_to_process[0].filename, "test_file.fasta");
  strcat(path_test, "/");
  strcat(path_test, "test_file.fasta");
  files_to_process[0].deleted = 0;

  communicate_max_fact_length(0);
  fact_choice = 1;

  pthread_attr_init(&att);
  pthread_create(&tid, &att, process_file, &index);
  end_of_processing = 1;
  pthread_join(tid, NULL);

  if (getcwd(path_test, 255) == NULL) {
    printf("process_fasta test cannot be run cause it hasn't been possible to find the current path\n");
    printf("error: %s\n", strerror(errno));
    exit(1);
  }

  test_file = fopen("test_file.fasta", "w");
  if (test_file == NULL) {
    printf("test couldn't be completed cause test_file cannot be opened in reading mode\n");
    remove("test_file");
    exit(1);
  }

  fprintf(test_file, "%s\n%s\n", ">header1", "CGTTGCGGAAAGGTC");
  fprintf(test_file, "%s\n%s\n", ">header2", "GTCCCCCAAAAGGGCTC");
  fprintf(test_file, "%s\n%s\n", ">header3", "GTCTCCCACCTCAG");

  fclose(test_file);

  file = opendir(path_test);
  if (path_test == NULL) {
    printf("process_fasta test could not be completed cause directory cannot be opened\n");
    exit(1);
  }

  inner_file = readdir(file);
  if (inner_file == NULL) {
    printf("process fasta could not be completed cause directory cannot be read\n");
    exit(1);
  }

  while (strcmp(inner_file->d_name, "test_file.fasta") != 0){
    inner_file = readdir(file);
    if (inner_file == NULL) {
      printf("process fasta could not be completed cause directory cannot be read\n");
      exit(1);
    }
  }
  closedir(file);

  fact_choice = 1;
  end_of_processing = 0;
  window_dimension = 4;
  strcpy(files_to_process[0].directory_path, path_test);
  strcpy(files_to_process[0].filename, "test_file.fasta");
  strcat(path_test, "/");
  strcat(path_test, "test_file.fasta");
  files_to_process[0].deleted = 0;

  communicate_max_fact_length(0);
  fact_choice = 1;

  pthread_attr_init(&att);
  pthread_create(&tid, &att, process_file, &index);
  end_of_processing = 1;
  pthread_join(tid, NULL);

  if (getcwd(path_test, 255) == NULL) {
    printf("process_fasta test cannot be run cause it hasn't been possible to find the current path\n");
    printf("error: %s\n", strerror(errno));
    exit(1);
  }

  test_file = fopen("test_file2.fasta", "w");
  if (test_file == NULL) {
    printf("test couldn't be completed cause test_file cannot be opened in reading mode\n");
    remove("test_file");
    exit(1);
  }

  fprintf(test_file, "%s\n%s\n", ">header1", "CGTTGCGGCA");
  fprintf(test_file, "%s\n%s\n", ">header2", "GTCCCC");
  fprintf(test_file, "%s\n%s\n", ">header3", "AAGTCA");

  fclose(test_file);

  file = opendir(path_test);
  if (path_test == NULL) {
    printf("process_fasta test could not be completed cause directory cannot be opened\n");
    exit(1);
  }

  inner_file = readdir(file);
  if (inner_file == NULL) {
    printf("process fasta could not be completed cause directory cannot be read\n");
    exit(1);
  }

  while (strcmp(inner_file->d_name, "test_file2.fasta") != 0){
    inner_file = readdir(file);
    if (inner_file == NULL) {
      printf("process fasta could not be completed cause directory cannot be read\n");
      exit(1);
    }
  }
  closedir(file);

  end_of_processing = 0;
  window_dimension = 4;
  strcpy(files_to_process[0].directory_path, path_test);
  strcpy(files_to_process[0].filename, "test_file2.fasta");
  strcat(path_test, "/");
  strcat(path_test, "test_file2.fasta");
  files_to_process[0].deleted = 0;

  communicate_max_fact_length(0);
  fact_choice = 1;

  pthread_attr_init(&att);
  pthread_create(&tid, &att, process_file, &index);
  end_of_processing = 1;
  pthread_join(tid, NULL);

  test_file = fopen("test_file2-factorization", "r");
  if (test_file == NULL) {
    printf("test couldn't be completed cause test_file2-factorization cannot be opened in reading mode\n");
    printf("error: %s\n", strerror(errno));
    remove("test_file2-factorization");
    remove("test_file2-fingerprint");
    remove("test_file2-kfingerprint");
    remove("test_file2-oneformat");
    remove("test_file2.fasta");
    exit(1);
  }

  fgets(header_test, 300, test_file);
  fgets(genom, 300, test_file);
  header_test[strlen(header_test) - 1] = '\0';
  genom[strlen(genom) - 1] = '\0';
  assert(strcmp(header_test, ">header1") == 0);
  assert(strcmp(genom, "[ \"CGTTG\" \"CGG\" \"C\" \"A\" ]") == 0);
  fgets(header_test, 300, test_file);
  fgets(genom, 300, test_file);
  header_test[strlen(header_test) - 1] = '\0';
  genom[strlen(genom) - 1] = '\0';
  assert(strcmp(header_test, ">header2") == 0);
  assert(strcmp(genom, "[ \"GT\" \"C\" \"C\" \"C\" \"C\" ]") == 0);
  fgets(header_test, 300, test_file);
  fgets(genom, 300, test_file);
  header_test[strlen(header_test) - 1] = '\0';
  genom[strlen(genom) - 1] = '\0';
  assert(strcmp(header_test, ">header3") == 0);
  assert(strcmp(genom, "[ \"AAGTC\" \"A\" ]") == 0);

  if(fgets(s, 300, test_file) != NULL)
    assert(0);

  fclose(test_file);

  test_file = fopen("test_file2-fingerprint", "r");
  if (test_file == NULL) {
    printf("test couldn't be completed cause test_file2-fingerprint cannot be opened in reading mode\n");
    printf("error: %s\n", strerror(errno));
    remove("test_file2-factorization");
    remove("test_file2-fingerprint");
    remove("test_file2-kfingerprint");
    remove("test_file2-oneformat");
    remove("test_file2.fasta");
    exit(1);
  }

  fscanf(test_file, "%s %s\n", header_test, genom);
  assert(strcmp(header_test, ">header1") == 0);
  assert(strcmp(genom, "5,3,1,1") == 0);
  fscanf(test_file, "%s\n%s\n", header_test, genom);
  assert(strcmp(header_test, ">header2") == 0);
  assert(strcmp(genom, "2,1,1,1,1") == 0);
  fscanf(test_file, "%s\n%s\n", header_test, genom);
  assert(strcmp(header_test, ">header3") == 0);
  assert(strcmp(genom, "5,1") == 0);

  if(fgets(s, 300, test_file) != NULL)
    assert(0);
  fclose(test_file);

  test_file = fopen("test_file2-kfingerprint", "r");
  if (test_file == NULL) {
    printf("test couldn't be completed cause test_file2-kfingerprint cannot be opened in reading mode\n");
    remove("test_file2-factorization");
    remove("test_file2-fingerprint");
    remove("test_file2-kfingerprint");
    remove("test_file2-oneformat");
    remove("test_file2.fasta");
    exit(1);
  }

 if(fgets(s, 300, test_file) == NULL)
    assert(0);
  s[strlen(s) - 1] = '\0';

  sscanf(s, "%d %d %d %d %c %d %d %d %d %c %s %d", &f1, &f2, &f3, &f4, &c1, &z1, &z2, &z3, &z4, &c2, header_test, &cont);

  assert(f1 == 5);
  assert(f2 == 3);
  assert(f3 == 1);
  assert(f4 == 1);
  assert(c1 == '$');
  assert(z1 == 0);
  assert(z2 == 0);
  assert(z3 == 0);
  assert(z4 == 0);
  assert(c2 == '$');
  assert(strcmp(header_test, ">header1") == 0);
  assert(cont == 0);

 if(fgets(s, 300, test_file) == NULL)
    assert(0);
  s[strlen(s) - 1] = '\0';

  sscanf(s, "%d %d %d %d %c %d %d %d %d %c %s %d", &f1, &f2, &f3, &f4, &c1, &z1, &z2, &z3, &z4, &c2, header_test, &cont);

  assert(f1 == 2);
  assert(f2 == 1);
  assert(f3 == 1);
  assert(f4 == 1);
  assert(c1 == '$');
  assert(z1 == 0);
  assert(z2 == 0);
  assert(z3 == 0);
  assert(z4 == 0);
  assert(c2 == '$');
  assert(strcmp(header_test, ">header2") == 0);
  assert(cont == 0);

 if(fgets(s, 300, test_file) == NULL)
    assert(0);
  s[strlen(s) - 1] = '\0';

  sscanf(s, "%d %d %d %d %c %d %d %d %d %c %s %d", &f1, &f2, &f3, &f4, &c1, &z1, &z2, &z3, &z4, &c2, header_test, &cont);

  assert(f1 == 1);
  assert(f2 == 1);
  assert(f3 == 1);
  assert(f4 == 1);
  assert(c1 == '$');
  assert(z1 == 0);
  assert(z2 == 0);
  assert(z3 == 0);
  assert(z4 == 0);
  assert(c2 == '$');
  assert(strcmp(header_test, ">header2") == 0);
  assert(cont == 2);

 if(fgets(s, 300, test_file) == NULL)
    assert(0);
  s[strlen(s) - 1] = '\0';

  sscanf(s, "%d %d %d %d %c %d %d %d %d %c %s %d", &f1, &f2, &f3, &f4, &c1, &z1, &z2, &z3, &z4, &c2, header_test, &cont);

  assert(f1 == 5);
  assert(f2 == 1);
  assert(f3 == 0);
  assert(f4 == 0);
  assert(c1 == '$');
  assert(z1 == 0);
  assert(z2 == 0);
  assert(z3 == 0);
  assert(z4 == 0);
  assert(c2 == '$');
  assert(strcmp(header_test, ">header3") == 0);
  assert(cont == 0);

 if(fgets(s, 300, test_file) != NULL)
    assert(0);

  fclose(test_file);

  test_file = fopen("test_file2-oneformat", "r");
  if (test_file == NULL) {
    printf("test couldn't be completed cause test_file2-oneformat cannot be opened in reading mode\n");
    printf("error: %s\n", strerror(errno));
    remove("test_file2-factorization");
    remove("test_file2-fingerprint");
    remove("test_file2-kfingerprint");
    remove("test_file2-oneformat");
    remove("test_file2.fasta");
    exit(1);
  }

  if(fgets(s, 300, test_file) == NULL)
    assert(0);
  s[strlen(s) - 1] = '\0';
  assert(strcmp(s, ">header1 $ 5,3,1,1 $ [ \"CGTTG\" \"CGG\" \"C\" \"A\" ]") == 0);

  if(fgets(s, 300, test_file) == NULL)
    assert(0);
  s[strlen(s) - 1] = '\0';
  assert(strcmp(s, ">header2 $ 2,1,1,1,1 $ [ \"GT\" \"C\" \"C\" \"C\" \"C\" ]") == 0);

  if(fgets(s, 300, test_file) == NULL)
    assert(0);
  s[strlen(s) - 1] = '\0';
  assert(strcmp(s, ">header3 $ 5,1 $ [ \"AAGTC\" \"A\" ]") == 0);

  fclose(test_file);

  remove("test_file-factorization");
  remove("test_file-fingerprint");
  remove("test_file-kfingerprint");
  remove("test_file-oneformat");
  remove("test_file.fasta");

  remove("test_file2-factorization");
  remove("test_file2-fingerprint");
  remove("test_file2-kfingerprint");
  remove("test_file2-oneformat");
  remove("test_file2.fasta");

  printf("test_case_passed\n");

  printf("\n\nStart of process_file test\n");
  printf("test conditions: test of state transitions\n");

   if (getcwd(path_test, 255) == NULL) {
    printf("process_fasta test cannot be run cause it hasn't been possible to find the current path\n");
    printf("error: %s\n", strerror(errno));
    exit(1);
  }

  test_file = fopen("test_file.fasta", "w");
  if (test_file == NULL) {
    printf("test couldn't be completed cause test_file cannot be opened in reading mode\n");
    remove("test_file");
    exit(1);
  }

  fclose(test_file);

  file = opendir(path_test);
  if (path_test == NULL) {
    printf("process_file test could not be completed cause directory cannot be opened\n");
    exit(1);
  }

  inner_file = readdir(file);
  if (inner_file == NULL) {
    printf("process file test could not be completed cause directory cannot be read\n");
    exit(1);
  }

  while (strcmp(inner_file->d_name, "test_file.fasta") != 0){
    inner_file = readdir(file);
    if (inner_file == NULL) {
      printf("process file test could not be completed cause directory cannot be read\n");
      exit(1);
    }
  }
  closedir(file);

  window_dimension = 4;
  strcpy(files_to_process[0].directory_path, path_test);
  strcpy(files_to_process[0].filename, "test_file.fasta");
  strcat(path_test, "/");
  strcat(path_test, "test_file.fasta");
  end_of_processing = 0;

  printf("%s\n", files_to_process[0].directory_path);
  printf("%s\n", files_to_process[0].filename);

  communicate_max_fact_length(0);
  fact_choice = 1;
  test_activate_process_file_log = 1;
  running = 0;

  pthread_attr_init(&att);
  printf("thread should be in wait mode now\n");
  pthread_create(&tid, &att, process_file, &index);
  sleep(1);  

  printf("thread should be running\n");
  files_to_process[0].deleted = 0;
  sleep(1);
  if (running)
    printf("running\n");
  printf("thread should be running at first an thad exiting\n");
  files_to_process[0].deleted = 0;
  end_of_processing = 1;  
  pthread_join(tid, NULL);


  remove("test_file-factorization");
  remove("test_file-fingerprint");
  remove("test_file-kfingerprint");
  remove("test_file-oneformat");
  remove("test_file.fasta");

  test_activate_process_file_log = 0;

  printf("test_case_passed\n");

  printf("process_fasta test passed\n");
}


void test_process_file_create_file(char* filename, char* path_test, FILE *test_file) {

  if (getcwd(path_test, 255) == NULL) {
    printf("process_fasta test cannot be run cause it hasn't been possible to find the current path\n");
    printf("error: %s\n", strerror(errno));
    exit(1);
  }

  test_file = fopen(filename, "w");
  if (test_file == NULL) {
    printf("test couldn't be completed cause test_file cannot be opened in writing mode\n");
    remove("test_file");
    exit(1);
  }

  fclose(test_file);
}

void test_set_parameters(char* path_test, struct dirent **inner_file, char *filename, DIR *file) {
  file = opendir(path_test);
  if (path_test == NULL) {
    printf("process_fasta test could not be completed cause directory cannot be opened\n");
    exit(1);
  }

  *inner_file = readdir(file);
  if (*inner_file == NULL) {
    printf("process fasta could not be completed cause directory cannot be read\n");
    printf("error: %s\n", strerror(errno));
    exit(1);
  }

  while (strcmp((*inner_file)->d_name, filename) != 0){
    *inner_file = readdir(file);
    if (inner_file == NULL) {
      printf("process fasta could not be completed cause directory cannot be read\n");
      printf("error: %s\n", strerror(errno));
      exit(1);
    }
  }
  closedir(file);

  strcat(path_test, "/");
  strcat(path_test, filename);

  window_dimension = 4;

  communicate_max_fact_length(0);
  fact_choice = 1;
  test_finish = 1;
  end_of_processing = 0;
}

void test_process_file_create_threads() {

  int i, *j;
  pthread_attr_t *attr;

  for(i = 0; i < PROCESSORS_NUMBER; i++) {
    available_threads[i] = 1;
    files_to_process[i].deleted = 1;
    available_threads[i] = 1;
    j = malloc(sizeof(int));
    *j = i;
    attr = malloc(sizeof(pthread_attr_t));
    pthread_attr_init(attr);
    pthread_create(&thread_ids[i], attr, process_file, j);
  }
}

int verify_file_exists(char* filename, char* kind, FILE *test_file) {
  char complete_filename[255];

  strcpy(complete_filename, filename);
  strcat(complete_filename, kind);
  test_file = fopen(complete_filename, "r");
  if (test_file == NULL) { 
    return 0;
  }

  fclose(test_file);
  remove(complete_filename);
  return 1;
}

int test_process_file_file_exist(char* filename, char* path_test, struct dirent *inner_file, FILE *test_file) {

  int result1 = 0, result2 = 0, result3 = 0, result4 = 0;
  sleep(2);  //creating files takes time. Waiting for a little while is enough to fix the problem
  result1 = verify_file_exists(filename, "-factorization", test_file);
  result2 = verify_file_exists(filename, "-fingerprint", test_file);
  result3 = verify_file_exists(filename, "-kfingerprint", test_file);
  result4 = verify_file_exists(filename, "-oneformat", test_file);

  remove(filename);

  return result1 && result2 && result3 && result4;
}

/*If the file descripted is a fasta file, factorization, fingerprint, kfingerprint and the one format of the first three
  will be created and saved in the same directory containing the file to be processed
  pre-condition file_description: must descript an existing file with reads that respect the correct format
  pre-condition path: must refer to the descripted file
  pre-condition: fact_choice >= 1 || fact_choice <= 4

  post-condition: given the name nam of the fasta file without ".fasta" at the end, nam-factorization, nam-fingerprint,
      nam-kfingerprint, nam-oneformat will be created in the same directory of the fasta file with the respective output inside

void process_fasta(struct dirent *file_description, char *path)
*/
void test_process_fasta() {

  FILE *test_file;
  char path_test[300];
  DIR *file;
  struct dirent *inner_file;
  char filename[255];
  struct dirent *beckup;

  printf("\n\nStart of process_fasta test\n");
  printf("test conditions: precessing of one fasta file with one available thread\n");
  
  strcpy(filename, "test_file.fasta");
  strcpy(path_test, "");
  test_file = NULL;
  inner_file = NULL;
  file = NULL;
  
  test_process_file_create_file(filename, path_test, test_file);
  test_set_parameters(path_test, &inner_file, filename, file);
  beckup = malloc(sizeof(struct dirent));
  strcpy(beckup->d_name, inner_file->d_name);
  inner_file = beckup;  //cause location pointed by inner_file changes after test_process_file_create_threads() function execution
  test_process_file_create_threads();

  printf("processing fasta\n");
  process_fasta(inner_file, path_test);
  while(test_finish > 0);
  printf("end of processing\n");
  end_of_processing = 1;
  assert(test_process_file_file_exist("test_file", path_test, inner_file, test_file));
  printf("test_case_passed\n");


  printf("test conditions: processing of a file whose name length is less then length of \".fasta\"\n");
  
  strcpy(filename, "min");
  strcpy(path_test, "");
  test_file = NULL;
  inner_file = NULL;
  file = NULL;
  
  test_process_file_create_file(filename, path_test, test_file);
  test_set_parameters(path_test, &inner_file, filename, file);
  beckup = malloc(sizeof(struct dirent));
  strcpy(beckup->d_name, inner_file->d_name);
  inner_file = beckup;  //cause location pointed by inner_file changes after test_process_file_create_threads() function execution
  test_process_file_create_threads();

  printf("processing fasta\n");
  process_fasta(inner_file, path_test);
  while(test_finish > 0);
  printf("end of processing\n");
 
  assert(!test_process_file_file_exist("min", path_test, inner_file, test_file));
  printf("test_case_passed\n");

  printf("test conditions: processing of a file whose name length is major then length of \".fasta\" but is not a fasta file\n");
  
  strcpy(filename, "longer_name");
  strcpy(path_test, "");
  test_file = NULL;
  inner_file = NULL;
  file = NULL;
  
  test_process_file_create_file(filename, path_test, test_file);
  test_set_parameters(path_test, &inner_file, filename, file);
  beckup = malloc(sizeof(struct dirent));
  strcpy(beckup->d_name, inner_file->d_name);
  inner_file = beckup;  //cause location pointed by inner_file changes after test_process_file_create_threads() function execution
  test_process_file_create_threads();

  printf("processing fasta\n");
  process_fasta(inner_file, path_test);
  while(test_finish > 0);
  printf("end of processing\n");
 
  assert(!test_process_file_file_exist("longer_name", path_test, inner_file, test_file));
  printf("test_case_passed\n");



}
