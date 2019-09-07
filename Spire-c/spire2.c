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


int fact_choice;

char *root_path;     //...of the directory to process
void set_max_fact_length(int fn_max_fact_length);


char *header_read;  //refers to the current read
char *genom_read;

int window_dimension;


FILE* get_log_file();
/*necessary to set the dimension of the string returned by list_to_string (efficiency reasons)*/
int get_number_of_factors();
int get_number_of_delimeters();
void set_number_of_elements(int num);
void set_read_dimension(int value);
void communicate_max_fact_length(int c);
void print_statistics();

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
void process_fasta(struct dirent *file_description, char *directory_path);
void process_all_fasta_files(struct dirent *subdirectory_description, char* current_path);

void test_safe_fgets();
void test_foreach_element_in();
void test_flush();
void test_open_towrite_file();
void do_unit_testing();

void  test_fill_k_finger();
void test_create_fingerprint();
void test_process_fasta();
void do_integration_testing();

void set_window_dimension(int d);

int main() {

  int max_fact_length = 0;

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
  set_max_fact_length(max_fact_length);
  communicate_max_fact_length(max_fact_length);
  printf("fornisca il numero di elementi per ciascuna finestra per le k-fingerprint\n");
  scanf("%d", &window_dimension);
  set_window_dimension(window_dimension);
  time_t m;
  time_t now = time(NULL);

  for_each_element_in(root_path, process_all_fasta_files);

  print_statistics();

  m = difftime(time(NULL), now);
  printf("tempo totale in secondi: %ld\n",m);

  if (get_log_file()!= NULL)
    fclose(get_log_file());
  return 0;
}

