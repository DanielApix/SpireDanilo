#ifndef _sequenzial_program_h
#define _sequenzial_program_h
#include <stdio.h>
#include <string.h>
#include <dirent.h>
#include <errno.h>
#include <stdlib.h>
#include <time.h>
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


extern struct stat st;

extern int fact_choice;

extern int max_fact_length; //arbitrary chosen and requested to the user


extern char *header_read;  //refers to the current read
extern char *genom_read;
extern int current_header_size;
extern int current_genom_size;

extern int window_dimension;

extern time_t time_spent_to_read_file;
extern time_t time_spent_to_write_in_file;
extern FILE *log_file;

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
char* safe_fgets(char* buffer, int *current_size, FILE *stream);

#endif
