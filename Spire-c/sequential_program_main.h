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

extern int fact_choice;
extern int max_fact_length;
extern int window_dimension;
extern FILE *log_file;

void communicate_max_fact_length(int c);
void for_each_element_in(char* directory_path,  void (*apply_function) (struct dirent *, char *));
void print_statistics();
void process_all_fasta_files(struct dirent *subdirectory_description, char* current_path);
char* safe_fgets(char* buffer, int *current_size, FILE *stream);

#endif

