#include "sequential_program_main.h"
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

char *inputString(FILE* fp, size_t size, char ending_character);

int main() {
   char *root_path;     //...of the directory to process

// do_unit_testing();
// do_integration_testing();
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

  for_each_element_in(root_path, process_all_fasta_files);

  print_statistics();

  m = difftime(time(NULL), now);
  printf("tempo totale in secondi: %ld\n",m);

  if (log_file != NULL)
    fclose(log_file);
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
