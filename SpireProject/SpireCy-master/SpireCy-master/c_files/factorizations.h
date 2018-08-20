#ifndef FACTORIZATIONS_H_
#define FACTORIZATIONS_H_

#include "utils.h"

int index_in_alphabet(char t, char typ_alphabet_list[]);
node_t *CFL(char word[]);
node_t *CFL_for_alphabet(char word[], char list_alphabet[]);
node_t *find_pre(char word[]);
node_t *find_pre_for_alphabet(char word[], char list_alphabet[]);
void border(char p[], int **pi);
node_t *find_bre(char *w, char *v);
node_t *find_bre_for_alphabet(char *w, char *v, char list_alphabet[]);
void compute_icfl_recursive(char word[], node_t **curr_pointer_icfl);
void compute_icfl_recursive_for_alphabet(char word[], node_t **curr_pointer_icfl, char list_alphabet[]);
node_t *ICFL_recursive(char word[]);
node_t *ICFL_recursive_for_alphabet(char word[], char list_alphabet[]);
node_t *CFL_icfl(char word[], int C);
node_t *CFL_icfl_for_alphabet(char word[], int C, char list_alphabet[]);
node_t *ICFL_cfl(char word[], int C);
node_t *ICFL_cfl_for_alphabet(char word[], int C, char list_alphabet[]);
/*
void *ICFL_thread(void *args);
void *CFL_thread(void *args);
node_t *CFL_icfl_multithread_pool(char word[], int C);
node_t *CFL_icfl_for_alphabet_multithread_pool(char word[], int C, char list_alphabet[]);
node_t *ICFL_cfl_multithread_pool(char word[], int C);
node_t *ICFL_cfl_for_alphabet_multithread_pool(char word[], int C, char list_alphabet[]);
*/
#endif /* FACTORIZATIONS_H_ */
