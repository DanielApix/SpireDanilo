#ifndef UTILS_H_
#define UTILS_H_

//#include <pthread.h>

typedef struct node {
    char *factor;
    struct node *next;
} node_t;

typedef struct param {
	char *w;
	char *list_alphabet;
	node_t *start_d;
	node_t *end_d;
} params;


int check_word_and_alphabet(char word[], char list_alphabet[]);
void free_list(node_t *head);
void print_list_reverse(node_t *node);
void print_list(node_t *node);
char *substring(char word[], int x, int y);
char *list_to_string(node_t *list, int reverse);

#endif /* UTILS_H_ */
