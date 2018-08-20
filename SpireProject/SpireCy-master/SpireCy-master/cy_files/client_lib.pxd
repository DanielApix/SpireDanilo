cdef extern from "../c_files/utils.h":

    struct node:
        char * factor
        node * next
    ctypedef node node_t

    void free_list(node_t *head)

    void print_list_reverse(node_t *node)

    char *list_to_string(node_t *list, int reverse)

    void print_statistics()

    void communicate_max_fact_length(int c)

cdef extern from "../c_files/factorizations.h":

    int index_in_alphabet(char t, char typ_alphabet_list[])

    node_t *CFL(char word[])

    node_t *CFL_icfl(char word[], int C)

    node_t *ICFL_recursive(char word[]);

    node_t *ICFL_cfl(char word[], int C);
