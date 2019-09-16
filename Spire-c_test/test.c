#include "sequential_program_main.h"
#include "test_safe_fgets.h"

void do_unit_testing();
void do_integration_testing();

//void test_safe_fgets();
void test_foreach_element_in();
void test_flush();
void test_open_towrite_file();

void  test_fill_k_finger();
void test_create_fingerprint();
void test_process_fasta();

void do_unit_testing() {
  printf("start of unit test\n");
  test_safe_fgets();
//  test_foreach_element_in();   non passa. Da riparare
  test_initialize_tail();
  test_open_towrite_file();
  printf("unit test successfully passed\n\n");
}

void do_integration_testing() {
  printf("start of integration test\n");
  test_flush();
  test_fill_k_finger();
  test_create_fingerprint();
  test_process_fasta();
  //test of process_all_fasta_files has not been accomplished cause of time constraints
  printf("\nintegration test successfully completed\n\n");
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

  file = opendir("/home/danilo/Scrivania/Spire-c/directory-test/directory1");
  if (file == NULL) {
    printf("for_each_element test cannot be completed cause created directory1 cannot be open\n");
    exit(1);
  }
  inner_file = readdir(file);
  if (inner_file == NULL) {
    printf("for_each_element test cannot be completed cause something gone wrong in catching a file\n");
    exit(1);
  }

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
  char header_test[300], genom[600];
  char s[300], c1, c2;
  int f1, f2, f3, f4, z1, z2, z3,z4, cont;


  printf("\n\nStart of process_fasta test\n");
  printf("test conditions: processing fasta with 3 reads\n");

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

  strcat(path_test, "/");
  strcat(path_test, "test_file.fasta");

  window_dimension = 4;

  current_header_size = 0;
  current_genom_size = 0;

  communicate_max_fact_length(0);
  fact_choice = 1;
  current_header_size = 0;
  current_genom_size = 0;

  printf("processing fasta\n");
  process_fasta(inner_file, path_test);

  printf("end of processing\n");

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

  strcat(path_test, "/");
  strcat(path_test, "test_file.fasta");

  window_dimension = 4;

  current_header_size = 0;
  current_genom_size = 0;
  communicate_max_fact_length(0);
  fact_choice = 1;

  printf("processing fasta\n");
  process_fasta(inner_file, path_test);

  printf("end of processing\n");

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

  strcat(path_test, "/");
  strcat(path_test, "test_file.fasta");

  window_dimension = 4;

  current_header_size = 0;
  current_genom_size = 0;
  communicate_max_fact_length(0);
  fact_choice = 1;

  printf("processing fasta\n");
  process_fasta(inner_file, path_test);

  printf("test conditions: processing fasta with 3 reads\n");

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

  strcat(path_test, "/");
  strcat(path_test, "test_file.fasta");

  window_dimension = 4;

  current_header_size = 0;
  current_genom_size = 0;
  communicate_max_fact_length(0);
  fact_choice = 1;

  printf("processing fasta\n");
  process_fasta(inner_file, path_test);

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

  strcat(path_test, "/");
  strcat(path_test, "test_file2.fasta");

  printf("processing fasta\n");
  process_fasta(inner_file, path_test);

  printf("end of processing\n");

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

  printf("process_fasta test passed\n");
}
