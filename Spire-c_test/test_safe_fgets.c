#include "sequential_program_main.h"

void test_buffer_expansion();
void dispose_resources();
void test_second_reading();
void test_first_reading();
void test_when_buffer_size_different_from_current_size();
void write_some_mock_data_in_file();
void test_when_buffer_size_is_current_size();
FILE *open_file(char *mode);

FILE *test_file;
char *buffer;
int current_size;

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
  printf("\n\nstart of test_safe_fgets test\n");
  write_some_mock_data_in_file();
  test_when_buffer_size_is_current_size();
  test_when_buffer_size_different_from_current_size();
  printf("safe_fgets test passed\n");

  remove("test_file.txt");
}

void test_when_buffer_size_different_from_current_size() {

  test_file = open_file("r");

  printf("test conditions: dimension of buffer < dimension of the string to be read\n");
  test_buffer_expansion();
  dispose_resources();
}

void test_buffer_expansion() {
  printf("test case: buffer is dimension current_size = 5 and the sentence to be read si of dimension 26\n");

  buffer = (char *) malloc(5 * sizeof(char));
  current_size = 5;
  buffer = safe_fgets(buffer, &current_size, test_file);
  assert(strcmp(buffer, "this is the first sentence") == 0);

  printf("test case passed\n");

}


void test_when_buffer_size_is_current_size() {

  test_file = open_file("r");

  printf("test conditions: buffer dimension and variable current_size are synchronized\n");
  test_first_reading();
  test_second_reading();
  dispose_resources();

}

void test_second_reading() {
  printf("test case: buffer is of dimension current_size = 26 and second line is going to be read\n");
  buffer = safe_fgets(buffer, &current_size, test_file);
  assert(strcmp(buffer, "this is the second one") == 0);

  printf("test case passed\n");
}

void test_first_reading() {

  printf("test case: buffer is NULL and current_size is 0\n");

  buffer = NULL;
  current_size = 0;

  buffer = safe_fgets(buffer, &current_size, test_file);
  assert(strcmp(buffer, "this is the first sentence") == 0);

  printf("test case passed\n");
}

FILE *open_file(char *mode) {
  test_file = fopen("test_file.txt", mode);
  if (test_file == NULL) {
    printf("test_safe_fgets: test couldn't be executed cause test file cannot be opened in writing mode\n");
    printf("error: %s\n", strerror(errno));
    exit(1);
  }
  return test_file;
}


void write_some_mock_data_in_file() {
  test_file = open_file("w");
  fprintf(test_file, "%s\n%s\n%s", "this is the first sentence", "this is the second one", "this is the third");
  fclose(test_file);
}

void dispose_resources() {
  free(buffer);
  fclose(test_file);
}
