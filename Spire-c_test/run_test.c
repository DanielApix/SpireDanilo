#include "test.h"
#include <stdio.h>

int main() {
  do_unit_testing();
  do_integration_testing();
  printf("All tests passed\n");
  return 0;
}
