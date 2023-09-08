#include<stdio.h>

static FILE *fp;
static int is_open = 0;


double PrintMe(double val) {
  if (!is_open) {
    fp = fopen("output.csv", "w");  // Open file for writing
    if (fp == NULL) {
      printf("Error opening file!\n");
      return 1/0;
    }
    is_open = 1;
  }
  fwrite(&val, sizeof(double), 1, fp);
  return 0.0;
}
