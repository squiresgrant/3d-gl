#include <stdio.h>
#include <stdlib.h>
#ifndef __str__
#define __str__

typedef struct {
  char*str;
  int len;
}str;
str* str_init();
int str_size(str*str);
int ca_size(char*str);
void str_pushc(str*str,char ch);
void str_pushca(str*str,char*ins);
int str_cmp_str(str*,str*);
int str_cmp_ca(str*,char*);
void str_free(str*str);
str* str_init_only_str(str*str);
void str_clear(str*str);
#endif
