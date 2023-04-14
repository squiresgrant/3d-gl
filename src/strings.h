#include <stdio.h>
#include <stdlib.h>
#ifndef __str__
#define __str__
//#define init_str(buffer) (buffer); \
    {buffer = (char*)malloc(sizeof(char));\
    buffer[0]='\0';}
//#define str_size(str,sizeb) (sizeb); \
        for(sizeb=0;;sizeb++)if(str[sizeb]=='\0')break;
//#define str_push(str,ch)\
    int size_y = str_size(str);\
    str = (char*)realloc(str,sizeof(char*)*(size_y+1));\
    str[size_y]=ch;\
    str[size_y+1]='\0';
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
