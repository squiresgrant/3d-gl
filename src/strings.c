#include "strings.h"
#include "util.h"
int str_size(str*str){
  for(int i = 0;;i++)
    if(str->str[i]=='\0')
      return i;
  return 0;
}
int ca_size(char*str){
  for(int i = 0;;i++)
    if(str[i]=='\0')
      return i;
  return 0;
}
void str_pushc(str*str,char ch){
  str->str = realloc(str->str,sizeof(str->str)*(str->len+1));
  if(str->str==NULL)
    err("failed to realloc string",pexit);
  str->str[str->len]=ch;
  str->str[str->len+1]='\0';
  str->len++;
}
void str_pushca(str*str,char*ins){
  int size = ca_size(ins);
  str->str = realloc(str->str,sizeof(str->str)*(str->len+1+size));
  if(str->str==NULL)
    err("failed to realloc string",pexit);
  for(int i = 0; i!=size; i++)
    str_pushc(str,ins[i]);
}
str* str_init(){
  str* str = malloc(sizeof(*str));
  str->str = malloc(sizeof(str->str));
  if(str->str==NULL||str==NULL)
    err("failed to alloc string",pexit);
  str->str[0] = '\0';
  str->len = 0;
  return str;
}
int str_cmp_str(str*str1,str*str2){
  if(str1->len!=str2->len)
    return 0;
  for(int i = 0; i!=str1->len;i++)
    if(str1->str[i]!=str2->str[i])
      return 0;
  return 1;
}
int str_cmp_ca(str*str,char*ca){
  int len = ca_size(ca);
  if(len!=str->len)
    return 0;
  for(int i = 0; i!=str->len;i++)
    if(str->str[i]!=ca[i])
      return 0;
  return 1;
}
void str_free(str*str){
  free(str->str);
  free(str);
}
str* str_init_only_str(str*str){
  //str* str = malloc(sizeof(*str));
  str->str = malloc(sizeof(str->str));
  if(str->str==NULL||str==NULL)
    err("failed to alloc string",pexit);
  str->str[0] = '\0';
  str->len = 0;
  return str;
}
void str_clear(str*str){
  for(int i = 0; i!=str->len; i++){
    str->str[i]='\0';
  }
  str->len=0;
}
