#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#ifndef __util__
#define __util__
#include "string.h"
#define greater(a,b) (a>b?a:b)
#define lesser(a,b) (a>b?b:a)
void err_m(char*,void (*)(int),char*,int);
void warn_m(char*,char*,int ,...);
void info_m(char*,char*,int ,...);
void log_m(char*ca,char*f,int l,...);
void flag_handle(int argc,char* argv[]);
void sig_handle(void);
unsigned int_len(const unsigned n);
char* force_ca_length(char*inp,int len);
void pexit(int s);

#define err(s,f,...) err_m(s,f,__FILE__,__LINE__,##__VA_ARGS__);
#define warn(s) warn_m(s,__FILE__,__LINE__);
#define info(s) info_m(s,__FILE__,__LINE__);
#define logm(s) log_m(s,__FILE__,__LINE__);

#endif
