#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include "string.h"
#ifndef __util__
#define __util__

static const double FL_DIS = 1e-7;
static const double DB_DIS = 1e-6;
static const double CLOSE_ENOUGH = 1e-4;
static const double NaN = 0.0f/0.0f;
static const int forced_length = 20;

#define greater(a,b) ((a)>(b)?(a):(b))
#define lesser(a,b) ((a)>(b)?(b):(a))
#define diff(a,b) ((a)>(b)?(a)-(b):(b)-(a))

#if defined(memory_trace) || defined(memory_count)
#ifdef memory_trace 
//maybe add warning here
#endif
#define malloc(X) mmalloc(X,(char*)__FILE__,(int)__LINE__,(char*)__FUNCTION__);
#define free(X) ffree(X,(char*)__FILE__,(int)__LINE__,(char*)__FUNCTION__); X=NULL;
#endif

#ifndef stfu

#ifdef __debug 
#define debug(s) debug_m(s,__FILE__,__LINE__);
#else 
#define debug(s){};
#endif

#define err(s,f,...) err_m(s,f,__FILE__,__LINE__,##__VA_ARGS__);
#define warn(s) warn_m(s,__FILE__,__LINE__);
#define info(s) info_m(s,__FILE__,__LINE__);
#define logm(s) log_m(s,__FILE__,__LINE__);
#else 
#define printf(...){};
#define err(s,f,...){};
#define warn(s){};
#define info(s){};
#define logm(s){};
#endif 

const double binomial(int n, int k);
void* mmalloc(size_t,char*,int,char*);
void ffree(void*,char*,int,char*);
void err_m(char*,void (*)(int),char*,int);
void debug_m(char*,char*,int,...);
void warn_m(char*,char*,int ,...);
void info_m(char*,char*,int ,...);
void log_m(char*ca,char*f,int l,...);
void flag_handle(int argc,char* argv[]);
void sig_handle(void);
unsigned int_len(const unsigned n);
char* force_ca_length(char*inp,int len);
void pexit(int s);

#endif
