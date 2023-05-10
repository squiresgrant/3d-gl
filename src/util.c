#include <stdio.h>
#include <stdlib.h>
#include "util.h"
#include <math.h>
#include "strings.h"
ulong allocs = 0;
int forced_length = 15;
typedef struct {
	void* addr;
	char* file;
	char* function;
	int line;
	ulong size;
} alloc;
alloc* allocations = NULL;
double binomial(int n, int k){
  if(n==k)
    return 1.0;
  double v = 1.0; 
  for(int i = 1; i<=k; i++){
    v=v*((float)(n+1-i)/i);
  } 
  return v;
}
void* mmalloc(ulong X,char*file,int line,char*func){
	if(allocations==NULL){
		allocations=(malloc)(sizeof(*allocations)*2);
	}
	allocations=(realloc)(allocations,sizeof(*allocations)*(allocs+1));
	void* mal = (malloc)(X);
	allocations[allocs].addr = mal;
	allocations[allocs].function = func;
	allocations[allocs].file = file;
	allocations[allocs].line = line;
	allocations[allocs].size = X;
	allocs++;
	return mal;
}

void ffree(void* X,char*file,int line,char*func){
	for(ulong i = 0; i<=allocs; i++){
		if(allocations[i].addr==X){
			allocations[i].addr = NULL;
			break;
		}
	}
	(free)(X);
  allocs--;
}
int log_level = 0;
int __signal = 0;
void pexit(int s){
  __signal = s;
  exit(s);
}
void sig_handle(void){
	if(log_level<=-1){
		#ifndef skip_memory_trace 
		(free)(allocations);
		#endif
		return;
	}
	
	#ifdef stfu
	#ifndef skip_memory_trace 
	(free)(allocations);
	#endif
	return;
	#endif

	#ifndef skip_memory_trace
  if(allocs>0){
    char ssa[45];
    sprintf(ssa,"%s | (found %i)","uneven allocations, memory leak(s)",(int)nearbyint(allocs));
    warn(ssa);
		for(ulong i = 0; i<=allocs; i++){
			if(allocations[i].addr!=NULL){
				
				char ad[50];
				sprintf(ad,"%p",allocations[i].addr);
				ad[0]='\20';
				ad[1]='\20';
				printf("   | - <\x1b[90m0x\x1b[0m%s> %s:%s:%i, %lu bytes initially allocated\n",ad,allocations[i].file,allocations[i].function,allocations[i].line,allocations[i].size);
			}
		}
		(free)(allocations);
  }
  if(allocs==0)
    info("even allocations, no internal leaks");
	#endif
	if(__signal==0){
    printf("\x1b[90mexited with \x1b[32m\x1b[1msignal [ %i ] \x1b[0m\x1b[90mgraceful exit\x1b[0m (meow)\n",__signal); 
  } else if(__signal>0){
    printf("\x1b[90mexited with \x1b[31m\x1b[1msignal [ %i ] \x1b[0m\x1b[90mgraceful exit\x1b[0m\n",__signal);
    //extra cleanup if needed
  } else {
    printf("\x1b[90mexited with \x1b[31m\x1b[1msignal [ %i ] \x1b[0m\x1b[90mnon-graceful exit\x1b[0m\n",__signal);
  }
}
unsigned int_len(const unsigned n) {
    if (n < 10) return 1;
    return 1 + int_len(n / 10);
}//https://stackoverflow.com/a/3068415
char* force_ca_length(char*inp,int len){
  char* nya = malloc(sizeof(*nya)*(len+1));
  int skip = 0;
  for(int i = 0;; i++){
    
    if((inp[i]=='\0'||skip)&&i>=len)
      break;
    if(inp[i]=='\0')
      skip=1;
    if(i>=len){
      for(int y = 1; y<=len; y++)
        nya[y-1] = nya[y];
      nya[len-1] = inp[i];
      continue;
    }
    if(skip)
      nya[i] = ' '; 
    else
      nya[i] = inp[i];
  }
  if(!skip){
    nya[2] = '.';
    nya[1] = '.';
    nya[0] = '.';
  }
  nya[len+1]='\0';
  return nya;
}
void err_m(char*ca,void (*cb)(int),char*f,int l){
  if(log_level>-1){ 
   int len = ca_size(f) + int_len(l); 
   char nn[len]; 
   sprintf(nn,"%s:%i",f,l);
   char* aa = force_ca_length(nn,forced_length);
   printf("\x1b[90m%s \x1b[31m[ !err ]\x1b[0m %s\n",aa,ca);
   free(aa); 
  } 
  cb(1);
}
void warn_m(char*ca,char*f,int l,...){
  if(log_level==-1)
    return;
  int len = ca_size(f) + int_len(l); 
  char nn[len]; 
  sprintf(nn,"%s:%i",f,l);
  char* aa = force_ca_length(nn,forced_length);
  printf("\x1b[90m%s \x1b[33m[ warn ]\x1b[0m %s\n",aa,ca);
  free(aa);
}
void info_m(char*ca,char*f,int l,...){
  if(log_level<1)
    return;
  int len = ca_size(f) + int_len(l); 
  char nn[len]; 
  sprintf(nn,"%s:%i",f,l);
  char* aa = force_ca_length(nn,forced_length);
  printf("\x1b[90m%s [ info ] %s\x1b[0m\n",aa,ca);
  free(aa);
}
void log_m(char*ca,char*f,int l,...){
  if(log_level<2)
    return;
  int len = ca_size(f) + int_len(l); 
  char nn[len]; 
  sprintf(nn,"%s:%i",f,l);
  char* aa = force_ca_length(nn,forced_length);
  printf("\x1b[35m%s [ log  ] \x1b[90m%s\x1b[0m\n",aa,ca);
  free(aa);
}
void flag_handle(int argc,char* argv[]){
  for(int i = 0; i<argc;i++){
    if(argv[i][0]=='-'){   
      //partial
      for(int y = 1;;y++){
        if(argv[i][y]=='\0')
          break;
        switch(argv[i][y]){
          case 'q':
            log_level-=1;
          break;
          case 'd':case 'v':
            log_level+=1;
          break;	
        }
      }
    }
  }
}
