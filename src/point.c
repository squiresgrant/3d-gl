#include <stdio.h>
#include <stdlib.h>

#include <math.h>
#include "glfww.h"
#include "util.h"
#include <unistd.h>
#include <time.h>

typedef struct {
  double x;
  double y;
  double z;
  int vertex;
} cord;
cord poi_d(double x1,double y1,double x2, double y2,int len,GLfloat* pixels,int shortest,int ign){
	double m1 = (y2-y1)/(x2-x1);	
	double b1 = y1 - m1 * x1;
	if(x2-x1==0){
		m1=0;
		b1=0;	
	}
	cord aa;
	aa.x = 0;
	aa.y = 0;
	aa.z = -1;
		
	for(int yyu = 0; yyu!=len-1; yyu++){ 
    //if(yyu==ign||yyu-1==ign||yyu+1==ign)continue;
  	double x3 = pixels[yyu*2];
    double x4 = pixels[(yyu+1)*2];
    double y3 = pixels[yyu*2+1];
    double y4 = pixels[(yyu+1)*2+1];
    double m2 = (y4-y3)/(x4-x3);
  	
    double b2 = y3 - m2 * x3;
		
		
    double nsx = (b2-b1)/(m1-m2);
    double nsy = m1*nsx+b1;
		if (x4-x3==0){
			nsx=x3;
			nsy=m1*nsx+b1;
		}
		if(!(nsx >= greater(lesser(x1, x2), lesser(x3, x4)) && nsx <= lesser(greater(x1, x2), greater(x3, x4)))
    	||!(nsy >= greater(lesser(y1, y2), lesser(y3, y4)) && nsy <= lesser(greater(y1, y2), greater(y3, y4))))
				continue;

		if((diff(nsx,x2)<FL_DIS&&diff(nsy,y2)<FL_DIS)
			||(diff(nsx,x1)<FL_DIS&&diff(nsy,y1)<FL_DIS))
				continue;

			if(aa.z==-1||pow(x1-nsx,2)+pow(y1-nsy,2)<pow(x1-aa.x,2)+pow(y1-aa.y,2)){	
				aa.x = nsx;
				aa.y = nsy;
				
			}
			aa.z++;		
			if(!shortest)
				break; 
  }
	return aa;
}

typedef struct {
  float r;
  float g;
  float b;
} color;
typedef struct {
  cord at;
  color color;
} point;
typedef struct {
  point* c;
  point* vert;
  double vlen;
  double len;
} point_arr;
typedef struct {
	point_arr* at;
	int len;
} point_m;

typedef struct {
  GLfloat* at;
  int len;
} glfl_a;
typedef struct {
  glfl_a* at;
  int len;
} glfl_m;
GLuint prog;
static const char* vshader_src = 
  "#version 330\n"
  "layout (location = 0) in vec3 pos;\n"
  "layout (location = 1) in vec3 color;\n"
  "layout (location = 2) in float trans;\n"
	"out vec3 ncolor;\n"
  "out float ntrans;\n"
	"void main(){\n"
  "ncolor = color;\n"
  "ntrans = trans;\n"
	"gl_Position = vec4(pos,1.0);\n" 
  "};";
static const char* fshader_src = 
  "#version 330\n"
  "in vec3 ncolor;\n"
  "in float ntrans;\n"
	"out vec3 color;\n"
  "void main(){\n"
  "gl_FragColor = vec4(ncolor[0],ncolor[1],ncolor[2],0.2);\n"
  "};";

point_arr* basier2d(double*xx,double*yy,int n,float rr, float gg, float bb){ 
  
  n-=1;

  int lle = diff(get_w(),greater(yy[0]+yy[1],xx[0]+xx[1]))/2;
  lle+=1;
	double aaar = (1.0/lle);
  
  point_arr* pa;
  pa = malloc(sizeof(*pa));
  pa->c = malloc(sizeof(*pa->c)*(lle*60));
  pa->vert = malloc(sizeof(*pa->vert)*(n*60));

  if(pa==NULL||pa->c==NULL||pa->vert==NULL)
    err("failed to allocate basier array",pexit);
  
  pa->len = lle;	
	for(int iy = 0; iy<=lle;iy+=1){
    double t = aaar*iy;
    double bcx = 0;
    double bcy = 0;
    for(int i = 0; i <=n;i++){
      double pp = binomial(n,i) * pow((1 - t),(n - i)) * pow(t,i);
      bcx +=  pp * xx[i];
      bcy +=  pp * yy[i]; 
    } 
    pa->c[iy].at.x = bcx;
    pa->c[iy].at.y = bcy;
    pa->c[iy].at.vertex = 0;
    for(int as = 0; as<=n; as++){
      if(xx[as]==bcx&&yy[as]==bcy){
        pa->c[iy].at.vertex = 1;
        break;
      }
    }
    pa->c[iy].color.r = rr;
    pa->c[iy].color.g = gg;
    pa->c[iy].color.b = bb;
  }
  for(int i = 0; i<=n; i++){
    pa->vert[i].at.x = xx[i];
    pa->vert[i].at.y = yy[i]; 
  }
  pa->vlen = n; 
  if(pa==NULL||pa->c==NULL||pa->vert==NULL)
    err("failed to allocate basier array",pexit);
	return pa;
}
typedef struct {
	GLfloat* pix;
	GLfloat* col;
	GLfloat* trans;
	int len;
} glfl_ar;
void render_p(glfl_ar* bba){
	//glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC1_ALPHA);
	GLuint verta;
  glGenVertexArrays(1,&verta);
  glBindVertexArray(verta);

  GLuint vetb;
  glGenBuffers(1,&vetb);
  glBindBuffer(GL_ARRAY_BUFFER,vetb);
  glBufferData(GL_ARRAY_BUFFER,sizeof(*bba->pix)*(bba->len*3),bba->pix,GL_STATIC_DRAW);
  
  glEnableVertexAttribArray(0);
  glBindBuffer(GL_ARRAY_BUFFER,vetb);
  glVertexAttribPointer(0,2,GL_FLOAT,GL_FALSE,0,(void*)0);
  
  GLuint colb;
  glGenBuffers(1,&colb);
  glBindBuffer(GL_ARRAY_BUFFER,colb);
  glBufferData(GL_ARRAY_BUFFER,sizeof(*bba->col)*(bba->len*4),bba->col,GL_STATIC_DRAW);

  glEnableVertexAttribArray(1);
  glBindBuffer(GL_ARRAY_BUFFER,colb);
  glVertexAttribPointer(1,3,GL_FLOAT,GL_FALSE,0,(void*)0);

  GLuint trab;
  glGenBuffers(1,&trab);
  glBindBuffer(GL_ARRAY_BUFFER,trab);
  glBufferData(GL_ARRAY_BUFFER,sizeof(*bba->trans)*(bba->len*2),bba->trans,GL_STATIC_DRAW);

  glEnableVertexAttribArray(2);
  glBindBuffer(GL_ARRAY_BUFFER,trab);
  glVertexAttribPointer(2,1,GL_FLOAT,GL_FALSE,0,(void*)0);

	glDrawArrays(GL_POINTS,0,bba->len);
  glDeleteBuffers(1,&trab);
	glDeleteBuffers(1,&vetb);
  glDeleteBuffers(1,&colb);
}
glfl_ar* perspective_proj(GLFWwindow* b,point_arr* c,double ctx,double cty,double ctz,double cx, double cy, double cz){ 

	GLfloat* pixels = malloc(sizeof(*pixels)*((1+c->len)*3));
  GLfloat* colors = malloc(sizeof(*colors)*((1+c->len)*4));
  GLfloat* trans = malloc(sizeof(*trans)*((1+c->len)*2));
	if(pixels==NULL||colors==NULL||trans==NULL)
    err("failed to allocate perspective array:(",pexit);

  double coy = cos(cty);
  double sinz = sin(ctz);
  double coz = cos(ctz);
  double six = sin(ctx);
  double siy = sin(cty);
  double cox = cos(ctx);
  double fov = 0.002;
  double ex = cx;
  double ey = cy; 
  refresh_size(b); 
  GLuint fb = 0;
  int c_len = 0;
  //double ez=1/tan(fov/2); //i dont get this at all
  double ez=get_w()*2;
  glEnableClientState(GL_VERTEX_ARRAY);
  
  for(int i = 0; i!=c->len; i++){
    double ax = c->c[i].at.x;
    double ay = c->c[i].at.y;
    double az = c->c[i].at.z;
    
    double eyz = (coz*(ay-cy)-sinz*ax-cx);
    double yzm = (coy*(az-cz) + siy*(sinz*(ay-cy) + coz*(ax-cx)));
    double dx = coy * (sinz*(ay-cy) + coz*(ax-cx)) - (siy*(az-cz));
    double dy = six * yzm + cox*eyz;
    double dz = cox * yzm - six*eyz; 
   
    double bx = ez/dz*dx+dx;
    double by = ez/dz*dy+dy; 
    if(dz>-1){
      ab_to_vp(xa,ya,get_w(),get_h(),bx,by);
      pixels[c_len*2] = xa+1;
      pixels[c_len*2+1] = ya;
      colors[c_len*3] = c->c[i].color.r;
      colors[c_len*3+1] = c->c[i].color.g;
      colors[c_len*3+2] = c->c[i].color.b;
      trans[c_len] = 0.5;
			c_len++;
    }
  }

  double fc_len = c_len;
  for(int i = 0; i<=fc_len-1; i++){

		if(isinf(pixels[i*2]))
			continue;
    double x22[3] = {pixels[i*2],pixels[(i+1)*2]};
    double y22[3] = {pixels[i*2+1],pixels[(i+1)*2+1]};
    point_arr* bas = basier2d(x22,y22,2,0.1f,0.1f,0.1f);
		
		pixels = realloc(pixels,sizeof(*pixels) *((c_len+2+bas->len)*3));

		colors = realloc(colors,sizeof(*colors) *((c_len+2+bas->len)*4));

		trans = realloc(trans,sizeof(*trans) *((c_len+2+bas->len)*2));
		if(pixels==NULL||colors==NULL||trans==NULL)
					err("can't reallocate AAA",pexit);
		for(int zaa=0; zaa<=bas->len; zaa++){
     	
				pixels[c_len*2] = bas->c[zaa].at.x;
        pixels[c_len*2+1] = bas->c[zaa].at.y; 
        colors[c_len*3] = colors[i*3];
        colors[c_len*3+1] = colors[i*3+1]; 
        colors[c_len*3+2] = colors[i*3+2]; 
        trans[c_len] = 0.5;
				c_len++;
    }
		
    free(bas->c);
    free(bas->vert);
    free(bas);
  }
  
	double dclen = c_len;
  int vvi = 0;
  //printf("---\n");
  glfl_m* trline = malloc(sizeof(*trline)*get_w()*30);
  trline->len = 0;
  trline->at = malloc(sizeof(*trline->at)*get_w()*40);
  if(trline==NULL||trline->at==NULL)
		pexit(54);
	for(int i = 0; i<=fc_len; i++){
    if(c->c[i].at.vertex==1||1){ 
      if(pixels==NULL||colors==NULL)
        pexit(55);
      if(isinf(pixels[i*2]))
				continue;
			
			vvi++;
      
      trline->at[trline->len].at = malloc(sizeof(*trline->at[trline->len].at)*(c_len+get_w()*2)*30);
      trline->at[trline->len].len = 0;
			//printf("%f\n",dclen);
			int le = 2222;
			cord aaa = poi_d(pixels[i*2],pixels[i*2+1],le,pixels[i*2+1],c->len*2,pixels,1,i);
			cord aab = poi_d(pixels[i*2],pixels[i*2+1],-le,pixels[i*2+1],c->len*2,pixels,1,i);
			cord aac = poi_d(pixels[i*2],pixels[i*2+1],pixels[i*2],-le,c->len*2,pixels,1,i);
			cord aad = poi_d(pixels[i*2],pixels[i*2+1],pixels[i*2],le,c->len*2,pixels,1,i);

			/*if(aac.z||aad.z==-1){
				free(trline->at[trline->len].at);
        continue;
			}*/
			//printf("%f %f %f %f\n",aaa.z,aab.z,aac.z,aad.z);
			if(fmod(aaa.z,2)==1||fmod(aab.z,2)==1||
					fmod(aac.z,2)==1||fmod(aad.z,2)==1){
				free(trline->at[trline->len].at);
				continue;
			}	
			if(aab.z!=-1)
			{double aax[] = {pixels[i*2],aab.x};
			double aay[] = {pixels[i*2+1],aab.y};	
			point_arr* frl = basier2d(aax,aay,2,0.1,0.1,0.1);
			if(trans==NULL||pixels==NULL||colors==NULL)
        	pexit(55);
			for(int cci = 0; cci<=frl->len; cci++){
				trline->at[trline->len].at[trline->at[trline->len].len*2] = frl->c[cci].at.x;
        trline->at[trline->len].at[trline->at[trline->len].len*2+1] = frl->c[cci].at.y;
				pixels = realloc(pixels,sizeof *pixels *((c_len+1)*3)); 
      	colors = realloc(colors,sizeof *colors *((c_len+1)*4));
				trans = realloc(trans,sizeof *trans *((c_len+1)*2));
				if(trans==NULL||pixels==NULL||colors==NULL)
        	pexit(55);
				trans[c_len] = 0.5;
				pixels[c_len*2] = frl->c[cci].at.x;
      	pixels[c_len*2+1] = frl->c[cci].at.y; 
      	colors[c_len*3] = 1.0f;
      	colors[c_len*3+1] = 0.3f; 
      	colors[c_len*3+2] = 1.0f; 
      	c_len++;
				//trline->len++;	
				trline->at[trline->len].len++;	
			}
			free(frl->c);
			free(frl->vert);
			free(frl);}
			if(aaa.z!=-1)
			{double aax[] = {pixels[i*2],aaa.x};
			double aay[] = {pixels[i*2+1],aaa.y};
			//printf("a %f %f | %f %f\n",pixels[i*2],aab.x,pixels[i*2+1],aab.y);
			point_arr* frl = basier2d(aax,aay,2,0.1,0.1,0.1);
			if(trans==NULL||pixels==NULL||colors==NULL)
        	pexit(55);
			for(int cci = 0; cci<=frl->len; cci++){
				trline->at[trline->len].at[trline->at[trline->len].len*2] = frl->c[cci].at.x;
        trline->at[trline->len].at[trline->at[trline->len].len*2+1] = frl->c[cci].at.y;
				pixels = realloc(pixels,sizeof *pixels *((c_len+1)*3)); 
      	colors = realloc(colors,sizeof *colors *((c_len+1)*4));
				trans = realloc(trans,sizeof *trans *((c_len+1)*2));
				
				trans[c_len] = 0.5;
				pixels[c_len*2] = frl->c[cci].at.x;
      	pixels[c_len*2+1] = frl->c[cci].at.y; 
      	colors[c_len*3] = 1.0f;
      	colors[c_len*3+1] = 0.3f; 
      	colors[c_len*3+2] = 1.0f; 
      	c_len++;
				trline->at[trline->len].len++;
			}
			free(frl->c);
			free(frl->vert);
			free(frl);}
		trline->len++;
		}
  }
  if(trline->len>1){
  for(int ii = 0; ii!=trline->len-1; ii++){ 
    if(trline->at[ii].at[3]<trline->at[ii+1].at[3]){
      glfl_a temp = trline->at[ii];
      trline->at[ii] = trline->at[ii+1];
      trline->at[ii+1] = temp;
      ii=-1;
    }
      
  }}
  double ffclen = c_len;

  double lmax_t = -2.0;
  double lmin_t = 2.0;
  double lmax2_t = -2.0;
  double lmin2_t = 2.0;
  for(int zzi = 0; zzi!=trline->len; zzi++){
    double max_t = -2.0;
    double min_t = 2.0;
    double max2_t = -2.0;
    double min2_t = 2.0;
    for(int zzi2 = 0; zzi2!=trline->at[zzi].len;zzi2++){
      if(max_t<trline->at[zzi].at[zzi2*2]){
      	max_t = trline->at[zzi].at[zzi2*2];
      	max2_t = trline->at[zzi].at[zzi2*2+1];
      }
      if(min_t>=trline->at[zzi].at[zzi2*2]){
      	min_t = trline->at[zzi].at[zzi2*2];
      	min2_t = trline->at[zzi].at[zzi2*2+1];
      }
       
    } 
    free(trline->at[zzi].at);
		if (trline->at[zzi].len == 0)
      continue; 
    if(lmax2_t!=max2_t&&lmin_t!=2.0&&lmax_t!=-2.0){

    float color = (float)zzi/trline->len;
    float color2 = 1.0f-((float)zzi/trline->len);
    double di = pow(lmin_t-max_t,2)+pow(lmin2_t-max2_t,2);//sqrt should be used here
    double di2 = pow(lmax_t-min_t,2)+pow(lmax2_t-min2_t,2);
    double ux1, uy1, ux2, uy2;
    int ma_to_mi = 0;
		if (di > di2) {
    	ux1 = max_t;
    	uy1 = max2_t;
    	ux2 = lmin_t;
    	uy2 = lmin2_t;
			ma_to_mi=1;
    } else {
    	ux1 = min_t;
    	uy1 = min2_t;
    	ux2 = lmax_t;
    	uy2 = lmax2_t;
    }
		

    color2 = 1.0f;
    
		double x1 = ma_to_mi?max_t:min_t;
    double x2 = ma_to_mi?lmin_t:lmax_t;
		double ite = ((float)get_w()/4);
	
		double ite1 = fabs(max_t/ite);
		double ite2 = fabs(lmax_t/ite);
	
		double y1 = uy1;
    double y2 = uy2; 	
		
		for(;(ma_to_mi?x1>=x2:x1<=x2);(ma_to_mi?(x1-=ite1):(x1+=ite2))){	
		
		cord aaaa = poi_d(x1,y1,x2,y2,fc_len,pixels,0,-1);	
    
		if(aaaa.z!=-1){
			/*pixels = realloc(pixels,sizeof *pixels *((c_len+1)*3)); 
      colors = realloc(colors,sizeof *colors *((c_len+1)*4));
      //printf("%f %f\n",aaa.x,aaa.y);
			pixels[c_len*2] = x1;
      pixels[c_len*2+1] = y1; 
      colors[c_len*3] = 1.0f;
      colors[c_len*3+1] = 0.3f; 
      colors[c_len*3+2] = 1.0f; 
      c_len++;
			
			pixels = realloc(pixels,sizeof *pixels *((c_len+1)*3)); 
      colors = realloc(colors,sizeof *colors *((c_len+1)*4));
      //printf("%f %f\n",aaa.x,aaa.y);
			pixels[c_len*2] = aaaa.x;
      pixels[c_len*2+1] = aaaa.y; 
      colors[c_len*3] = 1.0f;
      colors[c_len*3+1] = 0.0f; 
      colors[c_len*3+2] = 0.0f;
			c_len++;*/
			//printf("%f %f\n",aaaa.x,aaaa.y);
			continue;
		}
		double bb[] = {x1,x2};
    double bb2[] = {y1, y2};
    point_arr* asd = basier2d(bb,bb2,2,0.1,0.1,0.1);
		//printf("aa\n");
    for(int lli = 0; lli!=asd->len; lli++){
    	pixels = realloc(pixels,sizeof *pixels *((c_len+1)*4));
    	colors = realloc(colors,sizeof *colors *((c_len+1)*5)); 
    	trans = realloc(trans,sizeof *trans *((c_len+1)*2));
			trans[c_len] = 0.5;
			double dd = 10;
    	
    	pixels[c_len*2] = asd->c[lli].at.x;
    	pixels[c_len*2+1] = asd->c[lli].at.y; 
    	colors[c_len*3] = color2;
    	colors[c_len*3+1] = color;//vvi==3?0.1f:vvi==4?0.5f:1.0f; 
    	colors[c_len*3+2] = 0.0f; 
    	c_len++;
    }
    free(asd->c);
    free(asd->vert);
    free(asd);
		break;
		}/*
		pixels = realloc(pixels,sizeof *pixels *((c_len+1)*3)); 
      colors = realloc(colors,sizeof *colors *((c_len+1)*4));
      //printf("%f %f\n",aaa.x,aaa.y);
			pixels[c_len*2] = ux1;
      pixels[c_len*2+1] = uy1; 
      colors[c_len*3] = 1.0f;
      colors[c_len*3+1] = 0.3f; 
      colors[c_len*3+2] = 1.0f; 
      c_len++;
    pixels = realloc(pixels,sizeof *pixels *((c_len+1)*3)); 
      colors = realloc(colors,sizeof *colors *((c_len+1)*4));
      //printf("%f %f\n",aaa.x,aaa.y);
			pixels[c_len*2] = ux2;
      pixels[c_len*2+1] = uy2; 
      colors[c_len*3] = 1.0f;
      colors[c_len*3+1] = 1.0f; 
      colors[c_len*3+2] = 1.0f; 
      c_len++;*/
		}
    lmax_t = max_t;
    lmin_t = min_t;
    lmax2_t = max2_t;
    lmin2_t = min2_t; 
  }
  free(trline->at);
  free(trline);
  glPointSize(4.0f);
   
  glfl_ar* rea = malloc(sizeof(*rea));
	rea->col = colors;
	rea->pix = pixels;
	rea->len = c_len;
	rea->trans = trans;
	return rea; 
}
point_arr* basier3d(double*xx,double*yy,double*zz,int n,float rr, float gg, float bb){
  point_arr* pa;
  pa = malloc(sizeof(*pa));
  pa->c = malloc(sizeof(*pa->c)*(get_w()*60));
  pa->vert = malloc(sizeof(*pa->vert)*(n*60));

  if(pa->c==NULL||pa->vert==NULL)
    err("failed to allocate basier array",pexit); 
  n-=1;

  
  double aaar = (1.0/get_w());
  double am = 2;
  pa->len = am;
  for(int iy = 0; iy<=2;iy+=1){
    double t = (1.0/am)*iy;
    double bcx = 0;
    double bcy = 0;
    double bcz = 0;
    for(int i = 0; i <=n;i++){
      double pp = binomial(n,i) * pow((1 - t),(n - i)) * pow(t,i);
      bcx +=  pp * xx[i];
      bcy +=  pp * yy[i];
      bcz +=  pp * zz[i]; 
    } 
    pa->c[iy].at.x = bcx;
    pa->c[iy].at.y = bcy;
    pa->c[iy].at.z = bcz;
    pa->c[iy].at.vertex = 0;
    for(int as = 0; as<=n; as++){
      if(xx[as]==bcx&&yy[as]==bcy&&zz[as]==bcz){
        pa->c[iy].at.vertex = 1;
        break;
      }
    }
    pa->c[iy].color.r = rr;
    pa->c[iy].color.g = gg;
    pa->c[iy].color.b = bb;
  }
  for(int i = 0; i<=n; i++){
    pa->vert[i].at.x = xx[i];
    pa->vert[i].at.y = yy[i]; 
    pa->vert[i].at.z = zz[i]; 
  }
  pa->vlen = n; 
  return pa;
}
void join_cords(point_arr* a, point_arr* b){
  int a_len = a->len;
  int a_vlen = a->vlen;
  a->c = realloc(a->c,sizeof(*a->c)*(a->len+b->len)*60);
  a->vert = realloc(a->vert,sizeof(*a->vert)*(a->vlen+b->vlen)*60);
  a->len+=b->len;
  a->vlen+=b->vlen;
  if(a->c==NULL||a->vert==NULL)
    err("failed to reallocate cords",pexit);
  for(int i = 0; i<=b->len; i++){
    a->c[a_len+i].at = b->c[i].at; 
    a->c[a_len+i].color = b->c[i].color; 
  }
  for(int i = 0; i<=b->vlen; i++){
    a->vert[a_vlen+i].at = b->vert[i].at;  
  }
}
void join_glfl_a(glfl_ar* a, glfl_ar* b){ 
	int a_len = a->len;
  a->pix = realloc(a->pix,sizeof(*a->pix)*(a->len+b->len+1)*20);
  a->col = realloc(a->col,sizeof(*a->col)*(a->len+b->len+1)*20);
  a->trans = realloc(a->trans,sizeof(*a->trans)*(a->len+b->len+1)*20);
	a->len+=b->len;
  
  if(a->pix==NULL||a->col==NULL||a->trans==NULL)
    err("failed to reallocate float array",pexit);
  for(int i = 0; i<=b->len*2; i++){
    a->pix[a_len*2+i] = b->pix[i];  
  }
	for(int i = 0; i<=b->len*3; i++){
    a->col[a_len*3+i] = b->col[i];  
  }
	for(int i = 0; i<=b->len; i++){
    a->trans[a_len+i] = b->trans[i];  
  }
}
point_arr* polygon3d(double* vx, double*vy, double* vz, int n){
  /*double xx[2] = {vx[0],vx[1]};
  double yy[2] = {vy[0],vy[1]};
  double zz[2] = {vz[0],vz[1]};
  point_arr* y = basier3d(xx,yy,zz,2,1.0f,0.0f,0.5f);
  for(int i = 1; i<=n-2; i++){
    double xx1[2] = {vx[i],vx[i+1]};
    double yy1[2] = {vy[i],vy[i+1]};
    double zz1[2] = {vz[i],vz[i+1]};
    point_arr* aa = basier3d(xx1,yy1,zz1,2,1.0f,0.0f,0.5f);
    join_cords(y,aa);
    free(aa->c);
    free(aa->vert);
    free(aa); 
  }*/
	point_arr* pa;
  pa = malloc(sizeof(*pa));
  pa->c = malloc(sizeof(*pa->c)*(get_w()*60));
  pa->vert = malloc(sizeof(*pa->vert)*(n*60));

  if(pa->c==NULL||pa->vert==NULL)
    err("failed to allocate polygon array",pexit); 
  n-=1;
	for(int i = 0; i<=n; i++){
		pa->c[i].at.x = vx[i];
		pa->c[i].at.y = vy[i];	
		pa->c[i].at.z = vz[i];

		pa->vert[i].at.x = vx[i];
		pa->vert[i].at.y = vy[i];	
		pa->vert[i].at.z = vz[i];
		
		pa->c[i].color.r = 0.1f;
		pa->c[i].color.g = 0.1f;
		pa->c[i].color.b = 1.0f;
	}
	pa->len = n;
	pa->vlen= n;
  return pa;
}
point_arr* square_gen(double* tl, double* tr, double* bl, double*br,float rr, float gg, float bb){
  double xx[3] = {tl[0],tr[0]};
  double yy[3] = {tl[1],tr[1]};
  double zz[3] = {tl[2],tr[2]};
  point_arr* a = polygon3d(xx,yy,zz,2);
  
  double xx1[3] = {tl[0],bl[0]};
  double yy1[3] = {tl[1],bl[1]};
  double zz1[3] = {tl[2],bl[2]};
  point_arr* b = polygon3d(xx1,yy1,zz1,2);
  
  double xx2[3] = {tr[0],br[0]};
  double yy2[3] = {tr[1],br[1]};
  double zz2[3] = {tr[2],br[2]};
  point_arr* c = polygon3d(xx2,yy2,zz2,2);
  
  double xx3[3] = {bl[0],br[0]};
  double yy3[3] = {bl[1],br[1]};
  double zz3[3] = {bl[2],br[2]};
  point_arr* d = polygon3d(xx3,yy3,zz3,2);
  
  join_cords(a,b);
  join_cords(a,c);
  join_cords(a,d); 
	free(b->c);
  free(b->vert);
  free(b);
  free(c->c);
  free(c->vert);
  free(c);
  free(d->c);
  free(d->vert);
  free(d);
  return a;
}
void fill3d(point_arr* a){
  warn("please dont use this lol");
  double m_x = 0.0;
  double m_y = 0.0;
  double m_z = 0.0;
  double mi_x = 0.0;
  double mi_y = 0.0;
  double mi_z = 0.0;
  for(int i = 0; i<=a->len; i++){
    if(a->c[i].at.x>m_x)
      m_x = a->c[i].at.x;
    if(a->c[i].at.y>m_y)
      m_y = a->c[i].at.y; 
    if(a->c[i].at.z>m_z)
      m_z = a->c[i].at.z; 
    if(a->c[i].at.x<mi_x||mi_x==0)
      mi_x = a->c[i].at.x;
    if(a->c[i].at.y<mi_y||mi_y==0)
      mi_y = a->c[i].at.y; 
    if(a->c[i].at.z<mi_z||mi_z==0)
      mi_z = a->c[i].at.z; 
  } 
  int a_l = a->len;
  a->c = realloc(a->c,sizeof(*a->c)*(m_x*m_y+a->len)*60);
  for(double y = mi_y; y<=m_y; y+=.1){
    
    double* zz = malloc(sizeof(*zz));
    int zzl = 0;
    double* xx = malloc(sizeof(*xx));
    double* yy = malloc(sizeof(*yy));
    for(int i = 0; i<=a_l; i++){
      
      if(diff(a->c[i].at.y,y)<.1){
        zz = realloc(zz,sizeof(*zz)*(zzl+1));
        zz[zzl] = a->c[i].at.z;
        xx = realloc(xx,sizeof(*xx)*(zzl+1));
        xx[zzl] = a->c[i].at.x;
        yy = realloc(yy,sizeof(*yy)*(zzl+1));
        yy[zzl] = y;
        zzl++;
      }
    }
    point_arr* pp = basier3d(xx,yy,zz,zzl,1.0,0.0,0.0);
    join_cords(a,pp);
    free(pp->c);
    free(pp);
    free(xx);
    free(yy);
    free(zz);
  }
}
point_arr* cube_gen(double* tl, double* tr, double* bl, double*br,
               double* tl2, double* tr2, double* bl2, double*br2,
                    float rr, float gg, float bb){
   
  
	point_arr* a = square_gen(tl,tr,bl,br,rr,gg,bb); 
  point_arr* b = square_gen(tl2,tr2,bl2,br2,rr,gg,bb); 
  double s;
  join_cords(a,b);
  free(b->c);
  free(b->vert);
  free(b);
  point_arr* c = square_gen(tl2,tr2,tl,tr,rr,gg,bb);
  join_cords(a,c);
  
  free(c->c);
  free(c->vert);
  free(c);
  point_arr* d = square_gen(bl2,br2,bl,br,rr,gg,bb);
  join_cords(a,d);
  free(d->c);
  free(d->vert);
  free(d);
  point_arr* e = square_gen(tl,tl2,bl,bl2,rr,gg,bb);
  join_cords(a,e);
  free(e->c);
  free(e->vert);
  free(e);
  point_arr* f = square_gen(br2,br,tr2,tr,rr,gg,bb);
  join_cords(a,f);
  free(f->c);
  free(f->vert);
  free(f);
  return a;
}
point_m* rect3d_gen(double* tl, double* tr, double* bl, double*br,
               double* tl2, double* tr2, double* bl2, double*br2,
                    float rr, float gg, float bb){
	point_m* mm = malloc(sizeof * mm * 8);
	mm->len = 0;
	{
	double xx1[5]={tl[0],tr[0],br[0],bl[0]};
	double yy1[5]={tl[1],tr[1],br[1],bl[1]};
	double zz1[5]={tl[2],tr[2],br[2],bl[2]};
	mm[mm->len].at = polygon3d(xx1,yy1,zz1,5);
	mm->len++;
	}
	{
	double xx1[5]={tl2[0],tr2[0],br2[0],bl2[0]};
	double yy1[5]={tl2[1],tr2[1],br2[1],bl2[1]};
	double zz1[5]={tl2[2],tr2[2],br2[2],bl2[2]};
	mm[mm->len].at = polygon3d(xx1,yy1,zz1,5);
	mm->len++;
	}
	
	{
	double xx1[5]={tl2[0],tr2[0],tr[0],tl[0]};
	double yy1[5]={tl2[1],tr2[1],tr[1],tl[1]};
	double zz1[5]={tl2[2],tr2[2],tr[2],tl[2]};
	mm[mm->len].at = polygon3d(xx1,yy1,zz1,5);
	mm->len++;
	}
	{
	double xx1[5]={bl2[0],br2[0],br[0],bl[0]};
	double yy1[5]={bl2[1],br2[1],br[1],bl[1]};
	double zz1[5]={bl2[2],br2[2],br[2],bl[2]};
	mm[mm->len].at = polygon3d(xx1,yy1,zz1,5);
	mm->len++;
	}

	{
	double xx1[5]={tl2[0],bl2[0],bl[0],tl[0]};
	double yy1[5]={tl2[1],bl2[1],bl[1],tl[1]};
	double zz1[5]={tl2[2],bl2[2],bl[2],tl[2]};
	mm[mm->len].at = polygon3d(xx1,yy1,zz1,5);
	mm->len++;
	}
	{
	double xx1[5]={tr2[0],br2[0],br[0],tr[0]};
	double yy1[5]={tr2[1],br2[1],br[1],tr[1]};
	double zz1[5]={tr2[2],br2[2],br[2],tr[2]};
	mm[mm->len].at = polygon3d(xx1,yy1,zz1,5);
	mm->len++;
	}
	return mm;
}
int main(int argc,char*argv[]){
	
  flag_handle(argc,argv);
  atexit(sig_handle);
  GLFWwindow* w = glfw_init();
  refresh_size(w);
  GLenum err = glewInit();
  if (err != GLEW_OK)
    pexit(1); 
  if (!GLEW_VERSION_2_1)  
    pexit(1);
  GLuint vid = vshader_comp(vshader_src);
  GLuint fid = fshader_comp(fshader_src);
  prog = build_shader(vid,fid);
  glUseProgram(prog); 
  logm("built shaders");
  
  

  double tl[3] = {5.0,200.0,200.0};
  double tr[3] = {200.0,200.0,200.0};
  double bl[3] = {5.0,5.0,200.0};
  double br[3] = {200.0,5.0,200.0};
  
  double tl2[3] = {5.0,200.0,5.0};
  double tr2[3] = {200.0,200.0,5.0};
  double bl2[3] = {5.0,5.0,5.0};
  double br2[3] = {200.0,5.0,5.0};
	float rr = 0.0;
	float gg = 0.0;
	float bb = 1.0;
	point_m* aaaa = rect3d_gen(tl,tr,bl,br,tl2,tr2,bl2,br2,rr,gg,bb);
	
	int max_r = 630;
  double half_max_r = (double)max_r/2/2;
  double pl_x = 0;
  double pl_y = 0;
  double pl_z = 0;
  double plr_x = 0;
  double plr_y = 0;
  clock_t t;
  double frames = 0;
  t = clock();
  double frame_limit = 60;
	double single_frame = (float)1/frame_limit*1E+6;
	char sf_time_msg[60];
	sprintf(sf_time_msg,"%.3fµs per frame ≈ %.3ffps",single_frame,frame_limit);
	info(sf_time_msg);
	
	struct timespec start_t, end_t,tem_t2, tem_t;
	
	clock_gettime(CLOCK_REALTIME, &start_t);
	info("entering main loop");
	for(;;){	
	 	clock_gettime(CLOCK_REALTIME,&tem_t);
    double p1 = plr_x*0.01; 
    double p2 = plr_y*0.01;
    double p3 = 0;
    double p4 = pl_y;
    double p5 = -pl_y+pl_z;
    double p6 = pl_x;
    /*
    glfl_ar* bba = perspective_proj(w,e,p1,p2,p3,p4,p5,p6);
		glfl_ar* bbb = perspective_proj(w,b,p1,p2,p3,p4,p5,p6);	
		glfl_ar* bbc = perspective_proj(w,c,p1,p2,p3,p4,p5,p6);	
		glfl_ar* bbd = perspective_proj(w,d,p1,p2,p3,p4,p5,p6);
		glfl_ar* bbe = perspective_proj(w,e,p1,p2,p3,p4,p5,p6);	
		glfl_ar* bbf = perspective_proj(w,f,p1,p2,p3,p4,p5,p6);
		join_glfl_a(bba,bbb);
		join_glfl_a(bba,bbc);
		join_glfl_a(bba,bbd);
		join_glfl_a(bba,bbe);
		join_glfl_a(bba,bbf);
		*/

		if(aaaa->len>=0||1){
		glfl_ar* bba = perspective_proj(w,aaaa[0].at,p1,p2,p3,p4,p5,p6);
		for(int i = 1; i<=aaaa->len-1; i++){

			glfl_ar* bbb = perspective_proj(w,aaaa[i].at,p1,p2,p3,p4,p5,p6);

			join_glfl_a(bba,bbb);

			free(bbb->col);
			free(bbb->pix);
			free(bbb->trans);
			free(bbb);
		}

		render_p(bba);

		free(bba->trans);

		free(bba->col);

		free(bba->pix);

		free(bba);
		}

		glfw_load(w);

    int mod_move=2;
    double run_mul=2;
    glfwPollEvents(); 
    if(glfwGetKey(w,GLFW_KEY_R)){
      pl_x = 0;
      pl_y = 0;
      pl_z = 0;
      plr_x = 0;
      plr_y = 0;
    }
    if(glfwGetKey(w,GLFW_KEY_P))
      printf("(x:%f,y:%f,z:%f),rot(x:%f,y:%f,z:%f)||l(p1:%f,p2:%f,p3:%f,p4:%f,p5:%f,p6:%f)\n",pl_x,pl_y,pl_z,plr_x,plr_y,0.0,p1,p2,p3,p4,p5,p6);
    if(glfwGetKey(w,GLFW_KEY_I)){
      if(glfwGetKey(w,GLFW_KEY_LEFT_SHIFT))
        plr_x-=mod_move;
      else 
        plr_x--;
    plr_x = fmod(plr_x,max_r);
    }
    if(glfwGetKey(w,GLFW_KEY_K)){
      if(glfwGetKey(w,GLFW_KEY_LEFT_SHIFT))
        plr_x+=mod_move;
      else 
        plr_x++;
    plr_x = fmod(plr_x,max_r);
    } 
    if(glfwGetKey(w,GLFW_KEY_J)){
      if(glfwGetKey(w,GLFW_KEY_LEFT_SHIFT))
        plr_y-=mod_move;
      else 
        plr_y--;
    plr_y = fmod(plr_y,max_r);
    } 
    if(glfwGetKey(w,GLFW_KEY_L)){
      if(glfwGetKey(w,GLFW_KEY_LEFT_SHIFT))
        plr_y+=mod_move;
      else 
        plr_y++;
    plr_y = fmod(plr_y,max_r);
    }
    
    if(glfwGetKey(w,GLFW_KEY_W)){
      double mul = 1;
      if(glfwGetKey(w,GLFW_KEY_LEFT_SHIFT))
        mul = run_mul;
      pl_x+=cosf(plr_y*0.01)*mul;
      pl_y+=sinf(plr_y*0.01)*mul; 
      pl_z-=sinf(plr_x*0.01)*mul; 
    }
    if(glfwGetKey(w,GLFW_KEY_S)){
      double mul = 1;
      if(glfwGetKey(w,GLFW_KEY_LEFT_SHIFT))
        mul = run_mul;
      pl_x-=cosf(plr_y*0.01)*mul;
      pl_y-=sinf(plr_y*0.01)*mul; 
      pl_z+=sinf(plr_x*0.01)*mul; 
    } 
    if(glfwGetKey(w,GLFW_KEY_A)){
      double mul = 1;
      if(glfwGetKey(w,GLFW_KEY_LEFT_SHIFT))
        mul = run_mul;
      pl_x-=cosf((half_max_r+plr_y)*0.01)*mul;
      pl_y-=sinf((half_max_r+plr_y)*0.01)*mul; 
    }
    if(glfwGetKey(w,GLFW_KEY_D)){
      double mul = 1;
      if(glfwGetKey(w,GLFW_KEY_LEFT_SHIFT))
        mul = run_mul;
      pl_x+=cosf((half_max_r+plr_y)*0.01)*mul;
      pl_y+=sinf((half_max_r+plr_y)*0.01)*mul; 
    }
    if(glfwWindowShouldClose(w)||(glfwGetKey(w,GLFW_KEY_Q)))break; 
		glfw_clear(w); 
		if(1){
		if(single_frame>0){
	 		clock_gettime(CLOCK_REALTIME,&tem_t2);
			double tt2_diff = (tem_t2.tv_sec - tem_t.tv_sec) * 1000.0 +
                    (tem_t2.tv_nsec - tem_t.tv_nsec) / 1000000.0; 
			usleep(greater(single_frame - tt2_diff *1000,1));
		}
		
    
    
    frames+=1;
		clock_gettime(CLOCK_REALTIME, &end_t);
		double tt_diff = (end_t.tv_sec - start_t.tv_sec) * 1000.0 +
                    (end_t.tv_nsec - start_t.tv_nsec) / 1000000.0;
    if(tt_diff>1000){
			char fpsc_dis[40];
      sprintf(fpsc_dis,"%f fps",(double)(tt_diff/1000*frame_limit)*0.5+frames*0.5);
      logm(fpsc_dis);
			frames=0;
      clock_gettime(CLOCK_REALTIME,&start_t);
    }
	}
	}/*
  free(a->c);
  free(a->vert);
  free(a);
  */
	for(int i = 0; i!=aaaa->len; i++){
		free(aaaa[i].at->c);
		free(aaaa[i].at->vert);
		free(aaaa[i].at);
	}
	free(aaaa);
	glfwDestroyWindow(w);
  win_clean();
  glDeleteShader(vid);
  glDeleteShader(fid);
  glDeleteShader(prog);
  logm("killed window:p");
  return 0;	
}
