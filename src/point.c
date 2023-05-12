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
int point_on_line(double x,double y,double lx1,double ly1, double lx2,double ly2){
	double crossproduct = (y-ly1)*(lx2-lx1) - (x-lx1)*(ly2-ly1);
	if(fabs(crossproduct)>FL_DIS)
		return 0;

	double dotproduct = (x-lx1)*(lx2-lx1) + (y-ly1)*(ly2-ly1);
	if(dotproduct < 0)
		return 0;

	double squaredlength = (lx2-lx1)*(lx2-lx1) + (ly2-ly1)*(ly2-ly1);
	if(dotproduct > squaredlength)
		return 0;
	return 1;
}
cord poi_d(double x1,double y1,double x2, double y2,int len,GLfloat* pixels,int shortest,int ign){
	
	cord aa;
	aa.x = 0;
	aa.y = 0;
	aa.z = -1;
	
	for(int yyu = 0; yyu!=len-1; yyu++){ 
    //if(yyu==ign)continue;
  	double x3 = pixels[yyu*2];
    double x4 = pixels[(yyu+1)*2];
    double y3 = pixels[yyu*2+1];
    double y4 = pixels[(yyu+1)*2+1];

		//formulas from https://en.wikipedia.org/wiki/Intersection_(geometry)
		double t = ((y3-y1)*(x2-x1)-(y2-y1)*x3+(y2-y1)*x1)/
								((y2-y1)*(x4-x3)-(y4-y3)*(x2-x1));
		double s = (x3-x1+t*(x4-x3))/(x2-x1);
		double nsx = x1 + s*(x2-x1);
		double nsy = y1 + s*(y2-y1);
		
		double nsx2 = x3 + t*(x4-x3);
		double nsy2 = y3 + t*(y4-y3);

		if(isnan(nsx)||isnan(nsy)||!(0<=s&&1>=s&&0<=t&&1>=t))
			continue;
		//printf("%f %f, %f %f | (%f,%f),(%f,%f)\n",nsx,nsy,nsx2,nsy2,x1,y1,x2,y2);		
		if(!(nsx>=lesser(x1,x2)&&nsx<=greater(x2,x1)&&nsy>=lesser(y1,y2)&&nsy<=greater(y2,y1))||(diff(nsx,x1)<FL_DIS&&diff(nsy,y1)<FL_DIS))
			continue;
		//printf("aaa\n");

		if(aa.z==-1||pow(nsx-x1,2)+pow(nsy-y1,2)<pow(aa.x-x1,2)+pow(aa.y-y1,2)){
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
	"out vec4 color;\n"
  "void main(){\n"
  //"gl_FragColor = vec4(1.0,0.0,1.0,1.0);\n"
  "gl_FragColor = vec4(ncolor,1.0);\n"
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
	GLfloat* tricol;
	int len;
	GLfloat* tri;
	int tlen;
	//double depth;
	cord max;
	cord min;
} glfl_ar;
void render_p(glfl_ar* bba){
	//glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC1_ALPHA);
	glEnable(GL_BLEND);
	glBlendFunc(GL_ONE, GL_ONE_MINUS_SRC_ALPHA);
	
	GLuint verta;
  glGenVertexArrays(1,&verta);
  glBindVertexArray(verta);
/*	
  GLuint vetb;
  glGenBuffers(1,&vetb);
  glBindBuffer(GL_ARRAY_BUFFER,vetb);
  glBufferData(GL_ARRAY_BUFFER,sizeof(*bba->pix)*(bba->len*3),bba->pix,GL_STATIC_DRAW);
  
  glEnableVertexAttribArray(0);
  glBindBuffer(GL_ARRAY_BUFFER,vetb);
  glVertexAttribPointer(0,2,GL_FLOAT,GL_FALSE,0,(void*)0);
 */ 
  GLuint vetb;
  glGenBuffers(1,&vetb);
  glBindBuffer(GL_ARRAY_BUFFER,vetb);
  glBufferData(GL_ARRAY_BUFFER,sizeof(*bba->tri)*(bba->tlen*12),bba->tri,GL_STATIC_DRAW);
  
  glEnableVertexAttribArray(0);
  glBindBuffer(GL_ARRAY_BUFFER,vetb);
  glVertexAttribPointer(0,2,GL_FLOAT,GL_FALSE,0,(void*)0);
	/*
	GLuint colb;
  glGenBuffers(1,&colb);
  glBindBuffer(GL_ARRAY_BUFFER,colb);
  glBufferData(GL_ARRAY_BUFFER,sizeof(*bba->col)*(bba->len*4),bba->col,GL_STATIC_DRAW);

  glEnableVertexAttribArray(1);
  glBindBuffer(GL_ARRAY_BUFFER,colb);
  glVertexAttribPointer(1,3,GL_FLOAT,GL_FALSE,0,(void*)0);
	*/
  
	GLuint colb;
  glGenBuffers(1,&colb);
  glBindBuffer(GL_ARRAY_BUFFER,colb);
  glBufferData(GL_ARRAY_BUFFER,sizeof(*bba->tricol)*(bba->tlen*30),bba->tricol,GL_STATIC_DRAW);

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

	glDrawArrays(GL_TRIANGLES,0,bba->tlen*3);
  glDeleteBuffers(1,&trab);
	glDeleteBuffers(1,&vetb);
  glDeleteBuffers(1,&colb);
}
glfl_ar* perspective_proj(GLFWwindow* b,point_arr* c,double ctx,double cty,double ctz,double cx, double cy, double cz){ 

	GLfloat* pixels = malloc(sizeof(*pixels)*((1+c->len)*3));
  GLfloat* colors = malloc(sizeof(*colors)*((1+c->len)*4));
  GLfloat* trans 	= malloc(sizeof(*trans)*((1+c->len)*2));
	if(pixels==NULL||colors==NULL||trans==NULL)
    err("failed to allocate perspective array:(",pexit);
	double dep = 0.0;
	double depy = 0.0;
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
  double ez=get_w();
	//printf("---\n");
  glEnableClientState(GL_VERTEX_ARRAY);
  cord maxf;
	cord minf;
	maxf.x = -INFINITY;
  maxf.y = -INFINITY;
	maxf.z = -INFINITY;
	
	minf.x = INFINITY;
  minf.y = INFINITY;
	minf.z = INFINITY;
	
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
		maxf.z = greater(maxf.z,dz);
		minf.z = lesser(minf.z,dz);
		
		maxf.x = greater(maxf.x,bx);
		minf.x = lesser(minf.x,bx);

		maxf.y = greater(maxf.y,by);
		minf.y = lesser(minf.y,by);
		if(dz>=0){
			//printf("-- %f\n",dz);
			//if(dz>dep) dep = dz;
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
	//printf("%f\n",dep);
  int fc_len = c_len;
  /*
	for(int i = 0; 0&&i<=fc_len-1; i++){
		//printf("%f %f | %f %f\n",pixels[i*2],pixels[i*2+1],pixels[(i+1)*2],pixels[(i+1)*2+1]);
		if(isinf(pixels[i*2])||isinf(pixels[i*2+1])||
				isinf(pixels[(i+1)*2])||isinf(pixels[(i+1)*2+1]))
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
  */
	//double dclen = c_len;
  int vvi = 0;
	
  glfl_m* trline = malloc(sizeof(*trline));
  trline->len = 0;
  trline->at = malloc(sizeof(*trline->at)*fc_len);
  if(trline==NULL||trline->at==NULL)
		pexit(54);
	for(int i = 0; i<=fc_len-1; i++){
    if(c->c[i].at.vertex==1||1){ 
      if(pixels==NULL||colors==NULL)
        pexit(55);
      if(isinf(pixels[i*2]))
				continue;
			//printf("%f %f\n",pixels[i*2],pixels[i*2+1]);
			vvi++;
      
      
			//printf(" // %f %f\n",pixels[i*2],pixels[i*2+1]);
			int le = 2222;
			int lentt = fc_len;
			cord aaa = poi_d(pixels[i*2],pixels[i*2+1],le,pixels[i*2+1],lentt,pixels,1,i);
			cord aab = poi_d(pixels[i*2],pixels[i*2+1],-le,pixels[i*2+1],lentt,pixels,1,i);
			cord aac = poi_d(pixels[i*2],pixels[i*2+1],pixels[i*2],-le,lentt,pixels,1,i);
			cord aad = poi_d(pixels[i*2],pixels[i*2+1],pixels[i*2],le,lentt,pixels,1,i);
			trline->at[trline->len].at = malloc(sizeof(*trline->at[trline->len].at)*((1+c_len+get_w())*2)*20);
      trline->at[trline->len].len = 0;
			//printf("%f\n",dclen);
			if(trline->at[trline->len].at==NULL)
				abort();
			/*if(aac.z||aad.z==-1){
				free(trline->at[trline->len].at);
        continue;
			}*/	
			//printf("%.1f %.1f %.1f %.1f | %f %f, %f %f\n",aaa.z,aab.z,aac.z,aad.z,aaa.x,aaa.y,aab.x,aab.y);
			if(fmod(aaa.z,2)==1||fmod(aab.z,2)==1||
					fmod(aac.z,2)==1||fmod(aad.z,2)==1){
				//printf("%f %f %f %f\n",aaa.z,aab.z,aac.z,aad.z);
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
				//printf("aab %f %f\n",frl->c[cci].at.x,frl->c[cci].at.y);
				/*pixels = realloc(pixels,sizeof *pixels *((c_len+1)*3)); 
      	colors = realloc(colors,sizeof *colors *((c_len+1)*4));
				trans = realloc(trans,sizeof *trans *((c_len+1)*2));
				if(trans==NULL||pixels==NULL||colors==NULL)
        	pexit(55);
				trans[c_len] = 0.0;
				pixels[c_len*2] = frl->c[cci].at.x;
      	pixels[c_len*2+1] = frl->c[cci].at.y; 
      	colors[c_len*3] = 1.0f;
      	colors[c_len*3+1] = 0.3f; 
      	colors[c_len*3+2] = 1.0f; 
      	c_len++;
				*/
				//trline->len++;	
				trline->at[trline->len].len++;	
			}
			free(frl->c);
			free(frl->vert);
			free(frl);
			}
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
				//printf("aab %f %f\n",frl->c[cci].at.x,frl->c[cci].at.y);
				/*pixels = realloc(pixels,sizeof *pixels *((c_len+1)*3)); 
      	colors = realloc(colors,sizeof *colors *((c_len+1)*4));
				trans = realloc(trans,sizeof *trans *((c_len+1)*2));
				
				trans[c_len] = 0.0;
				pixels[c_len*2] = frl->c[cci].at.x;
      	pixels[c_len*2+1] = frl->c[cci].at.y; 
      	colors[c_len*3] = 1.0f;
      	colors[c_len*3+1] = 0.3f; 
      	colors[c_len*3+2] = 1.0f; 
      	c_len++;
				*/
				trline->at[trline->len].len++;
			}
			free(frl->c);
			free(frl->vert);
			free(frl);
			}
			if(aaa.z==-1&&aab.z==-1){
				double fx = pixels[i*2];
				double fy = pixels[i*2+1];
				trline->at[trline->len].at[trline->at[trline->len].len*2] = fx;
        trline->at[trline->len].at[trline->at[trline->len].len*2+1] = fy;
				//printf("aa\n");
				/*
				pixels = realloc(pixels,sizeof *pixels *((c_len+1)*3)); 
      	colors = realloc(colors,sizeof *colors *((c_len+1)*4));
				trans = realloc(trans,sizeof *trans *((c_len+1)*2));
				if(trans==NULL||pixels==NULL||colors==NULL)
        	pexit(55);
				trans[c_len] = 0.0;
				pixels[c_len*2] = fx;
      	pixels[c_len*2+1] = fy; 
      	colors[c_len*3] = 1.0f;
      	colors[c_len*3+1] = 0.3f; 
      	colors[c_len*3+2] = 1.0f; 
      	c_len++;*/
				
				trline->at[trline->len].len++;	
			}
		trline->len++;
		}
  }
	
  if(trline->len>1){
  for(int ii = 0; ii<trline->len-1; ii++){ 
    if(trline->at[ii].at[1]<trline->at[ii+1].at[1]){
      glfl_a temp = trline->at[ii];
      trline->at[ii] = trline->at[ii+1];
      trline->at[ii+1] = temp;
      ii=-1;
    }
      
  }}

	/*for(int ii = 0; ii<trline->len-1; ii++){ 
    printf("%i : %f %f\n",ii,trline->at[ii].at[0],trline->at[ii].at[1]);
      
  }*/
	//printf("%i\n",trline->len);
  GLfloat* tria = malloc(sizeof(*tria)*fc_len*2*6);
	GLfloat* tricol = malloc(sizeof(*tricol)*fc_len*3*9);
	int tric = 0;
	//double ffclen = c_len;
  double lmax_t = -INFINITY;
  double lmin_t = INFINITY;
  double lmax2_t = -INFINITY;
  double lmin2_t = INFINITY;
  for(int zzi = 0; zzi<trline->len; zzi++){
    double max_t = -INFINITY;
    double min_t = INFINITY;
    double max2_t = -INFINITY;
    double min2_t = INFINITY;
    for(int zzi2 = 0; zzi2<trline->at[zzi].len;zzi2++){
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
		
    if(lmax2_t!=max2_t&&lmin_t!=INFINITY&&lmax_t!=-INFINITY){
		
    float color = (float)zzi/trline->len;
    float color2 = 1.0f-((float)zzi/trline->len);
    double di = pow(lmin_t-max_t,2)+pow(lmin2_t-max2_t,2);//sqrt should be used here
    double di2 = pow(lmax_t-min_t,2)+pow(lmax2_t-min2_t,2);
    double ux1, uy1, ux2, uy2, ix1, iy1, ix2, iy2;
    int ma_to_mi = 0;
		if (di >= di2) {
    	ux1 = max_t;
    	uy1 = max2_t;
    	ux2 = lmin_t;
    	uy2 = lmin2_t;
			
			ix1 = min_t;
    	iy1 = min2_t;
    	ix2 = lmax_t;
    	iy2 = lmax2_t;
			ma_to_mi=1;
    } else {
    	ux1 = min_t;
    	uy1 = min2_t;
    	ux2 = lmax_t;
    	uy2 = lmax2_t;
   		
			ix1 = max_t;
    	iy1 = max2_t;
    	ix2 = lmin_t;
    	iy2 = lmin2_t;

		}
	
		//printf("%i %f,%f %f,%f\n",trline->at[zzi].len,ux1,uy1,ux2,uy2);
    color2 = 1.0f;
    
		double x1 = ux1;//ma_to_mi?max_t:min_t;
    double x2 = ux2;//ma_to_mi?lmin_t:lmax_t;
		
		//double itx = get_w();
		//double ity = get_w();
		double y1 = uy1;
    double y2 = uy2;

		//double chx = (x2-x1)/itx;
		//double chy = (y2-y1)/ity;
		int it = 0;
		//printf("%f->%f by %f\n%f->%f by %f\n",x1,x2,chx,y1,y2,chy);
		double xx1[] = {x1,x2};
		double yy1[] = {y1,y2};

		point_arr* aacc = basier2d(xx1,yy1,2,0.0,0.1,0.1);
		for(int paaa = 0; paaa<=aacc->len; paaa++){
			cord aaaa;
			aaaa = poi_d(x1,y1,aacc->c[paaa].at.x,aacc->c[paaa].at.y,fc_len,pixels,0,-1);
			//x1+=chx;
			//y1+=chy;
			
			if(trline->at[zzi].len<it)
				break;
			it++;
			//printf("%i\n",it);
			if(aaaa.z!=-1){	
				continue;
			}
			/*
			double bb[] = {x2,aacc->c[paaa].at.x};
    	double bb2[] = {y2, aacc->c[paaa].at.y};
			//printf("%f,%f -> %f,%f\n",bb[0],bb2[0],bb[1],bb2[1]);
    	point_arr* asd = basier2d(bb,bb2,2,0.1,0.1,0.1);
			//printf("aa\n");
    	for(int lli = 0; lli!=asd->len; lli++){
    	pixels = realloc(pixels,sizeof *pixels *((c_len+1)*4));
    	colors = realloc(colors,sizeof *colors *((c_len+1)*5)); 
    	trans = realloc(trans,sizeof *trans *((c_len+1)*2));
			trans[c_len] = 1.0;
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
		*/
			tria[tric*6] = x2;
			tria[tric*6+1] = y2;
			tria[tric*6+2] = aacc->c[paaa].at.x;
			tria[tric*6+3] = aacc->c[paaa].at.y;
			tria[tric*6+4] = ix2;
			tria[tric*6+5] = iy2;
			//printf("%i %f %f\n",zzi,c->len,pixels[zzi*3]);
			tricol[tric*9] = colors[zzi*3];
			tricol[tric*9+1]=colors[zzi*3+1];
			tricol[tric*9+2]=colors[zzi*3+2];

			tricol[tric*9+3] =colors[zzi*3];
			tricol[tric*9+4]=	colors[zzi*3+1];
			tricol[tric*9+5]=	colors[zzi*3+2];
			
			tricol[tric*9+6] =colors[zzi*3];
			tricol[tric*9+7]=	colors[zzi*3+1];
			tricol[tric*9+8]=	colors[zzi*3+2];
			tric++;
		
		//printf("%f %f, %f %f, %f %f, %f %f\n",x1,y1,x2,y2,ix1,iy1,ix2,iy2);
			tria[tric*6] = x2;
			tria[tric*6+1] = y2;
			tria[tric*6+2] = aacc->c[paaa].at.x;
			tria[tric*6+3] = aacc->c[paaa].at.y;
			tria[tric*6+4] = ix1;
			tria[tric*6+5] = iy1;
			
			tricol[tric*9] = colors[zzi*3];
			tricol[tric*9+1]=colors[zzi*3+1];
			tricol[tric*9+2]=colors[zzi*3+2];

			tricol[tric*9+3] =colors[zzi*3];
			tricol[tric*9+4]=	colors[zzi*3+1];
			tricol[tric*9+5]=	colors[zzi*3+2];
			
			tricol[tric*9+6] =colors[zzi*3];
			tricol[tric*9+7]=	colors[zzi*3+1];
			tricol[tric*9+8]=	colors[zzi*3+2];
			tric++;

			
			break;
		}
		free(aacc->c);
		free(aacc->vert);
		free(aacc);
		
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
	rea->tlen = tric;
	rea->tri = tria;
	//rea->depth = dep;
	rea->max = maxf;
	rea->min = minf;
	//printf("  %i\n",tric);
	rea->tricol = tricol;
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
  int a_tlen = a->tlen;
	a->pix = realloc(a->pix,sizeof(*a->pix)*(a->len+b->len+1)*20);
  a->col = realloc(a->col,sizeof(*a->col)*(a->len+b->len+1)*20);
  a->trans = realloc(a->trans,sizeof(*a->trans)*(a->len+b->len+1)*20);
	a->tri = realloc(a->tri,sizeof(*a->tri)*(a->tlen+b->tlen+1)*60);
	a->tricol = realloc(a->tricol,sizeof(*a->tricol)*(a->tlen+b->tlen+1)*30);
	a->len+=b->len;
  a->tlen+=b->tlen;
  if(a->tri==NULL||a->pix==NULL||a->col==NULL||a->trans==NULL)
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
	for(int i = 0; i<=b->tlen*6; i++){
    a->tri[a_tlen*6+i] = b->tri[i];  
  }
	for(int i = 0; i<=b->tlen*9; i++){
    a->tricol[a_tlen*9+i] = b->tricol[i];  
  }
}
point_arr* polygon3d(double* vx, double*vy, double* vz, int n,float r, float g, float b){
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
		pa->c[i].at.vertex = 1;

		pa->vert[i].at.x = vx[i];
		pa->vert[i].at.y = vy[i];	
		pa->vert[i].at.z = vz[i];
		pa->vert[i].at.vertex = 1;

		pa->vert[i].color.r = r;
		pa->vert[i].color.g = g;	
		pa->vert[i].color.b = b;
		
		pa->c[i].color.r = r;
		pa->c[i].color.g = g;
		pa->c[i].color.b = b;
	}
	pa->len = n;
	pa->vlen= n;
  return pa;
}
point_arr* square_gen(double* tl, double* tr, double* bl, double*br,float rr, float gg, float bb){
  double xx[3] = {tl[0],tr[0]};
  double yy[3] = {tl[1],tr[1]};
  double zz[3] = {tl[2],tr[2]};
  point_arr* a = polygon3d(xx,yy,zz,2,rr,gg,bb);
  
  double xx1[3] = {tl[0],bl[0]};
  double yy1[3] = {tl[1],bl[1]};
  double zz1[3] = {tl[2],bl[2]};
  point_arr* b = polygon3d(xx1,yy1,zz1,2,rr,gg,bb);
  
  double xx2[3] = {tr[0],br[0]};
  double yy2[3] = {tr[1],br[1]};
  double zz2[3] = {tr[2],br[2]};
  point_arr* c = polygon3d(xx2,yy2,zz2,2,rr,gg,bb);
  
  double xx3[3] = {bl[0],br[0]};
  double yy3[3] = {bl[1],br[1]};
  double zz3[3] = {bl[2],br[2]};
  point_arr* d = polygon3d(xx3,yy3,zz3,2,rr,gg,bb);
  
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
	double xx1[5]={tl[0],tr[0],br[0],bl[0],	tl[0]};
	double yy1[5]={tl[1],tr[1],br[1],bl[1], tl[1]};
	double zz1[5]={tl[2],tr[2],br[2],bl[2], tl[2]};
	mm[mm->len].at = polygon3d(xx1,yy1,zz1,6,1.0,0.0,0.0);
	mm->len++;
	}
	{
	double xx1[5]={tl2[0],tr2[0],br2[0],bl2[0],	tl2[0]};
	double yy1[5]={tl2[1],tr2[1],br2[1],bl2[1], tl2[1]};
	double zz1[5]={tl2[2],tr2[2],br2[2],bl2[2], tl2[2]};
	mm[mm->len].at = polygon3d(xx1,yy1,zz1,6,0.5,0.5,0.0);
	mm->len++;
	}
	
	{
	double xx1[5]={tl2[0],tr2[0],tr[0],tl[0],	tl2[0]};
	double yy1[5]={tl2[1],tr2[1],tr[1],tl[1],	tl2[1]};
	double zz1[5]={tl2[2],tr2[2],tr[2],tl[2], tl2[2]};
	mm[mm->len].at = polygon3d(xx1,yy1,zz1,6,0.0,1.0,0.0);
	mm->len++;
	}
	{
	double xx1[5]={bl2[0],br2[0],br[0],bl[0], bl2[0]};
	double yy1[5]={bl2[1],br2[1],br[1],bl[1], bl2[1]};
	double zz1[5]={bl2[2],br2[2],br[2],bl[2], bl2[2]};
	mm[mm->len].at = polygon3d(xx1,yy1,zz1,6,0.0,0.5,0.5);
	mm->len++;
	}

	{
	double xx1[5]={tl2[0],bl2[0],bl[0],tl[0], tl2[0]};
	double yy1[5]={tl2[1],bl2[1],bl[1],tl[1], tl2[1]};
	double zz1[5]={tl2[2],bl2[2],bl[2],tl[2], tl2[2]};
	mm[mm->len].at = polygon3d(xx1,yy1,zz1,6,0.0,0.0,0.1);
	mm->len++;
	}
	{
	double xx1[5]={tr2[0],br2[0],br[0],tr[0], tr2[0]};
	double yy1[5]={tr2[1],br2[1],br[1],tr[1], tr2[1]};
	double zz1[5]={tr2[2],br2[2],br[2],tr[2], tr2[2]};
	mm[mm->len].at = polygon3d(xx1,yy1,zz1,6,0.3,0.3,0.3);
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
  
  
	///*
  double tl[3] = {5.0,200.0,200.0};
  double tr[3] = {200.0,200.0,200.0};
  double bl[3] = {5.0,5.0,200.0};
  double br[3] = {200.0,5.0,200.0};
  
  double tl2[3] = {5.0,200.0,5.0};
  double tr2[3] = {200.0,200.0,5.0};
  double bl2[3] = {5.0,5.0,5.0};
  double br2[3] = {200.0,5.0,5.0};
	float rr = 0.0;
	float gg = 0.5;
	float bb = 1.0;
	point_m* aaaa = rect3d_gen(tl,tr,bl,br,tl2,tr2,bl2,br2,rr,gg,bb);
	//*/
	/*
	double xxx[4] = {2.0,100.0,50.0,2.0};
	double yyy[4] = {2.0,2.0,100.0,2.0};
	double zzz[4] = {2.0,2.0,2.0,2.0};	
	point_m* aaaa = malloc(sizeof(*aaaa)*5);
	aaaa->len = 0;
	aaaa[0].at = polygon3d(xxx,yyy,zzz,5,1.0,1.0,0.0);
	*/
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
	double single_frame = ((float)1/frame_limit)*1E+6;
	char sf_time_msg[60];
	sprintf(sf_time_msg,"%.3fµs per frame ≈ %.3ffps",single_frame,frame_limit);
	info(sf_time_msg);
	
	struct timespec start_t, end_t,tem_t2, tem_t;
	
	clock_gettime(CLOCK_REALTIME, &start_t);
	info("entering main loop");
	for(;;){	
	 	clock_gettime(CLOCK_REALTIME,&tem_t);
    glfw_clear(w);
		double p1 = plr_x*0.01; 
    double p2 = plr_y*0.01;
    double p3 = 0;
    double p4 = pl_y;
    double p5 = -pl_y+pl_z;
    double p6 = pl_x; 

		if(aaaa->len>=0){
		glfl_ar**con = malloc(sizeof **con * aaaa->len);
		int con_len = 0;
		glfl_ar* bba = perspective_proj(w,aaaa[0].at,p1,p2,p3,p4,p5,p6);
		con[con_len] = bba;
		con_len++;
		if(aaaa->len>0){
		for(int i = 1; i<=aaaa->len-1; i++){
			glfl_ar* bbb = perspective_proj(w,aaaa[i].at,p1,p2,p3,p4,p5,p6);
			//printf("%f\n",bbb->depth);
			con[con_len] = bbb;
			con_len++;
			//join_glfl_a(bba,bbb);
			//free(bbb->col);
			//free(bbb->pix);
			//free(bbb->trans);
			//free(bbb->tri);
			//free(bbb->tricol);
			//free(bbb);
			
			}
		}
		for(int i = 0; i<=aaaa->len-2; i++){
			//printf("%i %i\n",i,aaaa->len-1);
			if(con[i]->max.z<con[i+1]->min.z||
					con[i]->min.z<con[i+1]->min.z||
					con[i]->max.z<con[i+1]->max.z){
				
				glfl_ar* tempp = con[i];	
				con[i] = con[i+1];	
				con[i+1] = tempp;	

				i-=1;
			}
		}
		
		//char* absabd = malloc(1000);
		//absabd = "aaa\n";
		glfl_ar* push = con[0];
		for(int i = 1; i<=con_len-1;i++){
			//printf("%f %f\n",con[i]->max.x,con[i]->min.x);
			join_glfl_a(push,con[i]);
			free(con[i]->tri);
			free(con[i]->tricol);
			free(con[i]->pix);
			free(con[i]->col);
			free(con[i]->trans);
			free(con[i]);	
		}
		render_p(push);
		free(con[0]->tri);
		free(con[0]->tricol);
		free(con[0]->pix);
		free(con[0]->col);
		free(con[0]->trans);
		free(con[0]);
		free(con);
		//free(bba->tricol);
		//free(bba->tri);
		//free(bba->trans);
		//free(bba->col);
		//free(bba->pix);
		//free(bba);
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
	free(aaaa[0].at->c);
	free(aaaa[0].at->vert);
	free(aaaa[0].at);
	for(int i = 1; i<=aaaa->len-1; i++){
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
