#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "glfww.h"
#include <unistd.h>
float NUU = 1;
typedef struct {
  double x;
  double y;
  double z;
  int vertex;
} cord;
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
double binominal(int n, int k){
  if(n==k)
    return 1.0;
  double v = 1.0; 
  for(int i = 1; i<=k; i++){
    v=v*((float)(n+1-i)/i);
  } 
  return v;
}
int ma = 4;
void mul_matr(double a[][ma],double b[][ma], double res[][ma]){
  for(int i = 0; i<ma; i++){
    for(int j = 0; j<ma;j++){
      res[i][j] = 0;
      for(int k = 0; k<ma;k++){
        res[i][j]+=a[i][k]*b[k][j];
        
      }
    }
  }
}
GLuint prog;
const char* vshader_src = 
  "#version 330\n"
  "layout (location = 0) in vec3 pos;\n"
  "layout (location = 1) in vec3 color;\n"
  "out vec3 ncolor;\n"
  "void main(){\n"
  "ncolor = color;\n"
  "gl_Position = vec4(pos,1.0);\n" 
  "};";
const char* fshader_src = 
  "#version 330\n"
  "in vec3 ncolor;\n"
  "out vec3 color;\n"
  "void main(){\n"
  "gl_FragColor = vec4(ncolor,1.0);\n"
  "};";
GLuint vshader_comp(const char* shader_src){
  GLuint vertid = glCreateShader(GL_VERTEX_SHADER);
  glShaderSource(vertid,1,(const GLchar**)&shader_src, NULL);
  glCompileShader(vertid);
  return vertid;
}
GLuint fshader_comp(const char* shader_src){
  GLuint fragid = glCreateShader(GL_FRAGMENT_SHADER);
  glShaderSource(fragid,1,(const GLchar**)&shader_src, NULL);
  glCompileShader(fragid);
  return fragid;
}
GLuint build_shader(GLuint vertid, GLuint fragid){
  GLuint progid = glCreateProgram();
  glAttachShader(progid,vertid);
  glAttachShader(progid,fragid);
  glLinkProgram(progid);
  return progid;
}
void perspective_proj(GLFWwindow* b,point_arr* c,double ctx,double cty,double ctz,double cx, double cy, double cz){
  //double cx = -100;
  //double cy = 0;
  //double cz = 0;
  //double ctx = 0;
  //double cty = 0;
  //double ctz = 0;
  //cz=100-cz;
  //printf("%i\n",glGetError());
  //glColor3f(1.0f,0.0f,0.0f);
  //GLfloat* pixels = malloc(sizeof(*pixels)*c->len*5);
  //printf("s\n");
  //GLfloat pixels[(int)ceil(c->len*3)];
  //GLfloat colors[(int)ceil(c->len*3)];
  GLfloat* pixels = malloc(sizeof(*pixels)*(c->len*3));
  GLfloat* abpixels = malloc(sizeof(*abpixels)*(c->len*3));
  GLfloat* colors = malloc(sizeof(*colors)*(c->len*3));
  //GLfloat* colors = malloc(sizeof(*colors)*c->len*4);
  if(pixels==NULL||colors==NULL)
    err("failed to allocate perspective array:(",pexit);
  
  
  //glfwSwapBuffers(b);
  //return;
  double coy = cos(cty);
  double siz = sin(ctz);
  double coz = cos(ctz);
  double six = sin(ctx);
  double siy = sin(cty);
  double cox = cos(ctx);
  double fov = .01;
  double ex = cx;
  double ey = cy;
  //glfwGetFramebufferSize(b,&w,&h);
  refresh_size(b);
  //glColor3f(1.0f,0.0f,0.0f);
  //double ez = 1/tan(fov/2);
  //glBegin(GL_POINTS);
  //printf("e\n");
  GLuint fb = 0;
  double ez=get_h()*2; //i dont get this at all
  //glGenFramebuffers(1,&fb);
  glEnableClientState(GL_VERTEX_ARRAY);
  //glEnableClientState(GL_VERTEX_ARRAY);
  //glEnableClientState(GL_COLOR_ARRAY);
  for(int i = 0; i!=c->len; i++){
    double ax = c->c[i].at.x;
    double ay = c->c[i].at.y;
    double az = c->c[i].at.z;
    
    double eyz = (coz*(ay-cy)-siz*ax-cx);
    double yzm = (coy*(az-cz) + siy*(siz*(ay-cy) + coz*(ax-cx)));
    double dx = coy * (siz*(ay-cy) + coz*(ax-cx)) - (siy*(az-cz));
    double dy = six * yzm + cox*eyz;
    double dz = cox * yzm - six*eyz; 
    
    //printf("%f %f %f, %f %f\n",dx,dy,ez,ez/dz*dx+dx,dy*ez);
    double bx = ez/dz*dx+dx;
    double by = ez/dz*dy+dy;
    //printf("%f %f | %f %f %f\n",bx,by,ctx,cty,ctz);
    //int aaa = round((150-dz)/2);
    if(dz>=0){
      ab_to_vp(xa,ya,get_w(),get_h(),bx,by);
    //return;
    //printf("%i:%f %f | %f %f\n",i*2,xa,ya,bx,by);
    //return;
      pixels[i*2] = xa+1;
      pixels[i*2+1] = ya;
      abpixels[i*2] = bx;
      abpixels[i*2+1] = by;
      colors[i*3] = c->c[i].color.r;
      colors[i*3+1] = c->c[i].color.g;
      colors[i*3+2] = c->c[i].color.b;
    } else {
      pixels[i*2] = -20;
      pixels[i*2+1]= -20;
      abpixels[i*2] = -20;
      abpixels[i*2+1]= -20;
      colors[i*3] = 0.0f;
      colors[i*3+1] = 1.0f;
      colors[i*3+2] = 0.0f; 
    }
    //glfw_circle_partial(b,nearbyint(bx),nearbyint(by),/*(0>=aaa?1:*/1/*)*/); 
    //glfw_pixel_partial(b,round(bx),round(by)); 
  }
  //highlight vertic
  int c_len = c->len;
  for(int i = 0; i<=c->len; i++){
    if(c->c[i].at.vertex==1){ 
      if(pixels==NULL||colors==NULL)
        abort();
      double ttt = -1.0;
      int p_b = 0;
      int p_b2 = 0;
      double ttt2 = 1.0;
      for(int jj = 0; jj<=c->len; jj++){
        float sad = 1;
        float sad2 = 1;
        if(pixels[jj*2]<pixels[i*2]&&diff(pixels[jj*2],pixels[i*2])>(float)sad2/get_w()
          &&diff(pixels[jj*2+1],pixels[i*2+1])<(float)sad/get_w()&&i!=jj){
          //printf("%f %f\n",ttt,pixels[jj*2]);
          //if(ttt<pixels[jj*2])
          p_b=1;
          if(pixels[jj*2]>ttt)
          ttt=pixels[jj*2];
        }
        if(diff(pixels[jj*2],pixels[i*2])>(float)sad2/get_w()
          &&diff(pixels[jj*2+1],pixels[i*2+1])<(float)sad/get_w()&&i!=jj&&
          pixels[jj*2]>pixels[i*2]){
          //printf("%f %f\n",ttt,pixels[jj*2]);
          //if(ttt<pixels[jj*2])
          p_b2=1;
          if(pixels[jj*2]<ttt2)
          ttt2=pixels[jj*2];
        } 
      }
    if(p_b)
      for(double zz = pixels[i*2]; diff(zz,ttt)>(float)5/get_w()&&zz>-1;zz-=(float)1/get_w()){
        pixels = realloc(pixels,sizeof *pixels *((c_len+1)*3));
        colors = realloc(colors,sizeof *colors *((c_len+1)*4));
        int brea = 0;
        
        pixels[c_len*2] = zz;
        pixels[c_len*2+1] = pixels[i*2+1]; 
        colors[c_len*3] = 0.1f;
        colors[c_len*3+1] = 1.0f; 
        colors[c_len*3+2] = 1.0f; 
        c_len++;

       
      }
    if(p_b2)
      for(double zz = pixels[i*2]; diff(zz,ttt2)>(float)1/get_w()&&zz<1;zz+=(float)1/get_w()){
        pixels = realloc(pixels,sizeof *pixels *((c_len+1)*3));
        colors = realloc(colors,sizeof *colors *((c_len+1)*4));
        int brea = 0;
        
        pixels[c_len*2] = zz;
        pixels[c_len*2+1] = pixels[i*2+1]; 
        colors[c_len*3] = 0.1f;
        colors[c_len*3+1] = 1.0f; 
        colors[c_len*3+2] = 1.0f; 
        c_len++;

       
      }
    }
  }
  //return;
  glPointSize(2.0f);
  //glUseProgram(prog);
  //glVertexAttribPointer(0,3,GL_FLOAT,0, 7*4,0);
  //glVertexAttribPointer(3,4,GL_FLOAT,0, 7*4,3*4);
  //GLint posAttrib = glGetAttribLocation(0,"position");
  //GLuint vertbuff;
  //glGenBuffers(1,&vertexBufferObject);
  //for(int i = 0; i!=c->len*3; i++){
  //  colors[i] = 1.0f;
  //} 
  /*GLint uni = glGetUniformLocation(prog,"r");
  glUniform1f(uni,1.0);
  uni = glGetUniformLocation(prog,"g");
  glUniform1f(uni,0.2);
  uni = glGetUniformLocation(prog,"b");
  glUniform1f(uni,0.3);*/
  /*glVertexAttribPointer(0,2,GL_FLOAT,GL_FALSE, c->len/2*sizeof(float),(void*)0);
  glEnableVertexAttribArray(0);
  int aaa3 = c->len/2*sizeof(float);
  glVertexAttribPointer(2,5,GL_FLOAT,GL_FALSE, c->len/2*sizeof(float),(void*)(&aaa3));
  glEnableVertexAttribArray(1);*/
 /* GLuint colorb;
  glGenBuffers(1,&colorb);
  glBindBuffer(GL_ARRAY_BUFFER,colorb);
  glBufferData(GL_ARRAY_BUFFER,sizeof(colors),colors,GL_STATIC_DRAW);
  glEnableVertexAttribArray(1);
  glBindBuffer(GL_ARRAY_BUFFER,colorb);
  glVertexAttribPointer(1,3,GL_FLOAT,GL_FALSE,0,(void*)0);
*/// return; 
  
  const GLfloat tr[] = {
    1.0f, 1.0f,
    1.0f, -1.0f, 
    -1.0f, -1.0f,

    -1.0f, -1.0f,
    -1.0f, 1.0f, 
    1.0f, 1.0f, 
  };
  GLuint verta;
  glGenVertexArrays(1,&verta);
  glBindVertexArray(verta);

  GLuint vetb;
  glGenBuffers(1,&vetb);
  glBindBuffer(GL_ARRAY_BUFFER,vetb);
  glBufferData(GL_ARRAY_BUFFER,sizeof(*pixels)*(c_len*3),pixels,GL_STATIC_DRAW);

  glEnableVertexAttribArray(0);
  glBindBuffer(GL_ARRAY_BUFFER,vetb);
  glVertexAttribPointer(0,2,GL_FLOAT,GL_FALSE,0,(void*)0);
  
  GLuint colb;
  glGenBuffers(1,&colb);
  glBindBuffer(GL_ARRAY_BUFFER,colb);
  glBufferData(GL_ARRAY_BUFFER,sizeof(*colors)*(c_len*4),colors,GL_STATIC_DRAW);

  glEnableVertexAttribArray(1);
  glBindBuffer(GL_ARRAY_BUFFER,colb);
  glVertexAttribPointer(1,3,GL_FLOAT,GL_FALSE,0,(void*)0);

  glDrawArrays(GL_POINTS,0,c_len);
  glDeleteBuffers(1,&vetb);
  glDeleteBuffers(1,&colb);
  free(pixels);
  free(colors);
  free(abpixels);
  //glDisableVertexAttribArray(0);
  //glDisableVertexAttribArray(1);
  //glVertexPointer(2,GL_FLOAT,0,pixels);
  //glVertexPointer(3,GL_FLOAT,0,colors);
  //glDrawArrays(GL_POINTS,0,c->len);
  //glDisableClientState(GL_VERTEX_ARRAY);
  //glDisableClientState(GL_COLOR_ARRAY);
  //free(pixels);
  //free(colors);
  /*for(int i = 0; i!=get_w();i++){
    glfw_pixel(b,i,get_h()/2);
  }
  for(int i = 0; i!=get_h();i++){
    glfw_pixel(b,get_w()/2,i);
  }*/
  //glEnd();
}
point_arr* basier3d(double*xx,double*yy,double*zz,int n,float rr, float gg, float bb){
  point_arr* pa;
  pa = malloc(sizeof(*pa));
  pa->c = malloc(sizeof(*pa->c)*(get_w()*60));
  pa->vert = malloc(sizeof(*pa->vert)*(n*60));

  if(pa->c==NULL||pa->vert==NULL)
    err("failed to allocate basier array",pexit);
  
  //double xx[5] = {5.0,5.0,50.0,5.0,5.0};
  //double yy[5] = {5.0,5.0,5.0,5.0,10.0};
  //double zz[5] = {10.0,5.0,5.0,5.0,5.0};
  //int n = 5-1;
  pa->len = get_w();
  n-=1;

  
  double aaar = (1.0/get_w());
  for(int iy = 0; iy<=get_w();iy+=1){
    double t = aaar*iy;
    double bcx = 0;
    double bcy = 0;
    double bcz = 0;
    for(int i = 0; i <=n;i++){
      double pp = binominal(n,i) * pow((1 - t),(n - i)) * pow(t,i);
      bcx +=  pp * xx[i];
      bcy +=  pp * yy[i];
      bcz +=  pp * zz[i];
      //if(ii==3){
      //  printf("pp: %f bi : %f p1 : %f p2 : %f| %f, %f, %f\n",pp,binominal(n,i),bcx,bcy,bcz,pow((1-t),(n-i)),pow(t,i));
      //}
    }
    //int ii = floor(100*(double)(t/1.0));
    //printf("%i\n",iy);
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
    //printf("(%f,%f,%f)\n",bcx,bcy,bcz);
  }
  for(int i = 0; i<=n; i++){
    pa->vert[i].at.x = xx[i];
    pa->vert[i].at.y = yy[i]; 
    pa->vert[i].at.z = zz[i]; 
  }
  pa->vlen = n;
  //bw b = sdl_init();
  /*for(int i = 0; i!=100; i ++){
    SDL_SetRenderDrawColor(b.r,0,255-255*(pa[i].z/100),0,255);
    SDL_RenderDrawPoint(b.r,pa[i].x,pa[i].y);
    sdl_circle(b,round(pa[i].x),round(pa[i].y),round(100-pa[i].z));
    //printf("%i : (%f,%f,%f)\n",i,pa[i].x,pa[i].y,(100-pa[i].z));
  }*/
  //perspective_proj(b,pa,100);
  //sdl_loop(b);
  return pa;
}
void join_cords(point_arr* a, point_arr* b){
  //printf("%lu\n",sizeof(*a)*(a_len+b_len+2));
  int a_len = a->len;
  int a_vlen = a->vlen;
  a->c = realloc(a->c,sizeof(*a->c)*(a->len+b->len)*60);
  a->vert = realloc(a->vert,sizeof(*a->vert)*(a->vlen+b->vlen)*60);
  a->len+=b->len;
  a->vlen+=b->vlen;
  if(a->c==NULL)
    err("failed to reallocate cords",pexit);
  for(int i = 0; i<=b->len; i++){
    //printf("%i/%f\n",i,b->len);
    a->c[a_len+i].at = b->c[i].at; 
    a->c[a_len+i].color = b->c[i].color; 
  }
  for(int i = 0; i<=b->vlen; i++){
    //printf("%i/%f\n",i,b->len);
    a->vert[a_vlen+i].at = b->vert[i].at; 
    //a->vert[a_len+i].color = b->vert[i].color; 
  }
  //return a.len+b.len;
}
point_arr* square_gen(double* tl, double* tr, double* bl, double*br,float rr, float gg, float bb){
  double xx[3] = {tl[0],tr[0]};
  double yy[3] = {tl[1],tr[1]};
  double zz[3] = {tl[2],tr[2]};
  point_arr* a = basier3d(xx,yy,zz,2,rr,gg,bb);
  
  double xx1[3] = {tl[0],bl[0]};
  double yy1[3] = {tl[1],bl[1]};
  double zz1[3] = {tl[2],bl[2]};
  point_arr* b = basier3d(xx1,yy1,zz1,2,rr,gg,bb);
  
  double xx2[3] = {tr[0],br[0]};
  double yy2[3] = {tr[1],br[1]};
  double zz2[3] = {tr[2],br[2]};
  point_arr* c = basier3d(xx2,yy2,zz2,2,rr,gg,bb);
  
  double xx3[3] = {bl[0],br[0]};
  double yy3[3] = {bl[1],br[1]};
  double zz3[3] = {bl[2],br[2]};
  point_arr* d = basier3d(xx3,yy3,zz3,2,rr,gg,bb);
  
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
  //printf("AAA");
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
  //printf("%f %f %f | %f %f %f\n",m_x,m_y,m_z,mi_x,mi_y,mi_z);
  int a_l = a->len;
  a->c = realloc(a->c,sizeof(*a->c)*(m_x*m_y+a->len)*60);
  for(double y = mi_y; y<=m_y; y+=.1){
    //printf("%f/%f\n",y,m_y);
    double* zz = malloc(sizeof(*zz));
    int zzl = 0;
    double* xx = malloc(sizeof(*xx));
    double* yy = malloc(sizeof(*yy));
    for(int i = 0; i<=a_l; i++){
      //printf("%i/%f\n",i,a->len);
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
  /*for(int i = 0; i!=a->len+m_x*m_y; i++){
    a->c[i].at.x = ;
    a->c[i].at.y = fmod(i,m_y); 
    a->c[i].at.z = 20;   
    a->c[i].color.r = 0.0f;
    a->c[i].color.g = 1.0f; 
    a->c[i].color.b = 0.1f; 
  }*/
  /*
  for(int y = mi_y; y<=m_y; y++){
    for(int x = mi_x; x<=m_x; x++){
      a->len++;
      int i = a->len;
     // printf("%i\n",i);
      a->c[a_l+i].at.x = x;
      a->c[a_l+i].at.y = y; 
      a->c[a_l+i].at.z = mi_z;   
      a->c[a_l+i].color.r = 0.0f;
      a->c[a_l+i].color.g = 1.0f; 
      a->c[a_l+i].color.b = 0.1f; 
    }
  }
  printf("done\n");
  */ 
  //join_cords(a,bb);
}
point_arr* cube_gen(double* tl, double* tr, double* bl, double*br,
               double* tl2, double* tr2, double* bl2, double*br2,
                    float rr, float gg, float bb){
   
  point_arr* a = square_gen(tl,tr,bl,br,rr,gg,bb); 
  //fill3d(a);
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
point_arr* polygon3d(double* vx, double*vy, double* vz, int n){
  double xx[2] = {vx[0],vx[1]};
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
  }
  return y;
}
int main(int argc,char*argv[]){
  flag_handle(argc,argv);
  atexit(sig_handle);
  GLFWwindow* w = glfw_init();
  refresh_size(w);
  GLenum err = glewInit();
  if (err != GLEW_OK)
    exit(1); // or handle the error in a nicer way
  if (!GLEW_VERSION_2_1)  // check that the machine supports the 2.1 API.
    exit(1); // or handle the error in a nicer way
  GLuint vid = vshader_comp(vshader_src);
  GLuint fid = fshader_comp(fshader_src);
  prog = build_shader(vid,fid);
  glUseProgram(prog); 
  info("built shaders");
  
  
  
  double tl2[3] = {5.0,200.0,400.0};
  double tr2[3] = {200.0,200.0,400.0};
  double bl2[3] = {5.0,5.0,200.0};
  double br2[3] = {200.0,5.0,200.0};
  
  double tl1[3] = {5.0,200.0,200.0};
  double tr1[3] = {200.0,200.0,200.0};
  double bl1[3] = {5.0,5.0,5.0};
  double br1[3] = {200.0,5.0,5.0};
  point_arr* a = cube_gen(tl1,tr1,bl1,br1,tl2,tr2,bl2,br2,0.1f,0.1f,1.0f);
  /*double xx[8] = {0.0, 15.0, 50.0, 60.0,40.0,30.0, 0.0,0.0};
  double yy[8] = {5.0, 15.0, 30.0, 45.0,64.0, 45.0,55.0,5.0};
  double zz[8] = {50.0,50.0,50.0,50.0,50.0,50.0,50.0,50.0};
  point_arr* a = polygon3d(xx,yy,zz,8);*/
  
  
  /*for(int i = 0;i!=a->vlen;i++){
    printf("%f %f %f\n",a->vert[i].at.x,a->vert[i].at.y,a->vert[i].at.z);
  }*/
  //exit(0);
  int max_r = 630;
  double half_max_r = (double)max_r/2/2;
  //printf("%s\n",glGetString(GL_VERSION));
  //glfwSetKeyCallback(w,key_press);
  double pl_x = 0;
  double pl_y = 0;
  double pl_z = 0;
  double plr_x = 0;
  double plr_y = 0;
  for(;;){
    double p1 = plr_x*0.01; 
    double p2 = plr_y*0.01;
    double p3 = 0;
    double p4 = pl_y;
    double p5 = -pl_y+pl_z;
    double p6 = pl_x;
    
    perspective_proj(w,a,p1,p2,p3,p4,p5,p6);
    
    glfw_load(w);
    //break;
    int mod_move=1;
    double run_mul=2;
    //usleep(10000);
    glfwPollEvents(); 
    if(glfwGetKey(w,GLFW_KEY_I)){
      if(glfwGetKey(w,GLFW_KEY_LEFT_SHIFT))
        plr_x-=mod_move;
      else 
        plr_x--;
    plr_x = fmod(plr_x,max_r);
    //plr_y = fmod(plr_y,max_r);
    }
    if(glfwGetKey(w,GLFW_KEY_K)){
      if(glfwGetKey(w,GLFW_KEY_LEFT_SHIFT))
        plr_x+=mod_move;
      else 
        plr_x++;
    plr_x = fmod(plr_x,max_r);
    //plr_y = fmod(plr_y,max_r);
    } 
    if(glfwGetKey(w,GLFW_KEY_J)){
      if(glfwGetKey(w,GLFW_KEY_LEFT_SHIFT))
        plr_y-=mod_move;
      else 
        plr_y--;
    //plr_x = fmod(plr_x,max_r);
    plr_y = fmod(plr_y,max_r);
    } 
    if(glfwGetKey(w,GLFW_KEY_L)){
      if(glfwGetKey(w,GLFW_KEY_LEFT_SHIFT))
        plr_y+=mod_move;
      else 
        plr_y++;
    //plr_x = fmod(plr_x,max_r);
    plr_y = fmod(plr_y,max_r);
    }
    
    if(glfwGetKey(w,GLFW_KEY_W)){
      double mul = 1;
      if(glfwGetKey(w,GLFW_KEY_LEFT_SHIFT))
        mul = run_mul;
      pl_x+=cosf(plr_y*0.01)*mul;
      pl_y+=sinf(plr_y*0.01)*mul; 
      pl_z-=sinf(plr_x*0.01)*mul; 
      //printf("%f\n",pl_z);
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
    if(glfwGetKey(w,GLFW_KEY_O))
      NUU+=0.1;
    if(glfwGetKey(w,GLFW_KEY_P))
      NUU-=0.1;
    printf("%f\n",NUU);
    //printf("%f %f %f\n",plr_y,cosf(plr_y*0.01),sinf(plr_y*0.01));
    usleep(1000*1000/60);
    glfw_clear(w); 
    
    if(glfwWindowShouldClose(w)||(glfwGetKey(w,GLFW_KEY_Q)))break;
  }
  free(a->c);
  free(a->vert);
  free(a);
  glfwDestroyWindow(w);
  win_clean();
  glDeleteShader(vid);
  glDeleteShader(fid);
  glDeleteShader(prog);
  info("killed window:p");
  return 0;
}
