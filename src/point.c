#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "glfww.h"
#include <unistd.h>
typedef struct {
  double x;
  double y;
  double z;
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
  //"in vec4 s_vPosition;\n"
  //"in vec5 aPos;\n"
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
  GLfloat pixels[(int)ceil(c->len*15)];
  GLfloat colors[(int)ceil(c->len*15)];

  //GLfloat* colors = malloc(sizeof(*colors)*c->len*4);
  //if(pixels==NULL)
  //  err("failed to allocate perspective array:(",pexit);
  
  
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
      colors[i*3] = c->c[i].color.r;
      colors[i*3+1] = c->c[i].color.g;
      colors[i*3+2] = c->c[i].color.b;
    } else {
      pixels[i*2] = -20;
      pixels[i*2+1]= -20;
      colors[i*3] = 0.0f;
      colors[i*3+1] = 1.0f;
      colors[i*3+2] = 0.0f; 
    }
    //glfw_circle_partial(b,nearbyint(bx),nearbyint(by),/*(0>=aaa?1:*/1/*)*/); 
    //glfw_pixel_partial(b,round(bx),round(by)); 
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
  

  GLuint verta;
  glGenVertexArrays(1,&verta);
  glBindVertexArray(verta);GLuint vetb;
  glGenBuffers(1,&vetb);
  
  glBindBuffer(GL_ARRAY_BUFFER,vetb);
  glBufferData(GL_ARRAY_BUFFER,sizeof(pixels),pixels,GL_STATIC_DRAW);

  glEnableVertexAttribArray(0);
  glBindBuffer(GL_ARRAY_BUFFER,vetb);
  glVertexAttribPointer(0,2,GL_FLOAT,GL_FALSE,0,(void*)0);
  
  GLuint colb;
  glGenBuffers(1,&colb);
  glBindBuffer(GL_ARRAY_BUFFER,colb);
  glBufferData(GL_ARRAY_BUFFER,sizeof(colors),colors,GL_STATIC_DRAW);

  glEnableVertexAttribArray(1);
  glBindBuffer(GL_ARRAY_BUFFER,colb);
  glVertexAttribPointer(1,3,GL_FLOAT,GL_FALSE,0,(void*)0);

  glDrawArrays(GL_POINTS,0,c->len);
  glDeleteBuffers(1,&vetb);
  glDeleteBuffers(1,&colb);
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
  if(pa->c==NULL)
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
    pa->c[iy].color.r = rr;
    pa->c[iy].color.g = gg;
    pa->c[iy].color.b = bb;
    //printf("(%f,%f,%f)\n",bcx,bcy,bcz);
  }
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
  a->c = realloc(a->c,sizeof(*a->c)*(a->len+b->len)*60);
  a->len+=b->len;
  if(a->c==NULL)
    err("failed to reallocate cords",pexit);
  for(int i = 0; i<=a->len+b->len; i++){
    a->c[a_len+i].at = b->c[i].at; 
    a->c[a_len+i].color = b->c[i].color; 
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
  free(b);
  free(c->c);
  free(c);
  free(d->c);
  free(d);
  return a;
}
point_arr* cube_gen(double* tl, double* tr, double* bl, double*br,
               double* tl2, double* tr2, double* bl2, double*br2,
                    float rr, float gg, float bb){
   
  point_arr* a = square_gen(tl,tr,bl,br,rr,gg,bb); 
  
  point_arr* b = square_gen(tl2,tr2,bl2,br2,rr,gg,bb); 
  double s;
  join_cords(a,b);
  free(b->c);
  free(b);
  point_arr* c = square_gen(tl2,tr2,tl,tr,rr,gg,bb);
  join_cords(a,c);
  free(c->c);
  free(c);
  point_arr* d = square_gen(bl2,br2,bl,br,rr,gg,bb);
  join_cords(a,d);
  free(d->c);
  free(d);
  point_arr* e = square_gen(tl,tl2,bl,bl2,rr,gg,bb);
  join_cords(a,e);
  free(e->c);
  free(e);
  point_arr* f = square_gen(br2,br,tr2,tr,rr,gg,bb);
  join_cords(a,f);
  free(f->c);
  free(f);
  return a;
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
  point_arr* a = cube_gen(tl1,tr1,bl1,br1,tl2,tr2,bl2,br2,0.1f,0.1f,0.1f); 
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
    //printf("%f %f %f\n",plr_y,cosf(plr_y*0.01),sinf(plr_y*0.01));
    usleep(1000*1000/60);
    glfw_clear(w); 
    
    if(glfwWindowShouldClose(w))break;
  }
  free(a->c);
  free(a);
  glfwDestroyWindow(w);
  win_clean();
  glDeleteShader(vid);
  glDeleteShader(fid);
  glDeleteShader(prog);
  info("killed window:p");
  return 0;
}
