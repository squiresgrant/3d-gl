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
  int r;
  int g;
  int b;
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

void perspective_proj(GLFWwindow* b,cord* c, int len,double ctx,double cty,double ctz,double cx, double cy, double cz){
  //double cx = -100;
  //double cy = 0;
  //double cz = 0;
  //double ctx = 0;
  //double cty = 0;
  //double ctz = 0;
  //cz=100-cz;
  glColor3f(1.0f,0.0f,0.0f);
  GLfloat* pixels = malloc(sizeof(*pixels)*len*2);
  if(pixels==NULL)
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
  glColor3f(1.0f,0.0f,0.0f);
  //double ez = 1/tan(fov/2);
  //glBegin(GL_POINTS);
  GLuint fb = 0;
  double ez=get_h(); //i dont get this at all
  //glGenFramebuffers(1,&fb);
  
  for(int i = 0; i!=len; i++){
    double ax = c[i].x;
    double ay = c[i].y;
    double az = c[i].z;
    
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
    } else {
      pixels[i*2] = -20;
      pixels[i*2+1]= -20;
    }
    //glfw_circle_partial(b,nearbyint(bx),nearbyint(by),/*(0>=aaa?1:*/1/*)*/); 
    //glfw_pixel_partial(b,round(bx),round(by)); 
  }
  //return;
  glPointSize(2.0f);
  glEnableClientState(GL_VERTEX_ARRAY);
  glVertexPointer(2,GL_FLOAT,0,pixels);
  glDrawArrays(GL_POINTS,0,len);
  glDisableClientState(GL_VERTEX_ARRAY);
  free(pixels);
  /*for(int i = 0; i!=get_w();i++){
    glfw_pixel(b,i,get_h()/2);
  }
  for(int i = 0; i!=get_h();i++){
    glfw_pixel(b,get_w()/2,i);
  }*/
  //glEnd();
}
cord* basier3d(double*xx,double*yy,double*zz,int n){
  cord* pa = malloc(sizeof(*pa)*(get_w()*30));
  if(pa==NULL)
    err("failed to allocate basier array",pexit);
  //double xx[5] = {5.0,5.0,50.0,5.0,5.0};
  //double yy[5] = {5.0,5.0,5.0,5.0,10.0};
  //double zz[5] = {10.0,5.0,5.0,5.0,5.0};
  //int n = 5-1;
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
    pa[iy].x = bcx;
    pa[iy].y = bcy;
    pa[iy].z = bcz;
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
int join_cords(cord* a, cord* b, int a_len, int b_len){
  //printf("%lu\n",sizeof(*a)*(a_len+b_len+2));
  a = realloc(a,sizeof(*a)*(a_len+b_len+2));
  if(a==NULL)
    err("failed to reallocate cords",pexit);
  for(int i = 0; i<=a_len+b_len; i++){
    a[a_len+i] = b[i]; 
  }
  return a_len+b_len;
}
cord* square_gen(double* tl, double* tr, double* bl, double*br){
  double xx[3] = {tl[0],tr[0]};
  double yy[3] = {tl[1],tr[1]};
  double zz[3] = {tl[2],tr[2]};
  cord* a = basier3d(xx,yy,zz,2);
  
  double xx1[3] = {tl[0],bl[0]};
  double yy1[3] = {tl[1],bl[1]};
  double zz1[3] = {tl[2],bl[2]};
  cord* b = basier3d(xx1,yy1,zz1,2);
  
  double xx2[3] = {tr[0],br[0]};
  double yy2[3] = {tr[1],br[1]};
  double zz2[3] = {tr[2],br[2]};
  cord* c = basier3d(xx2,yy2,zz2,2);
  
  double xx3[3] = {bl[0],br[0]};
  double yy3[3] = {bl[1],br[1]};
  double zz3[3] = {bl[2],br[2]};
  cord* d = basier3d(xx3,yy3,zz3,2);
  
  int ss = join_cords(a,b,get_w(),get_w());
  ss = join_cords(a,c,ss,get_w());
  ss = join_cords(a,d,ss,get_w());
  free(b);
  free(c);
  free(d);
  return a;
}
cord* cube_gen(double* tl, double* tr, double* bl, double*br,
               double* tl2, double* tr2, double* bl2, double*br2){
  
  double n_len = get_w()*4;
  cord*a = square_gen(tl,tr,bl,br); 
  cord*b = square_gen(tl2,tr2,bl2,br2); 
  double s;
  s = join_cords(a,b,n_len,n_len);
  free(b);
  cord*c = square_gen(tl2,tr2,tl,tr);
  s = join_cords(a,c,s,n_len);
  free(c);
  cord*d = square_gen(bl2,br2,bl,br);
  s = join_cords(a,d,s,n_len);
  free(d); 
  cord*e = square_gen(tl,tl2,bl,bl2);
  s = join_cords(a,e,s,n_len);
  free(e); 
  cord*f = square_gen(br2,br,tr2,tr);
  s = join_cords(a,f,s,n_len);
  free(f); 
  return a;
}
int main(){
  GLFWwindow* w = glfw_init();
  refresh_size(w);
  double tl2[3] = {5.0,200.0,400.0};
  double tr2[3] = {200.0,200.0,400.0};
  double bl2[3] = {5.0,5.0,200.0};
  double br2[3] = {200.0,5.0,200.0};
  
  double tl1[3] = {5.0,200.0,200.0};
  double tr1[3] = {200.0,200.0,200.0};
  double bl1[3] = {5.0,5.0,5.0};
  double br1[3] = {200.0,5.0,5.0};
  cord*a = cube_gen(tl1,tr1,bl1,br1,tl2,tr2,bl2,br2);
  double s = get_w()*24;
  
  
  int max_r = 630;
  double half_max_r = (double)max_r/2/2;
  
  //glfwSetKeyCallback(w,key_press);
  double pl_x = 0;
  double pl_y = 0;
  double pl_z = 0;
  double plr_x = 0;
  double plr_y = 0;
  for(double rr = 0.01;rr!=0;rr+=0.01){
    double p1 = plr_x*0.01; 
    double p2 = plr_y*0.01;
    double p3 = 0;
    double p4 = pl_y;
    double p5 = -pl_y+pl_z;
    double p6 = pl_x;
    
    perspective_proj(w,a,s,p1,p2,p3,p4,p5,p6);

    glfw_load(w);
    int mod_move=1;
    double run_mul=2;
    usleep(10000);
    glfwPollEvents(); 
    if(glfwGetKey(w,GLFW_KEY_I)){
      if(glfwGetKey(w,GLFW_KEY_LEFT_SHIFT))
        plr_x-=mod_move;
      else 
        plr_x--;
    }
    if(glfwGetKey(w,GLFW_KEY_K)){
      if(glfwGetKey(w,GLFW_KEY_LEFT_SHIFT))
        plr_x+=mod_move;
      else 
        plr_x++;
    } 
    if(glfwGetKey(w,GLFW_KEY_J)){
      if(glfwGetKey(w,GLFW_KEY_LEFT_SHIFT))
        plr_y-=mod_move;
      else 
        plr_y--;
    } 
    if(glfwGetKey(w,GLFW_KEY_L)){
      if(glfwGetKey(w,GLFW_KEY_LEFT_SHIFT))
        plr_y+=mod_move;
      else 
        plr_y++;
    }
    plr_x = fmod(plr_x,max_r);
    plr_y = fmod(plr_y,max_r);
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
    glfw_clear(w); 
    
    if(glfwWindowShouldClose(w))break;
  }
  free(a); 
  glfwDestroyWindow(w);
  win_clean();
  return 0;
}
