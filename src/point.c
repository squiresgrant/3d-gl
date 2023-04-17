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
  int roll; //x axis
  int pitch; //y axis
  int yaw; //z axis
} rot;
typedef struct {
  float opacity;
} opts;
typedef struct {
  cord loc;
  rot rot;
  int uid;
  int id;
  opts opts;
} point;

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
  
  GLfloat* pixels = malloc(sizeof(*pixels)*len*800*500*2);

  
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
  double ez=500; //i dont get this at all
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
  cord* pa = malloc(sizeof(*pa)*(3000));
  //double xx[5] = {5.0,5.0,50.0,5.0,5.0};
  //double yy[5] = {5.0,5.0,5.0,5.0,10.0};
  //double zz[5] = {10.0,5.0,5.0,5.0,5.0};
  //int n = 5-1;
  n-=1;
  for(int iy = 0; iy<=100;iy+=1){
    double t = 0.01*iy;
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
void join_cords(cord* a, cord* b, int a_len, int b_len){
  printf("%lu\n",sizeof(*a)*(a_len+b_len+2));
  a = realloc(a,sizeof(*a)*(a_len+b_len+2));
  for(int i = 0; i<=a_len+b_len; i++){
    a[a_len+i] = b[i]; 
  }
}

int main(){
  double xx[5] = {5.0,5.0,5.0};
  double yy[5] = {5.0,100.0,200.0};
  double zz[5] = {95.0,95.0,95.0};
  cord* a = basier3d(xx,yy,zz,3);
 
  double xx2[5] = {200.0,200.0,200.0};
  double yy2[5] = {5.0,100.0,200.0};
  double zz2[5] = {95.0,95.0,95.0};
  cord* b = basier3d(xx2,yy2,zz2,3);
  
  double xx3[5] = {5.0,100.0,200.0};
  double yy3[5] = {5.0,5.0,5.0};
  double zz3[5] = {95.0,95.0,95.0};
  cord* c = basier3d(xx3,yy3,zz3,3);
  
  double xx4[5] = {5.0,100.0,200.0};
  double yy4[5] = {200.0,200.0,200.0};
  double zz4[5] = {95.0,95.0,95.0};
  cord* d = basier3d(xx4,yy4,zz4,3); 
  
  double xx5[5] = {5.0,100.0,200.0};
  double yy5[5] = {200.0,200.0,200.0};
  double zz5[5] = {290.0,290.0,290.0};
  cord* e = basier3d(xx5,yy5,zz5,3); 
  
  double xx6[5] = {5.0,100.0,200.0};
  double yy6[5] = {5.0,5.0,5.0};
  double zz6[5] = {290.0,290.0,290.0};
  cord* f = basier3d(xx6,yy6,zz6,3);
 
  double xx7[5] = {200.0,200.0,200.0};
  double yy7[5] = {5.0,100.0,200.0};
  double zz7[5] = {290.0,290.0,290.0};
  cord* g = basier3d(xx7,yy7,zz7,3);
  
  double xx8[5] = {5.0,5.0,5.0};
  double yy8[5] = {5.0,100.0,200.0};
  double zz8[5] = {290.0,290.0,290.0};
  cord* h = basier3d(xx8,yy8,zz8,3);
  //next
  double xx9[5] = {5.0,5.0,5.0};
  double yy9[5] = {200.0,200.0,200.0};
  double zz9[5] = {290.0,290.0,95.0};
  cord* ii = basier3d(xx9,yy9,zz9,3);
 
  double xx10[5] = {5.0,5.0,5.0};
  double yy10[5] = {5.0,5.0,5.0};
  double zz10[5] = {290.0,290.0,95.0};
  cord* j = basier3d(xx10,yy10,zz10,3);

  double xx11[5] = {200.0,200.0,200.0};
  double yy11[5] = {200.0,200.0,200.0};
  double zz11[5] = {290.0,290.0,95.0};
  cord* k = basier3d(xx11,yy11,zz11,3);
  
  double xx12[5] = {200.0,200.0,200.0};
  double yy12[5] = {5.0,5.0,5.0};
  double zz12[5] = {290.0,290.0,95.0};
  cord* l = basier3d(xx12,yy12,zz12,3);
  join_cords(a,b,100,100);
  free(b);
  join_cords(a,c,200,100);
  free(c);
  join_cords(a,d,300,100);
  free(d);
  join_cords(a,e,400,100);
  free(e);
  join_cords(a,f,500,100);
  free(f);
  join_cords(a,g,600,100);
  free(g);
  join_cords(a,h,700,100);
  free(h);
  join_cords(a,ii,800,100);
  free(ii);
  join_cords(a,j,900,100);
  free(j);
  join_cords(a,k,1000,100);
  free(k);
  join_cords(a,l,1100,100);
  free(l);
  int max_r = 630;
  double half_max_r = (double)max_r/2/2;
  GLFWwindow* w = glfw_init();
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
    perspective_proj(w,a,1200,p1,p2,p3,p4,p5,p6);

    glfw_load(w);
    
    usleep(10000);
    glfwPollEvents(); 
    if(glfwGetKey(w,GLFW_KEY_I))
      plr_x--;
    if(glfwGetKey(w,GLFW_KEY_K))
      plr_x++;
    if(glfwGetKey(w,GLFW_KEY_J))
      plr_y--;
    if(glfwGetKey(w,GLFW_KEY_L))
      plr_y++;
    plr_x = fmod(plr_x,max_r);
    plr_y = fmod(plr_y,max_r);
    if(glfwGetKey(w,GLFW_KEY_W)){
      pl_x+=cosf(plr_y*0.01);
      pl_y+=sinf(plr_y*0.01); 
      pl_z-=sinf(plr_x*0.01); 
      printf("%f\n",pl_z);
    }
    if(glfwGetKey(w,GLFW_KEY_S)){
      pl_x-=cosf(plr_y*0.01);
      pl_y-=sinf(plr_y*0.01); 
      pl_z+=sinf(plr_x*0.01); 
    } 
    if(glfwGetKey(w,GLFW_KEY_A)){
      pl_x-=cosf((half_max_r+plr_y)*0.01);
      pl_y-=sinf((half_max_r+plr_y)*0.01); 
    }
    if(glfwGetKey(w,GLFW_KEY_D)){
      pl_x+=cosf((half_max_r+plr_y)*0.01);
      pl_y+=sinf((half_max_r+plr_y)*0.01); 
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
