#include "glfww.h"
int w,h;

GLFWwindow* glfw_init(){
  GLFWwindow* window;

  if(!glfwInit())
    err("failed to init glfw",pexit);

  window = glfwCreateWindow(800,500,"nya",NULL,NULL);
  if(!window){
    glfwTerminate();
    err("failed to create window",pexit);
  }
  glEnableClientState(GL_VERTEX_ARRAY);
  //glEnableClientState(GL_COLOR_ARRAY);
  //glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
  //glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
  glfwMakeContextCurrent(window);
  int w,h;
  glfwGetFramebufferSize(window,&w,&h);
  glViewport(0,0,w,h);
  logm("yay:D i made a window uwu");
  return window;
  /*
  while(!glfwWindowShouldClose(window)){
    glClear(GL_COLOR_BUFFER_BIT);
    glColor3f(1.0f,0.0f,0.0f);
    glBegin(GL_POINTS);
    for(float x = 0.0; x!=200.0;x+=1.0){
    ab_to_vp(ax,ay,w,h,x,x);
    glVertex2f(ax,ay);
    }
    glEnd();
    glfwSwapBuffers(window);
    glfwPollEvents();
  }
  glfwTerminate();
  return 0;
  */
}
void refresh_size(GLFWwindow*wi){
  glfwGetFramebufferSize(wi,&w,&h);
}
int get_w(){
  return w;
}
int get_h(){
  return h;
}
#define glfw_load(w) glfwSwapBuffers(w);
void glfw_loop(GLFWwindow*window){
  while(!glfwWindowShouldClose(window)){
    
    //glfw_load(window);
    glfwPollEvents();
  }
  glfwTerminate();
}
void glfw_pixel_partial(GLFWwindow*wi,int x, int y){ 
  ab_to_vp(ax,ay,w,h,x,y);
  glVertex2f(ax,ay);
  
}

void glfw_clear(GLFWwindow*w){
  glClear(GL_COLOR_BUFFER_BIT);
}

void glfw_circle(GLFWwindow* wi,int x, int y, int r){
  //SDL_SetRenderDrawColor(w.r,255,255,0,0);
  glfwGetFramebufferSize(wi,&w,&h);
  glBegin(GL_POINTS);
  glColor3f(1.0f,0.0f,0.0f);
  for(int i = 1; i!=360; i++){
    float cf = cosf(i)*r;
    float sf = sinf(i)*r;
    //for(int z = 1; z<=r; z++){
    //int x2 = x + z * cf;
    //int y2 = y + z * sf; 
    
    glfw_pixel_partial(wi,x+cf,y+sf);
    //}
  } 
  glEnd();
}
void glfw_circle_partial(GLFWwindow* wi,int x, int y, int r){
  //SDL_SetRenderDrawColor(w.r,255,255,0,0); 
  if(r==0){
    glfw_pixel_partial(wi,x,y);
    return;
  }
  for(int i = 1; i!=360; i++){
    float cf = cosf(i)*r;
    float sf = sinf(i)*r;
    //for(int z = 1; z<=r; z++){
    //int x2 = x + z * cf;
    //int y2 = y + z * sf; 
    
    glfw_pixel_partial(wi,x+cf,y+sf);
    //}
  } 
}
void glfw_square(GLFWwindow* wi,int x, int y, int r){
  //SDL_SetRenderDrawColor(w.r,255,255,0,0);
  glfwGetFramebufferSize(wi,&w,&h);
  glBegin(GL_POINTS);
  glColor3f(1.0f,0.0f,0.0f); 
  for(int i = 0; i<=r; i++){
    for(int j = 0; j<=r; j++){
      glfw_pixel_partial(wi,i-r/2+x,j-r/2+y);
    }
  }
  glEnd();
}
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
int __glw_main(){
  GLFWwindow* w = glfw_init();
  glfw_circle(w,50,50,20);
  glfw_load(w);
  sleep(5);
  glfwTerminate();
  //glfw_loop(w);
  return 0;
}
