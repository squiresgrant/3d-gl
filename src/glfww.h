#include <stdio.h>
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <GL/gl.h>
#include "util.h"
#include <unistd.h>
#include <stdlib.h>
#include <math.h>
#ifndef __glfww_
#define __glfww_
#define win_clean() glfwTerminate();
#define ab_to_vp(x,y,w,h,x1,y1) float x = 2 * ((float)x1/(w)) -1;\
    float y = 2 * ((float)y1/(h)) -1;
#define vp_to_ab(w,x1) ((float)x1 +1.0/2)*w
#define glfw_load(w) glfwSwapBuffers(w);
#define glfw_pixel(wi,x,y)\
  glBegin(GL_POINTS);\
  glfw_pixel_partial(wi,x,y);\
  glEnd();
GLuint vshader_comp(const char*);
GLuint fshader_comp(const char*);
GLuint build_shader(GLuint, GLuint);
GLFWwindow* glfw_init();
void glfw_loop(GLFWwindow*);
int get_h();
int get_w();
void glfw_pixel_partial(GLFWwindow*,int, int);
void glfw_clear(GLFWwindow*);
void refresh_size(GLFWwindow*);

void glfw_circle(GLFWwindow*,int, int, int);
void glfw_circle_partial(GLFWwindow*,int, int, int);
#endif
