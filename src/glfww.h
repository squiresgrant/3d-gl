#include <stdio.h>
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
GLFWwindow* glfw_init();
#define glfw_load(w) glfwSwapBuffers(w);
void glfw_loop(GLFWwindow*window);
int get_h();
int get_w();
void glfw_pixel_partial(GLFWwindow*wi,int x, int y);
void glfw_clear(GLFWwindow*w);
void refresh_size(GLFWwindow*);
#define glfw_pixel(wi,x,y)\
  glBegin(GL_POINTS);\
  glfw_pixel_partial(wi,x,y);\
  glEnd();
void glfw_circle(GLFWwindow* w,int x, int y, int r);
void glfw_circle_partial(GLFWwindow* w,int x, int y, int r);
#endif
