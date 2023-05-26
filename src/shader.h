#include "glfww.h"
#ifndef __shader 
#define __shader
static const char* vshader_src = 
  "#version 330\n"
  "layout (location = 0) in vec3 pos;\n"
  "layout (location = 1) in vec3 color;\n" 
	"out vec3 ncolor;\n"
	"void main(){\n"
  "ncolor = color;\n"
	"gl_Position = vec4(pos,1.0);\n" 
  "};";
static const char* fshader_src = 
  "#version 330\n"
  "in vec3 ncolor;\n"
	"out vec4 color;\n"
  "void main(){\n"
  //"gl_FragColor = vec4(1.0,0.0,1.0,1.0);\n"
  "gl_FragColor = vec4(ncolor,1.0);\n"
	"};";

GLuint vshader_comp(const char*);
GLuint fshader_comp(const char*);
GLuint build_shader(GLuint, GLuint);
#endif
