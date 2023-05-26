#include "shader.h"

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
