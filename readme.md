
optimized build : `clang ./src/*.c -lm -lglfw -lGL -lGLEW -fno-trapping-math`

#known issues
- glfw & glew causes memory leak (not lost, related to video drivers)
