
optimized build : `clang ./src/*.c -lm -lglfw -lGL -lGLEW -fno-trapping-math`

optional flags : -Dskip_memory_trace -Dstfu

# todo

- find a better epsilon

- fix the fps counter, its probably wrong 

# known issues

- glfw & glew causes memory leak (not lost, related to video drivers)
