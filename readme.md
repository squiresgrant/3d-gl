
optimized build : `clang ./src/*.c -lm -lglfw -lGL -fno-trapping-math`

#known issues
- glfw causes memory leak (not lost, related to video drivers)
