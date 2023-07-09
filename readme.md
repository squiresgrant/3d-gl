wanted to learn about 3d rendering and how rendering engines work at the lower levels

optimized build : `clang ./src/*.c -lm -lglfw -lGL -lGLEW -fno-trapping-math -ffast-math`

build flags : -Dmemory_trace -Dstfu -D__debug -Dmemory_count

once these todo tasks are done, i will most likely drop this project, drawing everything with triangles was a fun challenge but,
it will simply make other things (ie, color, collisions, z-buffering) impossible or more difficult and resource demanding.

# todo

- z-buffer

- (not likely) moving parts 

- (maybe) auto separate colliding polygons 

- test transparency 

# known issues

- ordering is fucked (just for now hopefully)

- [limitation] no colliding polygons (todo)
