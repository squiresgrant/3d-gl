wanted to learn about 3d rendering and how rendering engines work at the lower levels

(mostly) optimized build : `clang ./src/*.c -lm -lglfw -lGL -lGLEW -fno-trapping-math -ffast-math -L/usr/local/lib -I/usr/local/include`

debug/other build flags : -Dmemory_trace -Dstfu -D__debug -Dmemory_count

windows? you're on your own (though it should work)

once these todo tasks are done, i will most likely drop this project, drawing everything with triangles was a fun challenge but,
it will simply make other things (ie, color, collisions, z-buffering) impossible or more difficult and resource demanding.

# todo

- sort stuff into files better:3

- z-buffer

- (not likely) moving parts 

- (maybe) auto separate colliding polygons 

- fix transparency 

# known issues

- ordering is now only barely fucked (WITH MEMORY LEAKS, for now)

- [limitation] no colliding polygons (todo)
