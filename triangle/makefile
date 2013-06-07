
CC=g++
OBJS =	new_mesh.o read_nodes.o read_edges.o read_elements.o compute_topology.o vector_operations.o geometric.o \
	submesh.o plot_mesh.o write_polyfile.o 

# Output: 
#  mesh_lib.o: an object file to be called by other programs.

all : mesh_lib.o

mesh_lib.o: $(OBJS)
	ld -r $(OBJS) -o $@


# compile each .cpp source file seperately
%.o: %.cpp mesh.hpp
	$(CC) -c $<


clean_objs:
	-mv *.o $(HOME)/.local/share/Trash/files/

complete_clean:
	-make clean_objs
