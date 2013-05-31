
CC=g++
OBJS =	new_mesh.o read_nodes.o read_edges.o read_elements.o compute_topology.o vector_operations.o geometric.o \
	submesh.o plot_mesh.o write_polyfile.o 

# Two outputs: 
# 1. text.exe is for a test run of "triangle"; 
# 2. mesh_lib.o: an object file to be called by other programs.

all : test.exe mesh_lib.o

test.exe : $(OBJS) test_mesh.o
	$(CC) $(OBJS) test_mesh.o -o test.exe

mesh_lib.o: $(OBJS)
	ld -r $(OBJS) -o $@


# compile each .cpp source file seperately
%.o: %.cpp mesh.hpp
	$(CC) -c $<


clean_objs:
	mv $(OBJS) test_mesh.o $(HOME)/.local/share/Trash/files/
clean_exe:
	mv *.exe $(HOME)/.local/share/Trash/files/
clean_data:
	mv *.node *.edge *.ele *.poly *.dat *.gnuplot *.png *.txt $(HOME)/.local/share/Trash/files/
complete_clean:
	make clean_objs clean_exe clean_data
