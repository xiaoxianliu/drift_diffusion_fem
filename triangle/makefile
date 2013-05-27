
CC=g++
OBJS=read_nodes.o read_edges.o read_elements.o compute_topology.o vector_operations.o geometric.o \
	write_gnuplot.o write_polyfile.o

# Two outputs: 
# 1. text.exe is for a test run of "triangle"; 
# 2. mesh_lib.o: an object file to be called by other programs.

all : test.exe mesh_lib.o

test.exe : $(OBJS) test_mesh.o
	$(CC) $(OBJS) test_mesh.o -o test.exe

mesh_lib.o: $(OBJS)
	ld -r $(OBJS) -o $@


# compile each .cpp source file seperately
test_mesh.o: test_mesh.cpp
	$(CC) -c test_mesh.cpp -o test_mesh.o


read_nodes.o: read_nodes.cpp mesh.hpp
	$(CC) -c read_nodes.cpp -o read_nodes.o

read_edges.o: read_edges.cpp mesh.hpp
	$(CC) -c read_edges.cpp -o read_edges.o

read_elements.o: read_elements.cpp mesh.hpp
	$(CC) -c read_elements.cpp -o read_elements.o

compute_topology.o: compute_topology.cpp 
	$(CC) -c compute_topology.cpp -o compute_topology.o

vector_operations.o: vector_operations.cpp
	$(CC) -c vector_operations.cpp -o vector_operations.o

geometric.o: geometric.cpp
	$(CC) -c geometric.cpp -o geometric.o
write_gnuplot.o: write_gnuplot.cpp
	$(CC) -c write_gnuplot.cpp -o write_gnuplot.o
write_polyfile.o: write_polyfile.cpp
	$(CC) -c write_polyfile.cpp -o $@


clean_objs:
	mv $(OBJS) test_mesh.o $(HOME)/.local/share/Trash/files/
clean_exe:
	mv *.exe $(HOME)/.local/share/Trash/files/
clean_data:
	mv *.1.* *.dat *.gnuplot *.png *.txt $(HOME)/.local/share/Trash/files/
complete_clean:
	make clean_objs clean_exe clean_data
