CC=g++

TRIANGLE_DIR = /home/xiaoxian/programs_cpp/triangle/
MESHLIB_OBJ = /home/xiaoxian/programs_cpp/triangle/mesh_lib.o

TRASH_DIR = /home/xiaoxian/.local/share/Trash/files/

OBJS = test.o mesh.o solve.o linear_system.o gnuplot_solution.o	functionals.o \
       	polyfile_smart_interface.o ### Only local ".o" files are included in OBJS!!!!

CFLAGS = -I$(TRIANGLE_DIR) -larmadillo

all: test.exe

test.exe: $(OBJS)
	$(CC) $(OBJS) $(CFLAGS) -o $@


test.o: test.cpp
	$(CC) $(CFLAGS) -c test.cpp
solve.o: solve.cpp
	$(CC) $(CFLAGS) -c solve.cpp
functionals.o: functionals.cpp
	$(CC) $(CFLAGS) -c functionals.cpp
gnuplot_solution.o: gnuplot_solution.cpp
	$(CC) $(CFLAGS) -c gnuplot_solution.cpp
linear_system.o: linear_system.cpp
	$(CC) $(CFLAGS) -c linear_system.cpp
mesh.o: triangle_mesh.o
	ld -r triangle_mesh.o $(MESHLIB_OBJ) -o $@
triangle_mesh.o: triangle_mesh.cpp
	$(CC) $(CFLAGS) -c triangle_mesh.cpp
polyfile_smart_interface.o: polyfile_smart_interface.cpp
	$(CC) $(CFLAGS) -c polyfile_smart_interface.cpp -o $@



clean_objs:
	mv $(OBJS) triangle_mesh.o $(TRASH_DIR)
clean_data:
	mv *.poly *.node *.ele *.edge *.png *.dat *.gnuplot $(TRASH_DIR)
clean_exe:
	mv *.exe $(TRASH_DIR)
complete_clean:	
	make clean_objs
	make clean_data
	make clean_exe

