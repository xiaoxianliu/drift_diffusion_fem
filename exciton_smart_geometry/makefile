
CC=g++
OBJS = main.o generate_mesh.o
TrashDIR = /home/xiaoxian/.local/share/Trash/files/
LOCAL_HPP = exciton_smart_geometry.hpp
CFLAGS = -larmadillo

all: main.exe

main.exe: 	mesh_libs.a fem_lib.a local_objs
	g++ *.o *.a -o main.exe $(CFLAGS)


local_objs: *.cpp
	g++ -c $^



fem_lib.a: ../my_fem/*.cpp
	g++ -c ../my_fem/*.cpp
	ar cr $@ *.o
	make clean_objs

mesh_libs.a: ../triangle/*.cpp
	g++ -c $^
	ar cr $@ *.o
	make clean_objs








clean_objs:
	- mv *.o $(TrashDIR)
clean_data:
	- mv *.node *.edge *.ele *.poly *.dat *.gnuplot *.png *.txt $(TrashDIR)
clean_libs:
	- mv *.a $(TrashDIR)
clean_exe:
	- mv *.exe $(TrashDIR)
clean:
	make clean_objs clean_data clean_exe clean_libs
