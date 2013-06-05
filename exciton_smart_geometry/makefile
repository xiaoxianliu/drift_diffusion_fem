
CC=g++
OBJS = main.o generate_mesh.o
TrashDIR = /home/xiaoxian/.local/share/Trash/files/
LOCAL_HPP = exciton_smart_geometry.hpp
CFLAGS = -larmadillo

all: main.exe

main.exe: 	fem_assemble.o mesh_lib.o \
		generate_mesh.o functions.o boundary_conditions.o plotting.o \
		main.o
	g++ *.o -o main.exe $(CFLAGS)
	make clean_objs



generate_mesh.o: generate_mesh.cpp $(LOCAL_HPP)
	g++ -c $<

functions.o: functions.cpp
	g++ -c $<

boundary_conditions.o: boundary_conditions.cpp $(LOCAL_HPP)
	g++ -c $<

plotting.o: plotting.cpp $(LOCAL_HPP)
	g++ -c $<

main.o: main.cpp $(LOCAL_HPP)
	g++ -c $<




fem_assemble.o: ../fem_assemble/*.cpp
	g++ -c ../fem_assemble/*.cpp

mesh_lib.o: ../triangle/*.cpp
	g++ -c ../triangle/*.cpp








clean_objs:
	- mv *.o $(TrashDIR)
clean_data:
	- mv *.node *.edge *.ele *.poly *.dat *.gnuplot *.png *.txt $(TrashDIR)
clean_exe:
	- mv *.exe $(TrashDIR)
complete_clean:
	make clean_objs clean_data clean_exe
