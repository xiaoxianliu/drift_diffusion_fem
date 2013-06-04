
CC=g++
OBJS = main.o generate_mesh.o
TrashDIR = /home/xiaoxian/.local/share/Trash/files/

all: main.exe

main.exe: fem_assemble.o mesh_lib.o generate_mesh.o functions.o main.o
	g++ *.o -o main.exe
	make clean_objs

main.o: main.cpp exciton_single_geometry.hpp
	g++ -c $<

functions.o: functions.cpp
	g++ -c $<

generate_mesh.o: generate_mesh.cpp exciton_single_geometry.hpp
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
