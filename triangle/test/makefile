TrashDIR = /home/xiaoxian/.local/share/Trash/files/

all: test.exe

test.exe: mesh_lib test.o
	g++ *.o -o $@
	make clean_objs

test.o: test_mesh.cpp
	g++ -c test_mesh.cpp

mesh_lib: ../*.cpp
	g++ -c ../*.cpp 



clean_data:
	- mv *.poly *.ele *.node *.edge *.png *.gnuplot *.dat *.txt $(TrashDIR)
clean_objs:
	- mv *.o $(TrashDIR)
clean_exe:
	- mv *.exe $(TrashDIR)
complete_clean:
	make clean_data
	make clean_objs
	make clean_exe
