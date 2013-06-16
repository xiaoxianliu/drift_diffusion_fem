TrashDIR = /home/xiaoxian/.local/share/Trash/files/

all: test.exe

test.exe: mesh_lib.a test.o
	g++ *.o *.a -o $@
#	make clean_objs

test.o: test_refine.cpp
	g++ -c test_refine.cpp

mesh_lib.a: ../*.cpp
	g++ -c $^
	ar cr $@ *.o
	make clean_objs



clean_data:
	- mv *.poly *.ele *.node *.edge *.area *.png *.gnuplot *.dat *.txt $(TrashDIR)
clean_objs:
	- mv *.o $(TrashDIR)
clean_libs:
	- mv *.a $(TrashDIR)
clean_exe:
	- mv *.exe $(TrashDIR)
clean:
	make clean_data
	make clean_objs
	make clean_libs
	make clean_exe
