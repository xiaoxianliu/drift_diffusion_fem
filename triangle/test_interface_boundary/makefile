TrashDIR = /home/xiaoxian/.local/share/Trash/files/

all: test.exe

test.exe: mesh_lib.a local_objs
	g++ *.o *.a -o $@
#	make clean_objs

local_objs: *.cpp
	g++ -c $^

mesh_lib.a: ../*.cpp
	g++ -c $^
	ar cr $@ *.o
	make clean_objs



clean_data:
	- mv *.poly *.ele *.node *.edge *.png *.gnuplot *.dat *.txt $(TrashDIR)
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
