TrashDIR = /home/xiaoxian/.local/share/Trash/files/
CFLAGS = -larmadillo 

all: test.exe

test.exe: mesh_lib.a fem_assemble.a test.o
	g++ *.o *.a -o test.exe $(CFLAGS)


test.o: test.cpp
	g++ -c test.cpp $(CFLAGS)

fem_assemble.a: ../*.cpp
	g++ -c ../*.cpp $(CFLAGS)
	ar cr $@ *.o
	make clean_objs

mesh_lib.a: ../../triangle/*.cpp
	g++ -c ../../triangle/*.cpp
	ar cr $@ *.o
	make clean_objs

clean_data:
	-mv *.poly *.ele *.node *.edge *.png *.gnuplot *.dat $(TrashDIR)
clean_objs:
	-mv *.o $(TrashDIR)
clean_libs:
	-mv *.a $(TrashDIR)
clean_exe:
	-mv *.exe $(TrashDIR)
clean:
	make clean_data
	make clean_objs
	make clean_libs
	make clean_exe
