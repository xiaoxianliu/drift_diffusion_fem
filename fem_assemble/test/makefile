TrashDIR = /home/xiaoxian/.local/share/Trash/files/
CFLAGS = -larmadillo 

all: test.exe

test.exe: mesh_lib.o fem_assemble.o test.o
	g++ *.o -o test.exe $(CFLAGS)
	make clean_objs

test.o: test.cpp
	g++ -c test.cpp $(CFLAGS)

fem_assemble.o: ../*.cpp
	g++ -c ../*.cpp $(CFLAGS)

mesh_lib.o: ../../triangle/*.cpp
	g++ -c ../../triangle/*.cpp


clean_data:
	-mv *.poly *.ele *.node *.edge *.png *.gnuplot *.dat $(TrashDIR)
clean_objs:
	-mv *.o $(TrashDIR)
clean_exe:
	-mv *.exe $(TrashDIR)
complete_clean:
	make clean_data
	make clean_objs
	make clean_exe
