TrashDIR = /home/xiaoxian/.local/share/Trash/files/
CFLAGS = -larmadillo 
LOCAL_HPP = test.hpp

all: main.exe

quick_build:
	g++ *.o *.a -o main.exe $(CFLAGS)

main.exe: mesh_lib.a fem_assemble.a \
		main.o solve.o functionals.o
	g++ *.o *.a -o $@ $(CFLAGS)
#	make clean_objs

main.o: main.cpp
	g++ -c main.cpp $(CFLAGS)

functionals.o: functionals.cpp $(LOCAL_HPP)
	g++ -c functionals.cpp
solve.o: solve.cpp $(LOCAL_HPP)
	g++ -c $< $(CFLAGS)


fem_assemble.a: ../*.cpp
	g++ -c ../*.cpp $(CFLAGS)
	ar cr $@ *.o
	make clean_objs

mesh_lib.a: ../../triangle/*.cpp
	g++ -c $^
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
