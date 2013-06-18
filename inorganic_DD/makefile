TrashDIR = /home/xiaoxian/.local/share/Trash/files/
CFLAGS = -larmadillo 


all: main.exe

main.exe: mesh_libs.a fem_assemble.a local_objs
	g++ *.o *.a -o $@ $(CFLAGS)

local_objs: *.cpp *.hpp
	g++ -c *.cpp $(CFLAGS)


fem_assemble.a: ../fem_assemble/*.cpp
	g++ -c $^ $(CFLAGS)
	ar cr $@ *.o
	make clean_objs

mesh_libs.a: ../triangle/*.cpp
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
