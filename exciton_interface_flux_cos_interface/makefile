TrashDIR = /home/xiaoxian/.local/share/Trash/files/
CFLAGS = -lm -larmadillo
MAIN_HPP = exciton_cos_interface.hpp
LOCAL_HPP = misc.hpp



##### Linking everything
all: main.exe

main.exe: 	mesh_lib.a fem_lib.a local_objs
	g++ *.o *.a -o $@ $(CFLAGS)


local_objs: *.cpp *.hpp
	g++ -c *.cpp $(CFLAGS)

##### Files from other folder
fem_lib.a: ../my_fem/*.cpp
	g++ -c $^
	ar cr $@ *.o
	make clean_objs

mesh_lib.a: ../triangle/*.cpp
	g++ -c ../triangle/*.cpp
	ar cr $@ *.o
	make clean_objs





clean_data:
	- mv *.poly *.edge *.node *.ele *.gnuplot *.png *.dat *.txt $(TrashDIR)

clean_objs:
	- mv *.o $(TrashDIR)
clean_libs:
	- mv *.a $(TrashDIR)

clean_exe:
	- mv *.exe $(TrashDIR)

clean:
	make clean_data clean_objs clean_libs clean_exe
