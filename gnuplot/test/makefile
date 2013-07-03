
TrashDir = /home/xiaoxian/.local/share/Trash/files/


all : test.exe

test.exe: *.cpp \
		mesh_lib.a fem_lib.a gnuplot_lib.a
	g++ *.cpp -c
	g++ *.o *.a -o $@


# libraries from other directories
gnuplot_lib.a: ../*.cpp ../my_gnuplot.hpp
	g++ -c $^
	ar -cr $@ *.o
	make clean_objs
	

fem_lib.a: ../../my_fem/*.cpp
	g++ -c $^
	ar -cr $@ *.o
	make clean_objs

mesh_lib.a: ../../triangle/*.cpp
	g++ -c $^
	ar -cr $@ *.o
	make clean_objs


#clean files
clean_data:
	- mv *.ele *.edge *.node *.poly *.dat *.gnuplot *.png *.txt $(TrashDir)
clean_objs:
	- mv *.o $(TrashDir)
clean_libs:
	- mv *.a $(TrashDir)

clean:
	make clean_data
	make clean_objs

complete_clean:
	make clean
	make clean_libs
