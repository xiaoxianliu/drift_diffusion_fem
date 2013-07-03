
TrashDir = /home/xiaoxian/.local/share/Trash/files/


all : gnuplot_lib.a

gnuplot_lib.a: *.cpp my_gnuplot.hpp
	g++ -c $^
	ar -cr $@ *.o



clean_objs:
	- mv *.o $(TrashDir)
clean_libs:
	- mv *.a $(TrashDir)

clean:
	make clean_objs
	make clean_libs
