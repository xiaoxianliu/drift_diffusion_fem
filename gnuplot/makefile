
TrashDir = /home/xiaoxian/.local/share/Trash/files/


all : gnuplot_libs

gnuplot_libs: *.cpp
	g++ -c $^



clean_objs:
	- mv *.o $(TrashDir)

clean:
	make clean_objs
