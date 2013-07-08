TrashDir = /home/xiaoxian/.local/share/Trash/files/

all:main.exe

main.exe: mesh_lib.a main.cpp
	g++ main.cpp -c
	g++ *.o *.a -o $@

mesh_lib.a: ../*.cpp
	g++ -c $^
	ar cr $@ *.o
	make clean_objs


clean_arc:
	-mv *.a $(TrashDir)
clean_objs:
	-mv *.o $(TrashDir)
clean_exes:
	-mv *.exe $(TrashDir)
clean_data:
	-mv *.poly *.ele *.node *.edge *.png *.gnuplot *.dat *.txt $(TrashDir)

clean:
	make clean_data clean_objs clean_exes clean_arc
