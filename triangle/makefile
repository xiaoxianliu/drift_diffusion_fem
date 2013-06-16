
TrashDir = /home/xiaoxian/.local/share/Trash/files/


all : mesh_libs

mesh_libs: *.cpp
	g++ -c $^



clean_objs:
	- mv *.o $(TrashDir)

clean:
	make clean_objs
