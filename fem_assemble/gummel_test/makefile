MeshLibDir = ../../triangle
FemLibDir = ..
TrashDir = /home/xiaoxian/.local/share/Trash/files




all: assemble_gummel_matrix.exe

assemble_gummel_matrix.exe: mesh_lib.a fem_lib.a local_objs
	g++ *.o *.a -o $@ -larmadillo

local_objs: *.cpp
	g++ *.cpp -c #-larmadillo

fem_lib.a: $(FemLibDir)/*.cpp
	g++ -c $(FemLibDir)/*.cpp
	ar cr $@ *.o
	make clean_objs

mesh_lib.a: $(MeshLibDir)/*.cpp
	g++ -c $(MeshLibDir)/*.cpp
	ar cr $@ *.o
	make clean_objs


clean_data:
	-mv *.poly *.ele *.node *.edge *.png *.gnuplot *.dat $(TrashDir)
clean_objs:
	-mv *.o $(TrashDir)
clean_libs:
	-mv *.a $(TrashDir)
clean_exes:
	-mv *.exe $(TrashDir)
clean:
	make clean_data
	make clean_objs
	make clean_libs
	make clean_exes
