TrashDIR = /home/xiaoxian/.local/share/Trash/files/
CFLAGS = -lm
LOCAL_HPP = exciton_cos_interface.hpp

all: main.exe

main.exe: 	mesh_libs \
		generate_mesh.o main.o
	g++ *.o -o $@
	make clean_objs



main.o: main.cpp $(LOCAL_HPP)
	g++ -c main.cpp

generate_mesh.o: generate_mesh.cpp $(LOCAL_HPP)
	g++ -c generate_mesh.cpp $(CFLAGS)




mesh_libs: ../triangle/*.cpp
	g++ -c ../triangle/*.cpp





clean_data:
	- mv *.poly *.edge *.node *.ele *.gnuplot *.png *.dat *.txt $(TrashDIR)

clean_objs:
	- mv *.o $(TrashDIR)

clean_exe:
	- mv *.exe $(TrashDIR)

complete_clean:
	make clean_data clean_objs clean_exe
