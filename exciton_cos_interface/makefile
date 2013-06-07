TrashDIR = /home/xiaoxian/.local/share/Trash/files/
CFLAGS = -lm -larmadillo
LOCAL_HPP = exciton_cos_interface.hpp



##### Linking everything
all: main.exe

main.exe: 	mesh_libs fem_assemble\
		generate_mesh.o state_eq.o adjoint_eq.o shape_gradient.o plotting.o main.o
	g++ *.o -o $@ $(CFLAGS)

quick_build: 	
	g++ *.o -o main.exe $(CFLAGS)





##### Compile each individual module
main.o: main.cpp $(LOCAL_HPP)
	g++ -c main.cpp

# mesh
generate_mesh.o: generate_mesh.cpp $(LOCAL_HPP)
	g++ -c generate_mesh.cpp $(CFLAGS)



# state equation
state_eq.o: state_eq.cpp $(LOCAL_HPP)
	g++ -c state_eq.cpp $(CFLAGS)
# adjoint equation
adjoint_eq.o: adjoint_eq.cpp $(LOCAL_HPP)
	g++ -c adjoint_eq.cpp $(CFLAGS)

shape_gradient.o: shape_gradient.cpp $(LOCAL_HPP)
	g++ -c shape_gradient.cpp $(CFLAGS)

plotting.o: plotting.cpp $(LOCAL_HPP)
	g++ -c plotting.cpp $(CFLAGS)





##### Files from other folder
fem_assemble: ../fem_assemble/*.cpp
	g++ -c ../fem_assemble/*.cpp
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
