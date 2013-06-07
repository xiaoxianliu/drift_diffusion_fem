TrashDIR = /home/xiaoxian/.local/share/Trash/files/
CFLAGS = -lm -larmadillo
MAIN_HPP = exciton_cos_interface.hpp
LOCAL_HPP = misc.hpp



##### Linking everything
all: main.exe

main.exe: 	mesh_libs fem_assemble\
		generate_mesh.o misc.o functionals.o\
		state_eq.o adjoint_eq.o shape_gradient.o plotting.o main.o
	g++ *.o -o $@ $(CFLAGS)

quick_build: 	
	g++ *.o -o main.exe $(CFLAGS)





##### Compile each individual module
main.o: main.cpp $(MAIN_HPP)
	g++ -c main.cpp

# mesh
generate_mesh.o: generate_mesh.cpp $(MAIN_HPP)
	g++ -c generate_mesh.cpp $(CFLAGS)


# functions to be called by other programs
misc.o:	misc.cpp misc.hpp
	g++ -c misc.cpp


# state equation
state_eq.o: state_eq.cpp $(MAIN_HPP) $(LOCAL_HPP)
	g++ -c state_eq.cpp $(CFLAGS)
# adjoint equation
adjoint_eq.o: adjoint_eq.cpp $(MAIN_HPP) $(LOCAL_HPP)
	g++ -c adjoint_eq.cpp $(CFLAGS)

# shape gradient
shape_gradient.o: shape_gradient.cpp $(MAIN_HPP) $(LOCAL_HPP)
	g++ -c shape_gradient.cpp $(CFLAGS)

functionals.o: functionals.cpp $(LOCAL_HPP)
	g++ -c functionals.cpp

# plotting
plotting.o: plotting.cpp $(MAIN_HPP)
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
