TrashDIR = /home/xiaoxian/.local/share/Trash/files/
CFLAGS = -larmadillo 

all: fem_assemble

fem_assemble : *.cpp fem_assemble.hpp
	g++ -c $(CFLAGS) *.cpp 


partial_clean:
	-mv matrix*.o $(TrashDIR)
clean:
	-mv *.o $(TrashDIR)

