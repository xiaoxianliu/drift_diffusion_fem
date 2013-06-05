TrashDIR = /home/xiaoxian/.local/share/Trash/files/
CFLAGS = -larmadillo 

all: fem_assemble

#fem_assemble_lib.o: fem_assemble
#	ld -r *.o -o $@

fem_assemble : *.cpp fem_assemble.hpp
	g++ -c $(CFLAGS) *.cpp 


partial_clean:
	-mv matrix*.o $(TrashDIR)
clean:
	-mv *.o $(TrashDIR)

