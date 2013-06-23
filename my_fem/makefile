TrashDIR = /home/xiaoxian/.local/share/Trash/files/
CFLAGS = -larmadillo 

all: my_fem

my_fem : *.cpp my_fem.hpp
	g++ -c $(CFLAGS) *.cpp 


partial_clean:
	-mv matrix*.o $(TrashDIR)
clean:
	-mv *.o $(TrashDIR)

