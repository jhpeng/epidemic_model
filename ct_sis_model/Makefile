# 'make'	build executable file
# 'make clean'	removes all *.o and executalbe file

# define the C compiler
CC	= gcc

# define any compile-time flags
#CFLAGS = -Wall -g -fPIC -O0 -std=c99
#CFLAGS = -Wall -g -pg -fPIC -O0 -std=c99
CFLAGS = -Wall -fPIC -O3 -std=c99

# define openmp flags
OPENMP  = -fopenmp
CUOPENMP  = -Xcompiler -fopenmp

# define the direction containing header file
INCLUDES= -I/usr/local/include -I./

# define the library path
LFLAGS	= -L/usr/local/lib

# define any libraries to link to executable
LIBS	= -lm -lgsl -lgslcblas

# define the C object files
OBJS	=  graphs.o kernel.o main.o


#define the directory for object
OBJSDIR = object

# define the executable file
MAIN	= exe

# define the main srouce file
#MAINSRC = projects/parallel/jq3_ladder_square/main.c


exe: $(OBJS)
	$(CC) $(CFLAGS) -o $(MAIN) $(OBJS) $(LIBS) $(LFLAGS) $(INCLUDES) 

graphs: graphs.o
	$(CC) $(CFLAGS) -o graphs graphs.o $(LIBS) $(LFLAGS) $(INCLUDES) 

%.o: %.c
	$(CC) $(CFLAGS) $(INCLUDES) -c $^



# clean the executable file and object files
clean:
	$(RM) $(METR_OBJS) $(OBJS) $(MAIN)
