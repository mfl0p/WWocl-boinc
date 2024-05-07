CC = g++
LD = $(CC)

.SUFFIXES:
.SUFFIXES: .o .c .h .cl .cpp

APP = WWocl-win64

SRC = main.cpp wwcpu.cpp simpleCL.c simpleCL.h kernels/clearn.cl kernels/clearresult.cl kernels/getsegprimes.cl kernels/wieferich.cl kernels/wallsunsun.cl
KERNEL_HEADERS = kernels/clearn.h kernels/clearresult.h kernels/getsegprimes.h kernels/wieferich.h kernels/wallsunsun.h
OBJ = main.o simpleCL.o wwcpu.o

OCL_INC = 
OCL_LIB = -L . -lOpenCL -lprimesieve

OCL_INC = 
OCL_LIB = OpenCL.dll libprimesievewin.a

BOINC_DIR = C:/mingwbuilds/boinc
BOINC_INC = -I$(BOINC_DIR)/lib -I$(BOINC_DIR)/api -I$(BOINC_DIR) -I$(BOINC_DIR)/win_build
BOINC_LIB = -L$(BOINC_DIR)/lib -L$(BOINC_DIR)/api -L$(BOINC_DIR) -lboinc_opencl -lboinc_api -lboinc -lstdc++

CFLAGS  = -I . -I kernels -O3 -m64 -Wall
LDFLAGS = $(CFLAGS) -static

all : clean $(APP)

$(APP) : $(OBJ)
	$(LD) $(LDFLAGS) $^ $(OCL_LIB) $(BOINC_LIB) -o $@

main.o : $(SRC) $(KERNEL_HEADERS)
	$(CC) $(CFLAGS) $(OCL_INC) $(BOINC_INC) -c -o $@ main.cpp

wwcpu.o : $(SRC)
	$(CC) $(CFLAGS) $(OCL_INC) $(BOINC_INC) -c -o $@ wwcpu.cpp

simpleCL.o : $(SRC)
	$(CC) $(CFLAGS) $(OCL_INC) $(BOINC_INC) -c -o $@ simpleCL.c

.cl.h:
	perl cltoh.pl $< > $@

clean :
	del *.o
	del kernels\*.h
	del $(APP).exe

