# do not use -O3 (inline)

CXX = g++
CFLAGS = -O3 -g -std=c++11 -c


INC  = -I/home/hac210/Dropbox/packages/01.M2DO_code/M2DO_LSM/include -I/home/hac210/Dropbox/packages/eigen3/Eigen
LIB = -lblas -llapack -L /usr/local/lib -lnlopt -lmetis -lma57 -lgfortran -larpack -lblas -lmetis

TARGET = optim_base_adj
DEBUG  = optim_debug
OBJECTS = optim_EICR_adjoint.o f_nlgeom.o FEA_hy.o CoreElement.o  EICR_shell.o f_DKT_OPT.o  m2do_lin_sensitivity.o # Sensitivity.o 

all : $(TARGET)	

debug: CFLAGS = -g -W -std=c++11 -c -D__DEBUGFLAG__
debug: $(DEBUG)

$(TARGET): $(OBJECTS) # $^: same as $(OBJECT)
		$(CXX) -o $@ $^ $(INC) $(LIB)

$(DEBUG): $(OBJECTS) # $^: same as $(OBJECT)
		$(CXX) -o $@ $^ $(INC) $(LIB)

%.o : %.cpp # rule of *o files
		$(CXX) $(CFLAGS) -o $@ $< $(INC)

clean: 
	rm -rf *.o *.txt *.vtk $(TARGET) $(DEBUG)

Reset:
	rm -rf *.txt *.vtk
