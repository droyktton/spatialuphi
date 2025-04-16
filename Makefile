CXX = nvcc


INCLUDES = -I/opt/nvidia/hpc_sdk/Linux_x86_64/23.7/math_libs/12.2/include 
FLAGS = --expt-extended-lambda -lcufft -std=c++14 -arch=sm_61 
PARAMS = -DCu=1.0 -DCphi=1.0 -DEpsilon=0.001 -DNOISESPECTRA -DDOUBLE -DMONITORCONFIGS #-DTILT=0.001

LDFLAGS = -L/opt/nvidia/hpc_sdk/Linux_x86_64/23.7/math_libs/12.2/lib64 


spatialuphi: main.cu
	$(CXX) $(FLAGS) $(PARAMS) main.cu -o spatialuphi $(LDFLAGS) $(INCLUDES) 


update_git:
	git add *.cu Makefile *.h *.sh README.md ; git commit -m "program update"; git push

clean:
	rm spatialuphi