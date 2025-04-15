#include<thrust/device_vector.h>
#include<thrust/for_each.h>
#include <thrust/iterator/zip_iterator.h>
#include <thrust/tuple.h>
#include<thrust/reduce.h>
#include<fstream>
#include<cstdlib>
#include <thrust/transform_reduce.h>
#include <thrust/transform.h>
#include <cufft.h>
#include "cutil.h"
#include <chrono>
#include <iomanip>

#ifndef C2
#define C2 1.0   // u elastic constant
#endif

#ifndef C4
#define C4 1.0   // u elastic constant
#endif

#ifdef DOUBLE
typedef double real;
typedef cufftDoubleComplex complex;
#else
typedef float real;
typedef cufftComplex complex;
#endif

// file to log parameters of the run
std::ofstream logout("logfile.dat");


// main class:
class cuerda{

    public:
    cuerda(unsigned long _L, real _dt):L(_L),dt(_dt),fourierCount(0)
    {
        // interface position
        u.resize(L);
        
        // interface forces
        force_u.resize(L);
 
        // flat initial condition
        thrust::fill(u.begin(),u.end(),real(0.0));
        
        // plans for the interface structure factor
        #ifdef DOUBLE
        CUFFT_SAFE_CALL(cufftPlan1d(&plan_r2c,L,CUFFT_D2Z,1));
        #else
        CUFFT_SAFE_CALL(cufftPlan1d(&plan_r2c,L,CUFFT_R2C,1));
        #endif

	    int Lcomp=L/2+1;
	    Fou_u.resize(Lcomp); // interface position in fourier space

        acum_Sofq_u.resize(L); // average structure factor
        inst_Sofq_u.resize(L); // instantaneous structure factor

        // initialization of structure factors   
        thrust::fill(acum_Sofq_u.begin(),acum_Sofq_u.end(),real(0.0));

        #ifdef DEBUG
        std::cout << "L=" << L << ", dt=" << dt << std::endl;
        #endif
    }

    void reset_acum_Sofq(){
        thrust::fill(acum_Sofq_u.begin(),acum_Sofq_u.end(),real(0.0));
    }

    // returns the center of mass position
    real center_of_mass()
    {
        //DANGER: large sum over large numbers
        real cmu = thrust::reduce(u.begin(),u.end(),real(0.0))/L;
        return cmu;
    }

    // returns the center of mass velocity
    real center_of_mass_velocity()
    {
        //SAFE: velocities are bounded
        real vcmu = thrust::reduce(force_u.begin(),force_u.end(),real(0.0))/L;
        return vcmu;
    }

    // computes the instantaneous and acumulated structure factor
    void fourier_transform(){

        real *raw_u = thrust::raw_pointer_cast(&u[0]); 
        complex *raw_fou_u = thrust::raw_pointer_cast(&Fou_u[0]); 

        #ifdef DOUBLE
        CUFFT_SAFE_CALL(cufftExecD2Z(plan_r2c, raw_u, raw_fou_u));
        #else
	    CUFFT_SAFE_CALL(cufftExecR2C(plan_r2c, raw_u, raw_fou_u));
        #endif

        thrust::for_each(
            thrust::make_zip_iterator(
                thrust::make_tuple(Fou_u.begin(),acum_Sofq_u.begin(),inst_Sofq_u.begin())
            ),
            thrust::make_zip_iterator(
                thrust::make_tuple(Fou_u.end(),acum_Sofq_u.end(),inst_Sofq_u.end())
            ),
            [=] __device__ (thrust::tuple<complex,real &,real &> t)
            {
                complex fu=thrust::get<0>(t);
                thrust::get<1>(t) += fu.x*fu.x + fu.y*fu.y; 
                thrust::get<2>(t) = fu.x*fu.x + fu.y*fu.y;
            }
        );
        fourierCount++;
    }

    // computes the center of mass, the variance (roughness), and the leading and receding points of the interface 
    thrust::tuple<real, real, real, real> roughness()
    {
        // CHECK large numbers
        real cmu = thrust::reduce(u.begin(),u.end(),real(0.f),thrust::plus<real>())/real(L);
	    
	    real u0=u[0]; 
        real maxu = thrust::reduce(u.begin(),u.end(),u0,thrust::maximum<real>());
        real minu = thrust::reduce(u.begin(),u.end(),u0,thrust::minimum<real>());

        real cmu2 = 
        thrust::transform_reduce(
            u.begin(),u.end(),
            [=] __device__ (real x){
                return (x-cmu)*(x-cmu);
            },
            real(0.f),
            thrust::plus<real>()
        )/real(L);

        return thrust::make_tuple(cmu,cmu2,maxu,minu);
    }

    // just compute and prints center of mass in out stream
    void print_center_of_mass(std::ofstream &out)
    {
        real cm=center_of_mass();    
        out << cm << std::endl;
    }

    // rescale all position in order to avoid large numbers
    void rescale()
    {
        real cmu=center_of_mass();

        thrust::transform(u.begin(),u.end(),u.begin(),
        [=] __device__ (real u){
            return u-cmu;
        }
        );
    };

    // print roughness results
    void print_roughness(std::ofstream &out, real t, real &acum)
    {
        thrust::tuple<real,real,real,real> cm = roughness();
        real vcm=center_of_mass_velocity();

        //get cmu,cmu2,maxu,minu
        real cmu = thrust::get<0>(cm);
        real cmu2 = thrust::get<1>(cm);
        real maxu = thrust::get<2>(cm);
        real minu = thrust::get<3>(cm);
        acum+=cmu2;

        out << t << " " << vcm << " " << cmu << " " << " " << cmu2 << " " << maxu << " " << minu << std::endl;
    }

    // Computes the forces and advance one time step using Euler method
    void update(){

        real *raw_u = thrust::raw_pointer_cast(&u[0]); 

        // not elegant solution for lambda...
        real dt_=dt;
        unsigned long L_ = L;

        // Forces
        thrust::for_each(
            thrust::make_zip_iterator(
                thrust::make_tuple(force_u.begin(),thrust::make_counting_iterator((unsigned long)0))        
            ),
            thrust::make_zip_iterator(
                thrust::make_tuple(force_u.end(),thrust::make_counting_iterator((unsigned long)L))        
            ),
            [=] __device__ (thrust::tuple<real &,unsigned long> t)
            {
                unsigned long i=thrust::get<1>(t);
                unsigned long ileft = (i-1+L_)%L_;
                unsigned long iright = (i+1)%L_;

                real uleft = raw_u[ileft];
                real uright = raw_u[iright];
                
                // to imposed tilted boundary conditions
                #ifdef TILT
                if(i==0) {
                    uleft -= L_*TILT;
                }  
                if(i==L_-1){
                    uright += L_*TILT;
                }  
                #endif
                
                real lap_u = uright + uleft - 2.0*raw_u[i];
                
                thrust::get<0>(t) = C2*lap_u;
            } 
        );

        #ifdef DEBUG
        std::cout << "updating" << std::endl;
        #endif

        // Euler step
        thrust::for_each(
            thrust::make_zip_iterator(
                thrust::make_tuple(u.begin(), force_u.begin())        
            ),
            thrust::make_zip_iterator(
                thrust::make_tuple(u.end(),force_u.end())        
            ),
            [=] __device__ (thrust::tuple<real &,real> t)
            {
                thrust::get<0>(t) = thrust::get<0>(t) + dt_*thrust::get<1>(t);
            } 
        );
    };

    // print the whole configuration to a file
    void print_config(std::ofstream &out){
        real cm = center_of_mass();

        for(int i=0;i<L;i++){
            out << u[i] << " " << cm << "\n";
        }
        out << "\n" << std::endl;
    };

    // prints the whole averaged structure factor to a file
    void print_sofq(std::ofstream &out){
        for(int i=0;i<L;i++){
            out << acum_Sofq_u[i]/fourierCount << "\n";
        }
        out << "\n" << std::endl;
    };

    // prints the instantaneous structure factor to a file
    void print_inst_sofq(std::ofstream &out){
        for(int i=0;i<L;i++){
            out << inst_Sofq_u[i] << "\n";
        }
        out << "\n" << std::endl;
    };

    // variables and arrays of the class
    private:
        real dt;
        unsigned long L;
        
        real f0;
        thrust::device_vector<real> u;
        thrust::device_vector<real> force_u;

        int fourierCount;
        cufftHandle plan_r2c;
        thrust::device_vector<complex> Fou_u;

        thrust::device_vector<real> acum_Sofq_u;
	    thrust::device_vector<real> inst_Sofq_u;
	    
};

int main(int argc, char **argv){
    // Get the current CUDA device
    int device;
    cudaGetDevice(&device);

    // Get the properties of the current CUDA device
    cudaDeviceProp deviceProp;
    cudaGetDeviceProperties(&deviceProp, device);

    std::ofstream confout("conf.dat");
    confout << "#u[i]" << " " << "cmu" << "\n";

    std::ofstream sofqout("sofq.dat");
    sofqout << "#av_Sofq_u[i]" << "\n";

    std::ofstream instsofqout("inst_sofq.dat");
    instsofqout << "#inst_Sofq_u[i]" << "\n";

    std::ofstream cmout("cm.dat");
    cmout << "#t" << " " << "velu" << " " << "cmu" << " " << "cmu2" << " " << "maxu" << " " << "minu" << std::endl;

    std::ofstream lastconfout("lastconf.dat");
    lastconfout << "#u[i]" << " " << "cmu" << "\n";

    unsigned int L=atoi(argv[1]);
    real f0=atof(argv[2]);
    unsigned long Nrun = atoi(argv[3]);
    unsigned long Neq = unsigned(Nrun/2.0);
    real dt=0.1;

    unsigned int seed=1234;
    if(argc==5) seed=(unsigned int)atoi(argv[5]); 
    srand(seed);
    cuerda C(L,dt);


    #ifdef DOUBLE
    logout << "double precision\n";
    #else
    logout << "simple precision\n";
    #endif
    #ifdef TILT
    logout << "TILT= " << TILT << "\n";
    #endif
    logout 
	<< "C2= " << C2 
	<< ", dt= " << dt
	<< ", L= " << L
	<< ", seed= " << seed;
    logout.flush();

    // Start the timer
    auto start = std::chrono::high_resolution_clock::now();

    for(int i=0;i<Nrun;i++){
        C.update();

        #ifdef MONITORCONFIGS
        if(i%Nmes==0){
            C.print_config(confout);
            C.fourier_transform();
            C.print_inst_sofq(instsofqout);
        }
        #endif
    }

    // Stop the timer
    auto end = std::chrono::high_resolution_clock::now();

    C.print_config(confout);
    C.print_sofq(sofqout);

    // Calculate the duration
    std::chrono::duration<double> duration = end - start;
    // Output the duration
       
    logout << ", Time taken: " << duration.count() << " seconds"
    << ", device= " << deviceProp.name << std::endl;

    return 0;
}

/*
nvcc --expt-extended-lambda -lcufft main.cu -DCu=0.0 -DCphi=0.0 -DEpsilon=0.001 -std=c++14 -arch=sm_61 -o a0.out
*/
