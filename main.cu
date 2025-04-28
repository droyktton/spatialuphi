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
#include <filesystem> // C++17

#ifndef Cu
#define Cu 1.0   // u elastic constant
#endif

#ifndef Cphi
#define Cphi 1.0 // phi elastic constant
#endif

#ifndef Cel
#define Cel 1.0 // phi elastic constant
#endif

#ifndef Alpha
#define Alpha 1.0
#endif

// used to perturb initial flat condition
#ifndef Epsilon
#define Epsilon 0.0
#endif

#ifdef DOUBLE
typedef double real;
typedef cufftDoubleComplex complex;
#else
typedef float real;
typedef cufftComplex complex;
#endif

//#define NOISESPECTRA
//#define RK4

class cuerda{

    public:
    
    void initial_condition(){
        thrust::fill(u.begin(),u.end(),real(0.0));
        thrust::fill(phi.begin(),phi.end(),real(0.0));

        // random initial condition
        for(int i=0;i<L;i++){
            //u[i]=Epsilon*sin(i*30.0*2*M_PI/L) ; 	
	        u[i]=Epsilon*(rand()*1.0/RAND_MAX-0.5);
            //phi[i]=Epsilon*sin(i*20.0*2*M_PI/L); 
	        phi[i]=Epsilon*(rand()*1.0/RAND_MAX-0.5) + 2*M_PI*rand()*1.0/RAND_MAX;
            #ifdef TILT
            u[i]+= i*TILT;      
            //phi[i]+= i*TILT;      
            #endif  
        }

        std::string filename = "initial_configuration.inp";
        if (std::filesystem::exists(filename)) {
            std::cout << "Reading Initial Configuration: " << filename << '\n';
            std::ifstream file(filename);
            double a, b, c, d;    
            int i=0;
            while (file >> a >> b >> c >> d) {
                if(i>=L){
                    std::cerr << "Wrong Initial Configuration" << std::endl;   
                }
                u[i]=a;
                phi[i]=b;
                i++;
            }
        } else {
            std::cerr << "Initial Configuration File not found: " << filename << '\n';
            std::cerr << "Starting flat + noise" << '\n';
        }
    }
    
    
    cuerda(unsigned long _L, real _dt):L(_L),dt(_dt),fourierCount(0)
    {
        u.resize(L);
        phi.resize(L);

        force_u.resize(L);
        force_phi.resize(L);

        //for RK4
        k1_u.resize(L); k1_phi.resize(L);
        k2_u.resize(L); k2_phi.resize(L);
        k3_u.resize(L); k3_phi.resize(L);
        k4_u.resize(L); k4_phi.resize(L);

        initial_condition();

        f0=0.0;
        
        //
        #ifdef DOUBLE
        CUFFT_SAFE_CALL(cufftPlan1d(&plan_r2c,L,CUFFT_D2Z,1));
        #else
        CUFFT_SAFE_CALL(cufftPlan1d(&plan_r2c,L,CUFFT_R2C,1));
        #endif

	    int Lcomp=L/2+1;
	    Fou_u.resize(Lcomp);
	    Fou_phi.resize(Lcomp);

        acum_Sofq_u.resize(L);
        acum_Sofq_phi.resize(L);
        inst_Sofq_u.resize(L);
        inst_Sofq_phi.resize(L);

        thrust::fill(acum_Sofq_u.begin(),acum_Sofq_u.end(),real(0.0));
        thrust::fill(acum_Sofq_phi.begin(),acum_Sofq_phi.end(),real(0.0));

        #ifdef DEBUG
        std::cout << "L=" << L << ", dt=" << dt << std::endl;
        #endif
    }

    void reset_acum_Sofq(){
        thrust::fill(acum_Sofq_u.begin(),acum_Sofq_u.end(),real(0.0));
        thrust::fill(acum_Sofq_phi.begin(),acum_Sofq_phi.end(),real(0.0));
    }

    void set_f0(real x){
        f0=x;
    };

    thrust::tuple<real, real> center_of_mass()
    {
        //DANGER: large sum over large numbers
        real cmu = thrust::reduce(u.begin(),u.end(),real(0.0))/L;
        real cmphi = thrust::reduce(phi.begin(),phi.end(),real(0.0))/L;
        // real maxvelu = thrust::reduce(u.begin(),u.end(),real(0.0),thrust::maximum<real>);
        // real maxvelphi = thrust::reduce(phi.begin(),phi.end(),real(0.0),thrust::maximum<real>);
        // real minvelu = thrust::reduce(u.begin(),u.end(),real(0.0),thrust::minimum<real>);
        // real minvelphi = thrust::reduce(phi.begin(),phi.end(),real(0.0),thrust::minimum<real>);
        return thrust::make_tuple(cmu,cmphi);
    }

    thrust::tuple<real, real> center_of_mass_velocity()
    {
        //SAFE: velocities are bounded
        // #ifdef RK4
        // real vcmu = thrust::reduce(k1_u.begin(),k1_u.end(),real(0.0))/L;
        // real vcmphi = thrust::reduce(k1_phi.begin(),k1_phi.end(),real(0.0))/L;
        // #else
        real vcmu = thrust::reduce(force_u.begin(),force_u.end(),real(0.0))/L;
        real vcmphi = thrust::reduce(force_phi.begin(),force_phi.end(),real(0.0))/L;
        //#endif
        // real maxvelu = thrust::reduce(u.begin(),u.end(),real(0.0),thrust::maximum<real>);
        // real maxvelphi = thrust::reduce(phi.begin(),phi.end(),real(0.0),thrust::maximum<real>);
        // real minvelu = thrust::reduce(u.begin(),u.end(),real(0.0),thrust::minimum<real>);
        // real minvelphi = thrust::reduce(phi.begin(),phi.end(),real(0.0),thrust::minimum<real>);
        return thrust::make_tuple(vcmu,vcmphi);
    }

    thrust::tuple<real, real> particle_velocity(int i){
        return thrust::make_tuple(force_u[i],force_phi[i]);
    }

    void fourier_transform(){

        real *raw_u = thrust::raw_pointer_cast(&u[0]); 
        real *raw_phi = thrust::raw_pointer_cast(&phi[0]); 

        complex *raw_fou_u = thrust::raw_pointer_cast(&Fou_u[0]); 
        complex *raw_fou_phi = thrust::raw_pointer_cast(&Fou_phi[0]); 

        #ifdef DOUBLE
        CUFFT_SAFE_CALL(cufftExecD2Z(plan_r2c, raw_u, raw_fou_u));
	    CUFFT_SAFE_CALL(cufftExecD2Z(plan_r2c, raw_phi, raw_fou_phi));
        #else
	    CUFFT_SAFE_CALL(cufftExecR2C(plan_r2c, raw_u, raw_fou_u));
	    CUFFT_SAFE_CALL(cufftExecR2C(plan_r2c, raw_phi, raw_fou_phi));
        #endif

        thrust::for_each(
            thrust::make_zip_iterator(
                thrust::make_tuple(Fou_u.begin(),Fou_phi.begin(),acum_Sofq_u.begin(),acum_Sofq_phi.begin(),inst_Sofq_u.begin(),inst_Sofq_phi.begin())
            ),
            thrust::make_zip_iterator(
                thrust::make_tuple(Fou_u.end(),Fou_phi.end(),acum_Sofq_u.end(),acum_Sofq_phi.end(),inst_Sofq_u.end(),inst_Sofq_phi.end())
            ),
            [=] __device__ (thrust::tuple<complex,complex,real &,real &, real &,real &> t)
            {
                complex fu=thrust::get<0>(t);
                complex fphi=thrust::get<1>(t);

                thrust::get<2>(t) += fu.x*fu.x + fu.y*fu.y; // S_u(q)
                thrust::get<3>(t) += fphi.x*fphi.x + fphi.y*fphi.y;; // S_phi(q)
                thrust::get<4>(t) = fu.x*fu.x + fu.y*fu.y;
                thrust::get<5>(t) = fphi.x*fphi.x + fphi.y*fphi.y;
            }
        );
        fourierCount++;
    }

    #ifdef NOISESPECTRA    
    void velocity_spectra(std::ofstream &out){

        cufftHandle plan;
        int Trun = vuspectrum.size();

        #ifdef DOUBLE
        //std::cout << vuspectrum.size() << std::endl;
        CUFFT_SAFE_CALL(cufftPlan1d(&plan,Trun,CUFFT_D2Z,1));
        #else
        CUFFT_SAFE_CALL(cufftPlan1d(&plan,Trun,CUFFT_R2C,1));
        #endif

        thrust::device_vector<real> vu(vuspectrum);
        thrust::device_vector<real> vphi(vphispectrum);


        thrust::device_vector<complex> specu(Trun);
        thrust::device_vector<complex> specphi(Trun);

        thrust::device_vector<real> Svu(Trun,real(0));
        thrust::device_vector<real> Svphi(Trun,real(0));

        real *raw_u = thrust::raw_pointer_cast(&vu[0]); 
        real *raw_phi = thrust::raw_pointer_cast(&vphi[0]); 
        complex *raw_fou_u = thrust::raw_pointer_cast(&specu[0]); 
        complex *raw_fou_phi = thrust::raw_pointer_cast(&specphi[0]); 
        
        #ifdef DOUBLE
        CUFFT_SAFE_CALL(cufftExecD2Z(plan, raw_u, raw_fou_u));
	    CUFFT_SAFE_CALL(cufftExecD2Z(plan, raw_phi, raw_fou_phi));
        #else
	    CUFFT_SAFE_CALL(cufftExecR2C(plan, raw_u, raw_fou_u));
	    CUFFT_SAFE_CALL(cufftExecR2C(plan, raw_phi, raw_fou_phi));
        #endif

        thrust::for_each(
            thrust::make_zip_iterator(
                thrust::make_tuple(specu.begin(),specphi.begin(),Svu.begin(),Svphi.begin())
            ),
            thrust::make_zip_iterator(
                thrust::make_tuple(specu.end(),specphi.end(),Svu.end(),Svphi.end())
            ),
            [=] __device__ (thrust::tuple<complex,complex,real &,real &> t){
                complex fu=thrust::get<0>(t);
                complex fphi=thrust::get<1>(t);

                thrust::get<2>(t) += fu.x*fu.x + fu.y*fu.y; // Svu(w)
                thrust::get<3>(t) += fphi.x*fphi.x + fphi.y*fphi.y;; // Svphi(w)
            }
        );

        for(int i=0;i<Svu.size();i++){
            out << Svu[i] << " " << Svphi[i] << "\n";
        }
        out << std::endl;
    }
    #endif


    thrust::tuple<real, real, real, real, real , real, real, real> roughness()
    {
        real cmu = thrust::reduce(u.begin(),u.end(),real(0.f),thrust::plus<real>())/real(L);
        real cmphi = thrust::reduce(phi.begin(),phi.end(),real(0.f),thrust::plus<real>())/real(L);

	    real u0=u[0]; real phi0=phi[0];
        real maxu = thrust::reduce(u.begin(),u.end(),u0,thrust::maximum<real>());
        real maxphi = thrust::reduce(phi.begin(),phi.end(),phi0,thrust::maximum<real>());
        real minu = thrust::reduce(u.begin(),u.end(),u0,thrust::minimum<real>());
        real minphi = thrust::reduce(phi.begin(),phi.end(),phi0,thrust::minimum<real>());

        real cmu2 = 
        thrust::transform_reduce(
            u.begin(),u.end(),
            [=] __device__ (real x){
                return (x-cmu)*(x-cmu);
            },
            real(0.f),
            thrust::plus<real>()
        )/real(L);

        real cmphi2 = 
        thrust::transform_reduce(
            phi.begin(),phi.end(),
            [=] __device__ (real x){
                return (x-cmphi)*(x-cmphi);
            },
            real(0.f),
            thrust::plus<real>()
        )/real(L);

        return thrust::make_tuple(cmu,cmphi,cmu2,cmphi2,maxu,maxphi,minu,minphi);
    }


    void print_center_of_mass(std::ofstream &out)
    {
        thrust::tuple<real,real> cm=center_of_mass();
    
        real cmu=thrust::get<0>(cm);
        real cmphi=thrust::get<1>(cm);

        out << cmu << " " << cmphi << std::endl;
    }

    void rescale()
    {
        thrust::tuple<real,real> cm=center_of_mass();
        real cmu=thrust::get<0>(cm);
        real cmphi=floor(thrust::get<1>(cm)/(2.0*M_PI))*(2.0*M_PI);

        thrust::transform(u.begin(),u.end(),u.begin(),
        [=] __device__ (real u){
            return u- 2*M_PI*int(cmu/(2*M_PI));
        }
        );

        thrust::transform(phi.begin(),phi.end(),phi.begin(),
        [=] __device__ (real phi){
            return phi - 2*M_PI*int(cmphi/(2*M_PI));
        }
        );
    };

    void print_roughness(std::ofstream &out, real t, thrust::tuple<real,real> &acum)
    {
        thrust::tuple<real,real,real,real,real,real,real,real> cm=roughness();
        thrust::tuple<real, real> vcm=center_of_mass_velocity();

        real velu=thrust::get<0>(vcm);
        real velphi=thrust::get<1>(vcm);

        real cmu=thrust::get<0>(cm);
        real cmphi=thrust::get<1>(cm);
        real cmu2=thrust::get<2>(cm);
        real cmphi2=thrust::get<3>(cm);
        real maxu=thrust::get<4>(cm);
        real maxphi=thrust::get<5>(cm);

        thrust::get<0>(acum)+=cmu2;
        thrust::get<1>(acum)+=cmphi2;

        out << t << " " << velu << " " << velphi << " " << cmu << " " 
            << cmphi << " " << cmu2 << " " << cmphi2 << " " << maxu << " " << maxphi << std::endl;
    }

    void update(){

        real *raw_u = thrust::raw_pointer_cast(&u[0]); 
        real *raw_phi = thrust::raw_pointer_cast(&phi[0]); 

        // not elegant solution for lambda...
        real dt_=dt;
        real f0_=f0;
        unsigned long L_ = L;

        // Forces
        thrust::for_each(
            thrust::make_zip_iterator(
                thrust::make_tuple(force_u.begin(),force_phi.begin(),thrust::make_counting_iterator((unsigned long)0))        
            ),
            thrust::make_zip_iterator(
                thrust::make_tuple(force_u.end(),force_phi.end(),thrust::make_counting_iterator((unsigned long)L))        
            ),
            [=] __device__ (thrust::tuple<real &,real &,unsigned long> t)
            {
                unsigned long i=thrust::get<2>(t);
                unsigned long ileft = (i-1+L_)%L_;
                unsigned long iright = (i+1)%L_;

                real uleft = raw_u[ileft];
                real uright = raw_u[iright];
                real phileft = raw_phi[ileft];
                real phiright = raw_phi[iright];
                
                #ifdef TILT
                if(i==0) {
                    uleft -= L_*TILT;
                    //phileft -= L_*TILT;
                }  
                if(i==L_-1){
                    uright += L_*TILT;
                    //phiright += L_*TILT;
                }  
                #endif

                real lap_u = uright + uleft - 2.0*raw_u[i];
                real lap_phi = phiright + phileft - 2.0*raw_phi[i];
		        real sin_2phi = sinf(2.0*raw_phi[i]);

                // forces on u and phi elements
                thrust::get<0>(t) = 0.5*(Alpha*Alpha*f0_+sin_2phi+(Cu*Alpha*lap_u-Cphi*lap_phi));
                thrust::get<1>(t) = 0.5*(Alpha*f0_-Alpha*sin_2phi+(Cu*lap_u+Alpha*Cphi*lap_phi)); 
            } 
        );

        #ifdef DEBUG
        std::cout << "updating" << std::endl;
        #endif

        // Euler step
        thrust::for_each(
            thrust::make_zip_iterator(
                thrust::make_tuple(u.begin(),phi.begin(), force_u.begin(),force_phi.begin())        
            ),
            thrust::make_zip_iterator(
                thrust::make_tuple(u.end(),phi.end(), force_u.end(),force_phi.end())        
            ),
            [=] __device__ (thrust::tuple<real &,real &,real,real> t)
            {
                thrust::get<0>(t) = thrust::get<0>(t) + dt_*thrust::get<2>(t);
                thrust::get<1>(t) = thrust::get<1>(t) + dt_*thrust::get<3>(t);
            } 
        );
    };

    void acum_velocity(thrust::tuple<real, real> vcm){
        #ifdef NOISESPECTRA
        //thrust::tuple<real, real> vcm=center_of_mass_velocity();
        vuspectrum.push_back(thrust::get<0>(vcm));
        vphispectrum.push_back(thrust::get<1>(vcm));
        #endif
    }

    void reset_acum_velocity(){
        #ifdef NOISESPECTRA
        vuspectrum.clear();
        vphispectrum.clear();
        #endif
    }


    void update_runge_kutta4(){

        real *raw_u = thrust::raw_pointer_cast(&u[0]); 
        real *raw_phi = thrust::raw_pointer_cast(&phi[0]); 

        real *raw_k1_u = thrust::raw_pointer_cast(&k1_u[0]); 
        real *raw_k1_phi = thrust::raw_pointer_cast(&k1_phi[0]); 
        real *raw_k2_u = thrust::raw_pointer_cast(&k2_u[0]); 
        real *raw_k2_phi = thrust::raw_pointer_cast(&k2_phi[0]); 
        real *raw_k3_u = thrust::raw_pointer_cast(&k3_u[0]); 
        real *raw_k3_phi = thrust::raw_pointer_cast(&k3_phi[0]); 
        real *raw_k4_u = thrust::raw_pointer_cast(&k4_u[0]); 
        real *raw_k4_phi = thrust::raw_pointer_cast(&k4_phi[0]); 

        // not elegant solution for lambda...
        real dt_=dt;
        real f0_=f0;
        unsigned long L_ = L;

        // Forces k1
        thrust::for_each(
            thrust::make_zip_iterator(
                thrust::make_tuple(k1_u.begin(),k1_phi.begin(),thrust::make_counting_iterator((unsigned long)0))        
            ),
            thrust::make_zip_iterator(
                thrust::make_tuple(k1_u.end(),k1_phi.end(),thrust::make_counting_iterator((unsigned long)L))        
            ),
            [=] __device__ (thrust::tuple<real &,real &,unsigned long> t)
            {
                unsigned long i=thrust::get<2>(t);
                unsigned long ileft = (i-1+L_)%L_;
                unsigned long iright = (i+1)%L_;

                real lap_u=raw_u[iright]+raw_u[ileft]-2.0*raw_u[i];
                real lap_phi=raw_phi[iright]+raw_phi[ileft]-2.0*raw_phi[i];
		        real sin_2phi = sinf(2.0*raw_phi[i]);

                thrust::get<0>(t) = 0.5*(Alpha*Alpha*f0_+sin_2phi+(Cu*Alpha*lap_u-Cphi*lap_phi));
                thrust::get<1>(t) = 0.5*(Alpha*f0_-Alpha*sin_2phi+(Cu*lap_u+Alpha*Cphi*lap_phi)); 
            } 
        );

        // Forces k2
        thrust::for_each(
            thrust::make_zip_iterator(
                thrust::make_tuple(k2_u.begin(),k2_phi.begin(),thrust::make_counting_iterator((unsigned long)0))        
            ),
            thrust::make_zip_iterator(
                thrust::make_tuple(k2_u.end(),k2_phi.end(),thrust::make_counting_iterator((unsigned long)L))        
            ),
            [=] __device__ (thrust::tuple<real &,real &,unsigned long> t)
            {
                unsigned long i=thrust::get<2>(t);
                unsigned long ileft = (i-1+L_)%L_;
                unsigned long iright = (i+1)%L_;

                real raw_u_i=raw_u[i]+dt_*raw_k1_u[i]*0.5;
                real raw_u_iright=raw_u[iright]+dt_*raw_k1_u[iright]*0.5;
                real raw_u_ileft=raw_u[ileft]+dt_*raw_k1_u[ileft]*0.5;

                real raw_phi_i=raw_phi[i]+dt_*raw_k1_phi[i]*0.5;
                real raw_phi_iright=raw_phi[iright]+dt_*raw_k1_phi[iright]*0.5;
                real raw_phi_ileft=raw_phi[ileft]+dt_*raw_k1_phi[ileft]*0.5;

                real lap_u=(raw_u_iright)+(raw_u_ileft)-2.0*(raw_u_i);
                real lap_phi=(raw_phi_iright)+(raw_phi_ileft)-2.0*(raw_phi_i);
		        real sin_2phi = sinf(2.0*(raw_phi_i));

                thrust::get<0>(t) = 0.5*(Alpha*Alpha*f0_+sin_2phi+(Cu*Alpha*lap_u-Cphi*lap_phi));
                thrust::get<1>(t) = 0.5*(Alpha*f0_-Alpha*sin_2phi+(Cu*lap_u+Alpha*Cphi*lap_phi)); 
            } 
        );

        // Forces k3
        thrust::for_each(
            thrust::make_zip_iterator(
                thrust::make_tuple(k3_u.begin(),k3_phi.begin(),thrust::make_counting_iterator((unsigned long)0))        
            ),
            thrust::make_zip_iterator(
                thrust::make_tuple(k3_u.end(),k3_phi.end(),thrust::make_counting_iterator((unsigned long)L))        
            ),
            [=] __device__ (thrust::tuple<real &,real &,unsigned long> t)
            {
                unsigned long i=thrust::get<2>(t);
                unsigned long ileft = (i-1+L_)%L_;
                unsigned long iright = (i+1)%L_;

                real raw_u_i=raw_u[i]+dt_*raw_k2_u[i]*0.5;
                real raw_u_iright=raw_u[iright]+dt_*raw_k2_u[iright]*0.5;
                real raw_u_ileft=raw_u[ileft]+dt_*raw_k2_u[ileft]*0.5;

                real raw_phi_i=raw_phi[i]+dt_*raw_k2_phi[i]*0.5;
                real raw_phi_iright=raw_phi[iright]+dt_*raw_k2_phi[iright]*0.5;
                real raw_phi_ileft=raw_phi[ileft]+dt_*raw_k2_phi[ileft]*0.5;

                real lap_u=(raw_u_iright)+(raw_u_ileft)-2.0*(raw_u_i);
                real lap_phi=(raw_phi_iright)+(raw_phi_ileft)-2.0*(raw_phi_i);
		        real sin_2phi = sinf(2.0*(raw_phi_i));

                thrust::get<0>(t) = 0.5*(Alpha*Alpha*f0_+sin_2phi+(Cu*Alpha*lap_u-Cphi*lap_phi));
                thrust::get<1>(t) = 0.5*(Alpha*f0_-Alpha*sin_2phi+(Cu*lap_u+Alpha*Cphi*lap_phi)); 
            } 
        );

        // Forces k4
        thrust::for_each(
            thrust::make_zip_iterator(
                thrust::make_tuple(k4_u.begin(),k4_phi.begin(),thrust::make_counting_iterator((unsigned long)0))        
            ),
            thrust::make_zip_iterator(
                thrust::make_tuple(k4_u.end(),k4_phi.end(),thrust::make_counting_iterator((unsigned long)L))        
            ),
            [=] __device__ (thrust::tuple<real &,real &,unsigned long> t)
            {
                unsigned long i=thrust::get<2>(t);
                unsigned long ileft = (i-1+L_)%L_;
                unsigned long iright = (i+1)%L_;

                real raw_u_i=raw_u[i]+dt_*raw_k3_u[i]*0.5;
                real raw_u_iright=raw_u[iright]+dt_*raw_k3_u[iright]*0.5;
                real raw_u_ileft=raw_u[ileft]+dt_*raw_k3_u[ileft]*0.5;

                real raw_phi_i=raw_phi[i]+dt_*raw_k3_phi[i]*0.5;
                real raw_phi_iright=raw_phi[iright]+dt_*raw_k3_phi[iright]*0.5;
                real raw_phi_ileft=raw_phi[ileft]+dt_*raw_k3_phi[ileft]*0.5;

                real lap_u=(raw_u_iright)+(raw_u_ileft)-2.0*(raw_u_i);
                real lap_phi=(raw_phi_iright)+(raw_phi_ileft)-2.0*(raw_phi_i);
		        real sin_2phi = sinf(2.0*(raw_phi_i));

                thrust::get<0>(t) = 0.5*(Alpha*Alpha*f0_+sin_2phi+(Cu*Alpha*lap_u-Cphi*lap_phi));
                thrust::get<1>(t) = 0.5*(Alpha*f0_-Alpha*sin_2phi+(Cu*lap_u+Alpha*Cphi*lap_phi)); 
            } 
        );

        // Euler step
        thrust::for_each(
            thrust::make_zip_iterator(
                thrust::make_tuple(u.begin(),phi.begin(),thrust::make_counting_iterator((unsigned long)0), force_u.begin(),force_phi.begin())        
            ),
            thrust::make_zip_iterator(
                thrust::make_tuple(u.end(),phi.end(),thrust::make_counting_iterator((unsigned long)L), force_u.end(),force_phi.end())        
            ),
            [=] __device__ (thrust::tuple<real &,real &,unsigned long,real &,real &> t)
            {
                unsigned long i=thrust::get<2>(t);
                thrust::get<0>(t) = thrust::get<0>(t) + dt_*(raw_k1_u[i]+2*raw_k2_u[i]+2*raw_k3_u[i]+raw_k4_u[i])/6.0;
                thrust::get<1>(t) = thrust::get<1>(t) + dt_*(raw_k1_phi[i]+2*raw_k2_phi[i]+2*raw_k3_phi[i]+raw_k4_phi[i])/6.0;
                thrust::get<3>(t) = (raw_k1_u[i]+2*raw_k2_u[i]+2*raw_k3_u[i]+raw_k4_u[i])/6.0;
                thrust::get<4>(t) = (raw_k1_phi[i]+2*raw_k2_phi[i]+2*raw_k3_phi[i]+raw_k4_phi[i])/6.0;
            } 
        );        
    }


    void update_leap_frog(){

        real *raw_u = thrust::raw_pointer_cast(&u[0]); 
        real *raw_phi = thrust::raw_pointer_cast(&phi[0]); 

        // not elegant solution for lambda...
        real dt_=dt;
        real f0_=f0;
        unsigned long L_ = L;

        // Forces u
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

                real lap_u=raw_u[iright]+raw_u[ileft]-2.0*raw_u[i];
                real lap_phi=raw_phi[iright]+raw_phi[ileft]-2.0*raw_phi[i];
		        real sin_2phi = sinf(2.0*raw_phi[i]);

                thrust::get<0>(t) = 0.5*(Alpha*Alpha*f0_+sin_2phi+(Cu*Alpha*lap_u-Cphi*lap_phi));            
            } 
        );
        // Euler step u
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
                //thrust::get<1>(t) = thrust::get<1>(t) + dt_*thrust::get<3>(t);
            } 
        );
        // Forces phi
        thrust::for_each(
            thrust::make_zip_iterator(
                thrust::make_tuple(force_phi.begin(),thrust::make_counting_iterator((unsigned long)0))        
            ),
            thrust::make_zip_iterator(
                thrust::make_tuple(force_phi.end(),thrust::make_counting_iterator((unsigned long)L))        
            ),
            [=] __device__ (thrust::tuple<real &,unsigned long> t)
            {
                unsigned long i=thrust::get<1>(t);
                unsigned long ileft = (i-1+L_)%L_;
                unsigned long iright = (i+1)%L_;

                real lap_u=raw_u[iright]+raw_u[ileft]-2.0*raw_u[i];
                real lap_phi=raw_phi[iright]+raw_phi[ileft]-2.0*raw_phi[i];
		        real sin_2phi = sinf(2.0*raw_phi[i]);

                thrust::get<0>(t) = 0.5*(Alpha*f0_-Alpha*sin_2phi+(Cu*lap_u+Alpha*Cphi*lap_phi)); 
            } 
        );
        // Euler step phi
        thrust::for_each(
            thrust::make_zip_iterator(
                thrust::make_tuple(phi.begin(), force_phi.begin())        
            ),
            thrust::make_zip_iterator(
                thrust::make_tuple(phi.end(),force_phi.end())        
            ),
            [=] __device__ (thrust::tuple<real &,real> t)
            {
                thrust::get<0>(t) = thrust::get<0>(t) + dt_*thrust::get<1>(t);
                //thrust::get<1>(t) = thrust::get<1>(t) + dt_*thrust::get<3>(t);
            } 
        );


    };

    void print_config(std::ofstream &out){
        thrust::tuple<real,real> cm = center_of_mass();

        for(int i=0;i<L;i++){
            out << u[i] << " " << phi[i] << " " << thrust::get<0>(cm) << " " << thrust::get<1>(cm) << "\n";
        }
        out << "\n" << std::endl;
    };

    void print_sofq(std::ofstream &out){
        for(int i=0;i<L;i++){
            out << acum_Sofq_u[i]/fourierCount << " " << acum_Sofq_phi[i]/fourierCount << "\n";
        }
        out << "\n" << std::endl;
    };

    void print_inst_sofq(std::ofstream &out){
        for(int i=0;i<L;i++){
            out << inst_Sofq_u[i] << " " << inst_Sofq_phi[i] << "\n";
        }
        out << "\n" << std::endl;
    };

    private:
        real dt;
        unsigned long L;
        //unsigned int seed;
        real f0;
        thrust::device_vector<real> u;
        thrust::device_vector<real> phi;
        thrust::device_vector<real> force_u;
        thrust::device_vector<real> force_phi;
        thrust::device_vector<real> k1_u, k1_phi, k2_u, k2_phi, k3_u, k3_phi, k4_u, k4_phi;

        int fourierCount;
        cufftHandle plan_r2c;

        thrust::device_vector<complex> Fou_u;
	    thrust::device_vector<complex> Fou_phi;

        thrust::device_vector<real> acum_Sofq_u;
	    thrust::device_vector<real> acum_Sofq_phi;
        thrust::device_vector<real> inst_Sofq_u;
	    thrust::device_vector<real> inst_Sofq_phi;

        #ifdef NOISESPECTRA
        std::vector<real> vuspectrum;
        std::vector<real> vphispectrum;
        #endif
};

int main(int argc, char **argv){

    if(argc<6){
        std::cout << "Error: Not enough arguments" << std::endl;
        std::cout << "Arguments: L f0 f1 df Nrun [seed]" << std::endl;
        std::cout << "L: length of the system" << std::endl;
        std::cout << "f0: starting force" << std::endl;
        std::cout << "f1: stoping force" << std::endl;
        std::cout << "df: force step" << std::endl;
        std::cout << "Usage: " << argv[0] << " L f0 f1 df Nrun [seed]" << std::endl;
        return 1;
    }
    std::ofstream logout("log.txt");


    unsigned int L=atoi(argv[1]);

    real f0=atof(argv[2]);
    real f1=atof(argv[3]);
    real df=atof(argv[4]);

    unsigned long Nrun = atoi(argv[5]);
    unsigned long Neq = unsigned(Nrun/2.0);

    real dt=0.1;

    int Nmes=1000000;
    int Nrescale=100000000;

    unsigned int seed=1234;
    if(argc==7) seed=(unsigned int)atoi(argv[6]); 
    srand(seed);

    
    cuerda C(L,dt);


    // Get the current CUDA device
    int device;
    cudaGetDevice(&device);

    // Get the properties of the current CUDA device
    cudaDeviceProp deviceProp;
    cudaGetDeviceProperties(&deviceProp, device);

    std::ofstream confout("conf.dat");
    confout << "#u[i]" << " " << "phi[i]" << " " << "cmu" << " " << "cmphi" << "\n";

    std::ofstream sofqout("sofq.dat");
    sofqout << "#av_Sofq_u[i]" << " " << "av_Sofq_phi[i]" << "\n";

    std::ofstream instsofqout("inst_sofq.dat");
    instsofqout << "#inst_Sofq_u[i]" << " " << "inst_Sofq_phi[i]" << "\n";

    std::ofstream cmout("cm.dat");
    cmout << "#t" << " " << "velu" << " " << "velphi" << " " << "cmu" << " " 
    << "cmphi" << " " << "cmu2" << " " << "cmphi2" << " " << "maxu" << " " << "maxphi" << std::endl;

    std::ofstream logcmout("logcm.dat");
    logcmout << "#t" << " " << "velu" << " " << "velphi" << " " << "cmu" << " " 
    << "cmphi" << " " << "cmu2" << " " << "cmphi2" << " " << "maxu" << " " << "maxphi" << std::endl;


    std::ofstream lastconfout("final_configuration.inp");
    //lastconfout << "#u[i]" << " " << "phi[i]" << " " << "cmu" << " " << "cmphi" << "\n";

    std::ofstream vhout("VvsH.dat");
    vhout << "#f0" << " " << "vcmu" << " " << "vcmphi" << " "
    << "(cmu1-cmu0)/time" << " " 
    << "(cphi1-cphi0)/time" << " " 
    << "roughu1" << " " << "roughphi1" << " "
    << "roughu0" << " " << "roughphi0" << " "
    << "maxu1" << " " << "maxphi1" << " "
    << "minu0" << " " << "minphi0" << " "
    << std::endl;

    #ifdef NOISESPECTRA
    std::ofstream sofwout("Sofw.dat");
    #endif

    #ifdef DOUBLE
    logout << "DOUBLE precision\n";
    #else
    logout << "SIMPLE precision\n";
    #endif
    #ifdef RK4
    logout << "RK4\n";
    #else 
    logout << "EULER\n";
    #endif
    #ifdef NOISESPECTRA
    logout << "NOISESPECTRA\n";
    #endif
    #ifdef TILT
    logout << "TILT= " << TILT << "\n";
    #endif
    logout 
	<< "Cu= " << Cu << "\n"
	<< "Cphi= " << Cphi << "\n"
	<< "Epsilon= " << Epsilon << "\n"
	<< "dt= " << dt << "\n"
	<< "L= " << L << "\n"
    << "Nmes= " << Nmes << "\n"
    << "Nrescale= " << Nrescale << "\n"
    << "Nrun= " << Nrun << "\n"
    << "Neq= " << Neq << "\n"
    << "df= " << df << "\n"
    << "f0= " << f0 << "\n"
    << "f1= " << f1 << "\n"
	<< "seed= " << seed << "\n ======= \n";
    logout.flush();


    //real q0,phi0,q1,phi1;
    thrust::tuple<real,real,real,real,real,real,real,real> r0 = C.roughness();
    thrust::tuple<real,real> acumroughness = thrust::make_tuple(real(0.), real(0.));

    thrust::tuple<real,real> vcm;

    real vcmu,vcmphi;

    // Start the timer
    auto start = std::chrono::high_resolution_clock::now();
    logout << "Starting simulation..." << std::endl;
    std::cout << "Starting simulation..." << std::endl;

    // histeresis loop
    for (real f = f0, dir = df; (dir > 0 && f <= f1) || (dir < 0 && f >= f0); f += dir)
    {
        logout << "f= " << f << std::endl;
        logout.flush();

        std::cout << "f= " << f << std::endl;
        std::cout.flush();

        C.set_f0(f);
        vcmu=vcmphi=0.0;
        C.reset_acum_Sofq();
        C.reset_acum_velocity(); 
        unsigned long ime=0;


        unsigned long jlog=1;


        for(int i=0;i<Nrun+1;i++)
        {
            #ifdef RK4
            C.update_runge_kutta4();
            #else
            C.update();
            #endif

            if(i>Neq){
                    vcm=C.center_of_mass_velocity();
                    vcmu+=thrust::get<0>(vcm);
                    vcmphi+=thrust::get<1>(vcm);
                    C.acum_velocity(vcm); // for velocity spectrum
                    ime++;
            }

            if(i%Nrescale==0) {
                C.rescale();
                logout << "rescale at " << i << std::endl;
            }    

            #ifdef MONITORCONFIGS
            if(i%MONITORCONFIGS==0){
                C.print_config(confout);
                C.fourier_transform();
                C.print_inst_sofq(instsofqout);
            }
            #endif
            // monitor every Nmes steps
            if(i%Nmes==0){
                //std::cout << "i= " << i << std::endl;
                C.print_roughness(cmout,i*dt,acumroughness);
                if(i>Neq){
                    C.fourier_transform();
                }
            } 
            // if steady state...
            if(i==Neq) 
            {
                r0 = C.roughness();
                C.reset_acum_Sofq();
            }
            
            // log monitoring
            if(i%jlog==0){
                C.fourier_transform();
                C.print_inst_sofq(instsofqout);
                C.print_roughness(logcmout,i*dt,acumroughness);
                jlog=jlog*10;
            }
        }
               
        //thrust::tuple<real,real> cm1 = C.center_of_mass();
        thrust::tuple<real,real,real,real,real,real,real,real> r1 = C.roughness();
        //q1=C.u[0];phi1=C.phi[0];
        
        #ifdef NOISESPECTRA
        C.velocity_spectra(sofwout);
        #endif

        // last config and last structure factor
        C.print_config(confout);
        C.print_sofq(sofqout);

        // Averages for a single field
        vhout << std::setprecision(4) << f << " " << vcmu/real(ime) << " " << vcmphi/real(ime) << " "
        << (thrust::get<0>(r1)-thrust::get<0>(r0))/(Nrun*dt*0.5) << " " 
        << (thrust::get<1>(r1)-thrust::get<1>(r0))/(Nrun*dt*0.5) << " " 
        << thrust::get<2>(r1) << " " << thrust::get<3>(r1) << " "
        << thrust::get<2>(r0) << " " << thrust::get<3>(r0) << " "
        << thrust::get<4>(r0) << " " << thrust::get<5>(r0) << " "
        << thrust::get<4>(r1) << " " << thrust::get<5>(r1) << " "
        << thrust::get<0>(acumroughness)/real(ime) << " " << thrust::get<1>(acumroughness)/real(ime) << " "
        << dir << std::endl;

        // change direction
        if (f + dir > f1 && dir > 0) dir = -df;
    }

    // prints last configuration
    C.print_config(lastconfout);

    // Stop the timer
    auto end = std::chrono::high_resolution_clock::now();

    // Calculate the duration
    std::chrono::duration<double> duration = end - start;

    // Output the duration      
    logout << "Time taken: " << duration.count() << " seconds"
    << ", device= " << deviceProp.name << std::endl;

    return 0;
}

/*
nvcc --expt-extended-lambda -lcufft main.cu -DCu=0.0 -DCphi=0.0 -DEpsilon=0.001 -std=c++14 -arch=sm_61 -o a0.out
*/
