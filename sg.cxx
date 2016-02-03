#include <iostream>
#include <fstream>
#include <complex>
#include <sstream>
//-----------------------------------
using namespace std;
//-----------------------------------
typedef complex<double> cmplx;
//-----------------------------------
void step(cmplx* const psi0, cmplx* const psi1, const double dt, const double dx, const double xmin, const int Nx, const double k);

void init( cmplx* const psi0, const double alpha, const double lambda,
           const double dx, const double dt,
           const int Nx, const double xmin);


void writeToFile(const cmplx* const v, const string s, const double dx,
         const int Nx, const double xmin, const double alpha,
         const double lambda, const double omega, const double t);

void step();


//-----------------------------------
int main(){

	const int Nx = 300;
	const double xmin = -40.0;
  const double xmax = 40.0;
	const double Tend = 10 * M_PI;
	const double dx = (xmax - xmin)/(Nx -1);
	const double dt = 0.1 * dx;
  double t = 0;
	const int Na = 10;
	int Nk = int(Tend / Na / dt + 0.5);

	const double lambda = 10;
  const double omega = 0.2;
  const double k = omega * omega;
  const double alpha = pow(k, 0.25);
  stringstream strm;

	cmplx* psi0 = new cmplx[Nx];
	cmplx* psi1 = new cmplx[Nx];
	cmplx* h;

	init(psi0, alpha, lambda, dx, dt, Nx, xmin);

	writeToFile(psi0,"psi_0", dx,Nx,xmin, alpha, lambda, omega,t);



	for (int i = 1; i <= Na; i++) {
		for (int j = 1; j <= Nk-1; j++) {
			 t+=dt;
			 
			 step(psi0, psi1, dt, dx, xmin, Nx, k);	
			
			 h = psi0;
			 psi0 = psi1;
			 psi1 = h;
			
		}
		strm.str("");
		strm << "psi_" << i;
		writeToFile(psi0,strm.str(), dx,Nx,xmin, alpha, lambda, omega,t);
	}
  cout << "t = " << t << endl;

  
  	delete[] psi0;
  	delete[] psi1;
  
	return 0;
}
//-----------------------------------
void step(cmplx* const psi0, cmplx* const psi1, const double dt, const double dx, const double xmin, const int Nx, const double k) {
	
	// help vars
	const cmplx alpha = cmplx(0, -dt/(4 * dx * dx));
	cmplx* d = new cmplx[Nx];
	cmplx* psih = new cmplx[Nx];
	double x= xmin;
	
	// implement substep 1 - fill the matrix
	for (int i = 0; i < Nx; i++) {
		x += i * dx;
		d[i] = cmplx(1.0, dt/(2 * dx * dx) + (dt * k * x * x)/4);
	}
	
	// calculate the right side of the equitation system
	psih[0] = conj(d[0]) * psi0[0] + conj(alpha) * psi0[1];
	for (int j = 1; j < Nx -1; j++) {
		psih[j] = conj(alpha) * psi0[j-1] + conj(d[j]) * psi0[j] + conj(alpha) * psi0[j+1];
	}
	psih[Nx-1] = conj(alpha) * psi0[Nx-2] + conj(d[Nx-1]) * psi0[Nx-1];
	
	// forward substitution like lab 11
	for (int i = 1; i < Nx; i++) {
		d[i] -= (alpha * alpha)/d[i-1];
		psih[i] -= psih[i-1] * alpha / d[i-1];
	}
	
	// backward substitution
	psi1[Nx-1] = psih[Nx-1] / d[Nx-1];
	for (int i = Nx-2; i >= 0; i--) {
		psi1[i] = (psih[i] - alpha *psi1[i+1])/d[i];
	}
	
	delete[] d;
	delete[] psih;
}

//-----------------------------------
void writeToFile(const cmplx* const v, const string s, const double dx,
                 const int Nx, const double xmin, const double alpha,
                 const double lambda, const double omega, const double t)
{
	ofstream out(s.c_str());
  double x, xi, xil;
  double h1, h2, h3;
  cmplx ana;
	for(int i=0; i<Nx; i++){
		x = xmin + i * dx;
    xi = alpha * x;
    xil = alpha * lambda;
    h1 = -0.5 * pow(xi - xil*cos(omega*t),2 );
    h2 = omega*t/2 + xi * xil * sin(omega*t);
    h3 =  - 0.25 * xil*xil* sin(2*omega*t);
    ana = cmplx( h1 , h2 + h3  );
    ana = sqrt(alpha / sqrt(M_PI) ) * exp(ana);
		out << x << "\t" << norm(v[i]) << "\t" << v[i].real() << "\t" << v[i].imag()
         << "\t" << norm(ana) << "\t" << ana.real() << "\t" << ana.imag() <<  endl;
	}
	out.close();
}
//-----------------------------------

void init( cmplx* const psi0, const double alpha, const double lambda,
           const double dx, const double dt,
           const int Nx, const double xmin)
{
	const double x0 = dx*Nx * 0.5;
	for(int i=0;i<Nx; i++){
		double x = xmin + i*dx ;
		psi0[i] = sqrt(alpha/sqrt(M_PI)) * exp(- pow(alpha*(x-lambda),2)/2 ); // Initial function -> given.
	}
}
