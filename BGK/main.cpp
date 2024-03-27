///----------------------------------------------------------------------------------------------------------------------------------
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <cassert>
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <ctime>
///----------------------------------------------------------------------------------------------------------------------------------
using namespace std;
const bool plot_vtk = true;
/// Flow quantities
const bool cms_hydro = true, cms_phase = true;
const int nx = 200, ny = 4*nx, np = 9;
const double cs2 = 1./3., cs4 = cs2*cs2, cs6 = cs4*cs2, cs8 = cs6*cs2;
const double  xi = 4., d = (double)(nx-1), U_ref = 0.04, gravity = pow(U_ref,2)/d, Reynolds = 256., nu = sqrt(d*gravity)*d/Reynolds, At = 0.5, T = sqrt(d/gravity/At), Pe = 50;
const int nsteps = (int)(3*T+1), n_out = (int)(T/2);
const double rhoL = 1., rhoH = rhoL*(1.+At)/(1.-At), niH = nu, tauH = niH/cs2, niL = nu, tauL = niL/cs2, sigma = 1E-5;
// HYDRO            0   1   2   3   4   5  6  7  8
const vector<int> cx = {0, 1, 0, -1, 0, 1, -1, -1, 1},
									cy = {0, 0, 1, 0, -1, 1, 1, -1, -1},
                  opp = {0, 3, 4, 1, 2, 7, 8, 5, 6};
const vector<double> wf = {4./9., 1./9., 1./9., 1./9., 1./9., 1./36., 1./36., 1./36., 1./36.};
double A, B, C, R, U, V, ftemp, CX, CY, U2, V2, U2V2, UV, U2V, UV2, U0, V0, Press0;
double Press, tau, omega, omega1, feq, ni;
vector<double> temp_pop(np, 0.), f1(nx*ny*np, 0.), f2(nx*ny*np, 0.), press(nx*ny, 0.), u(nx*ny, 0.), v(nx*ny, 0.), rho(nx*ny, 0.), u_old(nx*ny, 0.), v_old(nx*ny, 0.);
double Fpx, Fpy, Fsx, Fsy, Fm, Fmx, Fmy, mu;
double gradx_u, grady_u, gradx_v, grady_v;
int newx, newy, check;
double k0, k1, k2, k3, k4, k5, k6, k7, k8;
double r0, r1, r2, r3, r4, r5, r6, r7, r8;
double third_order, fourth_order;
// PHASE
const double PhiH = 1., PhiL = 0., Phi0 = 0.5*(PhiH+PhiL);
const double beta = 12.*sigma/xi, kappa = 3.*sigma*xi/2.;
const double M = 0.1, tau_phase = M/cs2+0.5, omega_phase = 1./tau_phase, omega_phase1 = 1.-omega_phase;
double grad_phix, grad_phiy, laplPhi;
double Phi, Phi_prev, geq;
vector<double> g1(nx*ny*np, 0.), g2(nx*ny*np, 0.), temp_pop_phase(np, 0.), phase(nx*ny, 0.), phase_old(nx*ny, 0.);
double Fx, Fy, gtemp, P, Fx_phase, Fy_phase;
double k1_g, k2_g, k3_g, k4_g, k5_g, k6_g, k7_g, k8_g, extra_term;
double Gamma, Nx, Ny;
double Phii, thisPhi;

int newxx, newyy, id, idn;
double energy, enstrophy, gradx_of_v, grady_of_u, gradx_of_u, grady_of_v, dissipation_rate;
const double physical_ref_velocity = sqrt(1./9.81), scale_velocity = physical_ref_velocity/U_ref,
scale_length = 1./d, scale_time = scale_length/scale_velocity, scale_rho = 1.,
scale_enstrophy = 1./pow(scale_time,2), scale_energy = scale_rho*pow(scale_velocity,2),
scale_dissipation = scale_enstrophy*scale_length*scale_length/scale_time;
///----------------------------------------------------------------------------------------------------------------------------------
///----------------------------------------------------------------------------------------------------------------------------------
void diagnostics(int l)
{
    if(l%1==0)
    {
    	printf("%lf of %lf.\n", l/T, nsteps/T);
    }
}
///----------------------------------------------------------------------------------------------------------------------------------
void write_fluid_vtk(int time)
{
	/// Create filename
	stringstream output_filename;
	output_filename << "vtk_fluid/fluid_t" << time << ".vtk";
	ofstream output_file;

	/// Open file
	output_file.open(output_filename.str().c_str());

	/// Write VTK header
	output_file << "# vtk DataFile Version 3.0\n";
	output_file << "fluid_state\n";
	output_file << "ASCII\n";
	output_file << "DATASET RECTILINEAR_GRID\n";
	output_file << "DIMENSIONS " << nx << " " << ny << " 1" << "\n";
	output_file << "X_COORDINATES " << nx << " double\n";
	for(int i = 0; i < nx; ++i)
		output_file << i << " ";
	output_file << "\n";
	output_file << "Y_COORDINATES " << ny  << " double\n";
	for(int j = 0; j < ny ; ++j)
		output_file << j  << " ";
	output_file << "\n";
	output_file << "Z_COORDINATES " << 1 << " double\n";
	output_file << 0 << "\n";
	output_file << "POINT_DATA " << (nx) * (ny) << "\n";

	output_file << "SCALARS phase double 1\n";
	output_file << "LOOKUP_TABLE default\n";
	for(int Y = 0; Y < ny ; ++Y)
		for(int X = 0; X < nx; ++X)
		{
			if(phase[X*ny+Y]>1E-10)
				output_file << phase[X*ny+Y]<< "\n";
			else
				output_file << 0 << "\n";
			}

  // output_file << "SCALARS density double 1\n";
  // output_file << "LOOKUP_TABLE default\n";
  // for(int Y = 0; Y < ny ; ++Y)
  //   for(int X = 0; X < nx; ++X)
  //   	output_file << rho[X*ny+Y]<< "\n";

  // Write velocity
  //output_file << "VECTORS velocity_vector double\n";
  //for(int Y = 0; Y < ny ; ++Y)
  	//for(int X = 0; X < nx; ++X)
  		//output_file << u[X*ny+Y] << " " << v[X*ny+Y] << " 0\n";

	/// Close file
	output_file.close();
}
///----------------------------------------------------------------------------------------------------------------------------------
void compute_grad_phi(int x, int y)
{
	grad_phix = grad_phiy = laplPhi = 0.;
	gradx_of_u = gradx_of_v = grady_of_u = grady_of_v = 0.;
  thisPhi = phase_old[x*ny+y];
	for(int k=1; k<np; k++)
	{
		newx = x+cx[k];
		newy = y+cy[k];
		if(x==0 || x==nx-1)
			newx = (newx+nx)%nx;
		if(y==0 || y==ny-1)
			newy = (newy+ny)%ny;
		Phii = phase_old[newx*ny+newy];
		grad_phix += Phii*wf[k]*cx[k];
		grad_phiy += Phii*wf[k]*cy[k];
    laplPhi += wf[k]*(Phii-thisPhi);

		gradx_of_u += u[newx*ny+newy]*wf[k]*cx[k];
		gradx_of_v += v[newx*ny+newy]*wf[k]*cx[k];
		grady_of_u += u[newx*ny+newy]*wf[k]*cy[k];
		grady_of_v += v[newx*ny+newy]*wf[k]*cy[k];
	}
	grad_phix /= cs2;
	grad_phiy /= cs2;
	gradx_of_u /= cs2;
	gradx_of_v /= cs2;
	grady_of_u /= cs2;
	grady_of_v /= cs2;
  laplPhi *= 2./cs2;
}
///----------------------------------------------------------------------------------------------------------------------------------
void initial_state()
{
  double extra_term, X;
	int max = 1, min = -1;
	double h, an, bn;
	for(int x=0; x<nx; x++)
    for(int y=0; y<ny; y++)
		{
			id = x*ny+y;
      /*X = (double)x / (5.*(double)nx-1);
			h = 0;
			for(int n=30; n<40; n++)
			{
				an = rand()%(max-min+1)+min;
				bn = rand()%(max-min+1)+min;
				h += an*cos(2.*M_PI*n*X)+bn*sin(2.*M_PI*n*X);
			}
			h = 0.5*ny+nx*0.002*h;
			phase[id] = 0.5+0.5*tanh(2.*(y-h)/xi);*/
			X = double(x) / double(nx-1);
			if(y>2.*nx+0.1*nx*(cos(2.*M_PI*X)))
				phase[id] = PhiH;
			else
				phase[id] = PhiL;
      rho[id] = rhoL+phase[id]*(rhoH-rhoL);
		}
	for(int x=0; x<nx; x++)
		for(int y=0; y<ny; y++)
		{
			id = x*ny+y;
			u[id] = U = 0.;
			v[id] = V = 0.;
      R = rho[id];
      Phi = phase[id];
		  compute_grad_phi(x,y);
      Nx  = grad_phix/(sqrt(pow(grad_phix,2)+pow(grad_phiy,2))+1E-12);
      Ny  = grad_phiy/(sqrt(pow(grad_phix,2)+pow(grad_phiy,2))+1E-12);

      Fpx = -Press*cs2*(rhoH-rhoL)/(PhiH-PhiL)*grad_phix;
      Fpy = -Press*cs2*(rhoH-rhoL)/(PhiH-PhiL)*grad_phiy;
      mu = 4.*beta*(Phi-PhiL)*(Phi-PhiH)*(Phi-Phi0)-kappa*laplPhi;
      Fsx = mu*grad_phix;
      Fsy = mu*grad_phiy;
      Fx = Fpx+Fsx+Fmx;
      Fy = Fpy+Fsy+Fmy;
      for(int k=0; k<np;k++)
			{
        A = U*cx[k]+V*cy[k];
        third_order = 1./(2.*cs6)*((cx[k]*cx[k]-cs2)*cy[k]*U*U*V+(cy[k]*cy[k]-cs2)*cx[k]*U*V*V);
        fourth_order = 1./(4.*cs8)*((cx[k]*cx[k]-cs2)*(cy[k]*cy[k]-cs2)*U*U*V*V);
        extra_term = 4.*Phi*(1.-Phi)/xi*wf[k]*(cx[k]*Nx+cy[k]*Ny);
        geq = Phi*wf[k]*(1.+3.*A+4.5*((cx[k]*cx[k]-cs2)*U*U+(cy[k]*cy[k]-cs2)*V*V+
                2.*cx[k]*cy[k]*U*V)+third_order+fourth_order) -
              0.5*extra_term;
        g1[id*np+k] = g2[id*np+k] = geq+extra_term;

        feq = wf[k]*(Press+3.*A+4.5*((cx[k]*cx[k]-cs2)*U*U+(cy[k]*cy[k]-cs2)*V*V+
                   2.*cx[k]*cy[k]*U*V)+third_order+fourth_order)-
                   0.5*wf[k]/(R*cs2)*(cx[k]*Fx+cy[k]*Fy);
        f1[id*np+k] = f2[id*np+k] = feq+wf[k]/(R*cs2)*(cx[k]*Fx+cy[k]*Fy);
			}
		}
}
///----------------------------------------------------------------------------------------------------------------------------------
int algorithm_lattice_boltzmann()
{
  check = 0;
	energy = enstrophy = 0.;
	u_old = u;
	v_old = v;
	phase_old = phase;
	for(int x=0; x<nx; x++)
		for(int y=0; y<ny; y++)
		{
			id = x*ny+y;
			compute_grad_phi(x,y);
			U0 = u[id];
			V0 = v[id];
			energy += 0.5*rho[id]*(U0*U0+V0*V0);
			Press0 = press[id];
			U = V = Press = Phi = 0.;
			for(int k=0; k<np; k++)
      {
        temp_pop[k] = ftemp = f1[id*np+k];
        temp_pop_phase[k] = g1[id*np+k];
        Press += ftemp;
        U += ftemp*cx[k];
        V += ftemp*cy[k];
				Phi += g1[id*np+k];
      }
			phase[id] = Phi;
      press[id] = Press;
      Fpx = -Press*cs2*(rhoH-rhoL)*grad_phix;
      Fpy = -Press*cs2*(rhoH-rhoL)*grad_phiy;
      //mu = 1.5*sigma*(32.*Phi*(Phi-1.)*(Phi-0.5)/xi-xi*laplPhi);
			mu = 4.*beta*(Phi-PhiL)*(Phi-PhiH)*(Phi-Phi0)-kappa*laplPhi;
      Fsx = mu*grad_phix;
      Fsy = mu*grad_phiy;
      Fmx = Fmy = 0.;
      tau = tauL+(Phi-PhiL)/(PhiH-PhiL)*(tauH-tauL);
      ni = tau*cs2;
      omega = 1./(tau+0.5);
      omega1 = 1.-omega;
			Fmx = Fmy = 0.;
			/*for(int k=0; k<np; k++)
			{
				A = U0*cx[k]+V0*cy[k];
			  feq = wf[k]*(Press0+3.*A+4.5*((cx[k]*cx[k]-cs2)*U0*U0+(cy[k]*cy[k]-cs2)*V0*V0+
			            2.*cx[k]*cy[k]*U0*V0));
			  Fmx += cx[k]*cy[k]*(ftemp-feq);
			  Fmy += cy[k]*cx[k]*(ftemp-feq);
		 	}
      Fmx *= -ni/((tau+0.5)*cs2)*(rhoH-rhoL)/(PhiH-PhiL)*grad_phix;
      Fmy *= -ni/((tau+0.5)*cs2)*(rhoH-rhoL)/(PhiH-PhiL)*grad_phiy;*/
			Fmx = ni* (gradx_of_u*(rhoH-rhoL)/(PhiH-PhiL)*grad_phix + gradx_of_v*(rhoH-rhoL)/(PhiH-PhiL)*grad_phiy);
			Fmy = ni* (grady_of_u*(rhoH-rhoL)/(PhiH-PhiL)*grad_phix + grady_of_v*(rhoH-rhoL)/(PhiH-PhiL)*grad_phiy);
			rho[id] = R = rhoL+(Phi-PhiL)/(PhiH-PhiL)*(rhoH-rhoL);
      Fx = Fpx+Fsx+Fmx;
      Fy = Fpy+Fsy+Fmy;
      Fy += -(R - 0.5*(rhoH+rhoL))*gravity;
			U += 0.5*Fx/R;
      V += 0.5*Fy/R;
			Nx  = grad_phix/(sqrt(pow(grad_phix,2)+pow(grad_phiy,2))+1E-12);
			Ny  = grad_phiy/(sqrt(pow(grad_phix,2)+pow(grad_phiy,2))+1E-12);
			if(y==0 || y==ny-1)
				U = V = Fx = Fy = Nx = Ny = 0.;
			u[id] = U;
			v[id] = V;
			U2 = U*U;
			V2 = V*V;
			UV = U*V;
			for(int k=0; k<np; k++)
			{
        A = U*cx[k]+V*cy[k];
        feq = wf[k]*(Press+3.*A+4.5*((cx[k]*cx[k]-cs2)*U2+(cy[k]*cy[k]-cs2)*V2+
                         2.*cx[k]*cy[k]*UV))-
                         0.5*wf[k]/(R*cs2)*(cx[k]*Fx+cy[k]*Fy);
        f1[id*np+k] = omega1*f1[id*np+k]+omega*feq+wf[k]/(R*cs2)*(cx[k]*Fx+cy[k]*Fy);
        extra_term = 4.*Phi*(1.-Phi)/xi*wf[k]*(cx[k]*Nx+cy[k]*Ny);
        geq = Phi*wf[k]*(1.+3.*A+4.5*((cx[k]*cx[k]-cs2)*U2+(cy[k]*cy[k]-cs2)*V2+
                          2.*cx[k]*cy[k]*UV)) - 0.5*extra_term;
        g1[id*np+k] = omega_phase1*g1[id*np+k] + omega_phase*geq + extra_term;
				newx = x+cx[k];
				newy = y+cy[k];
				if(x==0 || x==nx-1)
					newx = (newx+nx)%nx;
				if(y==0 || y==ny-1)
					newy = (newy+ny)%ny;
				idn = newx*ny+newy;
				f2[idn*np+k] = f1[id*np+k];
				g2[idn*np+k] = g1[id*np+k];
			}
			if(fabs(U)>1. || isnan(U))
				check = 1;
		}
		energy /= nx*ny;
		enstrophy /= nx*ny;
    return check;
}
///----------------------------------------------------------------------------------------------------------------------------------
void boundary_conditions()
{
  for(int x=0; x<nx; x++)
    for(int k=0; k<np; k++)
    {
			if(cy[k]>0)
			{
				id = x*ny+0;
	      f2[id*np+k] = f1[id*np+opp[k]];
	      g2[id*np+k] = g1[id*np+opp[k]];
			}
			if(cy[k]<0)
			{
				id = x*ny+ny-1;
				f2[id*np+k] = f1[id*np+opp[k]];
	      g2[id*np+k] = g1[id*np+opp[k]];
			}
    }
}
///----------------------------------------------------------------------------------------------------------------------------------
int main(int argc, char *argv[])
{
	FILE *data = fopen("data.txt","wt");
	system("mkdir vtk_fluid");
	int check_mach = 0, t;
	initial_state();
  printf("rhoH = %lf, rhoL = %lf\n", rhoH, rhoL);
  printf("gravity = %lf, mobility = %lf, T = %lf\n", gravity, M, T);
  for(t=0; t<nsteps; t++)
  {
    check_mach = algorithm_lattice_boltzmann();
    boundary_conditions();
		f1 = f2;
		g1 = g2;
		if(plot_vtk==true && t%n_out==0)
			write_fluid_vtk(t);
		diagnostics(t);
		dissipation_rate = nu*enstrophy;
		fprintf(data,"%lf    %e    %e    %e\n", (double)t/T, energy*scale_energy/pow(scale_length,2), enstrophy*scale_enstrophy/pow(scale_length,2), dissipation_rate*scale_dissipation/pow(scale_length,2));
		if(check_mach==1) // check the Mach number...if too high, it exits!
      goto labelA;
  }
  labelA:
	  fclose(data);
  return 0;
}
///----------------------------------------------------------------------------------------------------------------------------------
///----------------------------------------------------------------------------------------------------------------------------------
