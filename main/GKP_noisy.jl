# Specific script for GKP state for noisy environment
include("../modulos/fisher_functions.jl")
using PyPlot,Polynomials,Printf

N=60
phi=2.0
dl=1e-4
tipo="GKP"

T=1e-4;
kappa=300;
time=LinRange(0.0,T,2);

n_GKP=[1.84,2.11,2.33,2.55,2.76]
niter=length(n_GKP)

number= tensor(num(N),eye(N)) + tensor(eye(N),num(N))
QFI=zeros(niter)
photons=zeros(niter)

c_ops=[sqrt(kappa)*tensor(destroy(N),eye(N)), sqrt(kappa)*tensor(eye(N),destroy(N))]

for i=1:niter
	psi=estado(tipo,N,n_GKP[i],0,40)
	psi0=interferometro(tipo,N,psi,phi)
	psi1=interferometro(tipo,N,psi,phi+dl)
	psi2=interferometro(tipo,N,psi,phi-dl)
	
	rho1=mesolve(zero(tensor(eye(N),eye(N))),psi1,time,c_ops, e_ops=[ ],progress_bar=Val(false))
	rho1=rho1.states[2];
	rho2=mesolve(zero(tensor(eye(N),eye(N))),psi2,time,c_ops, e_ops=[ ],progress_bar=Val(false))
	rho2=rho2.states[2];
	drho= (rho1-rho2) / (dl*2) 
	rho1=nothing
	rho2=nothing
	rho0=mesolve(zero(tensor(eye(N),eye(N))),psi0,time,c_ops, e_ops=[ ],progress_bar=Val(false))
	rho0=rho0.states[2];
	
	photons[i]=abs(tr(rho0*number))
	QFI[i]=get_QFI_derivada(rho0,drho)
	drho=nothing
	rho0=nothing
	@printf("Iteración número: %i\n",i)
end

io=open(tipo*"_QFI_300_final.txt","w")
@printf(io,"#N=%i;	n=%i;	dl=%1.1e\n",N,niter,dl)
@printf(io,"fot, QFI\n")
for i=1:length(QFI)
	@printf(io,"%9.5f %9.5f \n",photons[i],QFI[i])
end
close(io)
