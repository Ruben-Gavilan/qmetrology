# Script to compute QFI for different states in noisy channel. To obtain QFI for GKP states use specific script.
using Printf

include("../modules/fisher_functions.jl")
N=80
phi=2.0
dl=1e-4
nmax=5

T=1e-4; # s time of the evolution
kappa=300; # s^-1 Loss rate

QFI=zeros(nmax)
photons=zeros(nmax)
number=tensor(qeye(N),num(N)) + tensor(num(N),qeye(N))

c_ops=[sqrt(kappa)*tensor(destroy(N),eye(N)), sqrt(kappa)*tensor(eye(N),destroy(N))]

print("Input state: \n")
tipo=readline()

for n=1:nmax
	psi=estado(tipo,N,n,0,40)
	psi0=interferometro(tipo,N,psi,phi)
	psi1=interferometro(tipo,N,psi,phi+dl)
	psi2=interferometro(tipo,N,psi,phi-dl)
	@printf("Applied interferometer \n")
	
	rho1=mesolve(zero(tensor(eye(N),eye(N))),psi1,T,c_ops, e_ops=[ ],progress_bar=Val(false))
	rho1=rho1.states[1];
	rho2=mesolve(zero(tensor(eye(N),eye(N))),psi2,T,c_ops, e_ops=[ ],progress_bar=Val(false))
	rho2=rho2.states[1];
	
	drho= (rho1 - rho2) / (dl*2)
	rho1=nothing # Erase density matrices as soon as possible to save memory
	rho2=nothing
	
	rho0=mesolve(zero(tensor(eye(N),eye(N))),psi0,time,c_ops, e_ops=[ ],progress_bar=Val(false))
	rho0=rho0.states[2];
	
	photons[n]=abs(tr(rho0*number))
	QFI[n]=get_QFI_derivada(rho0,drho)
	drho=nothing
	rho0=nothing
	@printf("Iteraction: %i\n",n)
end
io=open(tipo*"_QFI_noisy.txt","w")
@printf(io,"#N=%i;	n=%i;	dl=%1.1e\n",N,nmax,dl)
@printf(io,"phot, QFI\n")
for i=1:length(QFI)
	@printf(io,"%9.5f %9.5f \n",photons[i],QFI[i])
end
close(io)
