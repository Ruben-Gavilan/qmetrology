# Script específico para GKP, ajustando el número de fotones para que encaje con el resto.

include("../modulos/fisher_functions.jl")
using PyPlot,Polynomials,Printf

N=75 # Para un maximo de 2.76 está convergido a partir de ~30
phi=2.0
dl=1e-4
tipo="GKP2"

T=1e-4;
kappa=300; #Valor típìco en microondas (REVISAR)
time=LinRange(0.0,T,2);

n_GKP=[1.84,2.11,2.33,2.55,2.76]
niter=length(n_GKP)

numero= tensor(num(N),eye(N)) + tensor(eye(N),num(N))
QFI=zeros(niter)
fotones=zeros(niter)

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
	
	fotones[i]=abs(tr(rho0*numero))
	QFI[i]=get_QFI_derivada(rho0,drho)
	drho=nothing
	rho0=nothing
	@printf("Iteración número: %i\n",i)
end

io=open("../../resultados/final/"*tipo*"_QFI_300_final.txt","w")
@printf(io,"#N=%i;	n=%i;	dl=%1.1e\n",N,niter,dl)
@printf(io,"fot, QFI\n")
for i=1:length(QFI)
	@printf(io,"%9.5f %9.5f \n",fotones[i],QFI[i])
end
close(io)

#reg1=fit(log.(fotones),log.(1 ./QFI),1)
#label1=@sprintf("QFI, m=%f",reg1[1])
#plot(fotones,1 ./QFI,"sb",label=label1)
#plot(fotones,1 ./(fotones.^2),"^r",label="HL")
#plot(fotones, 1 ./fotones,"^g",label="SQL")
#xscale("log")
#yscale("log")
#legend()
#grid()
#titulito=@sprintf("GKP State, N=%i, dl=%1.1e",N,dl)
#title(titulito)
#xlabel("Número de fotones")
#ylabel("Delta phi^2")
