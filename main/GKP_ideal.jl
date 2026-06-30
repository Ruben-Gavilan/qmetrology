# Specific script for GKP state, since it needs specific parameter for number of photons

include("../modulos/fisher_functions.jl")
using PyPlot,Polynomials,Printf

N=60
phi=0.1
dl=1e-4
tipo="GKP"

n_GKP=[1.84,2.11,2.33,2.55,2.76]
niter=length(n_GKP)

number= tensor(num(N),eye(N)) + tensor(eye(N),num(N))
QFI=zeros(niter)
photons=zeros(niter)

for i=1:niter
	psi=estado(tipo,N,n_GKP[i],0,40)
	psi0=interferometro(tipo,N,psi,phi)
	psi1=interferometro(tipo,N,psi,phi+dl)
	photons[i]=abs(psi0' * number * psi0)
	QFI[i] = get_QFI_ket(psi0,psi1,dl)
end

reg1=fit(log.(photons),log.(1 ./QFI),1)
label1=@sprintf("GKP, m=%f",reg1[1])
plot(photons,1 ./QFI,"sb",label=label1)
plot(photons,1 ./(photons.^2),"--k",label="HL")
plot(photons, 1 ./photons,".-k",label="SQL")
xscale("log")
yscale("log")
legend()
grid()
plot_title=@sprintf("GKP State, N=%i, dl=%1.1e",N,dl)
title(plot_title)
xlabel("Number of photons")
ylabel("Delta phi^2")
