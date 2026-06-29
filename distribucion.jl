# Script para mirar la distribución en el espacio de Fock de un estado dado

using QuantumToolbox, SparseArrays, LinearAlgebra
using PyPlot, Printf

include("../modulos/fisher_functions.jl")

N=160;
psi1=GKP2(N,2.5,40);
psi2=normalize(squeeze(N,-asinh(sqrt(5)))*fock(N,0))
psi3=normalize(coherent(N,sqrt(5)))
#@printf("Número de fotones del estado: %f\n",abs(psi'*num(N)*psi))

#rho1=ket2dm(psi1);

Pn1=real(diag(ket2dm(psi1)))
Pn2=real(diag(ket2dm(psi2)))
Pn3=real(diag(ket2dm(psi3)))
n=1:1:N
n_2=1:2:N
#bar(n,Pn1,"g",width=0.8);
#bar(n,Pn2,"r");
#bar(n,Pn3,"b");
plot(n,Pn3,".g",label="Coherent state")
plot(n,Pn2,".r",label="Squeezed state")
plot(n_2,Pn1[1:2:end],".b",label="GKP state")
#plot(n,0.001*ones(N),"r")
#plot(n,0.0001*ones(N),"r")
xlabel("Number of photons", size=14)
ylabel(L"$P_n$",size=14)
yscale("log")
ylim([1e-16,5e-1])
legend()
