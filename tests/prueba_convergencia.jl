# Script to test the convergence of QFI as a function of dl or cutoff
using QuantumToolbox, LinearAlgebra, SparseArrays
using PyPlot, Printf, Polynomials

include("../modulos/fisher_functions.jl")

phi=2.0
N=10:10:150
dl=0.001;
#N=80
#dl=[0.1,0.09,0.08,0.07,0.06,0.05,0.04,0.03,0.02,0.01,0.009,0.008,0.007,0.006,0.005,0.004,0.003,0.002,0.001,0.0009,0.0008,0.0007,0.0006,0.0005,0.0004,0.0003,0.0002,0.0001]
QFI=zeros(length(N),3);
tipo=["coherent","squeezed","GKP2"]
n_GKP=2.76 # number of photons for GKP
n=5 # number of photons
for j=1:length(tipo)-1
	for i=1:length(N)
		psi=estado(tipo[j],N[i],n,0,0)
		psi1=interferometro(tipo,N[i],psi,phi)
		psi2=interferometro(tipo,N[i],psi,phi+dl)
		QFI[i,j]= get_QFI_ket(psi1,psi2,dl)
		@printf("Iteración número %i \n",i)
	end
end

for i=1:length(N)
		psi=estado(tipo[3],N[i],n_GKP,0,50)
		psi1=interferometro(tipo,N[i],psi,phi)
		psi2=interferometro(tipo,N[i],psi,phi+dl)
		QFI[i,3]= get_QFI_ket(psi1,psi2,dl)
		@printf("Iteración número %i \n",i)
end

plot(N,QFI[:,1],".g",label="Coherent state")
plot(N,QFI[:,2],".r",label="Squeezed state")
plot(N,QFI[:,3],".b",label="GKP state")
xlabel("Dimension cutoff N",size=14)
ylabel(L"$F_Q$",size=14)
