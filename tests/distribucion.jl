# Script to test Fock distribution of a given state

using QuantumToolbox, SparseArrays, LinearAlgebra
using PyPlot, Printf

include("../modules/fisher_functions.jl")

N=160;
psi1=GKP2(N,2.5,40);
psi2=normalize(squeeze(N,-asinh(sqrt(5)))*fock(N,0))
psi3=normalize(coherent(N,sqrt(5)))

Pn1=real(diag(ket2dm(psi1)))
Pn2=real(diag(ket2dm(psi2)))
Pn3=real(diag(ket2dm(psi3)))
n=1:1:N
n_2=1:2:N

plot(n,Pn2,".r",label="Squeezed state")
plot(n_2,Pn1[1:2:end],".b",label="GKP state")

xlabel("Number of photons", size=14)
ylabel(L"$P_n$",size=14)
yscale("log")
ylim([1e-16,5e-1])
legend()
