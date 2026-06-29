# Script específico para GKP, ajustando el número de fotones para que encaje con el resto.

include("../modulos/fisher_functions.jl")
using PyPlot,Polynomials,Printf

N=80 # Para un maximo de 2.76 está convergido a partir de ~30
phi=0.1
dl=1e-4
tipo="GKP2"

n_GKP=[1.84,2.11,2.33,2.55,2.76]
#n_GKP=1.8:0.1:5
niter=length(n_GKP)

numero= tensor(num(N),eye(N)) + tensor(eye(N),num(N))
QFI=zeros(2,niter)
fotones=zeros(2,niter)

for i=1:niter
	psi=estado(tipo,N,n_GKP[i],0,40)
	psi0=interferometro(tipo,N,psi,phi)
	psi1=interferometro(tipo,N,psi,phi+dl)
	fotones[1,i]=abs(psi0' * numero * psi0)
	QFI[1,i] = get_QFI_ket(psi0,psi1,dl)
end

tipo="squeezed"
for i=1:5
	psi=estado(tipo,N,i,0,40)
	psi0=interferometro(tipo,N,psi,phi)
	psi1=interferometro(tipo,N,psi,phi+dl)
	fotones[2,i]=abs(psi0' * numero * psi0)
	QFI[2,i] = get_QFI_ket(psi0,psi1,dl)
end

reg1=fit(log.(fotones[1,:]),log.(1 ./QFI[1,:]),1)
reg2=fit(log.(fotones[2,:]),log.(1 ./QFI[2,:]),1)
label1=@sprintf("GKP, m=%f",reg1[1])
label2=@sprintf("Squeezed, m=%f",reg2[1])
plot(fotones[1,:],1 ./QFI[1,:],"sb",label=label1)
plot(fotones[2,:],1 ./QFI[2,:],"sr",label=label2)
plot(fotones[1,:],1 ./(fotones[1,:].^2),"--k",label="HL")
plot(fotones[1,:], 1 ./fotones[1,:],".-k",label="SQL")
xscale("log")
yscale("log")
legend()
grid()
titulito=@sprintf("GKP State, N=%i, dl=%1.1e",N,dl)
title(titulito)
xlabel("Número de fotones")
ylabel("Delta phi^2")
