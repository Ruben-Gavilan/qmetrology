# Script para hacer el ploteo de los distintos estados sin ruido.
# Por ahora haré el cálculo y el ploteo a la vez, al ser sin ruido debería ser sencillo.

using QuantumToolbox, SparseArrays, LinearAlgebra
using PyPlot, Printf, Polynomials

include("../modulos/fisher_functions.jl")

N=80
phi=2.0
dl=0.0001
nmax=5
s=40
tipos=["fock","noon","coherent","squeezed","GKP2"]
colores=["sb","sr","sy","sm","sc"]
n_GKP=[1.84,2.11,2.33,2.55,2.76]
QFI=zeros(nmax,length(tipos));
fotones=zeros(nmax,length(tipos));
numero=tensor(qeye(N),num(N)) + tensor(num(N),qeye(N))
for i=1:length(tipos)-1
	@printf("Tipo:	%s \n\n",tipos[i])
	for n=1:nmax
		psi=estado(tipos[i],N,n,0,s)
		fotones[n,i] = abs(psi' * numero * psi)
		psi1=interferometro(tipos[i],N,psi,phi)
		psi2=interferometro(tipos[i],N,psi,phi+dl)
		QFI[n,i] = get_QFI_ket(psi1,psi2,dl)
		@printf("Iteración: %i \n",n)
	end
#	plot(fotones[:,i], 1 ./QFI[:,i], colores[i], label=tipos[i])
#	reg=fit(log.(fotones[:,i]), log.(1 ./QFI[:,i]),1)
#	@printf("Tipo: %s; m=%f\n",tipos[i],reg[1])
end

for n=1:length(n_GKP)
	psi=estado("GKP2",N,n_GKP[n],0,s)
	fotones[n,end] = abs(psi' * numero * psi)
	psi1=interferometro(tipos[end],N,psi,phi)
	psi2=interferometro(tipos[end],N,psi,phi+dl)
	QFI[n,end] = get_QFI_ket(psi1,psi2,dl)
	@printf("Iteración: %i \n",n)
end

io=open("../../resultados/final/QFI_ideal.txt","w")
@printf(io,"fot, Fock,	NOON,	Coh,	SV,	GKP\n")
for i=1:length(QFI[:,1])
	@printf(io,"%9.5f	%9.5f	%9.5f	%9.5f	%9.5f	%9.5f\n",fotones[i],QFI[i,1],QFI[i,2],QFI[i,3],QFI[i,4],QFI[i,5])
end
close(io)

#plot(fotones[:,end], 1 ./QFI[:,end], colores[end], label="GKP")


#fot_plot=[2,maximum(fotones)]

#plot(fot_plot, 1 ./fot_plot , "--k",label="SQL")
#plot(fot_plot, 1 ./fot_plot.^2, "-.k",label="HL")

#for i=1:length(tipos)
#	plot(fotones[:,i], 1 ./QFI[:,i], colores[i], label=tipos[i])
#end
#xlim([1.9,maximum(fotones)+1])
#ylim([6e-3,1])
#xscale("log")
#yscale("log")
#grid()
#xlabel("Number of photons",size=14)
#ylabel(L"1/$F_Q$",size=14)
#title("Noisless case",size=14)
#legend(bbox_to_anchor=[1.05,1],loc=2,borderaxespad=0)
