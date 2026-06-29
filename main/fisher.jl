using QuantumToolbox
using PyPlot
using Printf
using SparseArrays
using LinearAlgebra
using Polynomials

include("../modulos/fisher_functions.jl")

N=78
phi=2.0
dl=0.0001
nmax=5
delta_tilde=1/5.0
s=40
QFI=zeros(nmax)
fotones=zeros(nmax)
numero=tensor(qeye(N),num(N)) + tensor(num(N),qeye(N))
print("Introduce el tipo de estado: \n")
tipo=readline()

for n=1:nmax
	psi=estado(tipo,N,n,delta_tilde,s)
	fotones[n]=abs(psi'*numero*psi)
	#rho1=ket2dm(interferometro(tipo,N,psi,phi))
	#rho2=ket2dm(interferometro(tipo,N,psi,phi+dl))
	psi1=interferometro(tipo,N,psi,phi)
	psi2=interferometro(tipo,N,psi,phi+dl)
	QFI[n]=get_QFI_ket(psi1,psi2,dl)
	@printf("Iteración número: %i",n)
end

#io=open(tipo*"_QFI_dm.txt","w")
#@printf(io,"#N=%i;	n=%i\n",N,nmax)
#@printf(io,"fot, QFI\n")
#for i=1:length(QFI)
#	@printf(io,"%9.7f %9.7f \n",fotones[i],QFI[i])
#end
#close(io)

#run(Cmd(`shutdown now`));
reg1=fit(log.(fotones),log.(1 ./QFI),1)
reg2=fit(log.(fotones),log.(1 ./(fotones .*(2 .+fotones))),1)
label1=@sprintf("QFI, m=%f",reg1[1])
label2=@sprintf("Teórico sin ruido: m=%f",reg2[1])
plot(fotones,1 ./QFI,"sb",label=label1)
plot(fotones,1 ./(fotones.^2),"^r",label="HL")
plot(fotones, 1 ./fotones,"^g",label="SQL")
plot(fotones, 1 ./(fotones .*(2 .+fotones)),"y",label=label2)
xscale("log")
yscale("log")
legend()
grid()
titulito=@sprintf("Squeezed Vacuum State, N=%i, dl=%1.1e",N,dl)
title(titulito)
xlabel("Número de fotones (sinh(r²))")
ylabel("Delta phi^2")

#title!("Squeezed state")
