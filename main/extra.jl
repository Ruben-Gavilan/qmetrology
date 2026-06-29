# script para calcular la QFI de la forma extra mecagoenmiputamadre

#using QuantumToolbox, SparseArrays, LinearAlgebra
using Printf
using Polynomials, PyPlot

include("../modulos/fisher_functions.jl")
N=50
phi=2.0
dl=0.0001
nmax=5

T=1e-4;
kappa=300; #Valor típìco en microondas (REVISAR)
time=LinRange(0.0,T,2);

QFI=zeros(nmax)
#QFI_noisless=zeros(nmax)
fotones=zeros(nmax)
#fotones_noisless=zeros(nmax)
numero=tensor(qeye(N),num(N)) + tensor(num(N),qeye(N))

c_ops=[sqrt(kappa)*tensor(destroy(N),eye(N)), sqrt(kappa)*tensor(eye(N),destroy(N))]

#print("Introduce el tipo de estado: \n")
#tipo=readline()
tipo="cat"

for n=1:nmax
	psi=estado(tipo,N,n,0,40)
	psi0=interferometro(tipo,N,psi,phi)
	psi1=interferometro(tipo,N,psi,phi+dl)
	psi2=interferometro(tipo,N,psi,phi-dl)
	@printf("Hecho interferómetro \n")
	
	
#	fotones_noisless[n]=abs(psi0' * numero * psi0)
#	drho_noisless=(ket2dm(psi1) - ket2dm(psi2))/(2*dl)
#	QFI_noisless[n]=get_QFI_derivada(ket2dm(psi0),drho_noisless)
#	@printf("QFI sin ruido hecha \n")
#	drho_noisless=nothing
	
	rho1=mesolve(zero(tensor(eye(N),eye(N))),psi1,time,c_ops, e_ops=[ ],progress_bar=Val(false))
	rho1=rho1.states[2];
	rho2=mesolve(zero(tensor(eye(N),eye(N))),psi2,time,c_ops, e_ops=[ ],progress_bar=Val(false))
	rho2=rho2.states[2];
#	rho0=loss_channel_kraus_n(psi0,N,eta,20)
#	rho1=loss_channel_kraus_n(psi1,N,eta,20)
#	rho2=loss_channel_kraus_n(psi2,N,eta,20)
	
	drho= (rho1 - rho2) / (dl*2)
	rho1=nothing
	rho2=nothing
	rho0=mesolve(zero(tensor(eye(N),eye(N))),psi0,time,c_ops, e_ops=[ ],progress_bar=Val(false))
	rho0=rho0.states[2];
	fotones[n]=abs(tr(rho0*numero))
	QFI[n]=get_QFI_derivada(rho0,drho)
	drho=nothing
	rho0=nothing
	@printf("Iteración número: %i\n",n)
end
#io=open("../../resultados/final/"*tipo*"_QFI_300.txt","w")
#@printf(io,"#N=%i;	n=%i;	dl=%1.1e\n",N,nmax,dl)
#@printf(io,"fot, QFI\n")
#for i=1:length(QFI)
#	@printf(io,"%9.5f %9.5f \n",fotones[i],QFI[i])
#end
#close(io)

#run(Cmd(`shutdown now`));

reg1=fit(log.(fotones),log.(1 ./QFI),1)
#reg2=fit(log.(fotones_noisless),log.(1 ./QFI_noisless),1)
#reg3=fit(log.(fotones), log.(1 ./(fotones .* (2 .+ fotones))),1)
#reg4=fit(log.(fotones_noisless), log.(1 ./(fotones_noisless .* (2 .+ fotones_noisless))),1)
label1=@sprintf("QFI, m=%f",reg1[1])
#label2=@sprintf("Noisless QFI, m=%f",reg2[1])
#label3=@sprintf("Teórico, m=%f",reg3[1])
#label4=@sprintf("Noisless Teórico, m=%f",reg4[1])
plot(fotones, 1 ./QFI, "sb",label=label1)
#plot(fotones_noisless, 1 ./QFI_noisless, "sm",label=label2)
plot(fotones, 1 ./(fotones.^2), "^r",label="HL")
plot(fotones, 1 ./fotones, "^g",label="SQL")
#if tipo=="squeezed"
#	plot(fotones, 1 ./(fotones .* (2 .+fotones)),"y",label=label3)
#	plot(fotones_noisless, 1 ./(fotones_noisless .* (2 .+fotones_noisless)),"k",label=label4)
#end
xscale("log")
yscale("log")
legend()
grid()
#title("N00N State, noise by Krauss 20,QFI nueva")
xlabel("Número de fotones")
ylabel("1/QFI")
