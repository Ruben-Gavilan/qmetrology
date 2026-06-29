# Main script to obtain QFI of given state in noisless case
# For GKP state, use the specific script
using PyPlot, Printf
using Polynomials

include("../modules/fisher_functions.jl")

N=80
phi=2.0
dl=1e-4
nmax=5
QFI=zeros(nmax)
photons=zeros(nmax)
number=tensor(qeye(N),num(N)) + tensor(num(N),qeye(N)) # Total number of photons operator
print("Input state: \n")
tipo=readline()

for n=1:nmax
	psi=estado(tipo,N,n,0)
	photons[n]=abs(psi'*number*psi)
	#rho1=ket2dm(interferometro(tipo,N,psi,phi))
	#rho2=ket2dm(interferometro(tipo,N,psi,phi+dl))
	psi1=interferometro(tipo,N,psi,phi)
	psi2=interferometro(tipo,N,psi,phi+dl)
	QFI[n]=get_QFI_ket(psi1,psi2,dl)
	@printf("Iteration: %i",n)
end

print("Print to file (F) or plot (P)?: \n")
action=readline()

if action=="F"
	io=open(tipo*"_QFI_ideal.txt","w")
	@printf(io,"#N=%i;	n=%i\n",N,nmax)
	@printf(io,"fot, QFI\n")
	for i=1:length(QFI)
		@printf(io,"%9.7f %9.7f \n",photons[i],QFI[i])
	end
	close(io)
elseif action=="P"
	reg1=fit(log.(fotones),log.(1 ./QFI),1)
	label1=@sprintf("QFI, m=%f",reg1[1])
	plot(photons,1 ./QFI,"sb",label=label1)
	plot(photons,1 ./(photons.^2),"^r",label="HL")
	plot(photons, 1 ./photons,"^g",label="SQL")
	xscale("log")
	yscale("log")
	legend()
	grid()
	xlabel("Average number of photons")
	ylabel(L"$1/F_Q \propto \Delta \varphi^2$")
end
