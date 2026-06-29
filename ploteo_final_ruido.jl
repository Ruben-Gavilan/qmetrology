# Script para plotear la QFI con ruido

using PyPlot
using Printf
using DelimitedFiles,CSV,DataFrames
using Polynomials

tipos=["fock","noon","coherent","squeezed","GKP2"]
colores=["sb","sr","sy","sm","sc"]
fotones=zeros(5,length(tipos))
QFI=zeros(5,length(tipos))

for i=1:length(tipos)
	dlm=readdlm("../../resultados/final/"*tipos[i]*"_QFI_300.txt")
	fotones[:,i]=dlm[3:end,1]
	QFI[:,i]=dlm[3:end,2]
	reg=fit(log.(fotones[:,i]),log.(1 ./QFI[:,i]),1)
	@printf("Tipo %s, m=%f \n",tipos[i],reg[1])
end

for i=1:length(tipos)
	if tipos[i]=="GKP2"
		plot(fotones[:,i], 1 ./QFI[:,i],colores[i],label="GKP")
	else
		plot(fotones[:,i], 1 ./QFI[:,i],colores[i],label=tipos[i])
	end
end


fot_plot=[minimum(fotones),maximum(fotones)]

plot(fot_plot, 1 ./fot_plot, "--k",label="SQL")
plot(fot_plot, 1 ./fot_plot.^2, "-.k",label="HL")

xscale("log")
yscale("log")
xlim([1.9,maximum(fotones)+1])
ylim([6e-3,1])
legend()
grid()
#title(L"Noise channel, $\kappa_t = 0.01$",size=14)
xlabel("Number of photons",size=14)
ylabel(L"$1/F_Q$",size=14)

#savefig("../resultados/sv.png")
