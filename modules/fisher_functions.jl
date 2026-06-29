# Transcripcción de fisher.py a julia a ver si va más rapido jejeej

using QuantumToolbox
using SparseArrays
using LinearAlgebra

function estado(tipo,N,n,r,s)
# -------------------INPUT------------------
# tipo:		string, tipo de estado a usar
# N:		int, dimensión del espacio de Hilbert (el espacio final será N²)
# n:		int, número de iteración, relacionado por lo general con n_fotones
# r:		float, módulo del factor de squeezing constante para squeezed_coherent y GKP
# s:		int, límites del sumatorio para GKP
# ------------------------------------------
	if tipo=="fock"
		psi=basis(N,n)
		psi=tensor(psi,psi)
	elseif tipo=="coherent"
		psi=coherent(N,sqrt(n)*exp(0.25im*pi))
		psi=tensor(psi,psi)
	elseif tipo=="squeezed"
		psi=squeeze(N,-asinh(sqrt(n)))*basis(N,0)
		psi=tensor(psi,psi)
	elseif tipo=="squeezed_coherent"
		psi=displace(N,sqrt(n))*squeeze(N,r)*fock(N,0)
		psi=tensor(psi,psi)
	elseif tipo=="noon"
		psi=tensor(basis(N,2*n),basis(N,0)) + tensor(basis(N,0),basis(N,2*n))
	elseif tipo=="GKP"
		psi=GKP(N,n,r,s)
		psi=tensor(psi,psi)
	elseif tipo=="GKP2"
		psi=GKP2(N,n,s)
		psi=tensor(psi,psi)
	elseif tipo=="cat"
		psi= coherent(N,sqrt(n)*exp(0.25im*pi)) + coherent(N,-sqrt(n)*exp(0.25im*pi))
		psi=tensor(psi,psi)
	end
	return normalize(psi)
end

function GKP(N,n,r,s)
	psi0=normalize(squeeze(N,asinh(sqrt(n)))*fock(N,0))
	psi1=zero_ket(N)
	for i=-s:1:s
		psi1= psi1 + exp(-2*pi*(r^2) * (i^2))*displace(N,i*sqrt(2*pi))*psi0
	end
	return psi1
end

function GKP2(N,n,s)
	r=1/n
	psi0=normalize(squeeze(N,-log(r))*fock(N,0))
	psi1=zero_ket(N)
	for i=-s:1:s
		psi1 += exp(-2*pi*(r^2)*(i^2))*displace(N,2*i*sqrt(pi))*psi0
	end
	return psi1
end

function interferometro(tipo,N,psi,phi)
	U_BS=exp(-0.25im*pi*(tensor(destroy(N),create(N)) + tensor(create(N),destroy(N))))
	U_PS=tensor(qeye(N),exp(1im*phi*num(N)))
	if tipo=="noon"
		return U_BS*U_PS*psi
	else
		return U_BS*U_PS*U_BS*psi
	end
end

#function interferometro_dm(tipo,N,rho,phi)
#	U_BS=exp(0.25im*pi*(tensor(destroy(N),create(N)) + tensor(create(N),destroy(N)) ))
#	U_PS=tensor(qeye(N),exp(1im*phi*num(N)))
#	if tipo=="noon"
#		return U_BS*U_PS*rho*U_PS'*U_BS'
#	else
#		return U_BS*U_PS*U_BS*rho*U_BS'*U_PS'*U_BS'
#	end
#end

function get_QFI_dm(rho1,rho2,dl)
	D=1.0-(fidelity(rho1,rho2))
	return 8*D/(dl^2)
end

function get_QFI_ket(psi1,psi2,dl)
	D= 1.0 - abs(psi1'*psi2)
	return 8*D/(dl^2)
end

#function get_QFI_gpu(rho1,rho2,dl)
#	D= 1.0 - tr(sqrt(sqrt(rho1)*rho2*sqrt(rho1)))
#	return 8*D/(dl^2)
#end

#function get_QFI_21(rho1,rho2,dl)
#	a1,psi1,U=eigenstates(rho1);
#	a2,psi2,U=eigenstates(rho2);
#	F= sum((((a2.-a1)./dl).^2)./a1);
#	for i=1:length(a1)
#		for j=1:length(a2)
#				F += 2 * ((a1[i]-a1[j])^2) * (abs(psi1[i]' * (psi2[j]-psi1[j]) ./dl))^2 / (a1[i] + a1[j])
#		end
#	end
#	return F
#end

function get_QFI_derivada(rho,rho_derivada)
# -------------------INPUT------------------
# rho:		operador, matriz densidad del interferómetro con phi
# rho1:		operador, matriz densidad del interferómetro con phi+dl
# rho2:		operador, matriz densidad del interferómetro con phi-dl
# ------------------------------------------
	#rho_derivada= (rho1 - rho2)/(2*dl)
	lambda,psi,U=eigenstates(rho)
	rho_derivada= dag(U) * rho_derivada.data * U
	F=0
	lambda=real(lambda);
	for i=1:length(lambda)
		for j=1:length(lambda)
			if lambda[i]+lambda[j] > 1e-8
				F += (abs(rho_derivada[j,i])^2)/(lambda[i] + lambda[j])
			end
		end
	end
	return 2*F
end
