# Script para implementar las funciones necesarias para el canal ruidoso.
using QuantumToolbox
using SparseArrays
using LinearAlgebra
#using AMDGPU

function lindblad(N,rho,eta,mode)
# Función que devuelve el Lindbladiano correspondiente a photon loss
# ----------------------INPUT----------------------------------
# rho:		matriz densidad a la que calcularle el lindbladiano
# mode:		0 ó 1, modo óptico al que aplicar. 0=> a1; 1=> a2
# kappa:	loss rate. real

	if mode==0
		a=tensor(destroy(N),qeye(N))
	elseif mode==1
		a=tensor(qeye(N),destroy(N))
	end
	L= eta * (a*rho*a' - 0.5*commutator(a'*a,rho,anti=true))
end

function ruidoso(N,rho,eta)
	L=lindblad(N,rho,eta,0) + lindblad(N,rho,eta,1)
	return exp(L)*rho*exp(L')
end

function lindblad_vec(N,eta,mode)
	if mode==0
		a=tensor(destroy(N),qeye(N))
	elseif mode==1
		a=tensor(qeye(N),destroy(N))
	end
	L= eta * (tensor(a',a) - 0.5*tensor(tensor(qeye(N),qeye(N)),a'*a) - 0.5*tensor(trans(a'*a),tensor(qeye(N),qeye(N))))
end

function ruidoso_vec(N,psi,Tk,kappa)
	L=exp(lindblad_vec(N,kappa,0) + lindblad_vec(N,kappa,1))
	return L*psi
end

function ruidoso_kraus(N,rho,eta,n)
# Función que implementa el canal de photon loss mediante operadores de Krauss
# ------------------------------INPUT-------------------------------------
# N:	dimensión del espacio de Hilbert
# rho:	matriz densidad a la que aplicar el canal de photon loss
# eta:	exp(-kappa*T)
# n:	número del mayor operador de Krauss
	if n==1
		gamma=1-eta
		E0= tensor(qeye(N),qeye(N)) - 0.5*gamma*(tensor(num(N),qeye(N)) + tensor(qeye(N),num(N)))
		E1= sqrt(gamma)*tensor(destroy(N),destroy(N))	
		return normalize(E0*rho*E0' + E1*rho*E1')
	elseif n==2
		gamma=1-eta
		E0=tensor(qeye(N),qeye(N)) - 0.5*gamma*(tensor(num(N),qeye(N)) + tensor(qeye(N),num(N)))
		E1= sqrt(gamma)*tensor(destroy(N),destroy(N))
		E2= (gamma/sqrt(2))* tensor(destroy(N),destroy(N))^2
		return E0*rho*E0' + E1*rho*E1' + E2*rho*E2'
	elseif n==3
		gamma=1-eta
		E0=tensor(qeye(N),qeye(N)) - 0.5*gamma*(tensor(num(N),qeye(N)) + tensor(qeye(N),num(N)))
		E1= sqrt(gamma)*tensor(destroy(N),destroy(N))
		E2= (gamma/sqrt(2))* (tensor(destroy(N),destroy(N))^2)
		E3= (gamma^(3/2)/sqrt(6)) * (tensor(destroy(N),destroy(N))^3)
		return E0*rho*E0' + E1*rho*E1' + E2*rho*E2' + E3*rho*E3'
	end
end

function krauss_n(N,eta,l)
# Función que define el operador de Kraus de orden l, con transmisividad eta y dimensión NxN
	eta_n=exp(0.5*log(eta)*num(N)) # eta^(n/2)
	return (sqrt(1-eta)^l/factorial(l))*eta_n * (destroy(N)^l)
end

function loss_channel_kraus_n(rho_or_psi,N,eta,lmax)
	if lmax==0
		lmax=N-1
	end
	if isket(rho_or_psi)
		rho=ket2dm(rho_or_psi)
	else
		rho=rho_or_psi
	end
	rho_out=0*rho
	for i=0:lmax
		Ks=tensor(krauss_n(N,eta,i),krauss_n(N,eta,i))
		rho_out += Ks*rho*(Ks')
	end
	return rho_out
end
