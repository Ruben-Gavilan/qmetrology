# Module script to load all necesary functions

using QuantumToolbox, SparseArrays, LinearAlgebra

function estado(tipo,N,n,r,s)
# -------------------INPUT------------------
# tipo:		string, state to compute
# N:		int, Hilbert space dimension (final space will be N²)
# n:		int, number of photons for the state
# s:		int, sum limit of GKP
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
	elseif tipo=="noon"
		psi=tensor(basis(N,2*n),basis(N,0)) + tensor(basis(N,0),basis(N,2*n))
	elseif tipo=="GKP"
		psi=GKP2(N,n,s)
		psi=tensor(psi,psi)
	elseif tipo=="cat"
		psi= coherent(N,sqrt(n)*exp(0.25im*pi)) + coherent(N,-sqrt(n)*exp(0.25im*pi))
		psi=tensor(psi,psi)
	end
	return normalize(psi)
end

function GKP(N,n,s)
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

function get_QFI_ket(psi1,psi2,dl)
	D= 1.0 - abs(psi1'*psi2)
	return 8*D/(dl^2)
end

function get_QFI_derivada(rho,drho)
# -------------------INPUT------------------
# rho:		operator, density matrix of system with phi
# dro:		operator, derivative of density matrix
# ------------------------------------------
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
