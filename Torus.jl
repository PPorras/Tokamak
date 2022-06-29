###################################################
############### Parameterization ##################
###################################################

using TaylorIntegration
using Dierckx
using Plots, LaTeXStrings
using DelimitedFiles
include("Auxiliary-Functions.jl")

###################################################

titleFile = "Torus0.35597881128974557modulus.txt"
titleFileOne = "Torus0.35597881128974557.txt"

###################################################

##const ε = 0.004 #Perturbation
const ε = 0.0 #Perturbation
psi_0  = 0.35
w(x) = (2 - x)*(2 - 2*x + x^2)/4
dw(x) = ((-1)*(2 -2*x + x^2) + ( 2 - x )*(-2 + 2*x))/4
w_0 = w(psi_0)
dw_0 = dw(psi_0)
#const omega = 3.6430966211675684 #Def con inf
#const omega = 3.6430704771333495 ### calculado con 1000000 puntos
##const omega = 3.6430690936048755
#const omega = 3.643074351 ### calculado por Arturo
const alpha = 1.0 
params = [w_0, dw_0]
num = 64

###################################################
#################### Field ########################
###################################################

function Field!(dx, x, params, t)
 	w_0 = params[1]
 	dw_0 = params[2]
	dx[1] = (2 - x[2])*(2 - 2*x[2] + x[2]^2)/4
 	dx[2] = ε*(2*sin(2*x[1] - t) + 3*sin(3*x[1] - 2*t)) - (ε^2)*dw_0*(2*cos(2*x[1] - t)/(2*w_0 - 1) + 3*cos(3*x[1] - 2*t)/(3*w_0 - 2))*(4*sin(2*x[1] - t)/(2*w_0 - 1) + 9*sin(3*x[1]- 2*t)/(3*w_0 - 2))
 	nothing

 end
###################################################
#################### Filter #######################
###################################################

#θ, out = filter(titleFile, num)

#x = out[:, 1]
#y = out[:, 2]

#torus = Spline1D(x, y, periodic = true)

###################################################
############## Two dimensional Torus ##############
###################################################

x, omegaa = dataGraph(titleFileOne)
const omega = omegaa
println("numero de rotacion = ", omega, "\n")
xtaux = convertTheta(x, omega)
xtaux[:, 1] = mod.(xtaux[:, 1], 2π)
bubbleSort!(xtaux)

xθ =  zeros(length(x[:, 1]) + 1, 3)

xθ[1: end-1, 1] = xtaux[:, 1]
xθ[1: end-1, 2] = xtaux[:, 2]
xθ[1: end-1, 3] = xtaux[:, 3]
xθ[end, :] = [2π xtaux[1, :][2] xtaux[1, :][3]]

xtaux = nothing

println("longitud \n", length(xθ), "\n")

#println("theta \n", xθ[:, 1], "\n")
#println("x= \n", xθ[:, 2], "\n")
#println("y= \n", xθ[:, 3], "\n")
#### s =  collect(range(0.0 , stop = 2π ,length = 64))

xp = Spline1D(xθ[:, 1], xθ[:, 2]; w=ones(length(xθ[:, 1])), k=3, periodic = true, s=0.0)
yp = Spline1D(xθ[:, 1], xθ[:, 3]; w=ones(length(xθ[:, 1])), k=3, periodic = true, s=0.0)

###################################################
####   ####
###################################################

function checkVectorRot(x0)

tv = 0.0:π:2π
xvs = taylorinteg(Field!,  x0, tv, 32, 1.0E-20, params; maxsteps=7_000_000)
#println("punto en el toro  = ", x0, "\n integracion 2π = ", xvs[end, :], "\n punto inicial mas  numero de rotacion =  ", [omega, alpha] + x0 )

end

#aux = mod2pi(1.5 - (omega*0.7)/(2π))
#xi = 1.5 - (omega*0.7)/(2π) + xp(aux)
#yi = yp(aux)
#checkVectorRot([xi, yi])
xθ =  nothing
###################################################
#### Construction of the two-dimensional torus ####
###################################################
function torusTwoD(lenθ, lenφ)

	
	X, Y = zeros(lenθ, lenφ), zeros(lenθ, lenφ)
	θ = collect(range(0.0 , stop = 2π, length = lenθ))
	φ = collect(range(0.0 , stop = 2π, length = lenφ))
	X[:, 1], Y[:, 1]= xp.(θ), yp.(θ)

	### println("X=[:,1]  ", X[:, 1])
	
	j = 2

	while j <= lenφ 
		
		i =  1

		while i <= lenθ	

			aux = mod(θ[i] - (omega*φ[j])/2π, 2π) 
			xinit = θ[i] - (omega*φ[j])/2π + xp(aux)
			yinit = yp(aux)

			tv = 0.0:0.5*φ[j]:φ[j] 
			x0 = [xinit, yinit]
			xvs = taylorinteg(Field!,  x0, tv, 3, 1.0E-20,
			params; maxsteps=7_000_000)

			X[(j - 1)*lenθ + i] = xvs[end, :][1] - θ[i]
			Y[(j - 1)*lenθ + i] = xvs[end, :][2] 
			#X[:, j][i] = xvs[end, :][1] - θ[i]
			#Y[:, j][i] = xvs[end, :][2] 
			i = i + 1
		end

		#### println("X=[:,$j]  ", X[:, j])
		j = j + 1

	end

#plot(size=(800, 600))
#plot(Xθ[:, 1], Y[:, 1])
#surf(Xθ, Y , φ)
#savefig("Toro2Dphi.pdf")
	return θ, φ, X, Y

end

#######################################################
###### Derivatives of the components of the torus #####
#######################################################

function derivativesData(θ, φ, data)

	lenθ = length(θ)
	lenφ = length(φ)
	
	derθ = zeros(lenθ, lenφ)
	derφ = zeros(lenθ, lenφ)

	i = 1
	
	while i <= lenφ

		SD = Spline1D(θ, data[:, i]; w=ones(length(θ)), k=3, periodic = true, s=0.0)
		derθ[:, i] = Dierckx.derivative(SD, θ)

		i = i + 1

	end

	i = 1
	
	while i <= lenθ

		SD = Spline1D(φ, data[i, :]; w=ones(length(φ)), k=3, periodic = true, s=0.0)
		derφ[i, :] = Dierckx.derivative(SD, φ)

		i = i + 1

	end
	
return derθ, derφ

end

#######################################################
################### Invariance error ##################
#######################################################

function invarianceError(θ::Vector{Float64}, φ::Vector{Float64},
 Xp::Matrix{Float64}, Yp::Matrix{Float64})::Float64

	lenθ, lenφ = length(θ), length(φ)

	DXPθ, DXPφ = derivativesData(theta, varphi, XP)
	DYPθ, DYPφ = derivativesData(theta, varphi, YP)

	j, invError = 1, 0

	while j <= lenφ

		i = 1
		
		#println("Complete for K_{$(j - 1)}")
		while i <= lenθ
		
			 
			xComp =	(1.0 + DXPθ[:, j][i])*((omega*alpha)/(2π)) + DXPφ[:, j][i]*alpha		
			yComp =	(DYPθ[:, j][i]*omega*alpha)/(2π) + DYPφ[:, j][i]*alpha		
			DK = [xComp, yComp]
			K = [θ[i] + Xp[:, j][i], Yp[:, j][i]] 
			ZK = TokamakField(K, φ[j], params)

			#println("Vector Field = ", ZK, "\n  Directional derivative= ", DK)		

			absError = abs.(ZK - DK)

			if  absError[1] > invError || absError[2] > invError
				
				invError = max(absError[1], absError[2])
			end
				
			i = i + 1

		end

		j = j + 1

	end
	
	return invError

end

########################################################
############## Save components of torus ################
########################################################

function saveFileTorus(θ::Vector{Float64}, φ::Vector{Float64},
			Xp::Matrix{Float64}, Yp::Matrix{Float64})

	lenθ, lenφ = length(θ), length(φ) 
	#println("len theta= ", lenθ)
	nameFile = "ParamTorus$ε.txt"
	freqomega = omega/2π
	f = open(nameFile, "w")

		write(f,"2"*"\n"*"2"*"\n"*"20"*"\n")
		write(f,"1.0e-13"*"\n"*"1.0e-15"*"\n")
		write(f,"$(freqomega)"*" "*"$alpha"*"\n")
		write(f,"$(lenθ - 1)"*" "*"$(lenφ - 1)"*"\n")
		write(f,"$errorInv \n")

	i = 1
	while i < lenθ
		
		j = 1
		while j < lenφ 

			#auxX = Xp[:, j][i]
			#auxY = Yp[:, j][i]

			auxX = Xp[i, :][j]
			auxY = Yp[i, :][j]
			write(f, "$(i - 1)"*" "*"$(j - 1)"*" "*"$auxX"*" "*"$auxY"*"\n")
			j += 1		
		end
		i += 1
	end
	close(f)
end

########################################################
########################################################
########################################################
theta, varphi, XP, YP = torusTwoD(2^9 + 1, 2^8 + 1)
#println("X= ", XP)
#println("Y= ", YP)
#theta, varphi, XP, YP = torusTwoD(2^3 + 1, 2^2 + 1)
#theta, varphi, XP, YP = torusTwoD(5, 4)
println("Complete two-dimensional torus")

global errorInv = invarianceError(theta, varphi, XP, YP)

println("Complete invariance error", errorInv)
#println(XP)
saveFileTorus(theta, varphi, XP, YP)

println("Complete to write file with invariance Torus")
#=
#DXPθ, DXPφ = derivativesData(theta, varphi, XP)
#DYPθ, DYPφ = derivativesData(theta, varphi, YP)


########################################################
################### Theta Components ###################
########################################################

anim = @animate for i = 1:1:length(varphi)
	
	plot(theta, YP[:, i], ylim = (-0.4, 0.5) , xlabel = L"\theta" ,label = L"Y_{p}(\theta, \varphi)")
	plot!(theta, DYPθ[:, i], label = L"\partial_{\theta} Y_{p}(\theta, \varphi)")
end

gif(anim, "ypThetaComponents.gif", fps = 30)

anim = @animate for i = 1:1:length(varphi)
	
	plot(theta, XP[:, i], ylim = (-1.0, 2.0), xlabel = L"\theta" ,label = L"X_{p}(\theta, \varphi)")
	plot!(theta, DXPθ[:, i], label = L"\partial_{\theta} X_{p}(\theta, \varphi)")
end

gif(anim, "xpThetaComponents.gif", fps = 30)


########################################################
################## Varphi Components ###################
########################################################

anim = @animate for i = 1:1:length(theta)
	
	plot(varphi, YP[i, :], ylim = (-0.25, 0.6), xlabel = L"\varphi" ,label = L"Y_{p}(\theta, \varphi)")
	plot!(varphi, DYPφ[i, :], label = L"\partial_{\varphi} Y_{p}(\theta, \varphi)")
end

gif(anim, "ypVarphiComponents.gif", fps = 30)

anim = @animate for i = 1:1:length(theta)
	
	plot(varphi, XP[i, :], ylim = (-1.0, 0.5 ), xlabel = L"\varphi", label = L"X_{p}(\theta, \varphi)")
	plot!(varphi, DXPφ[i, :], label = L"\partial_{\varphi} X_{p}(\theta, \varphi)")
end

gif(anim, "xpVarphiComponents.gif", fps = 30)


=#
