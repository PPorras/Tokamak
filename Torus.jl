###################################################
############### Parameterization ##################
###################################################

using TaylorIntegration
using Dierckx
using Plots, LaTeXStrings
include("Auxiliary-Functions.jl")

###################################################

titleFile = "Torus0.35597881128974557modulus.txt"
titleFileOne = "Torus0.35597881128974557.txt"

###################################################

const eps = 0.004 #Perturbation
psi_0  = 0.35
w(x) = (2 - x)*(2 - 2*x + x^2)/4
dw(x) = ((-1)*(2 -2*x + x^2) + ( 2 - x )*(-2 + 2*x))/4
w_0 = w(psi_0)
dw_0 = dw(psi_0)
#const omega = 3.6430966211675684 #Def con inf
#const omega = 3.6430704771333495 ### calculado con 1000000 puntos
const omega = 3.643074351 ### calculado por Arturo
const alpha = 2π 
params = [w_0, dw_0]
num = 64

###################################################
#################### Field ########################
###################################################

function Field!(dx, x, params, t)
 	w_0 = params[1]
 	dw_0 = params[2]
	dx[1] = (2 - x[2])*(2 - 2*x[2] + x[2]^2)/4
 	dx[2] = eps*(2*sin(2*x[1] - t) + 3*sin(3*x[1] - 2*t)) - (eps^2)*dw_0*(2*cos(2*x[1] - t)/(2*w_0 - 1) + 3*cos(3*x[1] - 2*t)/(3*w_0 - 2))*(4*sin(2*x[1] - t)/(2*w_0 - 1) + 9*sin(3*x[1]- 2*t)/(3*w_0 - 2))
 	nothing

 end
###################################################
#################### Filter #######################
###################################################

θ, out = filter(titleFile, num)

x = out[:, 1]
y = out[:, 2]

torus = Spline1D(x, y, periodic = true)

###################################################
############## Two dimensional Torus ##############
###################################################

x, omegaa = dataGraph(titleFileOne)
xtaux = convertTheta(x, omega)
xtaux[:, 1] = mod2pi.(xtaux[:, 1])
bubbleSort!(xtaux)

xθ = closest(xtaux, θ)

s =  collect(range(0.0 , stop = 2π ,length = 64))

xp = Spline1D(xθ[:, 1], xθ[:, 2], k=3, periodic = true)
yp = Spline1D(xθ[:, 1], xθ[:, 3], k=3, periodic = true)


###################################################
#### Construction of the two-dimensional torus ####
###################################################

#########
#SP(x) = surfacePhi(x, 0)
###################################################

###################################################
function torusTwoD(lenθ, lenφ)
	
	X, Y = zeros(lenθ, lenφ), zeros(lenθ, lenφ)
	θ = collect(range(0.0 , stop = 2π, length = lenθ))
	φ = collect(range(0.0 , stop = 2π, length = lenφ))
	X[:, 1], Y[:, 1]= xp.(θ), yp.(θ)

	j = 2

	while j <= lenφ 
		
		i =  1

		while i <= lenθ	

			aux = mod2pi(θ[i] - (omega*φ[j])/2pi) 

			xinit = θ[i] - (omega*φ[j])/2pi + xp(aux)
			
			yinit = yp(aux)

			tv = 0.0:0.5*φ[j]:φ[j] 
			x0 = [xinit, yinit]
			xvs = taylorinteg(Field!,  x0, tv, 32, 1.0E-20,
			params; maxsteps=7_000_000)

			X[(j - 1)*lenθ + i] = xvs[end, :][1] - θ[i]
			Y[(j - 1)*lenθ + i] = xvs[end, :][2] 
		
			i = i + 1
		end

		j = j + 1

	end

	return θ, φ, X, Y

end

#######################################################
###### Derivatives of the components of the torus #####
#######################################################

function derivativesData(θ, φ, data)

	lenθ = length(data[:, 1])
	lenφ = length(data[1, :])
	
	derθ = zeros(lenθ, lenφ)
	derφ = zeros(lenθ, lenφ)

	i = 1
	
	while i <= lenφ

		SD = Spline1D(θ, data[:, i], k=3, periodic = true)
		derθ[:, i] = Dierckx.derivative(SD, θ)

		i = i + 1

	end

	i = 1
	
	while i <= lenθ

		SD = Spline1D(φ, data[i, :], k=3, periodic = true)
		derφ[i, :] = Dierckx.derivative(SD, φ)

		i = i + 1

	end
	
return derθ, derφ

end

#######################################################
################### Invariance error ##################
#######################################################


function invarianceError(θ::Vector{Float64}, φ::Vector{Float64},
 Xp::Matrix{Float64}, Yp::Matrix{Float64})

	lenθ, lenφ = length(θ), length(φ)

	j = 1

	while j <= lenφ

		i = 1

		while i <= lenθ
			
			KXPθ = Spline1D(θ, Xp[:, j])  
			KXPφ = Spline1D(φ, Xp[i, :])  
			KYPθ = Spline1D(θ, Yp[:, j])  
			KYPφ = Spline1D(φ, Yp[i, :])  

			DKXPθ = Dierckx.derivative(KXPθ, θ[i])
			DKXPφ = Dierckx.derivative(KXPφ, φ[j])
			DKYPθ = Dierckx.derivative(KYPθ, θ[i])
			DKYPφ = Dierckx.derivative(KYPφ, φ[j])
		
			
			xComp =	omega + DKXPθ*omega + DKXPφ*alpha		
			yComp =	DKYPθ*omega + DKYPφ*alpha		
			
			DK = [xComp, yComp]

			K = [θ[i] + Xp[:, j][i], Yp[:, j][i]] 
			ZK = TokamakField(K, φ[j], params)
			
			aux1[i, :] = [Yp[:, j][i],  DKYPθ, Xp[:, j][i], DKXPθ] 

			#println("Invariance error = ", abs.(ZK .- DK))
	
			i = i + 1

		end

		j = j + 1

	end

	println(lenθ)
end

########################################################
################# Derivative for data ##################
########################################################


########################################################
########################################################
########################################################
theta, varphi, XP, YP = torusTwoD(256, 250)
println("Tamaño XP = ", length(XP[:, 1]), length(XP[1, :]))
println("Tamaño YP = ", length(YP[:, 1]), length(YP[1, :]))
println("Complete two-dimensional torus")

DXPθ, DXPφ = derivativesData(theta, varphi, XP)
DYPθ, DYPφ = derivativesData(theta, varphi, YP)

########################################################
################### Theta Components ###################
########################################################

anim = @animate for i = 1:1:length(varphi)
	
	plot(theta, YP[:, i], label = L"Y_{p}(\theta, \varphi)")
	plot!(theta, DYPθ[:, i], label = L"\partial_{\theta} Y_{p}(\theta, \varphi)")
end

gif(anim, "ypThetaComponents.gif", fps = 30)

anim = @animate for i = 1:1:length(varphi)
	
	plot(theta, XP[:, i], label = L"X_{p}(\theta, \varphi)")
	plot!(theta, DXPθ[:, i], label = L"\partial_{\theta} X_{p}(\theta, \varphi)")
end

gif(anim, "xpThetaComponents.gif", fps = 30)


########################################################
################## Varphi Components ###################
########################################################

anim = @animate for i = 1:1:length(theta)
	
	plot(varphi, YP[i, :], label = L"Y_{p}(\theta, \varphi)")
	plot!(varphi, DYPφ[i, :], label = L"\partial_{\varphi} Y_{p}(\theta, \varphi)")
end

gif(anim, "ypVarphiComponents.gif", fps = 30)

anim = @animate for i = 1:1:length(theta)
	
	plot(varphi, XP[i, :], label = L"X_{p}(\theta, \varphi)")
	plot!(varphi, DXPφ[i, :], label = L"\partial_{\varphi} X_{p}(\theta, \varphi)")
end

gif(anim, "xpVarphiComponents.gif", fps = 30)


########################################################
########################################################
########################################################
#invarianceError(theta, varphi, XP, YP)
println("Calculate of the Invariance Error ")
#=
########################################################
#######################################################
#println(length(theta), " ", length(varphi)," " ,length(psi))
f = Spline2D(theta, varphi, psi) 

xs = collect(0.0:0.1:2π)
as = collect(0.0:0.1:2π)
x_grid = [x for x = xs for y = as]
a_grid = [y for x = xs for y = as]


plot(size = (800, 600))
plot!(x_grid, a_grid, f.(x_grid, a_grid), st = :surface, xlabel = L"\theta", ylabel = L"\varphi", zlabel = L"\Psi", camera = (45, 45))
savefig("Parameterized-torus.pdf")

########################################################
#######################################################

=#
