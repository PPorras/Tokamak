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

xp = Spline1D(xθ[:, 1], xθ[:, 2], k=3)
yp = Spline1D(xθ[:, 1], xθ[:, 3], k=3)


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
######### Graph of the components of the torus ########
#######################################################

theta, varphi, XP, YP = torusTwoD(32, 36)
println("Complete two-dimensional torus")

#######################################################
################### Invariance error ##################
#######################################################


"""
invarianceError(θ::Vector{Float64}, φ::Vector{Float64}, Xp::Matrix{Float64}, Yp::Matrix{Float64})

Compute the invariance error 

### Input
- `θ` -- list of 1D vector; 
- `φ` -- list of 1D vector; 

### Output

List of vectors containing the 2D coordinates of the corner points of the
convex hull.

# Examples
       ```jldoctest
       julia> cube(2)
       8
       ```
### Notes
### Algorithm
"""

function invarianceError(θ::Vector{Float64}, φ::Vector{Float64},
 Xp::Matrix{Float64}, Yp::Matrix{Float64})

	lenθ, lenφ = length(θ), length(φ)
	i, j = 1, 1
	
	while i <= lenφ

		while j <= lenθ
			
			KXPθ = Spline1D(θ, Xp[:, i])  
			KXPφ = Spline1D(φ, Xp[j, :])  
			KYPθ = Spline1D(θ, Yp[:, i])  
			KYPφ = Spline1D(φ, Yp[j, :])  

			DKXPθ = Dierckx.derivative(KXPθ, θ[j])
			DKXPφ = Dierckx.derivative(KXPφ, φ[i])
			DKYPθ = Dierckx.derivative(KYPθ, θ[j])
			DKYPφ = Dierckx.derivative(KYPφ, φ[i])
			
			xComp =	omega + DKXPθ*omega + DKXPφ*alpha		
			yComp =	DKYPθ*omega + DKYPφ*alpha		
			
			DK = [xComp, yComp]

			K = [θ[j] + Xp[:, i][j], Yp[:, i][j]] 
			ZK = TokamakField(K, φ[i], params)
		
			println("Invariance error = ", abs.(ZK .- DK))
	
			j = j + 1

		end 
		
		i = i + 1

	end

	println(lenθ)
end

invarianceError(theta, varphi, XP, YP)

#=
println(grato)

KXP = Spline2D(theta, varphi, XP) 
KYP = Spline2D(theta, varphi, YP) 
#println(KXP(pi, pi))

point1 = π
point2 = 0.0

KYP0 = Spline1D(theta, YP[:, 1])
KXP0 = Spline1D(theta, YP[:, 1])

println("D_θK_0(π)  = ", Dierckx.derivative(KYP0, point))
println("D_θK_0(π)ω  = ", Dierckx.derivative(KYP0, point)*omega)

initCond = [point + KXP(point, 0.0), KYP(point, 0.0) ]
a = TokamakField(initCond, 0.0, params)
println("Z(K_0(π), 0)  = ", a)


KYP = Spline1D(theta, YP[:, 2])
println("Derivative of the spline in K_1= ", Dierckx.derivative(KYP, π))
 
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
