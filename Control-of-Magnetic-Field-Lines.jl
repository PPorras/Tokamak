################################################################

using TaylorIntegration
include("Auxiliary-Functions.jl")
using Dierckx

################################################################
########################## Constats ############################
################################################################

#const ε = 0.004 #Perturbation is 0.004
const ε = 0.0 #Perturbation is 0.004
psi_0  = 0.35
w(x) = (2 - x)*(2 - 2*x + x^2)/4
dw(x) = ((-1)*(2 -2*x + x^2) + ( 2 - x )*(-2 + 2*x))/4
w_0 = w(psi_0)
dw_0 = dw(psi_0)
init = 0.35597881128974557
const omega = 2π*w(init) 
#const omega = 3.643074351  ## Rotation number
#const omega = 3.6430704771333495
#const omega = 3.6430690936048755 ##Este es el bueno para 0.004
num = 64

################################################################

println("\t Perturbation ", ε, "\n")
println("\t Rotation number ", omega, "\n")
println("\t Initial condition (0.0, 0.35397881128974557) \n")

################################################################
##################### Equations of motion #####################
################################################################

function Field!(dx, x, params, t)
   w_0 = params[1]
   dw_0 = params[2]
   dx[1] = (2 - x[2])*(2 - 2*x[2] + x[2]^2)/4
   dx[2] = ε*(2*sin(2*x[1] - t) + 3*sin(3*x[1] - 2*t)) - (ε^2)*dw_0*(2*cos(2*x[1] - t)/(2*w_0 - 1) + 3*cos(3*x[1] - 2*t)/(3*w_0 - 2))*(4*sin(2*x[1] - t)/(2*w_0 - 1) + 9*sin(3*x[1]- 2*t)/(3*w_0 - 2))
    nothing

end

################################################################
##################### Equations of motion #####################
################################################################

fin = 100_000
params = [w_0, dw_0]			# parameters 
x0 = [0, init]				# initial condition

tv = 0.0:2π:4π			# time for stroboscopic map
xvS = taylorinteg(Field!,  x0, tv, 25, 1.0E-20, params; maxsteps=7_000_000) 
println(xvS)

ω = WBA(xvS[:, 1])

println("End the first run, the numbre is  ",  ω)

tv = 0.0:2π:2π*fin			# time for stroboscopic map
xvS = taylorinteg(Field!,  x0, tv, 25, 1.0E-20, params; maxsteps=7_000_000) 

const rho = WBA(xvS[:, 1])
println("rotation number, omega is = ", omega)
println("frecuece number, omega/2π is = ", omega/2π)
println("frecuence number rho is = ", rho)
println("rotarion number is 2pi*rho = ", rho/2π)
#println("Arturo's calculus = ", omega)
#const omega = 3.6535115247876093

################################################################
####################### Invariant torus ########################
########### The torus is close to 0.35397881128974557 ##########
################################################################

function torus(fin = 100000::Int64)

	### init = 0.35597881128974557
	params = [w_0, dw_0]			# parameters 
	x0 = [0.0, init]			# initial condition
	tv = 0.0:2π:2π*fin			# time for stroboscopic map
	nameFile = "Torus"*"$init"*".txt"
	nameFileTwo = "Torus"*"$init"*"modulus.txt"

	f = open(nameFile, "w")
	write(f,"$(fin + 1.0)"*" $omega"*"\n")
	close(f)

	f = open(nameFileTwo, "w")
	write(f,"$(fin + 1.0)"*" $omega"*"\n")
	close(f)

	xvS = taylorinteg(Field!,  x0, tv, 32, 1.0E-20, params; maxsteps=7_000_000) 
	#println(xvS)	
	writeFile(nameFile, xvS) # Function of Auxiliary-Function
	xTwoPi = mod.(xvS, 2π)
	writeFile(nameFileTwo, xTwoPi) # Function of Auxiliary-Function

end

################################################################
####################### stroboscopic map ####################### 
################################################################

points = 64

torus(10)
println("Completo 10")

println("Stroboscopic map for $points points \n")
torus(points)

################################################################
######################## writing files ######################### 
################################################################
#=
titleFile = "Torus0.35597881128974557modulus.txt"

println("Filter with $num knots for the Spline  \n")

#θ, out = filter(titleFile, num)
xout, omegaa = dataGraph(titleFile)
bubbleSort!(xout)

θ = addTwoPi(xout)
#println(θ, "\n")

x = θ[:, 1]
y = θ[:, 2]

#println(x)

tori = Spline1D(x, y, periodic = true)
paraTorus = zeros(length(θ[:, 1]), 2)
 
paraTorus[:, 1] = θ[:, 1]
paraTorus[:, 2] = tori.(θ[:, 1])

nameFile = "uniformData.txt"
f = open(nameFile, "w")
write(f, "$num"*" "*"0.004"*"\n")
close(f)

writeFile(nameFile, paraTorus)
 
nameFileTwo = "nonUniformData.txt"
f = open(nameFileTwo, "w")
write(f, "$num"*" "*"0.004"*"\n")
close(f)

writeFile(nameFileTwo, θ)
=#
println("Complete \n")


