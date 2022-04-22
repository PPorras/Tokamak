################################################################

using TaylorIntegration
include("Auxiliary-Functions.jl")
using Dierckx
################################################################
########################## Constats ############################
################################################################

const eps = 0.004 #Perturbation
psi_0  = 0.35
w(x) = (2 - x)*(2 - 2*x + x^2)/4
dw(x) = ((-1)*(2 -2*x + x^2) + ( 2 - x )*(-2 + 2*x))/4
w_0 = w(psi_0)
dw_0 = dw(psi_0)
const omega = 3.643074351  ## Rotation number
#const omega 3.6430704771333495
num = 300

################################################################

println("\t Perturbation ", eps, "\n")
println("\t Rotation number ", omega, "\n")
println("\t Initial condition (0.0, 0.35397881128974557) \n")

################################################################
##################### Equations of motion #####################
################################################################

function Field!(dx, x, params, t)
   w_0 = params[1]
   dw_0 = params[2]
   dx[1] = (2 - x[2])*(2 - 2*x[2] + x[2]^2)/4
   dx[2] = eps*(2*sin(2*x[1] - t) + 3*sin(3*x[1] - 2*t)) - (eps^2)*dw_0*(2*cos(2*x[1] - t)/(2*w_0 - 1) + 3*cos(3*x[1] - 2*t)/(3*w_0 - 2))*(4*sin(2*x[1] - t)/(2*w_0 - 1) + 9*sin(3*x[1]- 2*t)/(3*w_0 - 2))
    nothing

end

################################################################
##################### Equations of motion #####################
################################################################
#=
fin = 1_000_000
init = 0.35597881128974557
params = [w_0, dw_0]			# parameters 
x0 = [0, init]				# initial condition

tv = 0.0:2π:4π			# time for stroboscopic map
xvS = taylorinteg(Field!,  x0, tv, 25, 1.0E-20, params; maxsteps=7_000_000) 
ω = rotation_number(xvS)
ω2 = rotationNumber(xvS)

tv = 0.0:2π:2π*fin			# time for stroboscopic map
xvS = taylorinteg(Field!,  x0, tv, 25, 1.0E-20, params; maxsteps=7_000_000) 

ω = rotation_number(xvS)
ω2 = rotationNumber(xvS)


=#
################################################################
####################### Invariant torus ########################
########### The torus is close to 0.35397881128974557 ##########
################################################################

function torus(fin = 100000::Int64)

	init = 0.35597881128974557
	params = [w_0, dw_0]			# parameters 
	x0 = [0.0, init]				# initial condition
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
	
	
	writeFile(nameFile, xvS) # Function of Auxiliary-Function

	xTwoPi = mod.(xvS, 2π)
	writeFile(nameFileTwo, xTwoPi) # Function of Auxiliary-Function


end

################################################################
####################### stroboscopic map ####################### 
################################################################


torus(10)
println("Completo 10")

println("Stroboscopic map for 300 points \n")
torus(1000)

################################################################
######################## writing files ######################### 
################################################################

titleFile = "Torus0.35597881128974557modulus.txt"

θ, out = filter(titleFile, num)
x = out[:, 1]
y = out[:, 2]

tori = Spline1D(x, y, periodic = true)
paraTorus = zeros(length(θ), 2)
 
paraTorus[:, 1] = θ
paraTorus[:, 2] = tori.(θ)

nameFile = "uniformData.txt"
f = open(nameFile, "w")
write(f, "$num"*" "*"0.004"*"\n")
close(f)

writeFile(nameFile, paraTorus)
 
nameFileTwo = "nonUniformData.txt"
f = open(nameFileTwo, "w")
write(f, "$num"*" "*"0.004"*"\n")
close(f)

writeFile(nameFileTwo, out)

println("Complete \n")


