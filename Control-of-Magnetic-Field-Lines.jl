################################################################

using TaylorIntegration
include("Auxiliary-Functions.jl")
using Dierckx

################################################################
########################## Constats ############################
################################################################

const ε = 0.004 #Perturbation is 0.004
psi_0  = 0.35
w(x) = (2 - x)*(2 - 2*x + x^2)/4
dw(x) = ((-1)*(2 -2*x + x^2) + ( 2 - x )*(-2 + 2*x))/4
w_0 = w(psi_0)
dw_0 = dw(psi_0)
init = 0.35597881128974557 ### Para epsilon 0.0 campo de Ugo
#### init = 0.3566248878338341665
#init = 0.35
#const omega = 2π*w(init) 
#const omega = 3.643074351  ## Rotation number
#const omega = 3.6430704771333495
#const omega = 3.6430690936048755 ## Numéro rotación 0.004 y campo de Cristel
#const omega = 3.6499335957005423950 # Numéro de rotación para 0.0 y campo de Ugo
const omega = 3.64306909360487549777562045188055480
num = 64

################################################################


################################################################
##################### Equations of motion #####################
################################################################

################### Chandre's Vector Field  ####################

function Field!(dx, x, params, t)
   w_0 = params[1]
   dw_0 = params[2]
   dx[1] = (2 - x[2])*(2 - 2*x[2] + x[2]^2)/4
   dx[2] = ε*(2*sin(2*x[1] - t) + 3*sin(3*x[1] - 2*t)) - (ε^2)*dw_0*(2*cos(2*x[1] - t)/(2*w_0 - 1) + 3*cos(3*x[1] - 2*t)/(3*w_0 - 2))*(4*sin(2*x[1] - t)/(2*w_0 - 1) + 9*sin(3*x[1]- 2*t)/(3*w_0 - 2))
    nothing

end

################### Ugo's Vector Field  ####################
#=
function Field!(dx, x, params, t)
   dx[1] = (2 - x[2])*(2 - 2*x[2] + x[2]^2)/4
   dx[2] = ε*(2*sin(2*x[1] - t) + 3*sin(3*x[1] - 2*t)) + (3/2)*(ε^2)*(27*sin(6*x[1] - 4*t) + 6*(5*sin(5*x[1] - 3*t) - sin(x[1] - t)) - 8*sin(4*x[1] - 2*t))
    nothing

end
=#
################################################################
##################### Equations of motion #####################
################################################################
#=
fin = 1_000_000
params = [w_0, dw_0]			# parameters 
# params = [0, 0]			# parameters 
x0 = [0.0, init]				# initial condition

tv = 0.0:2π:10π			# time for stroboscopic map
xvS = taylorinteg(Field!,  x0, tv, 25, 1.0E-20, params; maxsteps=7_000_000) 
ω = WBA(xvS[:, 1])
println("End the first run, the numbre is  ",  ω)


#tv = 0.0:2π:2000π			# time for stroboscopic map
#xvS = taylorinteg(Field!,  x0, tv, 25, 1.0E-20, params; maxsteps=7_000_000) 
#ω = WBA(xvS[:, 1])
#println("End the first run, the numbre is  ",  ω)

tv = 0.0:2π:2π*fin			# time for stroboscopic map
xvS = taylorinteg(Field!,  x0, tv, 25, 1.0E-20, params; maxsteps=7_000_000) 
ω = WBA(xvS[:, 1])
println("Rotation number, WBA = ", ω)
#const omega = WBA(xvS[:, 1])
#const omega = 3.6535115247876093
=#
################################################################

println("\t Perturbation ", ε, "\n")
println("\t Rotation number ", omega, "\n")
println("\t Frequecie ", (omega)/(2pi), "\n")
println("\t Initial condition", (0.0, init), "\n")
####################### Invariant torus ########################
########### The torus is close to 0.35397881128974557 ##########
################################################################

function torus(fin = 100000::Int64)

	### init = 0.35597881128974557
	params = [w_0, dw_0]			# parameters 
	x0 = [0.0, init]			# initial condition
	tv = 0.0:2π:2π*fin			# time for stroboscopic map
	nameFile = "Torus"*"$init"*"$ε"*".txt"
	nameFileTwo = "Torus"*"$init"*"$ε"*"modulus.txt"

	f = open(nameFile, "w")
	write(f,"$(fin + 1.0)"*" $omega"*" $ε"*"\n")
	close(f)

	f = open(nameFileTwo, "w")
	write(f,"$(fin + 1.0)"*" $omega"*"\n")
	close(f)

	xvS = taylorinteg(Field!,  x0, tv, 5, 1.0E-20, params; maxsteps=7_000_000) 
	#println(xvS)	
	writeFile(nameFile, xvS) # Function of Auxiliary-Function
	xTwoPi = mod.(xvS, 2π)
	writeFile(nameFileTwo, xTwoPi) # Function of Auxiliary-Function

end

################################################################
####################### stroboscopic map ####################### 
################################################################

points = 100

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


### 3 + 1/(1 + 1/(1 + 1/(1 + 1/(4 + 1/(23 + 1/(1 + 1/(6 + 1/(51 + 1/(1 + 1/(1 + 1/(1 + 1/(3 + 1/(1 + 1/(1 + 1/(6 + 1/(3 + 1/...))))))))))))))))
