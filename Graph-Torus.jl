using Plots, LaTeXStrings
using Dierckx

include("Auxiliary-Functions.jl")
include("Torus.jl")
titleFile = "Torus0.35597881128974557modulus.txt"
titleFileOne = "Torus0.35597881128974557.txt" 
titleFileTwo = "uniformData.txt"
titleFileThree = "nonUniformData.txt"


xm, omegam  = dataGraph(titleFile) # This a function of Auxiliary-Functions 
x, omega  = dataGraph(titleFileOne) # This a function of Auxiliary-Functions 
xu, epsilonu  = dataGraph(titleFileTwo) # This a function of Auxiliary-Functions 
xnu, epsilonnu  = dataGraph(titleFileThree) # This a function of Auxiliary-Function
xtheta = convertTheta(x, omega)

################################################################
xmaux =  xm
bubbleSort!(xmaux)
torus = Spline1D(xmaux[:, 1], xmaux[:, 2], periodic = true)
a =  collect(range(0.0 , stop = 2π ,length = 50))

################################################################
xtaux =  xtheta
xtaux[:, 1] = xtaux[:, 1].%2pi
bubbleSort!(xtaux)

xp = Spline1D(xtaux[:, 1], xtaux[:, 2], periodic = true)
yp = Spline1D(xtaux[:, 1], xtaux[:, 3], periodic = true)
################################################################
########################### Plot all ###########################
################################################################

plot(size = (800, 600))

p1 = scatter(xtheta[:, 1].%2π, xtheta[:, 2], markersize = 2.0, xlabel = L"\theta", ylabel = L"x_{p}(\theta, \varphi = 0)", legend=false) 
plot!(a, xp.(a))

p2 = scatter(xtheta[:, 1].%2π, xtheta[:, 3], markersize = 2.0, markerstrokestyle=:solid, xlabel = L"\theta", ylabel = L"y_{p}(\theta, \varphi = 0)", legend=false)
plot!(a, yp.(a))

p3 = scatter(xm[:,1], xm[:,2], color=:blue, markersize = 2.0 ,label="Integracion", xlabel = L"x(\theta) = \theta + x_{p}(\theta)", ylabel = L"y(\theta)", legendfontsize = 4, legend =:bottomleft)
scatter!(xu[:,1], xu[:,2], markersize = 4.0, markerstrokestyle=:solid, color=:green, label = "malla uniforme")
scatter!(xnu[:,1], xnu[:,2], markersize = 3.0, markerstrokestyle=:solid, color=:red, label= "malla no uniforme")

# p4 = scatter(xm[:, 1], xm[:, 2], markersize = 2.0, markerstrokestyle=:solid, legend = false, xlabel = L"x(\theta)", ylabel = L"y(\theta)")
#plot!(a, torus.(a))


p4 = scatter(xtheta[:, 2] .+ xtheta[:, 1].%2π , xtheta[:, 3], markersize = 2.0, markerstrokestyle=:solid, legend = false, xlabel = L"x(\theta)", ylabel = L"y(\theta)")
plot!(a, torus.(a))

plot(p1, p2, p3, p4, layout = (2, 2))

savefig("All.pdf")

################################################################
################## Plot two-dimensional torus ##################
################################################################

theta, varphi, XP, YP = torusTwoD(128, 256)

plot(size = (800, 600))

p1 = scatter(theta, YP[:, 1], markersize = 2.0, xlabel = L"\theta", ylabel = L"y_{p}    (\theta, \varphi)", legend = false)
scatter!(theta, YP[:, end], markersize = 2.0)

p2 = scatter(varphi, YP[1, :], markersize = 2.0, xlabel = L"\varphi", ylabel = L"y_{    p}(\theta, \varphi)", legend = false)
scatter!(varphi, YP[end, :], markersize = 2.0)

p3 = scatter(theta, XP[:, 1], markersize = 2.0, xlabel = L"\theta", ylabel = L"x_{p}    (\theta, \varphi)", legend =  false)
scatter!(theta, XP[:, end], markersize = 2.0)

p4 = scatter(varphi, XP[1, :], markersize = 1.0, xlabel = L"\varphi", ylabel = L"x_{    p}(\theta, \varphi)", legend =  false)
scatter!(varphi, XP[end, :], markersize = 2.0)

plot(p1, p2, p3, p4, layout =(2, 2))
savefig("torusComponents.pdf")

println("Complete")

########################################################
########################################################


anim = @animate for i = 1:2:256

	plot(theta, YP[:, i], ylim = (0.25, 0.455))
end
gif(anim, "ypComponent.gif", fps = 30)
