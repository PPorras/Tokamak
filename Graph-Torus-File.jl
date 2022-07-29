using Plots, LaTeXStrings
using Dierckx
#using PyPlot
include("Auxiliary-Functions.jl")
titleFile = "ParamTorus0.004.txt"

f = open(titleFile)
headed = parse.(Float64, split(readline(f), " "))
headed = parse.(Float64, split(readline(f), " "))
headed = parse.(Float64, split(readline(f), " "))
headed = parse.(Float64, split(readline(f), " "))
headed = parse.(Float64, split(readline(f), " "))
headed = parse.(Float64, split(readline(f), " "))
n, m = Int64(headed[1]), Int(headed[2])
len = n*m
println(len)

dataFile = zeros(len, 4)
i = 1

while i <= len

	aux = parse.(Float64, split(readline(f), " "))
	dataFile[i, :]  = aux
	global i = i + 1
end

close(f)
xp = ((2π*dataFile[:,1])/n) .+ dataFile[:,3]
yp = ((2π*dataFile[:,2])/m) 
zp = dataFile[:,4]

#println(zp[1:3])
########################################################
#=
global R = 1.0

function xt(s, p, r)
return (R + r*cos(s))*cos(p)
end

function yt(s, p, r)
return (R + r*cos(s))*sin(p)
end

function zt(s, p, r)
return R + r*sin(s)
end

########################################################
xx = xt.(xp, yp, zp)
yy = yt.(xp, yp, zp)
zz = zt.(xp, yp, zp)

#println(xx)
=#
########################################################
########################################################

#plot(size = (800, 600))

plot(size = (800, 600))
#scatter!(xp[513:1024], zp[513:1024])
#scatter!(xp[512*511:512*512], zp[512*511:512*512])
surface(xp,yp,zp)
#surf(xx, yy, zz, camera=(-30,30))
#plot(xx, yy, zz, st = :surface)
savefig("Toro-init-surface.pdf")
########################################################
########################################################
########################################################
