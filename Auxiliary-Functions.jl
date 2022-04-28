#using TaylorIntegration, Plots, LaTeXStrings
#pyplot()
using Dierckx

################################################################
#################### Functions for Torus.jl #################### 
################################################################
################################################################

###################################################
#################### Filter #######################
###################################################

function filter(file, len)

	a = collect(range(0.0 , stop = 2π ,length = len))
	i = 1	
	new = zeros(len, 2)
	
	while i < len
	
		x, y = closest(file, a[i])
		new[i], new[i + len] = x, y
		i = i + 1

	end
	
	new[i], new[i + len] = a[end], new[len + 1]

	return  a, new
end

################################################################
######################## BubbleSort ############################
################################################################

function bubbleSort!(array)
    
	i, len = 1, length(array[:, 1])

	while i <= len

		j = 1
		
		while j <= len - i 

			if array[j] > array[j + 1]
				
				aux = array[j, :]
				array[j, :] = array[j + 1, :]	
				array[j + 1, :] = aux
			end

			j = j + 1

		end
		
		i = i  + 1		
			
	end

end

################################################################
####################### Sorted data ############################
################################################################

################################################################
##################### Rotation number ##########################
################################################################

function rotation_number(input)

	a = input[:, 1]
	l = length(a)
	return (a[end] - a[1])/l

end

################################################################

function rotationNumber(input)

	a = input[:, 1]
	n, len =  2, length(a)
	an = a[len] - a[1]
	aux = (an + 2pi)/len
	println(aux)

	while n <= len

		an = a[n] - a[1]
		approx = (an + 2pi)/n

		if approx < aux

			#println(n)
			aux = approx

		end	
		
		n = n + 1
	end

	return aux
end

################################################################
######################## Vector Field ##########################
################################################################

function TokamakField(x, t, params)

	w_0 = params[1]	
	dw_0 = params[2]

	qdot = (2 - x[2])*(2 - 2*x[2] + x[2]^2)/4
	pdot = ε*(2*sin(2*x[1] - t) + 3*sin(3*x[1] - 2*t)) - (ε^2)*dw_0*(2*cos(2*x[1] - t)/(2*w_0 - 1) + 3*cos(3*x[1] - 2*t)/(3*w_0 - 2))*(4*sin(2*x[1] - t)/(2*w_0 - 1) + 9*sin(3*x[1]- 2*t)/(3*w_0 - 2))

	return [qdot, pdot]

end

################################################################

function surfacePhi(theta, phi)
 
 	aux = 0.35 - eps*(2*cos(2*theta - phi)/(2*w(0.35) - 1) + 3*cos(3*theta - 2*phi)/(3*w(0.35) - 2 ))
 
	return aux
 
end
################################################################
##################### Open File for Plot #######################
################################################################

function writeFile(name, line)
 
 	i = 1
	f = open(name, "a")

 	while i <= length(line[:, 1])
 
 		x, y  = line[i, :][1], line[i, :][2]
 		write(f,"$x"*" "*"$y"*"\n")
 		i = i + 1
 
	end
 
	close(f)
end

################################################################
function writeFileA(name, line)

	i = 1
	f = open(name, "a")

	while i <= length(line[:, 1])

		x, y, z  = line[i, :][1], line[i, :][2], line[i, :][3]
		write(f,"$x"*" "*"$y"*" "*"$z"*"\n")
		i = i + 1

	end

	close(f)
end


################################################################
######################### Module 2 pi ##########################
################################################################

function modulusTwoPi(x, len)
	i = 1
	while i <= len
	
		x[i] = x[i]%2π
		i =  i + 1		

	end

	return x

end

################################################################

function modulusPi(x, omega)
	#aux = (len, 2)
	len = length(x[:, 1])
	i = 1
	while i <= len
	
		x[i] = x[i] - (i - 1)*omega
		i =  i + 1		

	end

	return x
end

################################################################
########### Convert to theta, x(theta) and y(theta) ############
################################################################

function convertTheta(x, omega)

	i, len  = 1, length(x[:, 1])
	aux = zeros(len, 3)

	while i <= len
		
		aux[i] = omega*(i - 1)
		aux[len + i] = x[i] - aux[i]
		aux[2*len + i] = x[len + i]

		i = i + 1
	end
	
	return aux

end
################################################################
##################### Open File for Plot #######################
################################################################

function dataGraph(nameFile)

	f = open(nameFile)
	headed = parse.(Float64, split(readline(f), " "))
	len, omega = Int64(headed[1]), headed[2]
	i = 1
	dataFile = zeros(len, 2)
	
	while i <= len 
	
		aux = parse.(Float64, split(readline(f), " "))

		dataFile[i, :]  = aux
		i = i + 1
	end

	close(f)
	return dataFile, omega
 end

################################################################
###################### Closest points  #########################
################################################################

function closest(file::String, x::Float64)

	f = open(file, "r")

	headed  = parse.(Float64, split(readline(f), " "))
	old  = parse.(Float64, split(readline(f), " "))
#	aux = abs(x - old[1]%2π)
	aux = abs(x - old[1])
	i, len = 1, Int64(headed[1])

	while i < len

		new = parse.(Float64, split(readline(f), " "))

#		next = new[1]%2π
		next = new[1]

		if  abs(x - next) < aux
			old = new

#			aux = abs(x - old[1]%2π)
			aux = abs(x - old[1])

		end
 
		i = i + 1
 
 
	 end
 
	close(f)

#	return old[1]%2π, old[2]
	return old[1], old[2]
end

################################################################

function closest(data, x)

	ld, lx = length(data[:, 1]), length(x)	
	i, new = 2, zeros(lx, 3)
	new[1, :] = data[1, :]
	new[end, :] = data[1, :]
	new[end, :][1] = data[end]
	
	while i < lx

		j = 2 
		dis = abs(x[i] - data[j])
		aux = data[j, :]

		while j < ld

			if  abs(x[i] - data[j]) < dis

				aux = data[j, :]
				dis = abs(x[i] - data[j])

			end
			
			j = j + 1

		end
		
		new[i, :] = aux 
 
		i = i + 1
 	end

	new[lx] = data[ld]

	return new
end

################################################################
################# Parameterization with Splines ################
################################################################

function parameterizationK(titleFile)

	data, omega = dataGraph(titleFile)
	xtaux = convertTheta(data, omega)
	xtaux[:, 1] = mod.(xtaux[:, 1], 2pi)
	bubbleSort!(xtaux)
	xtaux[end, :] = xtaux[1, :]
	l = length(xtaux[:, 1])
	xtaux[l] = 2pi
	x = Spline1D(xtaux[:, 1], xtaux[:, 2], periodic = true)
	y = Spline1D(xtaux[:, 1], xtaux[:, 3], periodic = true)
	
	return x, y	

end

