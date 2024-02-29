using Pkg; Pkg.add("DynamicalSystems")
using Pkg; Pkg.add("StaticArrays")
using Pkg; Pkg.add("OrdinaryDiffEq")
using Pkg; Pkg.add("CairoMakie")
using DynamicalSystems
using StaticArrays
using OrdinaryDiffEq
using CairoMakie


function toda(u,p,t)
    x,y = u
    A, d, f = p
    #dx = y
    #dy = - d * y + x - x^3 + A * sin(2*pi*f * t)
    dx = y 
    dy = - d * y + exp(-x) -1  + A * sin(2*pi * f * t + pi)
    return SVector(dx, dy)

end

function toda_with_gaussian(u,p,t)
    x,y = u
    A,A0, d, f, r = p
    t0, sigma = 300/(2.5/(2*pi)) , 50
    #dx = y
    #dy = - d * y + x - x^3 + A * sin(2*pi*f * t)
    dx = y 
    dy = - d * y + exp(-x) - 1 + Gaussian_puls(t,t0,A,A0,sigma,r) * sin(2*pi * f * t)
    return SVector(dx, dy)

end 

function Gaussian_puls(t,t0,A,A0,sigma,r)

    
    return @. (A0*exp.(-(t-t0)^2/(2*sigma^2)) + A) * 1/2 * (sign(t-t0)+1) +  1/2 * (sign(-t + t0) + 1) #(A0*exp.(-(t-t0)^2/(2*(sigma*r)^2)) + A) 
    #return @. (sign(t-t0)+1) #+ (sign(-t+t0)+1) 
    #return @. (A0/ ( exp((-t+t0)/sigma) + 1 ) + A)
end

#t = range(0,1400,10000)
#scatter(t,Gaussian_puls(t,1000,1.9,0.5,500,.000001))

u01 = [-1.85,0.]
p01 = [1.85, .75 , .05, 2.5/(2*pi),0.000001]

u02 = [1.5,-2.]
p02 = [2., .05, 2.5/(2*pi)]


diffeq = (alg = Vern9(), abstol = 1e-6, reltol = 1e-6)
toda_system_gauss = ContinuousDynamicalSystem(toda_with_gaussian,u01,p01; diffeq)
toda_system = ContinuousDynamicalSystem(toda,u02,p02; diffeq)
#poincaremap_toda = DeterministicIteratedMap(poincare_toda,u02,p02)
T = 1/(p02[3])
#strob_toda = StroboscopicMap(toda_system,T)
strob_toda = StroboscopicMap(T, toda, u02, p02;diffeq,t0 = 0)

iterationen = 2000
total_time = iterationen / (2.5/(2*pi))
X1, t = trajectory(strob_toda, iterationen)
X2, t = trajectory(toda_system, total_time)

f = Figure()
ax = Axis(f[1,1], title= "A=1.85, d=0.05, f=2.5/2pi",ylabel="dx/dt",xlabel="x")
a,b= 10000,length(t)
scatter!(X2[a:b,1],X2[a:b,2],label= "trajectory")
scatter!(X1[length(X1)-30:length(X1),1],X1[length(X1)-30:length(X1),2],label = "poincare map")
#scatter!(last(X1[1:34,1]),last(X1[1:34,2]),label = "last")
axislegend()

f
#save("Attractorchange.pdf", f)


xg = yg = range(-2.,2; length = 50)

mapper = AttractorsViaRecurrences(strob_toda, (xg, yg); sparse = false)

basins, attractors = basins_of_attraction(mapper; show_progress = true)


f_heat, ax_heat, hm = heatmap(xg, yg, basins)
#scatter!(u01[1],u01[2],color = "red")
scatter!(u02[1],u02[2],color = "green",label="initial cond")
scatter!(X1[length(X1)-30:length(X1),1],X1[length(X1)-30:length(X1),2],color = "red",label = "poincare map")

axislegend()
Colorbar(f_heat[1,2],hm)
f_heat
#save( "basins_toda_poincare_A_185.pdf",f_heat)