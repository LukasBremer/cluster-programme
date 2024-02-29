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
    dy = - d * y + exp(-x) -1  + A * sin(2*pi * f * t)
    return SVector(dx, dy)

end


function Gaussian_puls(t,t0,A,A0,sigma)
    #return @. sin(t-t0)
    return @. A*exp.(-(t-t0)^2/(2*sigma^2)) + A0
    
end

function toda_with_gaussian(u,p,t)
    x,y = u
    A, d, f = p
    t0, A0, sigma = 2000, .05, 20 
    #dx = y
    #dy = - d * y + x - x^3 + A * sin(2*pi*f * t)
    dx = y 
    dy = - d * y + exp(-x) -1  + Gaussian_puls(t,t0,A,A0,sigma) * sin(2*pi * f * t)
    return SVector(dx, dy)

end 


u01 = [-2.1,0.6]
u02 = [-2.,0.]
p0 = [1.85,.05,2.5/(2*pi)]


diffeq = (alg = Vern9(), abstol = 1e-9, reltol = 1e-9)
toda_system1 = ContinuousDynamicalSystem(toda_with_gaussian,u01,p0; diffeq)
toda_system2 = ContinuousDynamicalSystem(toda,u02,p0; diffeq)
total_time = 5000

X1, t = trajectory(toda_system1, total_time)
X2, t = trajectory(toda_system2, total_time)

f = Figure()
ax = Axis(f[1,1], title= "A=1.85, d=0.05, f=2.5/2pi",ylabel="dx/dt",xlabel="x")
scatter!(X2[total_time-1000:total_time,1],X2[total_time-1000:total_time,2],label= string(u02))
scatter!(X1[total_time-1000:total_time,1],X1[total_time-1000:total_time,2],label = string(u01))
axislegend()
f
save("Coex_attr_A19.pdf", f)


# define a state space grid to compute the basins on:
xg = yg = range(-5, 5; length = 201)
# find attractors using recurrences in state space:
mapper = AttractorsViaRecurrences(toda_system1, (xg, yg); sparse = false, force_non_adaptive=false, horizon_limit = 1e6,mx_chk_safety = Int(1e6))
# compute the full basins of attraction:
basins, attractors = basins_of_attraction(mapper; show_progress = false)


f_heat, ax_heat, hm = heatmap(xg, yg, basins)
Colorbar(f_heat[1,2],hm)
f_heat
#save( "basins_toda.pdf",f_heat)



