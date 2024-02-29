using DynamicalSystems
#using CairoMakie
using OrdinaryDiffEq
using Plots



function duffing(u,p,t)
    x,y = u
    A,d,f = p
    dx = y
    dy = x - x^3 - d * y + A * sin(2*pi * f * t)
    return SVector(dx,dy)
end

function tipped_duffing_linear_phase(u,p,t)
    x, y = u
    A, d, f, A0, phi, t0 = p
    dx = y
    dy = x - x^3 - d * y + A * sin(2*pi * f * t + phi)
    return SVector(dx,dy)
end

function toda(u,p,t) 
    x, y = u
    A, d, f = p
    dx = y
    dy = - d * y + exp(-x) -1  + A * sin(2*pi * f * t + pi)
    return SVector(dx, dy)
end

function tipped_toda(u,p,t)
    x, y = u
    A, d, f, A0, phi = p
    dx = y
    dy = - d * y + exp(-x) -1  + A * sin(2*pi * f * t + Gaussian_puls(t,0,phi,1/f))
    return SVector(dx, dy)
end
 
function tipped_toda_frozen_phase(u,p,t)
    x, y = u
    A, d, f, A0, phi, t0 = p
    dx = y
    dy = - d * y + exp(-x) -1  + A * sin(2*pi * f * t + Gaussian_puls(t0,0,phi,1/f))
    return SVector(dx, dy)
end

function tipped_toda_linear_phase(u,p,t)
    x, y = u
    A, d, f,A0, phi, t0 = p
    dx = y
    dy = - d * y + exp(-x) - 1  + A * sin(2*pi * f * t + phi)
    return SVector(dx, dy)
end

function Gaussian_puls(t,A,A0,T)

    t0, sigma, r = 100*T, 20, .001

    return @. (A0*exp.(-(t-t0)^2/(2*sigma^2)) + A) * 1/2 * (sign(t-t0)+1) +  1/2 * (sign(-t + t0) + 1) * (A0*exp.(-(t-t0)^2/(2*(sigma*r)^2)) + A) 
    #return @. (sign(t-t0)+1) #+ (sign(-t+t0)+1) 
    #return @. (A0/ ( exp((-t+t0)/sigma) + 1 ) + A)
end

function current_state(x,y,p,t)
    T = 1/p[3]
    index = Integer(floor(t/T))

    return x[index],y[index],index
end

function new_basins(basins,attractors)
    px = Int(sqrt(length(basins)))
    merge!(attractors,Dict(-1=> StateSpaceSet([SVector(0.,0.),SVector(0.,0.),SVector(0.,0.),SVector(0.,0.),SVector(0.,0.),SVector(0.,0.),SVector(0.,0.),SVector(0.,0.)])))
    new_b = zeros(px,px)
    for i in 1:px
        for j in 1:px
            new_b[i,j] = length(attractors[basins[i,j]])
        end
    end
    return new_b
end 

u02 = [0.,0.]
u01 = [0.,0.]
p01 = [1.85, .05, 2.5/(2*pi)]
p02 = [1.85, .05, 2.5/(2*pi),1,pi,0]

diffeq = (alg = Tsit5(), adaptive = false, dt = 1/p02[3]/10)

toda_system = ContinuousDynamicalSystem(toda,u02,p02; diffeq)
tipped_toda_system = ContinuousDynamicalSystem(tipped_toda,u02,p02; diffeq)

T = 1/(p02[3])
strob_toda = StroboscopicMap(toda_system, T)
strob_tipped_toda = StroboscopicMap(tipped_toda_system, T)
iterationen = 1000
total_time = iterationen * T
X1,t = trajectory(strob_tipped_toda,iterationen)
X2,t2 = trajectory(toda_system, total_time)
x,y = X1[:,1],X1[:,2]
x2,y2 = X2[:,1],X2[:,2]


# anim = @animate for i in 1:length(t)-10
#     scatter(X1[i:i+10,1], X1[i:i+10,2], ms=5, lab="", alpha = 1 - .05*i/length(t), 
#         xlim=(-2,5), ylim=(-4, 4))
# end

#gif(anim,"tipped_toda.gif",fps=15)

f = Figure()

ax = Axis(f[1,1], title= "A=1.85, d=0.05, f=2.5/2pi",ylabel="dx/dt",xlabel="x")
#scatter(X1[50:150,1],X1[50:150,2])
scatter(X2[length(t2)-100:length(t2),1],X2[length(t2)-100:length(t2),2], color = "blue")
scatter!(X1[length(t)-100:length(t),1],X1[length(t)-100:length(t),2], color = "red")

f
# Basins of attraction

xg = yg = range(-4,10, length = 100)
mapper_grid = (xg, yg)

mapper = AttractorsViaRecurrences(strob_toda, mapper_grid;
    sparse = true, Ttr = 10, mx_chk_fnd_att = 500,
)

basins, attractors = basins_of_attraction(mapper, mapper_grid; show_progress = true)
attractors

new_b = new_basins(basins,attractors)
heatmap(xg,yg, transpose(basins))
heatmap(xg, yg, transpose(new_b),c = cgrad(:matter, [0.1,.3, 0.8], rev = true, categorical = true))
scatter!(SVector(-.5),SVector(-1.5))

# Attractors via proximity: only find trajectories that are _close_ to
# pre-defined attractors. We found these predefined attractors by incrementally
# increasing the size of the state space box; because there were consistently
# some initial conditions diverging, but there shouldn't be any divergence.

# pre_attractors = Dict(
#     1 => StateSpaceSet([SVector(0.0235595,  -0.878125)]),
#     2 => StateSpaceSet([
#     SVector(-1.18748,   -2.56398),
#     SVector(2.48184,    0.829278),
#     SVector(3.29707,   -1.63142)]),
#     3 => StateSpaceSet([
#         SVector(-1.84898,  -3.45932),
#         SVector( 5.20975,   2.06901),
#         SVector(8.81844,  -0.621984),
#         SVector(6.06727,  -2.99651),#ä#ä###

#     ])
# )

# If closeness is small enough, there are plenty of initial conditions
# that do not converge to any of the above 3 attractors
# hinting that there is a 4th attractor

# closeness = 0.01

# mapper = AttractorsViaProximity(strob_toda, pre_attractors, closeness; mx_chk_lost = 1000)

using Pkg
Pkg.add("JLD")
using JLD

anim = Plots.Animation()
iterationen = 100
ti = range(95*T,140*T,iterationen)
X1,t = trajectory(strob_tipped_toda,300*T)
x,y = X1[:,1],X1[:,2]
data = [x,y]
save("trajectory_tipped_toda.jld","data",data)

for i in 1:iterationen
    p_current_time = [1.85, .05, 2.5/(2*pi), .45, Gaussian_puls(ti[i],0,pi,T), ti[i]]
    tipped_toda_system = ContinuousDynamicalSystem(tipped_toda_frozen_phase,u02,p_current_time; diffeq)
    strob_tipped_toda = StroboscopicMap(tipped_toda_system, T)
    mapper = AttractorsViaRecurrences(strob_tipped_toda, mapper_grid;
    sparse = true, Ttr = 10, mx_chk_fnd_att = 500,
    )
    basins, attractors = basins_of_attraction(mapper, mapper_grid; show_progress = true)
    new_b_temp = new_basins(basins,attractors)
    
    new_b = cat(dims=3,new_b, new_b_temp)
    #p2 = scatter(ti,Gaussian_puls(ti,0,2*pi,T))
    #plot(p1,p2,layout= (2,1))
    #xc,yc,pp = current_state(x,y,p02,time[i])
    #scatter!(SVector(xc),SVector(yc))
end

save("basins_tipped_gaussian.jld","data",new_b)

basins = load("basins.jld")["data"]


anim = @animate for i in 1:iterationen
    new_b = basins[:,:,i]
    heatmap(xg, yg, transpose(new_b),colorbar = false,colormap = cgrad(:nipy_spectral, rev = true, categorical = true))
    #p2 = scatter(ti,Gaussian_puls(ti,0,2*pi,T))
    #plot(p1,p2,layout= (2,1))
    #xc,yc,pp = current_state(x,y,p02,time[i])
    #scatter!(SVector(xc),SVector(yc))
    
end

mp4(anim,"basins_tipped_gaussian.mp4",fps=15)

#scatter!(u02[1],u02[2],color = "green",label="initial cond")
#scatter!(X1[length(X1)-30:length(X1),1],X1[length(X1)-30:length(X1),2],color = "red",label = "poincare map")
#scatter!(0.0235595  ,-0.878125,label="starting position",color= "black")
#axislegend()
#f_heat
#save("toda_basins_3.3over5pi_23.pdf",f_heat)

