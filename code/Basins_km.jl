using Pkg
using JLD
using GLMakie
using Colors, ColorSchemes
using LaTeXStrings

function Gaussian_puls(t,A,A0,T)

    t0, sigma, r = 25*T, 10*T, .001

    return @. (-A0*exp.(-(t-t0)^2/(2*sigma^2)) + A) * 1/2 * (sign(t-t0)+1) +  1/2 * (sign(-t + t0) + 1) * (-A0*exp.(-(t-t0)^2/(2*(sigma*r)^2)) + A) 
    #return @. (sign(t-t0)+1) #+ (sign(-t+t0)+1) 
    #return @. (A0/ ( exp((-t+t0)/sigma) + 1 ) + A)
end

function current_state(x,y,T,t)
    index = Integer(floor(t/T))
    
    return x[index],y[index],index
end

function current_state_vector(x,y,T,t)
    iterationen = length(t)
    x_out = zeros(iterationen)
    y_out = zeros(iterationen)
    for i in 1:iterationen
        x_out[i],y_out[i],index = current_state(x,y,T,t[i])
    end
    return x_out,y_out
end

function basin_index(phi)
    index = Integer(floor(phi*100/pi))
    return index
end

function basins_under_puls(basins,time,T)
    basins_puls = zeros(500,500,100)
    for i in 1:100
        phi = Gaussian_puls(time[i],0,pi,T)
        basins_puls[:,:,i] = basins[:,:,basin_index(phi)+1]
    end
    return basins_puls
end

iterationen = 100
R_n,P_s = 10*10^(-6), 90*10^3 #Pa
_alpha,beta,gamma,nu,kappa = 1,2,0.001,197*10^3,4/3
nu_0 = _alpha * beta / (gamma * R_n)
delta_phi = .65*pi
T_t= nu_0/nu
delta_phi = 0.65*pi

ti = range(15*T_t,80*T_t,iterationen)
x,y = load("trajectory_tipped_KM.jld","data")
points = Observable(Point2f[(x[1],y[1])])
point_g = Observable(0.)
time = Observable(15*T_t)
t_plot = range(15*T_t,80*T_t,1000)
x_p,y_p = current_state_vector(x,y,T_t,ti)
x_g = Gaussian_puls(ti,delta_phi,delta_phi,T_t)
basins = load("basins_tipped_gaussian_KM.jld")["data"]
#basins = basins_under_puls(basins_temp,ti,T)
lower_bound = (1.2,-2)
upper_bound = (460,2)
#lower_bound = SVector(4*10^(-6),-30)
#upper_bound = SVector(1.8*10^(-5),30)
xg = range(lower_bound[1], upper_bound[1], length = 10)
yg = range(lower_bound[2], upper_bound[2], length = 10)
mapper_grid = (xg, yg)

f = Figure()
time
fig,ax,hm = heatmap(xg, yg, basins[:,:,2],colormap = cgrad(:nipy_spectral, rev = true, categorical = true))
heatmap(f[1:3,1],xg, yg, hm[3],colormap = cgrad(:nipy_spectral, rev = true, categorical = true),
axis = (title = "Tipped Toda Oscillator",xlabel = L"x", ylabel = L"\dot{x}"))
scatter!(f[1:3,1],points,color = "black")
#Axis(f[1:3,1],title = "as")
#lines(f[4,1], randn(20),randn(20))
lines(f[4,1],t_plot,Gaussian_puls(t_plot,delta_phi,delta_phi,T_t),axis = (title = L"\text{Phase Shift} ϕ",xlabel = L"t", ylabel = L"\phi (rad)"))
scatter!(f[4,1],time,point_g,color = "blue")

f


frames = 3:iterationen

record(f, "animation_tipped_KM.mp4", frames;
        framerate = 8) do frame
        time[] = ti[frame]
        point_g[] = x_g[frame]
        points[] = Point2f[(x_p[frame], y_p[frame])]
        hm[3] = basins[:,:,frame] # update data
end
