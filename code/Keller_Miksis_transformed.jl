using DynamicalSystems
using CairoMakie
using OrdinaryDiffEq
using JLD


function Parameters()

    
    P_stat = 100*10^3 #Pa
    P_v = 2.33*10^3 #Pa 
    sigma = 0.0725 #N/m 
    rho = 998. #kg/m^3
    mu = 0.001 #Ns/m^2
    c = 1500. #m/s
    delta = alpha/gamma
    

    P = zeros(12)
    P[1] = nu/nu_0
    P[2] = P_s/(delta^2*beta*rho)
    P[3] = -3/(2*beta)
    P[4] = delta/(2*beta*c)
    P[5] = 1/(rho*delta^2*beta)*(P_stat-P_v + 2*sigma/R_n)
    P[6] = 3*delta/c*P[5]
    P[7] = 2*sigma/(delta^2*beta*rho*R_n)
    P[8] = 4*mu/(delta*beta*rho*R_n)
    P[9] = delta/c
    P[10] = (P_stat-P_v)/(rho*delta^2*beta)
    P[11] = 2*pi*nu*P_s*R_n/(delta^2*beta*rho*c)
    P[12] = 4*mu/(rho*c)
   return P,nu,nu_0,R_n
end
#both systems (transformed and not transformed) use the slow time scale for better numerical integration
function KM(x,p,t)
    R = x[1]
    U = x[2]
    
    P_stat = 100*10^3 #Pa
    P_v = 2.33*10^3 #Pa 
    sigma = 0.0725 #N/m 
    rho = 998. #kg/m^3
    mu = 0.001 #Ns/m^2
    c = 1500. #m/s
    delta = alpha/gamma

    dR = 1/nu_0 * U
    dU = 1/nu_0 * (-U^2/2*(3-U/c)+(1+(1-3*kappa)*U/c)*((P_stat-P_v)/rho+2*sigma/(rho*R_n))*(R_n/R)^(3*kappa)-2*sigma/(rho*R)
    -4*mu*U/(rho*R)-(1+U/c)*(P_stat-P_v+P_s*sin(2*pi*nu/nu_0*t + delta_phi))/rho-R*2*pi*nu*P_s/(rho*c)*cos(2*pi*nu/nu_0*t + delta_phi))/
    ((1-U/c)*R+4*mu/(rho*c))

    return SVector(dR,dU)
end

function transformed_KM(x,P,t)
    p, nu, nu_0,R_n = P
    y = x[2]/x[1]
    z = 1/beta*log(x[1]/alpha)
    
    dx1 = x[2]
    dx2 = x[1]*(y^2 + ((p[3]+p[4]*y)*y^2 + (p[5]-p[6]*y)*z^(-4) - (p[7]+p[8]*y)*z^(-1) - (1+p[9]*y)*(p[10]+p[2]*sin(2*pi*nu/nu_0*t+ delta_phi)) - p[11]*z*cos(2*pi*nu/nu_0*t + delta_phi))
    /((1-p[9]*y)*z + p[12]))
    return SVector(dx1,dx2)
end

function KM_trajectory(u01,p01,iterationen,discrete_continuous,KM_s,T)
    
    diffeq = (alg = Vern9(), abstol = 1e-9, reltol = 1e-9)

    if discrete_continuous == "c"
        KM_system = ContinuousDynamicalSystem(KM_s,u01,p01; diffeq = diffeq)
        total_time = iterationen * T 
        X2,t2 = trajectory(KM_system,total_time)
    else
        KM_system_1 = ContinuousDynamicalSystem(KM_s,u01,p01; diffeq = diffeq)
        KM_system = StroboscopicMap(KM_system_1, T)
        X2,t2 = trajectory(KM_system,iterationen)
    end
    
    x2,y2 = X2[:,1],X2[:,2]
    return x2,y2,t2
end

function RU_to_x1x2(R,U)
    x1 = @.alpha*exp(beta*R/R_n)
    x2 = @.gamma*U*exp(beta*R/R_n)
    return x1,x2
end

function x1x2_to_RU(x1,x2)
    R = @.R_n/beta*log(x1/alpha)
    U = @.alpha/gamma*x2/x1
    return R,U    
end

function find_nearest_bassin_val(xg,yg,basins,x)
    i = findmin(@.abs(xg - x[1]))[2]
    j = findmin(@.abs(yg - x[2]))[2]
    bas_val = basins(i,j)
    return bas_val
end

function calculate_basins(mapper_grid,system)
    diffeq = (alg = Vern9(), abstol = 1e-9, reltol = 1e-9)
    KM_system = ContinuousDynamicalSystem(system,u01_t,p01; diffeq = diffeq)
    strob_KM = StroboscopicMap(KM_system, T_t)
    #change
    mapper = AttractorsViaRecurrences(strob_KM, mapper_grid;
    sparse = true, Ttr = 10, mx_chk_fnd_att = 500)

    basins, attractors = basins_of_attraction(mapper, mapper_grid; show_progress = true)
    return basins, attractors
end
    
function Angle_of_Attractor(u01_t,p01,iterationen)
    
    xt,yt,t = KM_trajectory(u01_t,p01,iterationen,"c",transformed_KM,T_t)
    i = findmin(@.abs(t - T_t * (iterationen-1)))[2]
    phi =  @.(t * 2*pi*nu/nu_0)
    return xt[i+1:length(t)-1],yt[i+1:length(t)-1],t[i+1:length(t)-1], phi[i+1:length(t)-1] .% (2*pi)
end

function current_position_angle_attractor(u01_t,p01,iterationen,phi_0)
    
    x,y,t,phi = Angle_of_Attractor(u01_t,p01,iterationen)
    i = findmin(@.abs(phi-phi_0))[2]
    return x[i],y[i]
end

function find_Attractor_of_current_position(u01_t,p01,iterationen,phi,xg,yg,basins)
    
    x,y,t,phi = Angle_of_Attractor(u01_t,p01,iterationen)
    bas_val = zeros(length(phi))

    for k in 1:length(phi)
        x_c,y_c,t_c,phi_c = x[k],y[k],t[k],phi[k]
        i = findmin(@.abs(xg - x_c))[2]
        j = findmin(@.abs(yg - y_c))[2]
        bas_val[k] = basins[i,j]
    end

    return phi,bas_val
end
#functions for basin plots 
function Gaussian_puls(t,A,A0,T)

    t0, sigma, r = 25*T, 10*T, .001

    return @.(-A0*exp.(-(t-t0)^2/(2*sigma^2)) + A) * 1/2 * (sign(t-t0)+1) +  1/2 * (sign(-t + t0) + 1) * (-A0*exp.(-(t-t0)^2/(2*(sigma*r)^2)) + A) 
end

function transformed_tipped_KM(x,P,t)
    p, nu, nu_0,R_n = P
    y = x[2]/x[1]
    z = 1/beta*log(x[1]/alpha)
    
    dx1 = x[2]
    dx2 = x[1]*(y^2 + ((p[3]+p[4]*y)*y^2 + (p[5]-p[6]*y)*z^(-4) - (p[7]+p[8]*y)*z^(-1) - (1+p[9]*y)*(p[10]+p[2]*sin(2*pi*nu/nu_0*t + Gaussian_puls(t,delta_phi,delta_phi,nu_0/nu))) - p[11]*z*cos(2*pi*nu/nu_0*t + Gaussian_puls(t,delta_phi,delta_phi,nu_0/nu)))
    /((1-p[9]*y)*z + p[12]))
    return SVector(dx1,dx2)
end

function transformed_tipped_fp_KM(x,P,t)
    P_temp,t_x = P
    p, nu, nu_0,R_n= P_temp
    y = x[2]/x[1]
    z = 1/beta*log(x[1]/alpha)
    
    dx1 = x[2]
    dx2 = x[1]*(y^2 + ((p[3]+p[4]*y)*y^2 + (p[5]-p[6]*y)*z^(-4) - (p[7]+p[8]*y)*z^(-1) - (1+p[9]*y)*(p[10]+p[2]*sin(2*pi*nu/nu_0*t + Gaussian_puls(t_x,delta_phi,delta_phi,nu_0/nu))) - p[11]*z*cos(2*pi*nu/nu_0*t + Gaussian_puls(t_x,delta_phi,delta_phi,nu_0/nu)))
    /((1-p[9]*y)*z + p[12]))
    return SVector(dx1,dx2)
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

function compute_frame(i)

    p_current_time = [Parameters(), ti[i]]
    print(Gaussian_puls(ti[i],delta_phi,delta_phi,nu_0/nu),ti[i])
    tipped_KM_system = ContinuousDynamicalSystem(transformed_tipped_fp_KM,u01_t,p_current_time; diffeq)
    strob_tipped_KM = StroboscopicMap(tipped_KM_system, T_t)
    mapper = AttractorsViaRecurrences(strob_tipped_KM, mapper_grid;
    sparse = true, Ttr = 10, mx_chk_fnd_att = 500,
    )
    basins, attractors = basins_of_attraction(mapper, mapper_grid; show_progress = true)
    new_b = new_basins(basins,attractors)
    #new_b = cat(dims=3,new_b, new_basins(basins,attractors)
    return new_b
end

#global Parameters
R_n,P_s = 10*10^(-6), 90*10^3 #Pa
alpha,beta,gamma,nu,kappa = 1,2,0.001,197*10^3,4/3
nu_0 = alpha * beta / (gamma * R_n)
delta_phi = 0.65*pi
T_t= nu_0/nu


#initial conditions and transformed initial conditions
R,U= R_n,20*10^(-4) #m, m/s 
u01 = SVector(R,U)
u01_t = RU_to_x1x2(R,U)
u01_t = SVector(100,2.)
p01 = Parameters()


# #trajectories
# iterationen = 100
# x,y,t = KM_trajectory(u01,p01,iterationen,"c",KM,T_t)
# xt,yt,t = KM_trajectory(u01_t,p01,iterationen,"c",transformed_tipped_KM,T_t)
# xt_strob, yt_strob, t_s = KM_trajectory(u01_t,p01,iterationen,"d",transformed_tipped_KM,T_t)
# save("trajectory_tipped_KM.jld","data",[xt_strob,yt_strob])
# a = 70000

# #Figure
# f = Figure()
# axis = Axis(f[1,1],xlabel="t dimensionless",ylabel="R in [m]",title="f=200kHz")
# #lines!([0,1.5*10^4],[3.4,3.4])
# lines!(xt[length(t)-a:length(t)],yt[length(t)-a:length(t)],label="transformed")
# #scatter!(xt_strob[30:40],yt_strob[30:40])
# #l2 = lines!(y,x,label="untransformed",color = "red")
# axislegend()
# # f


# minimum(xt[length(t)-a:length(t)])
# maximum(xt[length(t)-a:length(t)])

# minimum(yt[length(t)-a:length(t)])
# maximum(yt[length(t)-a:length(t)])
#save("comparisson_KM_trans_untransformed200kHz.pdf",f)


#Basins
lower_bound = SVector(1.2,-2)
upper_bound = SVector(460,2)
#lower_bound = SVector(4*10^(-6),-30)
#upper_bound = SVector(1.8*10^(-5),30)
xg = range(lower_bound[1], upper_bound[1], length = 50)
yg = range(lower_bound[2], upper_bound[2], length = 50)
mapper_grid = (xg, yg)

#Plot
iterationen = 3
#basins = load("basins_KM_0_pi.jld","data")
ti = range(15*T_t,80*T_t,iterationen)
#basins, attractors = calculate_basins(mapper_grid,transformed_KM)
diffeq = (alg = Vern9(), abstol = 1e-9, reltol = 1e-9)
new_b = Array{Matrix{Float64}}(undef, iterationen)

for i in 1:iterationen
    new_b[i] = compute_frame(i)
end

#save("basins_tipped_gaussian_KM_test.jld","data",new_b)

# x,y,t,phi = Angle_of_Attractor(u01_t,p01,iterationen)
# x_0, y_0 = current_position_angle_attractor(u01_t,p01,iterationen, .65*pi)
# phi,bas_val = find_Attractor_of_current_position(u01_t,p01,iterationen,phi,xg,yg,basins)

# f = Figure()
# axis = Axis(f[1,1], xlabel="x1", ylabel="x2", title="Basins of attraction; transformed system")
# heatmap!(xg, yg, basins)
# pl = lines!(x,y,color= phi)
# Colorbar(f[1, 2], pl)
# #lines!(xt[length(t)-a:length(t)],yt[length(t)-a:length(t)],label="transformed")
# scatter!(x_0,y_0)
# f

# f2 = Figure()
# axis = Axis(f2[1,1], xlabel="phi", ylabel="attractor", title="Attractor angle")
# lines!(phi,bas_val)
# f2