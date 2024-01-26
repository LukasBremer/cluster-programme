using DynamicalSystems
using CairoMakie
using OrdinaryDiffEq
using JLD
using Random


function Parameters_time_dependence()

    
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
    #t0 = 2*pi*T_t*rand()
    t0 = delta_t
   return P,nu,nu_0,R_n,t0
end

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
    t0 = 2*pi*T_t*rand()

   return P,nu,nu_0,R_n,t0
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
function Gaussian_puls(t,A,A0,T,t0)

    t0, sigma, r = 15*T+t0, 10*T, .001

    return @. (A0*exp.(-(t-t0)^2/(2*sigma^2)) + A) * 1/2 * (sign(t-t0)+1) +  1/2 * (sign(-t + t0) + 1) * (A0*exp.(-(t-t0)^2/(2*(sigma*r)^2)) + A) 
end

function transformed_tipped_KM(x,P,t)
    p,nu,nu_0,R_n,t0 = P
    y = x[2]/x[1]
    z = 1/beta*log(x[1]/alpha)
    
    dx1 = x[2]
    dx2 = x[1]*(y^2 + ((p[3]+p[4]*y)*y^2 + (p[5]-p[6]*y)*z^(-4) - (p[7]+p[8]*y)*z^(-1) - (1+p[9]*y)*(p[10]+p[2]*sin(2*pi*nu/nu_0*t + Gaussian_puls(t,0,delta_phi,nu_0/nu,t0))) - p[11]*z*cos(2*pi*nu/nu_0*t + Gaussian_puls(t,0,delta_phi,nu_0/nu,t0)))
    /((1-p[9]*y)*z + p[12]))
    return SVector(dx1,dx2)
end

function Success_rate(iterations_strob,iterations_loop)

    success_value = 0
    for i in 1:iterations_loop
        p01 = Parameters()
        #global delta_phi = 2*pi*T_t*rand()
        xt_strob, yt_strob, t_s = KM_trajectory(u01_t,p01,iterations_strob,"d",transformed_tipped_KM,T_t)
        attractor =  length(unique(round.(xt_strob[iterations_strob-20:iterations_strob],digits=2)))
        if attractor == 2
            success_value = success_value + 1
        end

    end
    return success_value

end

function Success_rate_for_different_phases()
    
    success_list = zeros(100)
    phase_shift_list = zeros(100)

    for i in 1:100
        global delta_phi = 2*pi*(i/100)
        success_list[i] = Success_rate(120,300)
        phase_shift_list[i] = delta_phi
        print(i)
    end

    return success_list,phase_shift_list
end

function Success_rate_for_different_time()
    
    success_list = zeros(100)
    time_list = zeros(100)

    for i in 1:100
        global delta_t = T_t * (i/100)
        success_list[i] = Success_rate(120,300)
        time_list[i] = delta_t
        print(i)
    end

    return success_list,time_list
end

t = range(0,T_t,length = 100)
d_phi = range(0,2*pi,length = 100)
grid = zeros(100,100)
iterations_strob = 120

for i in 1:100
    
    global  delta_t = t[i]
    print(i)
    for j in 1:100
        
        global  delta_phi = d_phi[j]
        p01 = Parameters_time_dependence()
        #global delta_phi = 2*pi*T_t*rand()
        xt_strob, yt_strob, t_s = KM_trajectory(u01_t,p01,iterations_strob,"d",transformed_tipped_KM,T_t)
        attractor =  length(unique(round.(xt_strob[iterations_strob-20:iterations_strob],digits=2)))
        grid[i,j] = attractor
    end
end

f1= Figure()
ax = Axis(f1[1,1],ylabel= "phase shift [rad]",xlabel="tipping time [t/T]")
hm1 = heatmap!(ax,0:1,0:2*pi,grid)
#Colorbar(f1[1,2],hm1)
f1


#global Parameters
grid[50,60]
R_n,P_s = 10*10^(-6), 90*10^3 #Pa
alpha,beta,gamma,nu,kappa = 1,2,0.001,197*10^3,4/3
nu_0 = alpha * beta / (gamma * R_n)
delta_phi = .65*pi
T_t= nu_0/nu
delta_t = 0

#initial conditions and transformed initial conditions
R,U= R_n,20*10^(-4) #m, m/s 
u01 = SVector(R,U)
u01_t = RU_to_x1x2(R,U)
u01_t = SVector(100,2.)
p01 = Parameters()


#trajectories
iterationen = 120
x,y,t = KM_trajectory(u01,p01,iterationen,"d",KM,T_t)
xt,yt,t = KM_trajectory(u01_t,p01,iterationen,"c",transformed_tipped_KM,T_t)


success_value = Success_rate(120,100)
a = 70000

success_list,phase_shift_list = Success_rate_for_different_phases()
success_list

f = Figure() 
axis = Axis(f[1,1],title = "KM random time at constant phase",xlabel = "Tipping phase", ylabel = "success rate in %")
lines!(axis,phase_shift_list,success_list/3)
f

scatter(success_list)
scatter(phase_shift_list)

#Figure
f = Figure()
axis = Axis(f[1,1],xlabel="t dimensionless",ylabel="R in [m]",title="f=200kHz")
#lines!([0,1.5*10^4],[3.4,3.4])
scatter!(xt[length(t)-a:length(t)],yt[length(t)-a:length(t)],label="transformed")
scatter!(xt_strob[90:100],yt_strob[90:100])
#l2 = lines!(y,x,label="untransformed",color = "red")
axislegend()
f