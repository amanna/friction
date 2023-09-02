using DifferentialEquations, Plots, StaticArrays, CUDA, Distributed, LinearAlgebra

function prob_func(prob, i, repeat)
    """
    This function creates new ODEProblem by changing the parameters
    where i is the index of the trajectories
    """
    # mass of the system
    m = 1.0
    
    # different values of stiffness (k)
    k_arr = range(10, 11, 10)
    k = k_arr[i]
    c = 0.5

    # friction coefficients
    mu_s = 0.3
    mu_k = 0.2

    # forcing function values
    A = 0
    ω = 1  
    p = (m, k, c, mu_s, mu_k, A, ω )
    remake(prob, u0 = prob.u0,p = p)
end


# Define the differential equation with forcing function
function spring_mass_damper_forcing!(du, u, p, t)
    """
    This function is a differential equation to spring mass damper System
    
    parameters
    ----------
    du: velocity(v), acceleration(a)
    u: position (x), velocity (v)
    p: parameters (m, k, c, mu_s, mu_k, A, ω)
    t: time

    Returns
    -------
    nothing
    """
    
    m, k, c, mu_s, mu_k, A, ω = p
    x, v = u
    
    # Calculate the forces
    @inbounds f_spring = -k * x
    @inbounds f_damping = -c * v
    @inbounds f_forcing = A*sin(ω *t)   # Sinusoidal forcing function
    @inbounds F_norm = m*9.81
    

    @inbounds tol_val = 1e-08
    @inbounds v_abs = abs.(v)
    @inbounds v_sign = sign.(v)

    # striebeck friction conditions
    @inbounds f_friction = min.(mu_s .* F_norm, f_spring + f_damping) .* (v_abs .< tol_val ) .* v_sign + mu_k .* F_norm .* (v_abs .> tol_val) .* v_sign
    # @inbounds f_friction = 0
    # @inbounds f_damping = 0
    # @inbounds f_forcing = 0
    @inbounds a =  (f_spring + f_damping + f_friction + f_forcing) / m
    
    du[1] = v
    du[2] = a
    return nothing
end

# Set up initial conditions and parameters
u0 = [-10, 0]  # initial position and velocity

# parameters
p = [1.0, 10, 0.5, 0.3, 0.2 , 0, 1 ]  # mass, spring constant, damping coefficient, static friction coefficient, kinetic friction coefficient, amplitude, angular frequency

# time of simulation
tspan = [0.0, 1000.0]

# total trajectories
traj_vals = 10

# Define the problem
prob = ODEProblem(spring_mass_damper_forcing!, u0, tspan, p)

# Solve the problem

# Create an EnsembleProblem to be solved for different trajectories
# we can either change the parameters values or initial conditions

monteprob = EnsembleProblem(prob, prob_func = prob_func, safetycopy = false)
sim = solve(monteprob, AutoTsit5(Rosenbrock23()), EnsembleThreads(), trajectories = traj_vals, saveat = 0.01, maxiters = 1e10)


# Extract values
t = sim[1].t


# plotting velocity for different stiffness 
plot1 = plot()
for i in 1:10
    plot!(t, sim[i][2,:], label = "$i")
end
A = 0
ω = 1

# plot forcing function
plot!(t, A.*sin.(ω .*t), label = "FF")
xlabel!("time (s)")
ylabel!(" Velocity")
title!("Vel vs time")
display(plot1)


# phase plots
plot2 = plot()

for i in 1:10
    plot!(sim[i][1,:], sim[i][2,:], label = "$i", legend = false)
end
xlabel!("Position")
ylabel!(" Velocity")
title!("Phase potraits")
display(plot2)
