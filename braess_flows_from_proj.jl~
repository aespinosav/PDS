using JuMP, Gurobi, PyPlot


"""
Polyhedron type
Defined by normals of planes of half-spaces and displacement from O
"""
type Polyhedron
    norms
    disps
end


"""
Orbit type
Trajectory array with time array
"""
type Orbit
    trajectory
    t
    steps
end
function Orbit(t, s)
    Orbit(t,s, length(t))
end


"""
Projection on polyhedron K
Using JuMP for minimising distance to K as quadratic problem
"""
function P_K(p, d)
    projection_model = Model(solver=GurobiSolver())
    
    #Flow variables
    n_dims = length(p)#dims of space
    @variable(projection_model, x[1:n_dims])
    
    #Constrain to positive orthant
    for i in 1:n_dims
        @constraint(projection_model, x[i] >= 0)    
    end
    #Constrain to simplex of route flows
    @constraint(projection_model, sum(x) == d)

    #Objective: minimise distance to set
    ex = @expression(projection_model, 0.5*sum(x.*x) - sum(p.*x))
    @objective(projection_model, Min, ex)
    
    solve(projection_model)
    getvalue(x)
end


"""
Euler single step
For projected system dx/dt = Π(x, -F(x))
"""
function euler_p(F, d, x₀, Δt)
    p = x₀ - Δt*F(x₀)
    x_new = P_K(p, d)
end


"""
Solve Euler
Euler metod for initial value x₀, for given time steps
"""
function orbit_proj_euler(F, d, x₀, Δt, t_steps; tol=0)
                 
    method_onestep =  (x_i,dt) -> euler_p(F, d, x_i, dt)

    t = zeros(t_steps+1)
    conv = zeros(t_steps+1)
    dim = length(x₀)
    r = [zeros(dim) for i in 1:t_steps+1]

    r[1] = x₀
    t[1] = 0.0

    x_old = x₀
    t_old = 0.0

    total_steps = 0
    for i in 1:t_steps
        x_new = method_onestep(x_old, Δt)
        t_new = Δt
        
        r[i+1] = x_new
        t[i+1] =  t_old + Δt
        
        Δx = norm(x_new - x_old)
        conv[i+1] = Δx

        x_old = x_new
        t_old = t_new
        
        total_steps += 1
        if Δx < tol
            break
        end
    end
    
    if total_steps < t_steps
        println("Converged in $total_steps time steps.\n")
    else
        println("After all t_steps ($t_steps), current x.\n")
    end
    
    r = r[1:total_steps]
    Orbit(hcat(r...)', t[1:total_steps])
end


"""
Make affine cost function vector from route-link incidence matrix R and parameter vectors a and b  
"""
function make_cost(a, b, R)
    B = diagm(b)
    y -> begin
             #link flow composition....
             x = R'*y 
             #transform back to route vars
             R*(a + B*x)
         end
end




####################################################




# Braess network (Route flows)
# ============================

# General
# -------

#Dimension of route-flow space
n_dims = 3 #number of routes...
n_simplex_dims = 2

#Time steps
t_steps = 100
Δt = 0.01
#Hmmm.... maybe Euler is to inaccurate... and also this is very slow.... but I guess there are faster ways of calculating the flows... snap... must plot aalytical flows here.... ooops haha...

# Set-up
# ------

#Cost parameters
a = [0.5, 1, 0.1, 1, 0.5]
b = [1, 0.5, 0.1, 0.5, 1]
#Route-link incidence matrix
R = [1 0 0 1 0;
     1 0 1 0 1;
     0 1 0 0 1]

#Cost vector
C = make_cost(a, b, R)
#negative sign in definition... hmmm...

#Demand
d_range = collect(0.001:0.01:2)

#Containers
ue_array = fill(zeros(n_dims), length(d_range))
flag_array = Array{String,1}(length(d_range))

#Iterate for demand values to get eq point for each
for (i,d) in enumerate(d_range)
    
    #Initial conditon
    x0 = zeros(n_dims)
    x0[1] = d
    
    #Integrte orbit
    orb = orbit_proj_euler(C, d, x0, Δt, t_steps, tol=1e-6)
    flag = length(orb.t) < t_steps ? "C" : "T"
    
    #Get UE pattern
    ue = orb.trajectory[end,:]
    ue_array[i] = copy(ue)
end


# Plotting
# --------

# Something is weird here cause the flows are not quite symmetric... I think it might just be that Euler is not the best method to use. Or that I am using JuMP in a dumb way...

clf()

xs = [a[1] for a in ue_array]
ys = [a[2] for a in ue_array]
zs = [a[3] for a in ue_array]

#Flows route 1
plot(d_range, xs, label="y₁")
#Flows route 2
plot(d_range, ys, label="y₂")
#Flows route 3
plot(d_range, zs, label="y₃")

legend(loc="top left", frameon=false)
xlabel("d", fontsize=18)
ylabel("yᵢ", fontsize=18)

#savefig("braess_flows_proj.pdf", bbox_inches="tight")
