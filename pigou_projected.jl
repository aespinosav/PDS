using JuMP, Gurobi, PyPlot

# Definitions

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
Using JuMP for minimising distance to K as convex problem

In this case:

x1 + x2 = 1
x1, x2 >= 0
"""
function P_K(p)
    projection_model = Model(solver=GurobiSolver())
    
    n_dims = length(p)#dims of space
    @variable(projection_model, x[1:n_dims])

    for i in 1:n_dims
        @constraint(projection_model, x[i] >= 0)    
    end
    @constraint(projection_model, sum(x) == 1)

    ex = @expression(projection_model, 0.5*sum(x.*x) - sum(p.*x))
    @objective(projection_model, Min, ex)
    
    solve(projection_model)
    sol = getvalue(x)
end

"""
Euler single step
For projected system dx/dt = Π(x, -F(x))
"""
function euler_p(F, x₀, Δt)
    p = x₀ - Δt*F(x₀)
    x_new = P_K(p)
end

"""
Solve Euler
Euler metod for initial value x₀, for given time steps
"""
function orbit_proj_euler(F, x₀, Δt, t_steps; tol=1e-3)
                 
    method_onestep =  (x_i,dt) -> euler_p(F, x_i, dt)

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

#Plotting help
meshgrid(x, y) = (repeat(x, outer=length(y)), repeat(y, inner=length(x)))

#######################################################
#######################################################
#######################################################

# Script 

n_dims = 2
k = n_dims + 1 
#Faces of polyhedron
#Positive orthant
normals = [zeros(n_dims) for i in 1:n_dims]
displacements = zeros(n_dims)
for i in 1:n_dims
            normals[i][i] = -1.0
end
#Simplex ( x1 + x2 < d)
d = 1
push!(displacements, 0.5*d*sqrt(n_dims))
push!(normals, (1/sqrt(n_dims)*ones(n_dims)))

#Flow
F = x -> [1, 2x[2]]

#Initial conditon
x₀ = [0.1, 0.9]

t_steps = 1000
Δt = 0.01

orb = orbit_proj_euler(F, x₀, Δt, t_steps, tol=1e-8)

#Plotting
clf()

# K
plot([0,1],[1,0],"-g", linewidth=3)
hlines(0, -0.2, 1.2, colors="gray")
vlines(0, -0.2, 1.2, colors="gray")

#Vector field
range = -0.5:0.1:1.5
X, Y = meshgrid(range, range)
f1 = (x,y) -> -F([x,y])[1]
f2 = (x,y) -> -F([x,y])[2]
U = f1.(X, Y)
V = f2.(X, Y)
quiver(X, Y, U, V)

#for i in 1:length(trajectories)
#    o = trajectories[i]
#    plot(o.trajectory[:,1], o.trajectory[:,2], "-k")
#end
plot(orb.trajectory[:,1], orb.trajectory[:,2], "-", linewidth=3.5, color="orange")
#plot(orb_neg.trajectory[:,1], orb_neg.trajectory[:,2], "-c", linewidth=3)
plot(x₀[1], x₀[2], "or", markersize=7)
text_pos = x₀ + 0.015*[1,1] 
text(text_pos..., "x₀", fontsize=20)
plot(0.5, 0.5, "ok", markersize=7)
text(0.515, 0.515, "xˢᵒ", fontsize=20)

xlim([-0.2,1.2])
ylim([-0.2,1.2])

xlabel("x₁", fontsize=18)
ylabel("x₂", fontsize=18)

ax = gca()
ax[:set_aspect]("equal")

# Save
savefig("pigou_pds.pdf")
