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
"""
function P_K(p, K::Polyhedron)
    projection_model = Model(solver=GurobiSolver())
    
    n_dims = length(p)#dims of space
    @variable(projection_model, x[1:n_dims])

    A = hcat(K.norms...)'
    d = K.disps
    for i in 1:length(d)
        @constraint(projection_model, sum(A[i,j]*x[j] for j in 1:n_dims) <= d[i])    
    end
    #@constraint(projection_model, A*x .<= d)

    ex = @expression(projection_model, 0.5*sum(x.*x) - sum(p.*x))
    @objective(projection_model, Min, ex)
    
    solve(projection_model)
    sol = getvalue(x)
end

"""
Euler single step
For projected system dx/dt = Π(x, -F(x))
"""
function euler_p(F, K, x₀, Δt)
    p = x₀ - Δt*F(x₀)
    x_new = P_K(p, K)
end

"""
Solve Euler
Euler metod for initial value x₀, for given time steps
"""
function orbit_proj_euler(F, K, x₀, Δt, t_steps; tol=1e-3)
                 
    method_onestep =  (x_i,dt) -> euler_p(F, K, x_i, dt)

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
    println(r)
    
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

K = Polyhedron(normals, displacements)

#Flow:elliptic trajectories
F1 = x -> [-x[2], 4x[1]]
#Off-centred
c1 = ones(2)/sqrt(2)
F = x -> F1(x - c1)
#F = x -> F1(x - 0.7*[sqrt(2), sqrt(2)])
#F = x -> -([2, 2] - x)

#Initial conditon
x₀ = [0.18, 0.18]

t_steps = 500
Δt = 0.01

orb = orbit_proj_euler(F, K, x₀, Δt, t_steps)
#orb_neg = orbit_proj_euler(F, K, x₀, -Δt, 100)

#Other orbits
trajectories = Orbit[]
v1 = [1,0.0001]
for γ in 0:0.2:1
    o_ω1 = orbit_proj_euler(F, K, γ*v1, Δt, t_steps)
    push!(trajectories, o_ω1)
end 

for γ in 0:0.05:1
    v2 = γ*[-1,1] + v1
    o_ω2 = orbit_proj_euler(F, K, v2, Δt, t_steps)
    push!(trajectories, o_ω2)
end


#Plotting
clf()

# K
plot([0,0],[0,1],"-g", linewidth=5)
plot([0,1],[1,0],"-g", linewidth=5)
plot([1,0],[0,0],"-g", linewidth=5)

#Vector field
range = -0.5:0.1:1.5
X, Y = meshgrid(range, range)
f1 = (x,y) -> -F([x,y])[1]
f2 = (x,y) -> -F([x,y])[2]
U = f1.(X, Y)
V = f2.(X, Y)
quiver(X, Y, U, V)

for i in 1:length(trajectories)
    o = trajectories[i]
    plot(o.trajectory[:,1], o.trajectory[:,2], "-k")
end
plot(orb.trajectory[:,1], orb.trajectory[:,2], "-m", linewidth=3)
#plot(orb_neg.trajectory[:,1], orb_neg.trajectory[:,2], "-c", linewidth=3)
plot(x₀[1], x₀[2], "or")
plot(c1[1], c1[2], "oy")
xlim([-0.2,1.2])
ylim([-0.2,1.2])

ax = gca()
ax[:set_aspect]("equal")
