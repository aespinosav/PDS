using JuMP, Gurobi
set_sover = GurobiSolver

type Polyhedron
    norms
    disps
end

type Orbit
    trajectory
    t
end

function P_K(r, K::Polyhedron)
    projection_prog = Model(solver=GurobiSolver())
    n_dims = length(r)
    
    @variable(projection_prog, x[1:n_dims])
    
    ex = @expression(
                     projection_prog,
                     dot(x - r, x - r)
                     )
    @objective(projection_prog, Min, ex)
    
    solve(projection_prog)
    sol = getvalue(x)
end

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

    ex = @expression(projection_model, sum(x.*x) - 2sum(p.*x))
    @objective(projection_model, Min, ex)
    
    solve(projection_model)
    sol = getvalue(x)
end

function euler_proj_onestep(x₀, F, K::Polyhedron, Δt)
    x₁ = project( x₀ - Δt*F(x₀) , K)
end

function orbit_proj(x₀, F, K::Polyhedron, Δt, t_steps; 
                            method_onestep=euler_proj_onestep, tol=1e-6)
                    
    method_onestep_kf =  (x,dt) -> method_onestep(x, F, K, dt)

    t = zeros(t_steps+1)
    conv = zeros(t_steps+1)
    dim = length(x₀)
    orbit = [zeros(dim) for i in 1:t_steps+1]

    orbit[1] = x₀
    t[1] = 0.0

    x_old = x₀
    t_old = 0.0

    total_steps = 0
    for i in 1:t_steps
        x_new = method_onestep_kf(x_old, Δt)
        t_new = Δt
        
        orbit[i+1] = x_new
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
    
    orbit = orbit[1:total_steps]
    orbit = hcat(orbit...)'
    orbit, t[1:total_steps]
end

#Make test polyhedron<
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
push!(displacements, d*sqrt(n_dims)/2)
push!(normals, (1/sqrt(n_dims)*ones(n_dims)))

K = Polyhedron(normals, displacements)

F = x -> -([0.3, 1.3] - x)

x₀ = [0.1, 0.1]
Δt = 0.01
time_steps = 100

orb = Orbit(orbit_proj(x₀, F, K, Δt, time_steps)...)
