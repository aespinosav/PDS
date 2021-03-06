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
Using JuMP for minimising distance to K as convex problem

In this case:

x1 + x2 = 1
x1, x2 >= 0
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
function orbit_proj_euler(F, d, x₀, Δt, t_steps; tol=1e-3)
                 
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
             x = R'*y
             R*(a + B*x)
         end
end


####################################################


# Braess network (Route flows)
# ============================

#Cost parameters
a = [0.5, 1, 0.1, 1, 0.5]
b = [1, 0.5, 0.1, 0.5, 1]
#Route-link incidence matrix
R = [1 0 0 1 0;
     1 0 1 0 1;
     0 1 0 0 1]
#Define cost function (vector valued function for route flows)
C = make_cost(a, b, R)
F = y -> C(y)
#negative sign in definition?

#Time steps
t_steps = 10000
Δt = 0.001


#Demand
d = 1

#Initial conditon
n_dims = 3 # set dimension by hand (3 routes...)
y₀ = zeros(n_dims)
y₀[1] = d

# Integrate orbits
orb = orbit_proj_euler(F, d, y₀, Δt, t_steps, tol=1e-6)
#orb.trajectory


#####################################

# Plotting
#---------

# Setup...

# barycentric coords: (a,b,c)
a = orb.trajectory[:,1]
b = orb.trajectory[:,2]
c = orb.trajectory[:,3]

# Trying to normalise size of simplex
a /= d
b /= d
c /= d

# Transform the data to cartesian...
#define transformation into barycentric
bary_trans = (a,b,c) -> [
                         0.5 * (( 2*b+c ) / ( a+b+c )), 
                         0.5*sqrt(3) * (c / (a+b+c))
                         ]

# transform coordinates -> barycentric coordinates
x = 0.5 * (( 2.*b+c ) ./ ( a+b+c ))
y = 0.5*sqrt(3) .* (c ./ (a+b+c))

# corners of triangle
corners = hcat([0, 0], [1, 0], [0.5,  sqrt(3)*0.5])'
# centre of triangle
κ = bary_trans(1,1,1)


# Plotting...
clf()

#Plot centre
plot(κ[1], κ[2], "o", color="grey", markersize=5)
text_pos = κ - [0.05, 0.05]
text(text_pos..., "O", fontsize=15)
#Plot route flow axes
for i in 1:3
    p2 = κ + 1.5*(corners[i,:] - κ)
    plot([κ[1], p2[1]], [κ[2], p2[2]], linewidth=0.8, color="grey")
end

# Plot simplex
plot([corners[3,1], corners[1,1]],[corners[3,2], corners[1,2]], "-k", linewidth=3)
for k in 2:3
    plot([corners[k-1,1], corners[k,1]],[corners[k-1,2], corners[k,2]], "-k", linewidth=3)
end
# Label corners
text_pos = bary_trans(1,0,0) - [0.12, 0.08]
text(text_pos..., "(1, 0, 0)", fontsize=14)

text_pos = bary_trans(0,1,0) - [0.08, 0.08]
text(text_pos..., "(0, 1, 0)", fontsize=14)

text_pos = bary_trans(0,0,1) + [-0.07, 0.03]
text(text_pos..., "(0, 0, 1)", fontsize=14)

# plot trajectory
plot(x, y, "orange", linewidth=2)

plot(x[1], y[1], "ok", markersize=8)

text_pos = [x[1], y[1]] - 0.02*[2.85, -1.25]
text(text_pos..., "y₀", fontsize=20)

plot(x[end], y[end], "or", markersize=8)

text_pos = [x[end], y[end]] + 0.015*[2, 0.5]
text(text_pos..., "yᵁᴱ", fontsize=20)

xlim([-0.2, 1.2])
ylim([-0.2, 1.0])

ax = gca()
ax[:set_aspect]("equal")

tick_params(
    axis="both",          # changes apply to the x-axis
    which="both",      # both major and minor ticks are affected
    bottom=false,      # ticks along the bottom edge are off
    top=false,
    left=false,
    right=false,         # ticks along the top edge are off
    labelbottom=false,
    labelleft=false)

box(false)
#savefig("braess_ternary_projected.pdf", bbox_inches="tight")
