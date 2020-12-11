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
    num_steps
end
Orbit(t, s) = Orbit(t, s, length(t))



"""
Projection on polyhedron K
Using JuMP for minimising distance to K as quadratic problem
"""
function P_K(p, d)
    projection_model = Model(solver=GurobiSolver())
    
    #Flow variables
    n_dims = length(p) #dims of space
    @variable(projection_model, x[1:n_dims])
    
    #Constrain to positive orthant
    for i in 1:n_dims
        @constraint(projection_model, x[i] >= 0)    
    end
    #Constrain to simplex of route flows
    @constraint(projection_model, sum(x) == d)‪

    #Objective: minimise distance to set
    ex = @ expression(projection_model, 0.5*sum(x.*x) - sum(p.*x))
    @objective(projection_model, Min, ex)
    
    solve(projection_model)
    getvalue(x)
end


"""
Euler single step for projected system,
dx/dt = Π(x, -F(x))
"""
function euler_p(F, d, x₀, Δt)
    p = x₀ - Δt*F(x₀)
    x_new = P_K(p, d)
end


"""
Euler single step
For projected system dx/dt = Π(x, -F(x))

using dafermos-sparrow "exact" algorithm by 
decompositon into subprpblems.

singel OD for now...
"""
function euler_dafermos_sparrow(F, d, x₀, Δt)
    # xₜ₊₁=P_K(xₜ - Δt*F(xₜ)) or something...
    # Becomes...
    # xₜ₊₁ = min (1/2) * x ⋅ x - (xₜ - Δt*F(xₜ)) ⋅  x
    # Exact equilibration 
    
    hₜ = aₜ * C(yₜ) - yₜ
    
    h_ind_sort = sortperm(hₜ)
    inv_ind_sort =  sortperm(h_ind_sort)
    
    y_sorted = copy(yₜ[h_ind_sort])
    h_sorted = copy(hₜ[h_ind_sort])
    
    λ_q = length(y) # isn't this just the cost of the used route?
    for i in 1:length(y)
        λ = (sum(h_sorted[1:i]) + d)/i
        if h_sorted[i] < λ
            q = i # index of stop path
            break
        end
    end
    
    y_sorted = zeros(length(yₜ))
    for i in 1:q
        y_sorted[i] = λ_q - h_sorted
    end
    y_ue = y_sorted[inv_ind_sort]
end


"""
Runge-Kutta Order 4 for projected system,
dx/dt = Π(x, -F(x))
"""
function rk4_p(F, d, x₀, Δt)
    h₁ = F( x₀ )
    h₂ = F( x₀ - (Δt*h₁)/2 )
    h₃ = F( x₀ - (Δt*h₂)/2 )
    h₄ = F( x₀ -  Δt*h₃ )
    
    Fᵣₖ = (1/6)*(h₁ + 2h₂ + 2h₃ + 4h₄ )
    p = x₀ - Δt*Fᵣₖ
    
    x_new = P_K(p, d)
end


"""
Numerically integrate projected dynamical system, call as,
    
    `orbit_proj(F, d, x₀, Δt, t_steps; tol=0, method=euler_p)`

F   -   (unconstrained) field
d   -   demand (displacement of flow-simplex)
x₀  -   initial value
Δt  -   timestep size
t_steps -   number of max timesteps to try

tol -   abs val diff threshold to converge
method  -   integration method (one step)

Default method is set to Euler metod.
"""
function orbit_proj(F, d, x₀, Δt, t_steps; tol=1e-7, method=euler_p)
                 
    integ_onestep =  (x_i,dt) -> euler_p(F, d, x_i, dt)

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
        x_new = integ_onestep(x_old, Δt)
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


"""
Transformation to barycentric coordinates 3D -> 2D
"""
bary_trans = (a,b,c) -> begin
                            x = 0.5 * (( 2*b+c ) / ( a+b+c ))
                            y = 0.5*sqrt(3) * (c / (a+b+c)) 
                            [x, y]
                        end
barry_trans(v::Array{Float64,1}) = barry_trans(v...)
#careful, only works for 3d


####################################################


# Braess network (Route flows)
# ============================

# General
# -------

#Dimension of route-flow space
n_dims = 3 #number of routes...
n_simplex_dims = 2

#Time steps
t_steps = 400
Δt = 0.01

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
d_range = collect(0.1:0.1:2)

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

#Transform to barycentric coords
Barry = map((x,d)->bary_trans(x*(d/sqrt(3))...),
            ue_array,
            d_range)
Barry = hcat(Barry...)'


# Plotting
# --------

clf()

#corners of triangle
corners = hcat([0, 0], [1, 0], [0.5,  sqrt(3)*0.5])'
#centre
κ = bary_trans(1,1,1)
#Plot route-flow axes (y1, y2, y3)
for i in 1:3
    p2 = κ + 1.5*(corners[i,:] - κ)
    plot([κ[1], p2[1]], [κ[2], p2[2]], linewidth=0.8, color="grey")
end

# Plot origin
plot(κ[1], κ[2], "o", color="grey", markersize=5, markeredgecolor="grey")
text_pos = κ - [0.05, 0.05]
text(text_pos..., "O", fontsize=15)

# Plot simplex
simplex_color = "black"
# Corners of triangle
corners = hcat([0, 0], [1, 0], [0.5,  sqrt(3)*0.5])'
plot([corners[3,1], corners[1,1]],[corners[3,2], corners[1,2]], "-", color=simplex_color, linewidth=3)
for k in 2:3
    xs = [corners[k-1,1], corners[k,1]]
    ys = [corners[k-1,2], corners[k,2]]
    plot(xs, ys, "-", color=simplex_color, linewidth=3)
end

#Label corners
#corner y1
text_pos = bary_trans(1,0,0) - [0.12, 0.08]
text(text_pos..., "(d, 0, 0)", fontsize=14)
#corner y2
text_pos = bary_trans(0,1,0) - [0.08, 0.08]
text(text_pos..., "(0, d, 0)", fontsize=14)
#corner y3
text_pos = bary_trans(0,0,1) + [-0.07, 0.03]
text(text_pos..., "(0, 0, d)", fontsize=14)

#Plot trace of eq point
plot(Barry[:,1], Barry[:,2], "-" , color="grey", linewidth=2.5)
#Stagger plotted points
n_steps = 20
stride_range = 1:Int(floor(length(d_range)/n_steps)):length(d_range)

ax = gca()
ax[:set_aspect]("equal")

scatter(Barry[stride_range,1], Barry[stride_range,2], s=80, c=d_range[stride_range], cmap="gnuplot")

#Trimmings
xlim([-0.2, 1.2])
ylim([-0.2, 1.0])

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

scatter(Barry[stride_range,1], Barry[stride_range,2], s=80, c=d_range[stride_range], cmap="gnuplot")

colorbar(shrink=0.75)
text(1.45, 0.42, "d", fontsize=18)

#savefig("braess_trace_eq_simplex.pdf", bbox_inches="tight")
