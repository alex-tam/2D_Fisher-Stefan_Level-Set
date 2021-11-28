# Control function for 2D Fisher-Stefan level-set solutions
# Alex Tam, 02/11/2021

# Load packages
using Parameters
using DifferentialEquations
using Plots
using LaTeXStrings
using DelimitedFiles

# Include external files
include("domain.jl")
include("fkpp.jl")
include("velocity.jl")
include("level-set.jl")
include("reinitialisation.jl")
include("draw.jl")

"Data structure for model parameters"
@with_kw struct Params
    D::Float64 = 1.0 # [-] Diffusion coefficient
    λ::Float64 = 1.0 # [-] Reaction rate
    κ::Float64 = 1.0 # [-] Stefan parameter
    u_f::Float64 = 0.0 # [-] Density at interface
    u_b::Float64 = 0.0 # [-] Density at computational boundary
    α::Float64 = 0.5 # [-] Maximum initial density
    β::Float64 = 10.0 # [-] Initial interface radius
    θb::Float64 = 0.01 # [-] Threshold for whether a grid point is close to interface
    θ::Float64 = 1.99 # [-] Parameter for minmod flux-limiter
    Lx::Float64 = 20.0 # [-] Spatial domain limit (x)
    Ly::Float64 = 20.0 # [-] Spatial domain limit (y)
    T::Float64 = 100.0 # [-] End time
    Nx::Int = 201 # [-] Number of grid points (x)
    Ny::Int = 201 # [-] Number of grid points (y)
    Nt::Int = 50001 # [-] Number of time steps
    V_Iterations::Int = 20 # [-] Number of iterations for velocity extrapolation PDE
    ϕ_Iterations::Int = 20 # [-] Number of iterations for reinitialisation PDE
end

"Generate initial density and level-set function"
function ic(par::Params, x, y)
    U = Array{Float64}(undef, par.Nx, par.Ny) # Pre-allocate 2D array of U
    ϕ = Array{Float64}(undef, par.Nx, par.Ny) # Pre-allocate 2D array of ϕ
    # Loop over grid points
    for i = 1:length(x)
        for j = 1:length(y)
            # Circle
            if sqrt((x[i]-par.Lx/2)^2 + (y[j]-par.Ly/2)^2) <= par.β
                U[i,j] = par.α # Constant density (inside Ω)
            else
                U[i,j] = par.u_f # Initial density (outside Ω)
            end
            ϕ[i,j] = (sqrt((x[i]-par.Lx/2)^2 + (y[j]-par.Ly/2)^2) - par.β) # Initial signed-distance
            # # Rectangle
            # if (x[i] >= par.Lx/2 - 10.0/2) && (x[i] <= par.Lx/2 + 10.0/2) && (y[j] >= par.Ly/2 - 4.0/2) && (y[j] <= par.Ly/2 + 4.0/2)
            #     U[i,j] = par.α
            #     ϕ[i,j] = -0.5
            # else
            #     U[i,j] = par.u_f
            #     ϕ[i,j] = 0.5
            # end
        end
    end
    return U, ϕ
end

"Build a vector of 2D data, ordered by entries in D"
function build_vector(U::Array{Float64}, D)
    u = Vector{Float64}() # Pre-allocate empty vector
    for gp in D
        push!(u, U[gp.xInd, gp.yInd])
    end
    return u
end

"Build 2D array of data from vector"
function build_u_matrix(u::Vector, par, D)
    U = par.u_b*ones(par.Nx, par.Ny) # Pre-allocate (incorporate a Dirichlet condition on computational boundary)
    for i = 1:length(D)
        U[D[i].xInd, D[i].yInd] = u[i]
    end
    return U
end

"Build 2D array of data from vector"
function build_v_matrix(v::Vector, par, D)
    V = zeros(par.Nx, par.Ny) # Pre-allocate (incorporate a Dirichlet condition on computational boundary)
    for i = 1:length(D)
        V[D[i].xInd, D[i].yInd] = v[i]
    end
    return V
end

"Compute a solution"
function fisher_stefan_2d()
    # Parameters and domain
    par = Params() # Initialise data structure of model parameters
    x = range(0, par.Lx, length = par.Nx); dx = x[2] - x[1] # Computational domain (x)
    y = range(0, par.Ly, length = par.Ny); dy = y[2] - y[1] # Computational domain (y)
    t = range(0, par.T, length = par.Nt); dt = t[2] - t[1] # Time domain
    # Initial condition
    U, ϕ = ic(par, x, y) # Obtain initial density and ϕ
    ϕ = reinitialisation(ϕ, par, dx, dy, 100)
    draw_heat(par, x, y, U, ϕ, 0)
    draw_slices(par, x, y, U, ϕ, 0, 101, 101)
    Lx = Vector{Float64}() # Preallocate empty vector of interface position
    Ly = Vector{Float64}() # Preallocate empty vector of interface position
    push!(Lx, par.Lx/2 + par.β); push!(Ly, par.Ly/2 + par.β)
    # Time stepping
    for i = 1:par.Nt-1
        # 1. Find Ω, dΩ, and irregular grid points
        D = find_domain(par, ϕ)
        dΩ = find_interface(D, ϕ)
        # 2. Solve FKPP equation on Ω
        U = fkpp(D, U, ϕ, par, dx, dy, dt, i)
        # 3. Compute extension velocity field
        V = extend_velocity(D, dΩ, U, ϕ, par, dx, dy)
        # 4. Solve level-set equation
        ϕ = level_set(V, ϕ, par, dx, dy, dt)
        # 5. Re-initialise level-set function as a signed-distance
        if mod(i, 10) == 0
            ϕ = reinitialisation(ϕ, par, dx, dy, par.ϕ_Iterations)
        end
        # Optional: Post-processing
        if mod(i, 10) == 0
            draw_heat(par, x, y, U, V, ϕ, i)
            draw_slices(par, x, y, U, V, ϕ, i, 101, 101)
        end
        px, py = front_position(x, y, ϕ, 101, 101, dx, dy)
        push!(Lx, px); push!(Ly, py)
    end
    writedlm("Lx.csv", Lx)
    writedlm("Ly.csv", Ly)
    writedlm("t.csv", t)
end

@time fisher_stefan_2d()