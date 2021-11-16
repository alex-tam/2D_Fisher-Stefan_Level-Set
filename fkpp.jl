# Solve 2D Fisher-KPP equation with Dirichlet boundary conditions
# Alex Tam, 03/11/2021

"Solve Fisher-KPP equation"
function fkpp(D, U, ϕ, par, dx, dy, dt, i)
    # Manually assign value at grid points close to interface
    for gp in D
        if (gp.Ω == true) && (gp.dΩ == true)
            U[gp.xInd, gp.yInd] = par.u_f
        end
    end
    # Create vector from matrix data
    u = build_vector(U, D)
    # Construct and solve ODE problem using DifferentialEquations.jl
    prob = ODEProblem((du, u, p, t) -> fkpp_rhs!(du, u, p, t, D, ϕ, par, dx, dy), u, ((i-1)*dt, i*dt))
    sol = solve(prob, Tsit5(), reltol = 1e-3, abstol = 1e-6, saveat = i*dt)
    # Reshape solution to matrix
    U = build_u_matrix(sol[:,end], par, D)
    return U
end

"Construct right-hand vector for use in DifferentialEquations.jl"
function fkpp_rhs!(du, u, p, t, D, ϕ, par, dx, dy)
    U = build_u_matrix(u, par, D) # Generate matrix
    for i = 1:length(D) # Loop over interior grid points
        if (D[i].Ω == true) && (D[i].dΩ == false) # Grid points inside Ω, away from dΩ
            # Obtain density at stencil points
            u_mid = U[D[i].xInd, D[i].yInd]
            if ϕ[D[i].xInd-1, D[i].yInd] >= 0
                u_left = par.u_f
            else
                u_left = U[D[i].xInd-1, D[i].yInd]
            end
            if ϕ[D[i].xInd+1, D[i].yInd] >= 0
                u_right = par.u_f
            else
                u_right = U[D[i].xInd+1, D[i].yInd]
            end
            if ϕ[D[i].xInd, D[i].yInd-1] >= 0
                u_bottom = par.u_f
            else
                u_bottom = U[D[i].xInd, D[i].yInd-1]
            end
            if ϕ[D[i].xInd, D[i].yInd+1] >= 0
                u_top = par.u_f
            else
                u_top = U[D[i].xInd, D[i].yInd+1]
            end
            # Compute Laplacian and source term
            uxx = 2.0*( u_left/(D[i].θxm*(D[i].θxm+D[i].θxp)) - u_mid/(D[i].θxm*D[i].θxp) + u_right/(D[i].θxp*(D[i].θxm+D[i].θxp)) )/(dx^2)
            uyy = 2.0*( u_bottom/(D[i].θym*(D[i].θym+D[i].θyp)) - u_mid/(D[i].θym*D[i].θyp) + u_top/(D[i].θyp*(D[i].θym+D[i].θyp)) )/(dy^2)
            du[i] = par.D*(uxx + uyy) + par.λ*u_mid*(1-u_mid)
        else # Grid points outside Ω or close to dΩ
            du[i] = 0.0
        end
    end
end