# Compute extension velocity field for Fisher-Stefan
# Alex Tam, 09/11/2021

"Compute velocity extension by orthogonal extrapolation"
function extend_velocity(D, dΩ, U, ϕ, par, dx, dy)
    dτ = 0.5*(dx/5 + dy/5) # "Time" step for extrapolation PDE
    vi = interface_velocity(dΩ, U, ϕ, par, dx, dy) # Find speed at interface using Stefan condition
    V = velocity_init(D, dΩ, par, vi) # Obtain "initial" condition for V(x,y,t) using speed at interface
    v = build_vector(V, D) # Create vector using data at interior grid points
    # Solve PDE to extrapolate V outwards from interface
    prob = ODEProblem((du, u, p, t) -> vel_outward_rhs!(du, u, p, t, D, dΩ, vi, par, ϕ, dx, dy), v, (0, par.V_Iterations*dτ))
    sol = solve(prob, Tsit5(), reltol = 1e-3, abstol = 1e-6, saveat = par.V_Iterations*dτ)
    v = sol[:,end] # Update vector of V
    # # Solve PDE to extrapolate V inwards from interface
    prob = ODEProblem((du, u, p, t) -> vel_inward_rhs!(du, u, p, t, D, dΩ, vi, par, ϕ, dx, dy), v, (0, par.V_Iterations*dτ))
    sol = solve(prob, Tsit5(), reltol = 1e-3, abstol = 1e-6, saveat = par.V_Iterations*dτ)
    # Reshape solution to matrix
    V = build_v_matrix(sol[:,end], par, D)
    return V
end

"Apply Stefan condition to find speed at interface points"
function interface_velocity(dΩ, U, ϕ, par, dx, dy)
    vi = Vector{Float64}() # Pre-allocate empty vector of interface speeds
    for ip in dΩ # Loop over interface points
        # Compute du/dx, dϕ/dx, du/dy and dϕ/dy
        if ip.m.xInd != ip.p.xInd # If interface lies in x-direction
            # du/dx
            if (ip.m.Ω == true) && (ip.m.dΩ == false) # Left point lies inside Ω and not close to interface
                if ϕ[ip.m.xInd-1,ip.m.yInd] < 0 # Check adjacent grid point
                    um = U[ip.m.xInd-1,ip.m.yInd] # Set to grid value if inside Ω
                else
                    um = par.u_f # Set to interface value if outisde Ω
                end
                ux = (1/dx)*(ip.m.θxp*um/(ip.m.θxm*(ip.m.θxm+ip.m.θxp)) - (ip.m.θxm+ip.m.θxp)*U[ip.m.xInd,ip.m.yInd]/(ip.m.θxm*ip.m.θxp) + par.u_f*(2*ip.m.θxp+ip.m.θxm)/(ip.m.θxp*(ip.m.θxp+ip.m.θxm))) # Biased stencil
            elseif (ip.m.Ω == true) && (ip.m.dΩ == true) # Left point lies inside Ω and close to interface
                ϕm = ϕ[ip.m.xInd-1,ip.m.yInd]
                ϕmm = ϕ[ip.m.xInd-2,ip.m.yInd]
                if (ϕm < 0) && (ϕmm < 0) # Check two adjacent grid points
                    ux = 1/(2*dx)*(3*U[ip.m.xInd,ip.m.yInd] - 4*U[ip.m.xInd-1,ip.m.yInd] + U[ip.m.xInd-2,ip.m.yInd]) # Second-order
                elseif (ϕm < 0) && (ϕmm >= 0)
                    θ = ϕm/(ϕm - ϕmm)
                    if θ < par.θb # If opposite ghost node is close to interface
                        ux = 0.0
                    else
                        ux = (1/dx)*(par.u_f/(θ*(1+θ)) - (U[ip.m.xInd-1,ip.m.yInd]*(1+θ))/θ + (par.u_f*(2+θ))/(1+θ)) # Second-order difference
                    end
                else
                    ux = 0.0
                end
            elseif (ip.p.Ω == true) && (ip.p.dΩ == false) # Right point lies inside Ω and not close to interface
                if ϕ[ip.p.xInd+1,ip.p.yInd] < 0 # Check adjacent grid point
                    up = U[ip.p.xInd+1,ip.p.yInd] # Set to grid value if inside Ω
                else
                    up = par.u_f # Set to grid value if outside Ω
                end
                ux = (1/dx)*(-par.u_f*(2*ip.p.θxm + ip.p.θxp)/(ip.p.θxm*(ip.p.θxm+ip.p.θxp)) + (ip.p.θxm+ip.p.θxp)*U[ip.p.xInd,ip.p.yInd]/(ip.p.θxm*ip.p.θxp) - ip.p.θxm*up/(ip.p.θxp*(ip.p.θxm+ip.p.θxp))) # Biased stencil
            else # Right point lies inside Ω and close to interface
                ϕp = ϕ[ip.p.xInd+1,ip.p.yInd]
                ϕpp = ϕ[ip.p.xInd+2,ip.p.yInd]
                if (ϕp < 0) && (ϕpp < 0) # Check two adjacent grid points
                    ux = 1/(2*dx)*(-U[ip.p.xInd+2,ip.p.yInd] + 4*U[ip.p.xInd+1,ip.p.yInd] - 3*U[ip.p.xInd,ip.p.yInd]) # Second-order
                elseif (ϕp < 0) && (ϕpp >= 0)
                    θ = ϕp/(ϕp - ϕpp)
                    if θ < par.θb # If opposite ghost node is close to interface
                        ux = 0.0
                    else
                        ux = (1/dx)*(-(par.u_f*(2+θ))/(1+θ) + (U[ip.m.xInd+1,ip.m.yInd]*(1+θ))/θ - par.u_f/(θ*(1+θ))) # Second-order difference
                    end
                else
                    ux = 0.0
                end
            end
            # dϕ/dx
            ϕxl = (ϕ[ip.m.xInd+1, ip.m.yInd]-ϕ[ip.m.xInd-1, ip.m.yInd])/(2*dx)
            ϕxr = (ϕ[ip.p.xInd+1, ip.p.yInd]-ϕ[ip.p.xInd-1, ip.p.yInd])/(2*dx)
            ϕx = ip.m.θxp*ϕxr + (1-ip.m.θxp)*ϕxl
            # dϕ/dy
            ϕp = ip.m.θxp*ϕ[ip.p.xInd, ip.p.yInd+1] + (1-ip.m.θxp)*ϕ[ip.m.xInd, ip.m.yInd+1]
            ϕm = ip.m.θxp*ϕ[ip.p.xInd, ip.p.yInd-1] + (1-ip.m.θxp)*ϕ[ip.m.xInd, ip.m.yInd-1]
            ϕy = (ϕp-ϕm)/(2*dy) # Use ghost point method
            # du/dy
            if ϕp < 0 # Positive ghost point inside Ω
                up =  ip.m.θxp*U[ip.p.xInd, ip.p.yInd+1] + (1-ip.m.θxp)*U[ip.m.xInd, ip.m.yInd+1]
                ϕpp =  ip.m.θxp*ϕ[ip.p.xInd, ip.p.yInd+2] + (1-ip.m.θxp)*ϕ[ip.m.xInd, ip.m.yInd+2] # Second ghost point
                if ϕpp < 0
                    upp =  ip.m.θxp*U[ip.p.xInd, ip.p.yInd+2] + (1-ip.m.θxp)*U[ip.m.xInd, ip.m.yInd+2]
                    uy = (-upp + 4*up - 3*par.u_f)/(2*dy) # Second-order difference
                else
                    θ = ϕp/(ϕp - ϕpp) # Distance between middle ghost point and interface
                    if θ < par.θb
                        uy = 0.0 # Set to zero if ghost point close to interface
                    else
                        uy = ((-(2+θ)*par.u_f)/(1+θ) + (1+θ)*up/θ - par.u_f/(θ*(1+θ)))/dy # Second-order difference
                    end
                end
            elseif ϕm < 0 # Negative ghost point inside Ω
                um =  ip.m.θxp*U[ip.p.xInd, ip.p.yInd-1] + (1-ip.m.θxp)*U[ip.m.xInd, ip.m.yInd-1]
                ϕmm = ip.m.θxp*ϕ[ip.p.xInd, ip.p.yInd-2] + (1-ip.m.θxp)*ϕ[ip.m.xInd, ip.m.yInd-2]
                if ϕmm < 0
                    umm = ip.m.θxp*U[ip.p.xInd, ip.p.yInd-2] + (1-ip.m.θxp)*U[ip.m.xInd, ip.m.yInd-2]
                    uy = (3par.u_f - 4*um + umm)/(2*dy) # Second-order difference
                else
                    θ = ϕm/(ϕm - ϕmm) # Distance between middle ghost point and interface
                    if θ < par.θb
                        uy = 0.0 # Set to zero if ghost point close to interface
                    else
                        uy = ((2+θ)*par.u_f/(1+θ) - (1+θ)*um/θ + par.u_f/(θ*(θ+1)))/dy # Second-order difference
                    end
                end
            else # Both ghost points outside Ω
                uy = 0.0
            end
        else # If interface lies in y-direction
            # du/dx and du/dy
            if (ip.m.Ω == true) && (ip.m.dΩ == false) # Bottom point lies inside Ω and not close to interface
                if ϕ[ip.m.xInd,ip.m.yInd-1] < 0 # Check adjacent grid point
                    ub = U[ip.m.xInd,ip.m.yInd-1] # Set to grid value if inside Ω
                else
                    ub = par.u_f # Set to interface value if outside Ω
                end
                uy = (1/dy)*(ip.m.θyp*ub/(ip.m.θym*(ip.m.θym+ip.m.θyp)) - (ip.m.θym+ip.m.θyp)*U[ip.m.xInd,ip.m.yInd]/(ip.m.θym*ip.m.θyp) + par.u_f*(2*ip.m.θyp+ip.m.θym)/(ip.m.θyp*(ip.m.θyp+ip.m.θym)) ) # Biased stencil
            elseif (ip.m.Ω == true) && (ip.m.dΩ == true) # Bottom point lies inside Ω and close to interface
                ϕm = ϕ[ip.m.xInd,ip.m.yInd-1]
                ϕmm = ϕ[ip.m.xInd,ip.m.yInd-2]
                if (ϕm < 0) && (ϕmm < 0) # Check two adjacent points
                    uy = 1/(2*dy)*(3*U[ip.m.xInd,ip.m.yInd] - 4*U[ip.m.xInd,ip.m.yInd-1] + U[ip.m.xInd,ip.m.yInd-2]) # Second-order
                elseif (ϕm < 0) && (ϕmm >= 0)
                    θ = ϕm/(ϕm - ϕmm)
                    if θ < par.θb # If opposite ghost node is close to interface
                        uy = 0.0
                    else
                        uy = (1/dy)*(par.u_f/(θ*(1+θ)) - (U[ip.m.xInd,ip.m.yInd-1]*(1+θ))/θ + (par.u_f*(2+θ))/(1+θ)) # Second-order difference
                    end
                else
                    uy = 0.0
                end
            elseif (ip.p.Ω == true) && (ip.p.dΩ == false) # Top point lies inside Ω and not close to interface
                if ϕ[ip.p.xInd,ip.p.yInd+1] < 0 # Check adjacent grid point
                    ut = U[ip.p.xInd,ip.p.yInd+1] # Set to grid value if inside Ω
                else
                    ut = par.u_f # Set to interface value if outside Ω
                end
                uy = (1/dy)*(-par.u_f*(2*ip.p.θym + ip.p.θyp)/(ip.p.θym*(ip.p.θym+ip.p.θyp)) + (ip.p.θym+ip.p.θyp)*U[ip.p.xInd, ip.p.yInd]/(ip.p.θym*ip.p.θyp) - ip.p.θym*ut/(ip.p.θyp*(ip.p.θym+ip.p.θyp)) ) # Biased stencil
            else
                ϕp = ϕ[ip.m.xInd,ip.m.yInd+1]
                ϕpp = ϕ[ip.m.xInd,ip.m.yInd+2]
                if (ϕp < 0) && (ϕpp < 0) # Check two adjacent nodes
                    uy = 1/(2*dy)*(-U[ip.p.xInd,ip.p.yInd+2] + 4*U[ip.p.xInd,ip.p.yInd+1] - 3*U[ip.p.xInd,ip.p.yInd]) # Second-order
                elseif (ϕp < 0) && (ϕpp >= 0)
                    θ = ϕp/(ϕp - ϕpp)
                    if θ < par.θb # If opposite ghost node is close to interface
                        uy = 0.0
                    else
                        uy = (1/dy)*(-(par.u_f*(2+θ))/(1+θ) + (U[ip.m.xInd,ip.m.yInd+1]*(1+θ))/θ - par.u_f/(θ*(1+θ))) # Second-order difference
                    end 
                else
                    uy = 0.0
                end
            end
            # dϕ/dx
            ϕp = ip.m.θyp*ϕ[ip.p.xInd+1, ip.p.yInd] + (1-ip.m.θyp)*ϕ[ip.m.xInd+1, ip.m.yInd]
            ϕm = ip.m.θyp*ϕ[ip.p.xInd-1, ip.p.yInd] + (1-ip.m.θyp)*ϕ[ip.m.xInd-1, ip.m.yInd]
            ϕx = (ϕp-ϕm)/(2*dx) # Use ghost point method         
            # dϕ/dy
            ϕyb = (ϕ[ip.m.xInd, ip.m.yInd+1]-ϕ[ip.m.xInd, ip.m.yInd-1])/(2*dy)
            ϕyt = (ϕ[ip.p.xInd, ip.p.yInd+1]-ϕ[ip.p.xInd, ip.p.yInd-1])/(2*dy)
            ϕy = ip.m.θyp*ϕyt + (1-ip.m.θyp)*ϕyb
            # du/dx
            if ϕp < 0 # Positive ghost point inside Ω
                up = ip.m.θyp*U[ip.p.xInd+1, ip.p.yInd] + (1-ip.m.θyp)*U[ip.m.xInd+1, ip.m.yInd]
                ϕpp = ip.m.θyp*ϕ[ip.p.xInd+2, ip.p.yInd] + (1-ip.m.θyp)*ϕ[ip.m.xInd+2, ip.m.yInd]
                if ϕpp < 0
                    upp = ip.m.θyp*U[ip.m.xInd+2, ip.p.yInd] + (1-ip.m.θyp)*U[ip.m.xInd+2, ip.m.yInd]
                    ux = (-upp + 4*up -3*par.u_f)/(2*dx) # Second-order
                else
                    θ = ϕp/(ϕp - ϕpp) # Distance between middle ghost point and interface
                    if θ < par.θb
                        ux = 0.0 # Set to zero if ghost point close to interface
                    else
                        ux = ((-(2+θ)*par.u_f)/(1+θ) + (1+θ)*up/θ - par.u_f/(θ*(1+θ)))/dx # Second-order
                    end
                end
            elseif ϕm < 0 # Negative ghost point inside Ω
                um = ip.m.θyp*U[ip.m.xInd-1, ip.p.yInd] + (1-ip.m.θyp)*U[ip.m.xInd-1, ip.m.yInd]
                ϕmm = ip.m.θyp*ϕ[ip.p.xInd-2, ip.p.yInd] + (1-ip.m.θyp)*ϕ[ip.m.xInd-2, ip.m.yInd]
                if ϕmm < 0
                    umm = ip.m.θyp*U[ip.p.xInd-2, ip.p.yInd] + (1-ip.m.θyp)*U[ip.m.xInd-2, ip.m.yInd]
                    ux = (3*par.u_f - 4um + umm)/(2*dx) # Second-order
                else
                    θ = ϕm/(ϕm - ϕmm) # Distance between middle ghost point and interface
                    if θ < par.θb
                        ux = 0.0 # Set to zero if ghost point close to interface
                    else
                        ux = ((2+θ)*par.u_f/(1+θ) - (1+θ)*um/θ + par.u_f/(θ*(θ+1)))/dx # Second-order
                    end
                end
            else # Both ghost points outside Ω
                ux = 0.0
            end
        end
        push!(vi, -par.κ*(ux*ϕx + uy*ϕy)) # Append speed at current point to vector
    end
    return vi
end

"Obtain initial condition for V(x,y,t) using speed at interface"
function velocity_init(D, dΩ, par, vi)
    V = zeros(par.Nx, par.Ny) # Pre-allocate matrix of V(x,y,t)
    # Manually overwrite speed close to interface
    for gp in D # Loop over interior grid points
        if (gp.dΩ == true) # If grid point is close to interface
            # Locate closest interface point to current grid point
            if (gp.θxp < gp.θxm) && (gp.θxp < gp.θyp) && (gp.θxp < gp.θym) # Closest point in positive x-direction
                xm = gp.xInd; xp = gp.xInd+1; ym = gp.yInd; yp = gp.yInd
            elseif (gp.θxm < gp.θxp) && (gp.θxm < gp.θyp) && (gp.θxm < gp.θym) # Closest point in negative x-direction
                xm = gp.xInd-1; xp = gp.xInd; ym = gp.yInd; yp = gp.yInd
            elseif (gp.θyp < gp.θxm) && (gp.θyp < gp.θxp) && (gp.θyp < gp.θym) # Closest point in positive y-direction
                xm = gp.xInd; xp = gp.xInd; ym = gp.yInd; yp = gp.yInd+1
            else # Closest point in negative y-direction
                xm = gp.xInd; xp = gp.xInd; ym = gp.yInd-1; yp = gp.yInd
            end
            # Assign V to interface speed at closest point
            for i = 1:length(dΩ)
                if (xm == dΩ[i].m.xInd) && (ym == dΩ[i].m.yInd) && (xp == dΩ[i].p.xInd) && (yp == dΩ[i].p.yInd)
                    V[gp.xInd, gp.yInd] = vi[i]
                end
            end
        end
    end
    return V
end

"Obtain velocity at a specified interface point"
function obtain_interface_velocity(xm, ym, xp, yp, dΩ, vi)
    for i = 1:length(dΩ) # Loop over candidate points
        if (xm == dΩ[i].m.xInd) && (ym == dΩ[i].m.yInd) && (xp == dΩ[i].p.xInd) && (yp == dΩ[i].p.yInd) # If all indices are correct
            return vi[i]
        end
    end
end

"Construct right-hand vector of outward velocity extension PDE, for use in DifferentialEquations.jl"
function vel_outward_rhs!(du, u, p, t, D, dΩ, vi, par, ϕ, dx, dy)
    V = build_v_matrix(u, par, D) # Generate matrix
    for i = 1:length(D) # Loop over interior grid points
        if (D[i].Ω == true) || (D[i].dΩ == true) # Grid points in Ω or close to interface
            du[i] = 0.0 # Set RHS to zero
        else
            # Compute advection coefficients
            ϕx = 1/(2*dx)*(ϕ[D[i].xInd+1, D[i].yInd]-ϕ[D[i].xInd-1, D[i].yInd]) # dϕ/dx
            ϕy = 1/(2*dy)*(ϕ[D[i].xInd, D[i].yInd+1]-ϕ[D[i].xInd, D[i].yInd-1]) # dϕ/dy
            n = sqrt(ϕx^2 + ϕy^2) # Normalisation
            if n == 0 # Prevent division by zero
                a = 0; b = 0
            else
                a = ϕx/n # x-component
                b = ϕy/n # y-component
            end
            # Discretise dV/dx (first-order upwind)
            if a < 0 # Forward difference
                if D[i].θxp == 1.0 # Regular discretisation
                    Vx = (V[D[i].xInd+1, D[i].yInd] - V[D[i].xInd, D[i].yInd])/dx
                else # Irregular discretisation
                    v = obtain_interface_velocity(D[i].xInd, D[i].yInd, D[i].xInd+1, D[i].yInd, dΩ, vi)
                    Vx = (v - V[D[i].xInd, D[i].yInd])/(D[i].θxp*dx)
                end
            else # Backward difference
                if D[i].θxm == 1.0 # Regular discretisation
                    Vx = (V[D[i].xInd, D[i].yInd] - V[D[i].xInd-1, D[i].yInd])/dx
                else # Irregular discretisation
                    v = obtain_interface_velocity(D[i].xInd-1, D[i].yInd, D[i].xInd, D[i].yInd, dΩ, vi)
                    Vx = (V[D[i].xInd, D[i].yInd]-v)/(D[i].θxm*dx)
                end
            end
            # Discretise dV/dy (first-order upwind)
            if b < 0 # Forward difference
                if D[i].θyp == 1.0 # Regular discretisation
                    Vy = (V[D[i].xInd, D[i].yInd+1] - V[D[i].xInd, D[i].yInd])/dy
                else # Irregular discretisation
                    v = obtain_interface_velocity(D[i].xInd, D[i].yInd, D[i].xInd, D[i].yInd+1, dΩ, vi)
                    Vy = (v - V[D[i].xInd, D[i].yInd])/(D[i].θyp*dy)
                end
            else # Backward difference
                if D[i].θym == 1.0
                    Vy = (V[D[i].xInd, D[i].yInd] - V[D[i].xInd, D[i].yInd-1])/dy
                else
                    # Obtain interface velocity
                    v = obtain_interface_velocity(D[i].xInd, D[i].yInd-1, D[i].xInd, D[i].yInd, dΩ, vi)
                    Vy = (V[D[i].xInd, D[i].yInd]-v)/(D[i].θym*dy)
                end
            end
            du[i] = -a*Vx - b*Vy # Update RHS term
        end
    end
end

"Construct right-hand vector of inward velocity extension PDE, for use in DifferentialEquations.jl"
function vel_inward_rhs!(du, u, p, t, D, dΩ, vi, par, ϕ, dx, dy)
    V = build_v_matrix(u, par, D) # Generate matrix
    for i = 1:length(D) # Loop over interior grid points
        if (D[i].Ω == false) || (D[i].dΩ == true) # Grid points outside Ω or close to interface
            du[i] = 0.0 # Set RHS to zero
        else
            # Compute advection coefficients
            ϕx = 1/(2*dx)*(ϕ[D[i].xInd+1, D[i].yInd]-ϕ[D[i].xInd-1, D[i].yInd]) # dϕ/dx
            ϕy = 1/(2*dy)*(ϕ[D[i].xInd, D[i].yInd+1]-ϕ[D[i].xInd, D[i].yInd-1]) # dϕ/dy
            n = sqrt(ϕx^2 + ϕy^2) # Normalisation
            if n == 0 # Prevent division by zero
                a = 0; b = 0
            else
                a = -ϕx/n # x-component
                b = -ϕy/n # y-component
            end
            # Discretise dV/dx (first-order upwind)
            if a < 0 # Forward difference
                if D[i].θxp == 1.0 # Regular discretisation
                    Vx = (V[D[i].xInd+1, D[i].yInd] - V[D[i].xInd, D[i].yInd])/dx
                else # Irregular discretisation
                    v = obtain_interface_velocity(D[i].xInd, D[i].yInd, D[i].xInd+1, D[i].yInd, dΩ, vi)
                    Vx = (v - V[D[i].xInd, D[i].yInd])/(D[i].θxp*dx)
                end
            else # Backward difference
                if D[i].θxm == 1.0 # Regular discretisation
                    Vx = (V[D[i].xInd, D[i].yInd] - V[D[i].xInd-1, D[i].yInd])/dx
                else # Irregular discretisation
                    v = obtain_interface_velocity(D[i].xInd-1, D[i].yInd, D[i].xInd, D[i].yInd, dΩ, vi)
                    Vx = (V[D[i].xInd, D[i].yInd]-v)/(D[i].θxm*dx)
                end
            end
            # Discretise dV/dy (first-order upwind)
            if b < 0 # Forward difference
                if D[i].θyp == 1.0 # Regular discretisation
                    Vy = (V[D[i].xInd, D[i].yInd+1] - V[D[i].xInd, D[i].yInd])/dy
                else # Irregular discretisation
                    v = obtain_interface_velocity(D[i].xInd, D[i].yInd, D[i].xInd, D[i].yInd+1, dΩ, vi)
                    Vy = (v - V[D[i].xInd, D[i].yInd])/(D[i].θyp*dy)
                end
            else # Backward difference
                if D[i].θym == 1.0
                    Vy = (V[D[i].xInd, D[i].yInd] - V[D[i].xInd, D[i].yInd-1])/dy
                else
                    # Obtain interface velocity
                    v = obtain_interface_velocity(D[i].xInd, D[i].yInd-1, D[i].xInd, D[i].yInd, dΩ, vi)
                    Vy = (V[D[i].xInd, D[i].yInd]-v)/(D[i].θym*dy)
                end
            end
            du[i] = -a*Vx - b*Vy # Update RHS term
        end
    end
end