# Plot solutions to 2D Fisher-Stefan
# Alex Tam, 09/11/2021

"Plot solutions as 2D heat maps"
function draw_heat(par, x, y, U, V, ϕ, i)
    gr(); plot() # Load GR plotting backend and clear previous plots
    default(titlefont = (18, "Computer Modern"), guidefont = (26, "Computer Modern"), tickfont = (18, "Computer Modern"))
    heatmap(x,y,transpose(U), xlabel = L"$x$", ylabel = L"$y$", aspect_ratio=:equal, xlims=(0,par.Lx), ylims =(0,par.Ly), tick_direction=:out, c=:plasma, clims=(0.0, 1.0))
    savefig("u-$i.pdf")
    heatmap(x,y,transpose(V), xlabel = L"$x$", ylabel = L"$y$", aspect_ratio=:equal, xlims=(0,par.Lx), ylims =(0,par.Ly), tick_direction=:out, c=:plasma)
    savefig("V-$i.pdf")
    heatmap(x,y,transpose(ϕ), xlabel = L"$x$", ylabel = L"$y$", aspect_ratio=:equal, xlims=(0,par.Lx), ylims =(0,par.Ly), tick_direction=:out, c=:plasma)
    savefig("phi-$i.pdf")
end

function draw_heat(par, x, y, U, ϕ, i)
    gr(); plot() # Load GR plotting backend and clear previous plots
    default(titlefont = (18, "Computer Modern"), guidefont = (26, "Computer Modern"), tickfont = (18, "Computer Modern"))
    heatmap(x,y,transpose(U), xlabel = L"$x$", ylabel = L"$y$", aspect_ratio=:equal, xlims=(0,par.Lx), ylims =(0,par.Ly), tick_direction=:out, c=:plasma, clims=(0.0,1.0))
    savefig("u-$i.pdf")
    heatmap(x,y,transpose(ϕ), xlabel = L"$x$", ylabel = L"$y$", aspect_ratio=:equal, xlims=(0,par.Lx), ylims =(0,par.Ly), tick_direction=:out, c=:plasma)
    savefig("phi-$i.pdf")
end

# Edit slice plots to correctly identify the contact line
"Plot solutions as slices"
function draw_slices(par, x, y, U, V, ϕ, i, nx, ny)
    gr(); plot() # Load GR plotting backend and clear previous plots
    default(titlefont = (18, "Computer Modern"), guidefont = (26, "Computer Modern"), tickfont = (18, "Computer Modern"))
    plot(x, U[:,ny], xlabel = L"$x$", ylabel = L"$u(x,25,t)$", legend = false, xlims=(0,par.Lx))
    savefig("u_slice_x-$i.pdf")
    plot(x, V[:,ny], xlabel = L"$x$", ylabel = L"$V(x,25,t)$", legend = false, xlims=(0,par.Lx))
    savefig("V_slice_x-$i.pdf")
    plot(x, ϕ[:,ny], xlabel = L"$x$", ylabel = L"$\phi(x,25,t)$", legend = false, xlims=(0,par.Lx))
    savefig("phi_slice_x-$i.pdf")
    plot(y, U[nx,:], xlabel = L"$y$", ylabel = L"$u(25,y,t)$", legend = false, xlims=(0,par.Ly))
    savefig("u_slice_y-$i.pdf")
    plot(y, V[nx,:], xlabel = L"$y$", ylabel = L"$V(25,y,t)$", legend = false, xlims=(0,par.Ly))
    savefig("V_slice_y-$i.pdf")
    plot(y, ϕ[nx,:], xlabel = L"$y$", ylabel = L"$\phi(25,y,t)$", legend = false, xlims=(0,par.Ly))
    savefig("phi_slice_y-$i.pdf")
    writedlm("x.csv", x)
    writedlm("y.csv", y)
    writedlm("ux-$i.csv", U[:,ny])
    writedlm("uy-$i.csv", U[nx,:])
end

function draw_slices(par, x, y, U, ϕ, i, nx, ny)
    gr(); plot() # Load GR plotting backend and clear previous plots
    default(titlefont = (18, "Computer Modern"), guidefont = (26, "Computer Modern"), tickfont = (18, "Computer Modern"))
    plot(x, U[:,ny], xlabel = L"$x$", ylabel = L"$u(x,25,t)$", legend = false, xlims=(0,par.Lx))
    savefig("u_slice_x-$i.pdf")
    plot(x, ϕ[:,ny], xlabel = L"$x$", ylabel = L"$\phi(x,25,t)$", legend = false, xlims=(0,par.Lx))
    savefig("phi_slice_x-$i.pdf")
    plot(y, U[nx,:], xlabel = L"$y$", ylabel = L"$u(25,y,t)$", legend = false, xlims=(0,par.Ly))
    savefig("u_slice_y-$i.pdf")
    plot(y, ϕ[nx,:], xlabel = L"$y$", ylabel = L"$\phi(25,y,t)$", legend = false, xlims=(0,par.Ly))
    savefig("phi_slice_y-$i.pdf")
end

function front_position(x, y, ϕ, nx, ny, dx, dy)
    L = Vector{Float64}() # Initialise
    # Find front position using x-direction slice
    ϕv = ϕ[:,ny]
    for i = 1:length(ϕv)
        if (ϕv[i] < 0) && (ϕv[i+1] >= 0)
            θ = ϕv[i]/(ϕv[i] - ϕv[i+1])
            push!(L, x[i] + θ*dx)
        end
    end
    # Find front position using y-direction slice
    ϕv = ϕ[nx,:]
    for i = 1:length(ϕv)
        if (ϕv[i] < 0) && (ϕv[i+1] >= 0)
            θ = ϕv[i]/(ϕv[i] - ϕv[i+1])
            push!(L, y[i] + θ*dy)
        end
    end
    return L[1], L[2]
end