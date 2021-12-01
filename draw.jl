# Plot solutions to 2D Fisher-Stefan
# Alex Tam, 09/11/2021

using Plots
using LaTeXStrings
using DelimitedFiles

"Plot solutions as 2D heat maps"
function draw_heat(x, y, U, V, ϕ, i, Lx, Ly)
    gr(); plot() # Load GR plotting backend and clear previous plots
    default(titlefont = (18, "Computer Modern"), guidefont = (26, "Computer Modern"), tickfont = (18, "Computer Modern"))
    heatmap(x,y,transpose(U), xlabel = L"$x$", ylabel = L"$y$", aspect_ratio=:equal, xlims=(0,Lx), ylims =(0,Ly), tick_direction=:out, c=:plasma, clims=(0.0, 0.5))
    savefig("u-$i.pdf")
    heatmap(x,y,transpose(V), xlabel = L"$x$", ylabel = L"$y$", aspect_ratio=:equal, xlims=(0,Lx), ylims =(0,Ly), tick_direction=:out, c=:plasma)
    savefig("V-$i.pdf")
    heatmap(x,y,transpose(ϕ), xlabel = L"$x$", ylabel = L"$y$", aspect_ratio=:equal, xlims=(0,Lx), ylims =(0,Ly), tick_direction=:out, c=:plasma)
    savefig("phi-$i.pdf")
end

function draw_heat(x, y, U, ϕ, i, Lx, Ly)
    gr(); plot() # Load GR plotting backend and clear previous plots
    default(titlefont = (18, "Computer Modern"), guidefont = (26, "Computer Modern"), tickfont = (18, "Computer Modern"))
    heatmap(x,y,transpose(U), xlabel = L"$x$", ylabel = L"$y$", aspect_ratio=:equal, xlims=(0,Lx), ylims =(0,Ly), tick_direction=:out, c=:plasma, clims=(0.0, 0.5))
    savefig("u-$i.pdf")
    heatmap(x,y,transpose(ϕ), xlabel = L"$x$", ylabel = L"$y$", aspect_ratio=:equal, xlims=(0,Lx), ylims =(0,Ly), tick_direction=:out, c=:plasma)
    savefig("phi-$i.pdf")
end

"Plot solutions as slices"
function draw_slices(x, y, U, V, ϕ, i, nx, ny, Lx, Ly)
    gr(); plot() # Load GR plotting backend and clear previous plots
    default(titlefont = (18, "Computer Modern"), guidefont = (26, "Computer Modern"), tickfont = (18, "Computer Modern"))
    plot(x, U[:,ny], xlabel = L"$x$", ylabel = L"$u(x,25,t)$", legend = false, xlims=(0,Lx))
    savefig("u_slice_x-$i.pdf")
    plot(x, V[:,ny], xlabel = L"$x$", ylabel = L"$V(x,25,t)$", legend = false, xlims=(0,Lx))
    savefig("V_slice_x-$i.pdf")
    plot(x, ϕ[:,ny], xlabel = L"$x$", ylabel = L"$\phi(x,25,t)$", legend = false, xlims=(0,Lx))
    savefig("phi_slice_x-$i.pdf")
    plot(y, U[nx,:], xlabel = L"$y$", ylabel = L"$u(25,y,t)$", legend = false, xlims=(0,Ly))
    savefig("u_slice_y-$i.pdf")
    plot(y, V[nx,:], xlabel = L"$y$", ylabel = L"$V(25,y,t)$", legend = false, xlims=(0,Ly))
    savefig("V_slice_y-$i.pdf")
    plot(y, ϕ[nx,:], xlabel = L"$y$", ylabel = L"$\phi(25,y,t)$", legend = false, xlims=(0,Ly))
    savefig("phi_slice_y-$i.pdf")
    writedlm("x.csv", x)
    writedlm("y.csv", y)
    writedlm("ux-$i.csv", U[:,ny])
    writedlm("uy-$i.csv", U[nx,:])
end

function draw_slices(x, y, U, ϕ, i, nx, ny, Lx, Ly)
    gr(); plot() # Load GR plotting backend and clear previous plots
    default(titlefont = (18, "Computer Modern"), guidefont = (26, "Computer Modern"), tickfont = (18, "Computer Modern"))
    plot(x, U[:,ny], xlabel = L"$x$", ylabel = L"$u(x,25,t)$", legend = false, xlims=(0,Lx))
    savefig("u_slice_x-$i.pdf")
    plot(x, ϕ[:,ny], xlabel = L"$x$", ylabel = L"$\phi(x,25,t)$", legend = false, xlims=(0,Lx))
    savefig("phi_slice_x-$i.pdf")
    plot(y, U[nx,:], xlabel = L"$y$", ylabel = L"$u(25,y,t)$", legend = false, xlims=(0,Ly))
    savefig("u_slice_y-$i.pdf")
    plot(y, ϕ[nx,:], xlabel = L"$y$", ylabel = L"$\phi(25,y,t)$", legend = false, xlims=(0,Ly))
    savefig("phi_slice_y-$i.pdf")
end

# Control function for plotting
function draw()
    plot_times = convert(Vector{Int}, vec(readdlm("plot_times.csv")))
    x = vec(readdlm("x.csv"))
    y = vec(readdlm("y.csv"))
    nx::Int = (length(x)-1)/2; ny::Int = (length(y)-1)/2
    Lx = maximum(x); Ly = maximum(y)
    U = readdlm("U-0.csv")
    ϕ = readdlm("Phi-0.csv")
    draw_heat(x, y, U, ϕ, 0, Lx, Ly)
    draw_slices(x, y, U, ϕ, 0, nx, ny, Lx, Ly)
    for i in plot_times
        U = readdlm("U-$i.csv")
        V = readdlm("V-$i.csv")
        ϕ = readdlm("Phi-$i.csv")
        draw_heat(x, y, U, V, ϕ, i, Lx, Ly)
        draw_slices(x, y, U, V, ϕ, i, nx, ny, Lx, Ly)
    end
end

@time draw()