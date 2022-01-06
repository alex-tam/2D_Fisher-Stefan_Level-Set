# Plot solutions to 2D Fisher-Stefan
# Alex Tam, 09/11/2021

using Plots
using Measures
using LaTeXStrings
using DelimitedFiles

"Plot solutions as 2D heat maps"
function draw_heat(x, y, U, V, ϕ, i, Lx, Ly)
    gr(); plot() # Load GR plotting backend and clear previous plots
    default(titlefont = (18, "Computer Modern"), guidefont = (26, "Computer Modern"), tickfont = (18, "Computer Modern"))
    heatmap(x,y,transpose(U), xlabel = L"$x$", ylabel = L"$y$", margin=3mm, aspect_ratio=:equal, xlims=(0,Lx), ylims =(0,Ly), tick_direction=:out, c=:plasma, clims=(0.0, 1.0))
    savefig("u-$i.pdf")
    heatmap(x,y,transpose(V), xlabel = L"$x$", ylabel = L"$y$", margin=3mm, aspect_ratio=:equal, xlims=(0,Lx), ylims =(0,Ly), tick_direction=:out, c=:plasma)
    savefig("V-$i.pdf")
    heatmap(x,y,transpose(ϕ), xlabel = L"$x$", ylabel = L"$y$", margin=3mm, aspect_ratio=:equal, xlims=(0,Lx), ylims =(0,Ly), tick_direction=:out, c=:plasma)
    savefig("phi-$i.pdf")
end

function draw_heat(x, y, U, ϕ, i, Lx, Ly)
    gr(); plot() # Load GR plotting backend and clear previous plots
    default(titlefont = (18, "Computer Modern"), guidefont = (26, "Computer Modern"), tickfont = (18, "Computer Modern"))
    heatmap(x,y,transpose(U), xlabel = L"$x$", ylabel = L"$y$", margin=3mm, aspect_ratio=:equal, xlims=(0,Lx), ylims =(0,Ly), tick_direction=:out, c=:plasma, clims=(0.0, 1.0))
    savefig("u-$i.pdf")
    heatmap(x,y,transpose(ϕ), xlabel = L"$x$", ylabel = L"$y$", margin=3mm, aspect_ratio=:equal, xlims=(0,Lx), ylims =(0,Ly), tick_direction=:out, c=:plasma)
    savefig("phi-$i.pdf")
end

"Density slice plots"
function draw_slices(x, y, nx, ny, Lx, Ly, plot_times)
    gr(); plot() # Load GR plotting backend and clear previous plots
    default(titlefont = (18, "Computer Modern"), guidefont = (26, "Computer Modern"), tickfont = (18, "Computer Modern"))
    U = readdlm("U-0.csv")
    plot(x, U[:,ny], xlabel = L"$x$", ylabel = L"$u(x,10,t)$", linecolor = :black, grid = false, margin=3mm, legend = false, xlims=(0,Lx), ylims=(0,1))
    for i in plot_times[2:2:end]
        u = readdlm("ux-$i.csv")
        plot!(x, u, linecolor = :black, linestyle = :dash)
    end
    # plot!([7.595174442,7.595174442], [0,1], linecolor = :red)
    # plot!([12.404825558,12.404825558], [0,1], linecolor = :red)
    savefig("u_slice_x.pdf")
    plot(y, U[nx,:], xlabel = L"$y$", ylabel = L"$u(10,y,t)$", linecolor = :black, grid = false, margin=3mm, legend = false, xlims=(0,Ly), ylims=(0,1))
    for i in plot_times[2:2:end]
        u = readdlm("uy-$i.csv")
        plot!(y, u, linecolor = :black, linestyle = :dash)
    end
    # plot!([7.595174442,7.595174442], [0,1], linecolor = :red)
    # plot!([12.404825558,12.404825558], [0,1], linecolor = :red)
    savefig("u_slice_y.pdf")
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
    draw_slices(x, y, nx, ny, Lx, Ly, plot_times)
    for i in plot_times
        U = readdlm("U-$i.csv")
        V = readdlm("V-$i.csv")
        ϕ = readdlm("Phi-$i.csv")
        draw_heat(x, y, U, V, ϕ, i, Lx, Ly)
    end
end

@time draw()
