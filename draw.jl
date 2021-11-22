# Plot solutions to 2D Fisher-Stefan
# Alex Tam, 09/11/2021

"Plot solutions as 2D heat maps"
function draw_heat(par, x, y, U, V, ϕ, i)
    gr(); plot() # Load GR plotting backend and clear previous plots
    heatmap(x,y,transpose(U), xlabel = L"$x$", ylabel = L"$y$", aspect_ratio=:equal, fontfamily="Computer Modern", xlims=(0,par.Lx), ylims =(0,par.Ly), tick_direction=:out)
    savefig("u-$i.pdf")
    heatmap(x,y,transpose(V), xlabel = L"$x$", ylabel = L"$y$", aspect_ratio=:equal, fontfamily="Computer Modern", xlims=(0,par.Lx), ylims =(0,par.Ly), tick_direction=:out)
    savefig("V-$i.pdf")
    heatmap(x,y,transpose(ϕ), xlabel = L"$x$", ylabel = L"$y$", aspect_ratio=:equal, fontfamily="Computer Modern", xlims=(0,par.Lx), ylims =(0,par.Ly), tick_direction=:out)
    savefig("phi-$i.pdf")
end

function draw_heat(par, x, y, U, ϕ, i)
    gr(); plot() # Load GR plotting backend and clear previous plots
    heatmap(x,y,transpose(U), xlabel = L"$x$", ylabel = L"$y$", aspect_ratio=:equal, fontfamily="Computer Modern", xlims=(0,par.Lx), ylims =(0,par.Ly), tick_direction=:out)
    savefig("u-$i.pdf")
    heatmap(x,y,transpose(ϕ), xlabel = L"$x$", ylabel = L"$y$", aspect_ratio=:equal, fontfamily="Computer Modern", xlims=(0,par.Lx), ylims =(0,par.Ly), tick_direction=:out)
    savefig("phi-$i.pdf")
end

# Edit slice plots to correctly identify the contact line
"Plot solutions as slices"
function draw_slice(par, x, y, U, V, ϕ, i, n)
    gr(); plot() # Load GR plotting backend and clear previous plots
    plot(x, U[:,n], xlabel = L"$x$", ylabel = L"$u(x,25,t)$", legend = false, fontfamily="Computer Modern", xlims=(0,par.Lx))
    savefig("u_slice-$i.pdf")
    plot(x, V[:,n], xlabel = L"$x$", ylabel = L"$V(x,25,t)$", legend = false, fontfamily="Computer Modern", xlims=(0,par.Lx))
    savefig("V_slice-$i.pdf")
    plot(x, ϕ[:,n], xlabel = L"$x$", ylabel = L"$\phi(x,25,t)$", legend = false, fontfamily="Computer Modern", xlims=(0,par.Lx))
    savefig("phi_slice-$i.pdf")
    writedlm("x.csv", x)
    writedlm("u-$i.csv", U[:,n])
end

function draw_slice(par, x, y, U, ϕ, i, n)
    gr(); plot() # Load GR plotting backend and clear previous plots
    plot(x, U[:,n], xlabel = L"$x$", ylabel = L"$u(x,25,t)$", legend = false, fontfamily="Computer Modern", xlims=(0,par.Lx))
    savefig("u_slice-$i.pdf")
    plot(x, ϕ[:,n], xlabel = L"$x$", ylabel = L"$\phi(x,25,t)$", legend = false, fontfamily="Computer Modern", xlims=(0,par.Lx))
    savefig("phi_slice-$i.pdf")
end

function front_position(x, ϕ, n, dx)
    ϕv = ϕ[:,n]
    for i = 1:length(ϕv)
        if (ϕv[i] < 0) && (ϕv[i+1] >= 0)
            θ = ϕv[i]/(ϕv[i] - ϕv[i+1])
            return x[i] + θ*dx
        end
    end
end