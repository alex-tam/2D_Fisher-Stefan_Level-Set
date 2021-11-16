# Plot solutions to 2D Fisher-Stefan
# Alex Tam, 09/11/2021

"Plot solutions as 2D heat maps"
function draw_heat(x, y, U, V, ϕ, i)
    gr(); plot() # Load GR plotting backend and clear previous plots
    heatmap(x,y,transpose(U), xlabel = L"$x$", ylabel = L"$y$", aspect_ratio=:equal)
    savefig("u-$i.png")
    heatmap(x,y,transpose(V), xlabel = L"$x$", ylabel = L"$y$", aspect_ratio=:equal)
    savefig("V-$i.png")
    heatmap(x,y,transpose(ϕ), xlabel = L"$x$", ylabel = L"$y$", aspect_ratio=:equal)
    savefig("phi-$i.png")
end

function draw_heat(x, y, U, ϕ, i)
    gr(); plot() # Load GR plotting backend and clear previous plots
    heatmap(x,y,transpose(U), xlabel = L"$x$", ylabel = L"$y$", aspect_ratio=:equal)
    savefig("u-$i.png")
    heatmap(x,y,transpose(ϕ), xlabel = L"$x$", ylabel = L"$y$", aspect_ratio=:equal)
    savefig("phi-$i.png")
end

# Edit slice plots to correctly identify the contact line
"Plot solutions as slices"
function draw_slice(x, y, U, V, ϕ, i, n)
    gr(); plot() # Load GR plotting backend and clear previous plots
    plot(x, U[:,n], xlabel = L"$x$", ylabel = L"$u(x,25,t)$", legend = false)
    savefig("u_slice-$i.png")
    plot(x, V[:,n], xlabel = L"$x$", ylabel = L"$V(x,25,t)$", legend = false)
    savefig("V_slice-$i.png")
    plot(x, ϕ[:,n], xlabel = L"$x$", ylabel = L"$\phi(x,25,t)$", legend = false)
    savefig("phi_slice-$i.png")
    writedlm("x.csv", x)
    writedlm("u-$i.csv", U[:,n])
end

function draw_slice(x, y, U, ϕ, i, n)
    gr(); plot() # Load GR plotting backend and clear previous plots
    plot(x, U[:,n], xlabel = L"$x$", ylabel = L"$u(x,25,t)$", legend = false)
    savefig("u_slice-$i.png")
    plot(x, ϕ[:,n], xlabel = L"$x$", ylabel = L"$\phi(x,25,t)$", legend = false)
    savefig("phi_slice-$i.png")
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