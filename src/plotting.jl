#=
using Gadfly

function ticks(min_val, max_val)
    return collect(linspace(min_val, max_val, 5))
end

function imshow(x, y, z; x_label="", y_label="", z_label="", title="")

    xx = repeat(x, inner=length(y))
    yy = repeat(y, outer=length(x))

    p = plot(
        layer(color=vec(z), x=xx, y=yy, Geom.rectbin),
        Scale.x_continuous(minvalue=x[1], maxvalue=x[end]),
        Scale.y_continuous(minvalue=y[1], maxvalue=y[end]),
        Guide.xticks(ticks=ticks(x[1], x[end])),
        Guide.yticks(ticks=ticks(y[1], y[end])),
        Guide.XLabel(x_label),
        Guide.YLabel(y_label),
        Guide.title(title),
        Guide.colorkey(z_label),
        Theme(panel_fill=colorant"black", grid_line_width = 0inch))
end

p = imshow(x, y, arr)
draw(PDF(string(key) * ".pdf", 4inch, 3inch), p)
=#

using PyPlot


function plot_results(results)

    for (key, res) in results

        field = res[:mean]

        if length(field.dimensions) == 1
            p = plot(
                layer(y=errs, Geom.line),
                layer(y=w, Geom.line))
        elseif length(field.dimensions) == 2
            arr = as_array(field)

            x = field.dimensions[1]
            y = field.dimensions[2]

            # plotting the propagation dimension horizontally regardless of the original order
            if (y == res[:pdimension])
                x, y = y, x
                arr = arr'
            end

            fig = figure()
            s = fig[:add_subplot](1, 1, 1)
            im = s[:imshow](
                arr',
                extent=[x.grid.grid[1], x.grid.grid[end], y.grid.grid[1], y.grid.grid[end]],
                origin="lower",
                aspect="auto")
            s[:set_xlabel](x.name)
            s[:set_ylabel](y.name)
            fig[:colorbar](im, label=string(key))
        end

        fig[:savefig](string(key) * ".pdf")

    end
end
