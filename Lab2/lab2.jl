using LinearAlgebra
using Random
using Plots
using Printf

n = 100
k = 0.01
ratio_interval = (0.01, 2.0)
points = 80

function make_matrix(n, ratio)
    Random.seed!(1234)
    A = rand(1.0:100.0, n, n)
    for i in 1:n
        S = sum(abs, A[i, :]) - abs(A[i, i])
        A[i, i] = round(ratio * S, digits=2)
    end
    return A
end

function has_diagonal_dominance(A)
    return all(
        abs(A[i, i]) >= sum(abs, A[i, :]) - abs(A[i, i])
        for i in 1:size(A, 1)
    )
end

function gauss(A0, b0)
    A = copy(A0)
    b = copy(b0)
    n = length(b)
    for k in 1:n-1
        A[k, k] == 0 && error("Нулевой ведущий элемент на шаге $k")
        for i in k+1:n
            m = A[i, k] / A[k, k]
            for j in k:n
                A[i, j] -= m * A[k, j]
            end
            b[i] -= m * b[k]
        end
    end
    A[n, n] == 0 && error("Нулевой ведущий элемент")
    x = zeros(n)
    for i in n:-1:1
        s = 0.0
        for j in i+1:n
            s += A[i, j] * x[j]
        end
        x[i] = (b[i] - s) / A[i, i]
    end
    return x
end

function relative_error(x, x_t)
    return 100 * norm(x_t - x) / norm(x)
end

function vector_string(x, digits)
    format = Printf.Format("%.$(digits)f")
    return "[" * join((Printf.format(format, value) for value in x), ", ") * "]"
end

function log_ticks(values)
    positive = values[values .> 0]
    min_exp = floor(Int, log10(minimum(positive)))
    max_exp = ceil(Int, log10(maximum(positive)))

    if max_exp - min_exp < 5
        center = (min_exp + max_exp) / 2
        min_exp = floor(Int, center - 2.5)
        max_exp = min_exp + 5
    end

    step = max(1, ceil(Int, (max_exp - min_exp) / 5))
    exponents = collect(min_exp:step:(min_exp + 5 * step))
    ticks = 10.0 .^ exponents
    labels = ["10^$(e)" for e in exponents]

    return (ticks, labels), (minimum(ticks), maximum(ticks))
end

A = make_matrix(n, k)
x_t = ones(n)
b = A * x_t

@printf("n = %d\n", n)
@printf("k = %.2f\n", k)
println(
    "Диагональное преобладание: ",
    has_diagonal_dominance(A) ? "есть" : "нет"
)

println()
println("Точное решение")
println("x_t = ", vector_string(x_t, 2))

println()
println("b = ", vector_string(b, 15))

x_r = gauss(A, b)
error_r = relative_error(x_r, x_t)

println()
println("Метод Гаусса")
println("x_r = ", vector_string(x_r, 15))

println()
@printf("Относительная ошибка = %.10e %%\n", error_r)

x_b = A \ b
error_b = relative_error(x_b, x_t)

println()
println("Библиотека LAPACK")
println("x_b = ", vector_string(x_b, 15))

println()
@printf("Относительная ошибка = %.10e %%\n", error_b)

ratios = range(ratio_interval[1], ratio_interval[2], length=points)

error_r_graph = zeros(length(ratios))
error_b_graph = zeros(length(ratios))

for i in eachindex(ratios)
    A_current = make_matrix(n, ratios[i])
    x_t_current = ones(n)
    b_current = A_current * x_t_current
    x_r_current = gauss(A_current, b_current)
    x_b_current = A_current \ b_current
    error_r_graph[i] = relative_error(x_r_current, x_t_current)
    error_b_graph[i] = relative_error(x_b_current, x_t_current)
end

error_r_plot = max.(error_r_graph, eps(Float64))
error_b_plot = max.(error_b_graph, eps(Float64))

yticks_r, ylims_r = log_ticks(error_r_plot)
yticks_b, ylims_b = log_ticks(error_b_plot)

p1 = plot(
    ratios,
    error_r_plot,
    xlabel="|a_ii| / Σ|a_ij|",
    ylabel="Относительная ошибка, % (у.е.)",
    title="Метод Гаусса",
    xlims=(0, 2),
    xticks=0:0.2:2,
    ylims=ylims_r,
    yticks=yticks_r,
    yscale=:log10,
    marker=:circle,
    markersize=3,
    linewidth=2,
    legend=false,
    grid=true,
    size=(900, 500)
)

display(p1)

p2 = plot(
    ratios,
    error_b_plot,
    xlabel="|a_ii| / Σ|a_ij|",
    ylabel="Относительная ошибка, % (у.е.)",
    title="Библиотека LAPACK",
    xlims=(0, 2),
    xticks=0:0.2:2,
    ylims=ylims_b,
    yticks=yticks_b,
    yscale=:log10,
    marker=:circle,
    markersize=3,
    linewidth=2,
    legend=false,
    grid=true,
    size=(900, 500)
)

display(p2)
