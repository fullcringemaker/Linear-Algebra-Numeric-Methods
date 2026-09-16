using LinearAlgebra
using Random
using Plots

plotly()

n = 4
Random.seed!(42)

A = Float64.(rand(-100:100, n, n)) ./ 100.00

# Делаем матрицу диагонально преобладающей
# for i in 1:n
#     A[i, i] += n + 1.00
# end

b = [
    2.00,
    -1.00,
    3.00,
    4.00
]

delta_A = [
    0.01  0.00  0.00  0.00;
    0.00  0.01  0.00  0.00;
    0.00  0.00  0.01  0.00;
    0.00  0.00  0.00  0.01
]

delta_b = [
    0.01,
    -0.01,
    0.02,
    -0.02
]

# РЕШЕНИЕ СЛАУ
function solveSystem(A, b)
    Acopy = copy(A)
    bcopy = copy(b)
    x, factors, pivots = LAPACK.gesv!(Acopy, bcopy)
    return x
end

# НАХОЖДЕНИЕ ОБРАТНОЙ МАТРИЦЫ
function inverseMatrix(A)
    n = size(A, 1)
    E = Matrix{Float64}(I, n, n)
    Ainv, factors, pivots = LAPACK.gesv!(copy(A), E)
    return Ainv
end

# ЧИСЛО ОБУСЛОВЛЕННОСТИ nu(A)
function condNumber(A, p)
    Ainv = inverseMatrix(A)
    return opnorm(A, p) * opnorm(Ainv, p)
end

# АНАЛИЗ ОШИБКИ
function analyzeError(A, b, delta_A, delta_b, p)
    x = solveSystem(A, b)
    xNew = solveSystem(A + delta_A, b + delta_b)
    delta_x = xNew - x
    nu = condNumber(A, p)
    rel_b = norm(delta_b, p) / norm(b, p)
    rel_A = opnorm(delta_A, p) / opnorm(A, p)

    # Прямой расчет относительной ошибки
    errDirect = norm(delta_x, p) / norm(x, p)
    # Теоретическая оценка
    errEstimate = nu * (rel_b + rel_A)

    return (x, xNew, delta_x, nu, rel_b, rel_A, errDirect, errEstimate)
end

# ВЫВОД РЕЗУЛЬТАТОВ
function printAnalysis(A, b, delta_A, delta_b, p)
    x, xNew, delta_x, nu, rel_b, rel_A, errDirect, errEstimate = analyzeError(A, b, delta_A, delta_b, p)

    if p == Inf
        println()
        println("============== p = Inf ==============")
    else
        println()
        println("=============== p = 2 ===============")
    end

    println()
    println("Решение системы x")
    println(round.(x, digits = 6))
    println()
    println("Решение системы x~")
    println(round.(xNew, digits = 6))
    println()
    println("delta_x = x~ - x")
    println(round.(delta_x, digits = 6))
    println()
    println("nu(A)")
    println(round(nu, digits = 6))
    println()
    println("||delta_b||p / ||b||p")
    println(round(rel_b, digits = 6))
    println()
    println("||delta_A||p / ||A||p")
    println(round(rel_A, digits = 6))
    println()
    println("Прямая относительная ошибка")
    println(round(errDirect * 100.00, digits = 6)," %")
    println()
    println("Теоретическая оценка")
    println(round(errEstimate * 100.00, digits = 6)," %")
    println()
    println("Прямая ошибка <= оценка: ", errDirect <= errEstimate)
end

# СОЗДАНИЕ ДАННЫХ ДЛЯ 3D-ГРАФИКА
function buildSurfaces(A, b, delta_A, delta_b, p)
    x = solveSystem(A, b)
    nu = condNumber(A, p)
    
    # Изменение возмущений от 0 до полного значения
    scale = collect(0.00:0.05:1.00)
    count = length(scale)

    xAxis = zeros(count)
    yAxis = zeros(count)

    directSurf = zeros(count, count)
    estimateSurf = zeros(count, count)
    
    # Значения для оси X
    for j in eachindex(scale)
        cur_b = scale[j] * delta_b
        xAxis[j] = norm(cur_b, p)
    end
    
    # Значения для оси Y
    for i in eachindex(scale)
        cur_A = scale[i] * delta_A
        yAxis[i] = opnorm(cur_A, p)
    end

    # Вычисление значений двух поверхностей
    for i in eachindex(scale)
        cur_A = scale[i] * delta_A
        for j in eachindex(scale)
            cur_b = scale[j] * delta_b

            xNew = solveSystem(A + cur_A,b + cur_b)
            delta_x = xNew - x

            directSurf[i, j] = (norm(delta_x, p) / norm(x, p)) * 100.00
            estimateSurf[i, j] = nu * (norm(cur_b, p) / norm(b, p) + opnorm(cur_A, p) / opnorm(A, p)) * 100.00
        end
    end

    return (xAxis, yAxis, directSurf, estimateSurf)
end

# ПОСТРОЕНИЕ 3D-ГРАФИКА
function makeGraph(A, b, delta_A, delta_b, p)
    xAxis, yAxis, directSurf, estimateSurf = buildSurfaces(A, b, delta_A, delta_b, p)

    if p == Inf
        pName = "Inf"
    else
        pName = "2"
    end

    graph = surface(
        xAxis, 
        yAxis, 
        directSurf,
        xlabel = "||delta_b||p",
        ylabel = "||delta_A||p",
        zlabel = "delta(x), %",
        title = "p = " * pName,
        label = "Прямой расчет",
        c = :blues,
        seriesalpha = 0.75,
        colorbar = false,
        camera = (35, 25),
        size = (900, 650)
    )

    surface!(
        graph,
        xAxis,
        yAxis,
        estimateSurf,
        label = "Оценка",
        c = :reds,
        seriesalpha = 0.55,
        colorbar = false
    )

    return graph
end

# ВЫВОД ИСХОДНЫХ ДАННЫХ
println("Матрица A")
display(round.(A, digits = 2))

println()
println("Вектор b")
println(b)

println()
println("Матрица delta_A")
display(delta_A)

println()
println("Вектор delta_b")
println(delta_b)

# РАСЧЕТ ДЛЯ p = 2
printAnalysis(A, b, delta_A, delta_b, 2)

# РАСЧЕТ ДЛЯ p = Inf
printAnalysis(A, b, delta_A, delta_b, Inf)

# ГРАФИК ДЛЯ p = 2
graph2 = makeGraph(A, b, delta_A, delta_b, 2)

# ГРАФИК ДЛЯ p = Inf
graphInf = makeGraph(A, b, delta_A, delta_b, Inf)

# ВЫВОД ГРАФИКОВ
display(graph2)
display(graphInf)
