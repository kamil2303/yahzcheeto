using LinearAlgebra
using BenchmarkTools
using Random

# 1. Алгоритм обратного хода Гаусса для треугольных матриц
function backward_substitution(U::AbstractMatrix, b::AbstractVector)
    """
    Обратный ход метода Гаусса для верхней треугольной матрицы U.
    Решает систему Ux = b.
    """
    n = size(U, 1)
    x = similar(b, eltype(U))
    
    for i in n:-1:1
        # Проверяем, что диагональный элемент не нулевой
        if iszero(U[i, i])
            throw(SingularException(i))
        end
        
        # Вычисляем i-ю компоненту решения
        x[i] = b[i]
        for j in i+1:n
            x[i] -= U[i, j] * x[j]
        end
        x[i] /= U[i, i]
    end
    
    return x
end

# 2. Приведение матрицы к ступенчатому виду (метод Гаусса)
function gaussian_elimination!(A::AbstractMatrix{T}) where T
    """
    Приводит матрицу A к ступенчатому виду методом Гаусса.
    Модифицирует исходную матрицу.
    """
    m, n = size(A)
    pivot_row = 1
    
    # Преобразуем матрицу к типу с поддержкой деления, если это целые числа
    if T <: Integer
        A_float = convert(Matrix{float(T)}, A)
        return gaussian_elimination!(A_float)
    end
    
    for col in 1:min(m, n)
        # Находим опорный элемент (максимальный по модулю в столбце)
        max_val = zero(eltype(A))
        pivot_index = pivot_row
        for i in pivot_row:m
            if abs(A[i, col]) > abs(max_val)
                max_val = A[i, col]
                pivot_index = i
            end
        end
        
        if !iszero(A[pivot_index, col])
            # Меняем строки местами, если нужно
            if pivot_index != pivot_row
                A[[pivot_row, pivot_index], :] = A[[pivot_index, pivot_row], :]
            end
            
            # Обнуляем элементы под опорным
            for i in pivot_row+1:m
                factor = A[i, col] / A[pivot_row, col]
                A[i, col:end] .-= factor .* A[pivot_row, col:end]
            end
            
            pivot_row += 1
        end
    end
    
    return A
end

function gaussian_elimination(A::AbstractMatrix)
    """
    Приводит матрицу A к ступенчатому виду (без модификации исходной).
    """
    # Для целочисленных матриц преобразуем к Float64
    if eltype(A) <: Integer
        return gaussian_elimination!(convert(Matrix{Float64}, A))
    else
        return gaussian_elimination!(copy(A))
    end
end

# 3. Полный метод Гаусса для решения СЛАУ
function gauss_solve(A::AbstractMatrix, b::AbstractVector)
    """
    Решает систему Ax = b методом Гаусса.
    """
    n = size(A, 1)
    
    # Преобразуем к типу с поддержкой деления, если это целые числа
    if eltype(A) <: Integer
        A_float = convert(Matrix{Float64}, A)
        b_float = convert(Vector{Float64}, b)
        return gauss_solve(A_float, b_float)
    end
    
    # Создаем расширенную матрицу
    Ab = hcat(A, b)
    
    # Прямой ход: приведение к треугольному виду
    for col in 1:n
        # Выбор главного элемента по столбцу (без использования argmax)
        max_val = zero(eltype(Ab))
        pivot_row = col
        for i in col:n
            if abs(Ab[i, col]) > abs(max_val)
                max_val = Ab[i, col]
                pivot_row = i
            end
        end
        
        if iszero(Ab[pivot_row, col])
            throw(SingularException(col))
        end
        
        # Перестановка строк
        if pivot_row != col
            Ab[[col, pivot_row], :] = Ab[[pivot_row, col], :]
        end
        
        # Исключение переменной
        for i in col+1:n
            factor = Ab[i, col] / Ab[col, col]
            Ab[i, col:end] .-= factor .* Ab[col, col:end]
        end
    end
    
    # Обратный ход
    x = backward_substitution(Ab[:, 1:n], Ab[:, end])
    
    return x
end

# 4. Оптимизированная версия метода Гаусса (только для Float64)
function optimized_gauss_solve(A::AbstractMatrix{Float64}, b::AbstractVector{Float64})
    """
    Оптимизированная версия метода Гаусса для больших систем (только Float64).
    """
    n = size(A, 1)
    A_copy = copy(A)
    b_copy = copy(b)
    
    # Прямой ход с выбором главного элемента
    for k in 1:n-1
        # Выбор главного элемента (только для Float64 используем argmax)
        pivot_index = argmax(abs.(A_copy[k:n, k])) + k - 1
        
        if pivot_index != k
            # Перестановка строк в матрице A
            A_copy[[k, pivot_index], :] = A_copy[[pivot_index, k], :]
            # Перестановка элементов в векторе b
            b_copy[[k, pivot_index]] = b_copy[[pivot_index, k]]
        end
        
        # Исключение переменной (векторизованная операция)
        for i in k+1:n
            factor = A_copy[i, k] / A_copy[k, k]
            A_copy[i, k+1:end] .-= factor .* A_copy[k, k+1:end]
            b_copy[i] -= factor * b_copy[k]
        end
    end
    
    # Обратный ход
    x = similar(b_copy)
    for i in n:-1:1
        x[i] = (b_copy[i] - dot(A_copy[i, i+1:end], x[i+1:end])) / A_copy[i, i]
    end
    
    return x
end

# 5. Функция для вычисления ранга матрицы
function matrix_rank(A::AbstractMatrix)
    """
    Вычисляет ранг матрицы путем приведения к ступенчатому виду.
    """
    # Для целочисленных матриц преобразуем к Float64
    if eltype(A) <: Integer
        A_float = convert(Matrix{Float64}, A)
        return matrix_rank(A_float)
    end
    
    # Приводим матрицу к ступенчатому виду
    U = gaussian_elimination(A)
    
    # Считаем ненулевые строки
    rank = 0
    tolerance = 1e-10
    for i in 1:min(size(U)...)
        if !all(x -> abs(x) < tolerance, U[i, i:end])
            rank += 1
        end
    end
    
    return rank
end

# 6. Функция для вычисления определителя
function matrix_determinant(A::AbstractMatrix)
    """
    Вычисляет определитель квадратной матрицы методом Гаусса.
    """
    n = size(A, 1)
    if size(A, 2) != n
        throw(DimensionMismatch("Матрица должна быть квадратной"))
    end
    
    # Для целочисленных матриц преобразуем к Float64
    if eltype(A) <: Integer
        A_float = convert(Matrix{Float64}, A)
        return matrix_determinant(A_float)
    end
    
    A_copy = copy(A)
    det_val = one(eltype(A))
    sign_val = 1
    
    for col in 1:n-1
        # Выбор главного элемента (без использования argmax)
        max_val = zero(eltype(A_copy))
        pivot_row = col
        for i in col:n
            if abs(A_copy[i, col]) > abs(max_val)
                max_val = A_copy[i, col]
                pivot_row = i
            end
        end
        
        if iszero(A_copy[pivot_row, col])
            return zero(eltype(A))  # Вырожденная матрица
        end
        
        if pivot_row != col
            # Меняем строки местами
            A_copy[[col, pivot_row], :] = A_copy[[pivot_row, col], :]
            sign_val *= -1  # Меняем знак определителя
        end
        
        # Исключение переменной
        for i in col+1:n
            factor = A_copy[i, col] / A_copy[col, col]
            A_copy[i, col+1:end] .-= factor .* A_copy[col, col+1:end]
        end
        
        # Умножаем определитель на диагональный элемент
        det_val *= A_copy[col, col]
    end
    
    det_val *= A_copy[n, n] * sign_val
    return det_val
end

# Функция для генерации случайных рациональных чисел
function generate_rational_matrix(n::Int, m::Int)
    """
    Генерирует матрицу случайных рациональных чисел.
    """
    A = Matrix{Rational{Int64}}(undef, n, m)
    for i in 1:n
        for j in 1:m
            # Генерируем рациональное число с числителем и знаменателем от -9 до 9
            num = rand(-9:9)
            den = rand(1:9)  # Знаменатель не может быть 0
            A[i, j] = num // den
        end
    end
    return A
end

function generate_rational_vector(n::Int)
    """
    Генерирует вектор случайных рациональных чисел.
    """
    v = Vector{Rational{Int64}}(undef, n)
    for i in 1:n
        num = rand(-9:9)
        den = rand(1:9)
        v[i] = num // den
    end
    return v
end

# Специальная функция для создания треугольной матрицы из рациональных чисел
function generate_triangular_rational_matrix(n::Int)
    """
    Генерирует верхнюю треугольную матрицу из рациональных чисел.
    Гарантирует, что диагональные элементы не нулевые.
    """
    U = Matrix{Rational{Int64}}(undef, n, n)
    for i in 1:n
        for j in 1:n
            if j < i
                U[i, j] = 0 // 1  # Нижний треугольник - нули
            else
                # Для диагонали гарантируем ненулевые значения
                if i == j
                    num = rand([-9:-1; 1:9])  # Исключаем 0
                    den = rand(1:9)
                else
                    num = rand(-9:9)
                    den = rand(1:9)
                end
                U[i, j] = num // den
            end
        end
    end
    return U
end

# Функция для создания треугольной матрицы Float64
function generate_triangular_float64(n::Int)
    """
    Генерирует верхнюю треугольную матрицу из Float64.
    """
    U = zeros(Float64, n, n)
    for i in 1:n
        for j in i:n
            U[i, j] = randn()
        end
    end
    # Гарантируем ненулевые диагональные элементы
    for i in 1:n
        if U[i, i] == 0
            U[i, i] = rand([-1.0, 1.0]) * (rand() + 0.1)
        end
    end
    return U
end

# Функция для создания треугольной матрицы BigFloat
function generate_triangular_bigfloat(n::Int)
    """
    Генерирует верхнюю треугольную матрицу из BigFloat.
    """
    U = zeros(BigFloat, n, n)
    for i in 1:n
        for j in i:n
            U[i, j] = randn(BigFloat)
        end
    end
    # Гарантируем ненулевые диагональные элементы
    for i in 1:n
        if U[i, i] == 0
            U[i, i] = rand([-1.0, 1.0]) * (rand(BigFloat) + 0.1)
        end
    end
    return U
end

# Отдельная функция для тестирования рациональных чисел
function test_rational_backward_substitution()
    println("   Rational тест:")
    n = 5
    
    # Создаем верхнюю треугольную матрицу с рациональными числами
    U_rational = generate_triangular_rational_matrix(n)
    x_true_rational = generate_rational_vector(n)
    b_rational = U_rational * x_true_rational
    
    # Решаем систему
    x_rational = backward_substitution(U_rational, b_rational)
    
    # Проверяем решение
    residual_rational = U_rational * x_rational - b_rational
    println("   Rational невязка: ", norm(float.(residual_rational)))
    println("   Точное решение: ", all(iszero.(residual_rational)))
    
    return true
end

# Функция для тестирования
function test_algorithms()
    println("Тестирование алгоритмов:")
    println("="^50)
    
    # Тест 1: Обратный ход для треугольной матрицы
    println("1. Тест обратного хода:")
    n = 5
    
    # Float64
    U_float = generate_triangular_float64(n)
    x_true_float = randn(n)
    b_float = U_float * x_true_float
    x_calculated_float = backward_substitution(U_float, b_float)
    residual_float = norm(x_calculated_float - x_true_float)
    println("   Float64 невязка: ", residual_float)
    
    # Тест 2: Разные численные типы
    println("\n2. Тест разных численных типов:")
    
    # BigFloat
    U_big = generate_triangular_bigfloat(n)
    b_big = U_big * randn(BigFloat, n)
    x_big = backward_substitution(U_big, b_big)
    println("   BigFloat невязка: ", norm(U_big * x_big - b_big))
    
    # Rational{Int64} - отдельный тест
    test_rational_backward_substitution()
    
    # Тест 3: Производительность для больших систем (только Float64)
    println("\n3. Тест производительности для n=1000:")
    
    n_large = 1000
    A_large = randn(n_large, n_large)
    b_large = randn(n_large)
    
    println("   Встроенная функция \\:")
    @time x_builtin = A_large \ b_large
    
    println("   Наш метод Гаусса (только для Float64):")
    @time x_gauss = optimized_gauss_solve(A_large, b_large)
    
    println("   Невязка встроенного метода: ", norm(A_large * x_builtin - b_large))
    println("   Невязка нашего метода: ", norm(A_large * x_gauss - b_large))
    
    # Тест 4: Ранг и определитель (используем Float64 вместо Int64)
    println("\n4. Тест ранга и определителя:")
    
    A_test = [1.0 2.0 3.0; 4.0 5.0 6.0; 7.0 8.0 9.0]  # Вырожденная матрица
    println("   Ранг: наш = ", matrix_rank(A_test), ", встроенный = ", rank(A_test))
    
    A_test2 = randn(4, 4)
    det_our = matrix_determinant(A_test2)
    det_builtin = det(A_test2)
    println("   Определитель: наш = ", det_our, ", встроенный = ", det_builtin)
    println("   Относительная ошибка: ", abs(det_our - det_builtin) / abs(det_builtin))
end

# Запуск тестов
if abspath(PROGRAM_FILE) == @__FILE__
    test_algorithms()
end