using Base: NamedTuple

# 1. Реализация типа Residue для кольца вычетов
"""
    Residue{T, M}

Пользовательский тип данных, реализующий кольцо вычетов по модулю M.
Параметры:
- T: тип элементов (Int, Polynom и т.д.)
- M: модуль (может быть числом или кортежем коэффициентов для многочлена)
"""
struct Residue{T, M}
    value::T
    
    # Конструктор с проверкой корректности
    function Residue{T, M}(val) where {T, M}
        # Приводим значение к диапазону [0, M-1] для числовых типов
        if T <: Number
            new(mod(val, M))
        else
            # Для нечисловых типов (например, многочленов) используем как есть
            new(val)
        end
    end
end

# Базовые операции для Residue

# Сложение
function add_residue(a::Residue{T, M}, b::Residue{T, M}) where {T, M} 
    return Residue{T, M}(a.value + b.value)
end

# Вычитание
function sub_residue(a::Residue{T, M}, b::Residue{T, M}) where {T, M} 
    return Residue{T, M}(a.value - b.value)
end

# Унарный минус
function neg_residue(a::Residue{T, M}) where {T, M} 
    return Residue{T, M}(-a.value)
end

# Умножение
function mul_residue(a::Residue{T, M}, b::Residue{T, M}) where {T, M} 
    return Residue{T, M}(a.value * b.value)
end

# Возведение в степень
function pow_residue(a::Residue{T, M}, n::Integer) where {T, M}
    n < 0 && return inv_residue(a)^(-n)  # Обработка отрицательных степеней
    result = Residue{T, M}(1)     # Единичный элемент
    base = a
    # Быстрое возведение в степень
    while n > 0
        if n % 2 == 1
            result = mul_residue(result, base)
        end
        base = mul_residue(base, base)
        n ÷= 2
    end
    return result
end

# Переопределение операторов для удобства
Base.:+(a::Residue, b::Residue) = add_residue(a, b)
Base.:-(a::Residue, b::Residue) = sub_residue(a, b)
Base.:-(a::Residue) = neg_residue(a)
Base.:*(a::Residue, b::Residue) = mul_residue(a, b)
Base.:^(a::Residue, n::Integer) = pow_residue(a, n)

# Вывод в консоль
Base.show(io::IO, a::Residue{T, M}) where {T, M} = print(io, "Residue{$T, $M}(", a.value, ")")

# 2. Реализация типа Polynom для работы с многочленами
"""
    Polynom{T}

Тип данных для представления многочленов с коэффициентами типа T.
"""
struct Polynom{T}
    coeffs::Vector{T}  # Коэффициенты многочлена (от младшей степени к старшей)
    
    # Конструктор, удаляющий ведущие нули
    function Polynom{T}(coeffs::Vector{T}) where T
        # Удаляем ведущие нули
        i = length(coeffs)
        while i > 1 && iszero(coeffs[i])
            i -= 1
        end
        new(coeffs[1:i])
    end
end

# Конструктор по умолчанию
Polynom(coeffs::Vector{T}) where T = Polynom{T}(coeffs)

# Получение кортежа коэффициентов (для использования в качестве параметра M)
function get_tuple(p::Polynom)
    return tuple(p.coeffs...)
end

# Базовые операции для Polynom

# Сложение многочленов
function add_polynom(p1::Polynom{T}, p2::Polynom{T}) where T
    n = max(length(p1.coeffs), length(p2.coeffs))
    result = zeros(T, n)
    result[1:length(p1.coeffs)] .= p1.coeffs
    result[1:length(p2.coeffs)] .+= p2.coeffs
    return Polynom(result)
end

# Вычитание многочленов
function sub_polynom(p1::Polynom{T}, p2::Polynom{T}) where T
    n = max(length(p1.coeffs), length(p2.coeffs))
    result = zeros(T, n)
    result[1:length(p1.coeffs)] .= p1.coeffs
    result[1:length(p2.coeffs)] .-= p2.coeffs
    return Polynom(result)
end

# Умножение многочленов
function mul_polynom(p1::Polynom{T}, p2::Polynom{T}) where T
    n = length(p1.coeffs) + length(p2.coeffs) - 1
    result = zeros(T, n)
    for i in 1:length(p1.coeffs)
        for j in 1:length(p2.coeffs)
            result[i+j-1] += p1.coeffs[i] * p2.coeffs[j]
        end
    end
    return Polynom(result)
end

# Переопределение операторов для удобства
Base.:+(p1::Polynom, p2::Polynom) = add_polynom(p1, p2)
Base.:-(p1::Polynom, p2::Polynom) = sub_polynom(p1, p2)
Base.:*(p1::Polynom, p2::Polynom) = mul_polynom(p1, p2)

# Деление многочленов по модулю (для Residue)
function mod_polynom(p::Polynom{T}, mod_poly::Polynom{T}) where T
    # Реализация деления многочленов с остатком
    dividend = copy(p.coeffs)
    divisor = mod_poly.coeffs
    n = length(dividend)
    m = length(divisor)
    
    if n < m
        return Polynom(dividend)
    end
    
    for i in n:-1:m
        if !iszero(dividend[i])
            factor = dividend[i] ÷ divisor[m]
            for j in 1:m
                dividend[i-m+j] -= factor * divisor[j]
            end
        end
    end
    
    # Находим первую ненулевую позицию
    first_nonzero = findfirst(!iszero, dividend)
    if first_nonzero === nothing
        return Polynom([zero(T)])
    else
        return Polynom(dividend[first_nonzero:end])
    end
end

# Вывод многочлена
function Base.show(io::IO, p::Polynom)
    if isempty(p.coeffs) || all(iszero, p.coeffs)
        print(io, "0")
    else
        terms = []
        for (i, c) in enumerate(p.coeffs)
            if !iszero(c)
                if i == 1
                    push!(terms, "$c")
                else
                    degree = i - 1
                    if degree == 1
                        push!(terms, "$c*x")
                    else
                        push!(terms, "$c*x^$degree")
                    end
                end
            end
        end
        print(io, join(reverse(terms), " + "))
    end
end

# 3. Реализация расширенного алгоритма Евклида
"""
    extended_gcd(a, b)

Расширенный алгоритм Евклида. Возвращает кортеж (d, u, v), где:
- d = НОД(a, b)
- u, v такие, что d = u*a + v*b
"""
function extended_gcd(a::T, b::T) where T
    if b == zero(T)
        return (a, one(T), zero(T))
    else
        d, u1, v1 = extended_gcd(b, a % b)
        u = v1
        v = u1 - (a ÷ b) * v1
        return (d, u, v)
    end
end

# 4. Функция для нахождения обратного элемента в кольце вычетов
"""
    inv_residue(a::Residue)

Находит обратный элемент для a в кольце вычетов, если он существует.
Возвращает обратный элемент или nothing, если элемент необратим.
"""
function inv_residue(a::Residue{T, M}) where {T, M}
    d, u, v = extended_gcd(a.value, M)
    if d == one(T)
        return Residue{T, M}(u)
    else
        return nothing
    end
end

# 5. Решение диофантова уравнения
"""
    solve_diophantine(a::T, b::T, c::T) where T

Решает диофантово уравнение a*x + b*y = c.
Возвращает кортеж (x, y) с решением или nothing, если решения нет.
"""
function solve_diophantine(a::T, b::T, c::T) where T
    d, u, v = extended_gcd(a, b)
    if c % d != zero(T)
        return nothing  # Уравнение не имеет решений
    else
        k = c ÷ d
        return (u * k, v * k)
    end
end

# 6. Проверка числа на простоту
"""
    check_prime(n::Integer)

Проверяет, является ли число n простым.
Использует тест простоты методом пробного деления.
"""
function check_prime(n::Integer)
    n <= 1 && return false
    n == 2 && return true
    n % 2 == 0 && return false
    
    # Проверяем делители до квадратного корня из n
    i = 3
    while i * i <= n
        if n % i == 0
            return false
        end
        i += 2
    end
    return true
end

# 7. Решето Эратосфена
"""
    sieve_eratosthenes(n::Integer)

Возвращает массив всех простых чисел, не превосходящих n.
Реализует алгоритм "решето Эратосфена".
"""
function sieve_eratosthenes(n::Integer)
    n < 2 && return Int[]
    
    # Создаем массив и отмечаем составные числа
    is_prime_arr = trues(n)
    is_prime_arr[1] = false
    
    i = 2
    while i * i <= n
        if is_prime_arr[i]
            # Отмечаем все кратные i как составные
            j = i * i
            while j <= n
                is_prime_arr[j] = false
                j += i
            end
        end
        i += 1
    end
    
    # Собираем все простые числа
    primes = Int[i for i in 2:n if is_prime_arr[i]]
    return primes
end

# 8. Факторизация числа
"""
    factorize_number(n::T) where T <: Integer

Возвращает вектор именованных кортежей (div, deg) с простыми делителями
числа n и их кратностями.
"""
function factorize_number(n::T) where T <: Integer
    n <= 1 && return NamedTuple{(:div, :deg), Tuple{T, T}}[]
    
    list = NamedTuple{(:div, :deg), Tuple{T, T}}[]
    remaining = abs(n)
    
    # Проверяем делимость на 2
    if remaining % 2 == 0
        k = get_degree(remaining, T(2))
        push!(list, (div=T(2), deg=k))
        remaining ÷= 2^k
    end
    
    # Проверяем нечетные делители
    i = T(3)
    while i * i <= remaining
        if remaining % i == 0
            k = get_degree(remaining, i)
            push!(list, (div=i, deg=k))
            remaining ÷= i^k
        end
        i += 2
    end
    
    # Если осталось простое число больше 1
    if remaining > 1
        push!(list, (div=remaining, deg=T(1)))
    end
    
    return list
end

"""
    get_degree(n, p)

Вычисляет кратность делителя p в числе n.
"""
function get_degree(n, p)
    k = 0
    temp = n
    while temp % p == 0
        k += 1
        temp ÷= p
    end
    return k
end

# Демонстрация работы всех функций
function demo()
    println("Демонстрация работы:")
    println("="^50)
    
    # 1. Демонстрация Residue с числами
    println("1. Кольцо вычетов по модулю 7:")
    a = Residue{Int, 7}(5)
    b = Residue{Int, 7}(3)
    println("a = ", a)
    println("b = ", b)
    println("a + b = ", a + b)
    println("a * b = ", a * b)
    println("a^3 = ", a^3)
    println()
    
    # 2. Демонстрация Polynom
    println("2. Многочлены:")
    p1 = Polynom([1, 2, 3])  # 3x^2 + 2x + 1
    p2 = Polynom([0, 1, 1])  # x^2 + x
    println("p1 = ", p1)
    println("p2 = ", p2)
    println("p1 + p2 = ", p1 + p2)
    println("p1 * p2 = ", p1 * p2)
    println()
    
    # 3. Демонстрация extended_gcd
    println("3. Расширенный алгоритм Евклида:")
    d, u, v = extended_gcd(48, 18)
    println("extended_gcd(48, 18) = ($d, $u, $v)")
    println("Проверка: $u*48 + $v*18 = $(u*48 + v*18)")
    println()
    
    # 4. Демонстрация inv_residue
    println("4. Обратный элемент в кольце вычетов:")
    a = Residue{Int, 7}(3)
    inv_a = inv_residue(a)
    println("Обратный к $a: ", inv_a)
    if inv_a !== nothing
        println("Проверка: $a * $inv_a = ", a * inv_a)
    end
    println()
    
    # 5. Демонстрация solve_diophantine
    println("5. Решение диофантова уравнения:")
    solution = solve_diophantine(3, 5, 7)
    if solution !== nothing
        x, y = solution
        println("Решение 3x + 5y = 7: x=$x, y=$y")
        println("Проверка: 3*$x + 5*$y = $(3*x + 5*y)")
    else
        println("Уравнение 3x + 5y = 7 не имеет решений")
    end
    println()
    
    # 6. Демонстрация check_prime
    println("6. Проверка чисел на простоту:")
    for num in [2, 3, 4, 17, 25, 29]
        println("$num - простое? ", check_prime(num))
    end
    println()
    
    # 7. Демонстрация sieve_eratosthenes
    println("7. Решето Эратосфена для n=30:")
    primes = sieve_eratosthenes(30)
    println("Простые числа ≤ 30: ", primes)
    println()
    
    # 8. Демонстрация factorize_number
    println("8. Факторизация чисел:")
    for num in [12, 17, 60, 100]
        factors = factorize_number(num)
        println("$num = ", join(["$(f.div)^$(f.deg)" for f in factors], " * "))
    end
end

# Запуск демонстрации
if abspath(PROGRAM_FILE) == @__FILE__
    demo()
end

dsds = readline()

