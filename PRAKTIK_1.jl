# Подключаем модуль LinearAlgebra для использования некоторых функций, если потребуется.
using LinearAlgebra

# --- ЗАДАНИЕ 1: РЕАЛИЗАЦИЯ ТИПА ДАННЫХ POLYNOMIAL ---

# Определяем параметрический композитный тип `Polynomial`.
struct Polynomial{T} <: Number
    coeff::Vector{T}

    function Polynomial{T}(coeff::Vector{T}) where T
        # Удаляем ведущие нули
        first_nonzero = findfirst(!iszero, coeff)
        if first_nonzero === nothing
            return new{T}(T[zero(T)])
        end
        return new{T}(coeff[first_nonzero:end])
    end
end

# Внешний конструктор
Polynomial(coeff::Vector{T}) where T = Polynomial{T}(coeff)

# Степень многочлена
ord(p::Polynomial) = length(p.coeff) - 1

# --- АРИФМЕТИЧЕСКИЕ ОПЕРАЦИИ ---

# Сложение многочленов
function Base.:+(p::Polynomial{T}, q::Polynomial{T}) where T
    n, m = length(p.coeff), length(q.coeff)
    if n > m
        r = copy(p.coeff)
        r[n-m+1:end] .+= q.coeff
    else
        r = copy(q.coeff)
        r[m-n+1:end] .+= p.coeff
    end
    return Polynomial(r)
end

# Сложение с скаляром
Base.:+(p::Polynomial{T}, c::T) where T = p + Polynomial([c])
Base.:+(c::T, p::Polynomial{T}) where T = p + c

# Унарный минус
Base.:-(p::Polynomial{T}) where T = Polynomial(-p.coeff)

# Вычитание
Base.:-(p::Polynomial{T}, q::Polynomial{T}) where T = p + (-q)
Base.:-(p::Polynomial{T}, c::T) where T = p + (-c)
Base.:-(c::T, p::Polynomial{T}) where T = c + (-p)

# Умножение многочленов
function Base.:*(p::Polynomial{T}, q::Polynomial{T}) where T
    n, m = ord(p), ord(q)
    result = zeros(T, n + m + 1)
    for i in 1:length(p.coeff)
        for j in 1:length(q.coeff)
            result[i+j-1] += p.coeff[i] * q.coeff[j]
        end
    end
    return Polynomial(result)
end

# Умножение на скаляр
Base.:*(p::Polynomial{T}, c::T) where T = Polynomial(p.coeff .* c)
Base.:*(c::T, p::Polynomial{T}) where T = p * c

# Деление многочленов - ИСПРАВЛЕНИЕ: используем Float64 для коэффициентов при делении
function Base.divrem(p::Polynomial{Int}, q::Polynomial{Int})
    # Преобразуем Int коэффициенты в Float64 для избежания InexactError
    p_float = Polynomial(Float64.(p.coeff))
    q_float = Polynomial(Float64.(q.coeff))
    
    if ord(q_float) > ord(p_float)
        return (Polynomial([0.0]), p_float)
    end
    
    r = deepcopy(p_float)
    quotient_coeff = zeros(Float64, ord(p_float) - ord(q_float) + 1)

    while ord(r) >= ord(q_float) && !iszero(r.coeff[1])
        c = r.coeff[1] / q_float.coeff[1]
        deg_diff = ord(r) - ord(q_float)
        quotient_coeff[end-deg_diff] = c
        
        # Создаем многочлен для вычитания
        subtrahend_coeff = zeros(Float64, deg_diff + 1)
        subtrahend_coeff[1] = c
        subtrahend = Polynomial(subtrahend_coeff)
        
        r -= q_float * subtrahend
    end
    
    return (Polynomial(quotient_coeff), r)
end

# Для других типов оставляем оригинальную реализацию
function Base.divrem(p::Polynomial{T}, q::Polynomial{T}) where T
    if ord(q) > ord(p)
        return (Polynomial([zero(T)]), p)
    end
    
    r = deepcopy(p)
    quotient_coeff = zeros(T, ord(p) - ord(q) + 1)

    while ord(r) >= ord(q) && !iszero(r.coeff[1])
        c = r.coeff[1] / q.coeff[1]
        deg_diff = ord(r) - ord(q)
        quotient_coeff[end-deg_diff] = c
        
        # Создаем многочлен для вычитания
        subtrahend_coeff = zeros(T, deg_diff + 1)
        subtrahend_coeff[1] = c
        subtrahend = Polynomial(subtrahend_coeff)
        
        r -= q * subtrahend
    end
    
    return (Polynomial(quotient_coeff), r)
end

Base.:/(p::Polynomial{T}, q::Polynomial{T}) where T = divrem(p, q)[1]
Base.:%(p::Polynomial{T}, q::Polynomial{T}) where T = divrem(p, q)[2]

# Вычисление значения
function (p::Polynomial)(x)
    result = p.coeff[1]
    for i in 2:length(p.coeff)
        result = result * x + p.coeff[i]
    end
    return result
end

# Вычисление значения и производной
function valdiff_poly(p::Polynomial, x)
    val = p.coeff[1]
    der = zero(val)
    for i in 2:length(p.coeff)
        der = der * x + val
        val = val * x + p.coeff[i]
    end
    return (val, der)
end

# Вывод многочлена
function Base.show(io::IO, p::Polynomial{T}) where T
    if isempty(p.coeff) || all(iszero, p.coeff)
        print(io, "0")
        return
    end
    
    terms = String[]
    n = ord(p)
    
    for (i, c) in enumerate(p.coeff)
        if iszero(c)
            continue
        end
        
        power = n - i + 1
        abs_c = abs(c)
        sign_str = c > 0 ? "+" : "-"
        
        if i == 1
            sign_str = c < 0 ? "-" : ""
        end
        
        coeff_str = (abs_c == 1 && power > 0) ? "" : string(abs_c)
        power_str = power == 0 ? "" : power == 1 ? "x" : "x^$power"
        
        term = sign_str * " " * coeff_str * power_str
        push!(terms, term)
    end
    
    result = join(terms, " ")
    if startswith(result, "+ ")
        result = result[3:end]
    end
    print(io, result)
end

# --- ЗАДАНИЕ 2: РЕАЛИЗАЦИЯ ТИПА ДАННЫХ DUAL ---

struct Dual{T} <: Number
    a::T
    b::T
end

# Правила промоутинга для Dual
Base.promote_rule(::Type{Dual{T}}, ::Type{S}) where {T,S} = Dual{promote_type(T, S)}

# Конвертация чисел в Dual
Base.convert(::Type{Dual{T}}, x::Number) where T = Dual{T}(convert(T, x), zero(T))
Base.convert(::Type{Dual{T}}, d::Dual{S}) where {T,S} = Dual{T}(convert(T, d.a), convert(T, d.b))

# Базовые операции
Base.real(d::Dual) = d.a

Base.:+(x::Dual{T}, y::Dual{T}) where T = Dual{T}(x.a + y.a, x.b + y.b)
Base.:+(x::Dual{T}, y::Number) where T = x + Dual{T}(y, zero(T))
Base.:+(y::Number, x::Dual{T}) where T = Dual{T}(y, zero(T)) + x

Base.:-(x::Dual) = Dual(-x.a, -x.b)
Base.:-(x::Dual{T}, y::Dual{T}) where T = Dual{T}(x.a - y.a, x.b - y.b)
Base.:-(x::Dual{T}, y::Number) where T = x - Dual{T}(y, zero(T))
Base.:-(y::Number, x::Dual{T}) where T = Dual{T}(y, zero(T)) - x

Base.:*(x::Dual{T}, y::Dual{T}) where T = Dual{T}(x.a * y.a, x.a * y.b + x.b * y.a)
Base.:*(x::Dual{T}, y::Number) where T = x * Dual{T}(y, zero(T))
Base.:*(y::Number, x::Dual{T}) where T = Dual{T}(y, zero(T)) * x

Base.:/(x::Dual{T}, y::Dual{T}) where T = Dual{T}(x.a / y.a, (x.b * y.a - x.a * y.b) / (y.a * y.a))
Base.:/(x::Dual{T}, y::Number) where T = x / Dual{T}(y, zero(T))
Base.:/(y::Number, x::Dual{T}) where T = Dual{T}(y, zero(T)) / x

# Математические функции
Base.sqrt(x::Dual{T}) where T = Dual{T}(sqrt(x.a), x.b / (2 * sqrt(x.a)))
Base.exp(x::Dual{T}) where T = Dual{T}(exp(x.a), x.b * exp(x.a))
Base.log(x::Dual{T}) where T = Dual{T}(log(x.a), x.b / x.a)
Base.sin(x::Dual{T}) where T = Dual{T}(sin(x.a), x.b * cos(x.a))
Base.cos(x::Dual{T}) where T = Dual{T}(cos(x.a), -x.b * sin(x.a))

# Константы
Base.zero(::Type{Dual{T}}) where T = Dual{T}(zero(T), zero(T))
Base.one(::Type{Dual{T}}) where T = Dual{T}(one(T), zero(T))

# --- ЗАДАНИЕ 3: ФУНКЦИЯ VALDIFF ---

function valdiff_func(f::Function, x::T) where T
    d = Dual{T}(x, one(T))
    result = f(d)
    return (real(result), result.b)
end

# --- ЗАДАНИЕ 4: VALDIFF ДЛЯ POLYNOMIAL ЧЕРЕЗ DUAL ---

function valdiff_poly_dual(p::Polynomial{T}, x::T) where T
    d = Dual{T}(x, one(T))
    result = p(d)
    return (real(result), result.b)
end

# --- ТЕСТИРОВАНИЕ ---

println("=== ТЕСТИРОВАНИЕ POLYNOMIAL ===")

# Тест 1: Базовые операции - используем Float64 вместо Int
p1 = Polynomial([1.0, -2.0, -3.0])  # x² - 2x - 3
p2 = Polynomial([2.0, 5.0])         # 2x + 5

println("p1 = ", p1)
println("p2 = ", p2)
println("p1 + p2 = ", p1 + p2)
println("p1 * p2 = ", p1 * p2)
println("p1(3) = ", p1(3))

# Тест 2: Деление
q, r = divrem(p1, p2)
println("p1 / p2 = ", q)
println("p1 % p2 = ", r)

# Тест 3: Производная через схему Горнера
val1, der1 = valdiff_poly(p1, 3.0)
println("Значение p1(3) = ", val1)
println("Производная p1'(3) = ", der1)

println("\n=== ТЕСТИРОВАНИЕ DUAL NUMBERS ===")

# Тест 4: Дуальные числа
d1 = Dual(3.0, 1.0)  # 3 + ε
d2 = Dual(2.0, 0.0)  # 2
println("d1 = ", d1.a, " + ε", d1.b)
println("d2 = ", d2.a)
println("d1 + d2 = ", (d1 + d2).a, " + ε", (d1 + d2).b)
println("d1 * d2 = ", (d1 * d2).a, " + ε", (d1 * d2).b)

# Тест 5: Операции с разными типами
println("3.0 + d1 = ", (3.0 + d1).a, " + ε", (3.0 + d1).b)
println("d1 * 2.0 = ", (d1 * 2.0).a, " + ε", (d1 * 2.0).b)

# Тест 6: Конвертация
d3 = convert(Dual{Float64}, 5)
println("convert(Dual{Float64}, 5) = ", d3.a, " + ε", d3.b)

# Тест 7: Функции с дуальными числами
f(x) = x^2 - 2x - 3
val2, der2 = valdiff_func(f, 3.0)
println("f(3) = ", val2)
println("f'(3) = ", der2)

println("\n=== ТЕСТИРОВАНИЕ POLYNOMIAL С DUAL ===")

# Тест 8: Многочлен с дуальными числами
val3, der3 = valdiff_poly_dual(p1, 3.0)
println("p1(3) через dual = ", val3)
println("p1'(3) через dual = ", der3)

# Тест 9: Сравнение методов
println("\n=== СРАВНЕНИЕ РЕЗУЛЬТАТОВ ===")
println("Горнер: значение = ", val1, ", производная = ", der1)
println("Dual:   значение = ", val3, ", производная = ", der3)
println("Разница: ", abs(val1 - val3), ", ", abs(der1 - der3))

# Дополнительный тест с целыми числами
println("\n=== ТЕСТ С ЦЕЛЫМИ ЧИСЛАМИ ===")
p_int1 = Polynomial([1, -2, -3])
p_int2 = Polynomial([2, 5])
q_int, r_int = divrem(p_int1, p_int2)
println("Целочисленное деление:")
println("p_int1 / p_int2 = ", q_int)
println("p_int1 % p_int2 = ", r_int)

sdsd = readline()
