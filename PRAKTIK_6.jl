using StaticArrays
using LinearAlgebra

# 1. Базовая структура Vector_2D
struct Vector_2D{T<:Real}
    x::T
    y::T
end

# Базовые операции для Vector_2D
Base.:+(a::Vector_2D, b::Vector_2D) = Vector_2D(a.x + b.x, a.y + b.y)
Base.:-(a::Vector_2D, b::Vector_2D) = Vector_2D(a.x - b.x, a.y - b.y)
Base.:*(a::Vector_2D, k::Real) = Vector_2D(a.x * k, a.y * k)
Base.:*(k::Real, a::Vector_2D) = a * k
Base.:/(a::Vector_2D, k::Real) = Vector_2D(a.x / k, a.y / k)
Base.:-(a::Vector_2D) = Vector_2D(-a.x, -a.y)

# Скалярное произведение
LinearAlgebra.dot(a::Vector_2D, b::Vector_2D) = a.x * b.x + a.y * b.y

# Векторное произведение (псевдоскаляр)
cross(a::Vector_2D, b::Vector_2D) = a.x * b.y - a.y * b.x

# Норма вектора
LinearAlgebra.norm(a::Vector_2D) = sqrt(a.x^2 + a.y^2)

# Нормализация
LinearAlgebra.normalize(a::Vector_2D) = a / norm(a)

# Проверка на равенство
Base.:(==)(a::Vector_2D, b::Vector_2D) = a.x == b.x && a.y == b.y

# Вывод в консоль
Base.show(io::IO, v::Vector_2D) = print(io, "Vector_2D($(v.x), $(v.y))")

# 2. Структура с наследованием от FieldVector
struct Vector2D{T<:Real} <: FieldVector{2, T}
    x::T
    y::T
end

# Конструктор из кортежа
Vector2D((x, y)::Tuple{T,T}) where T<:Real = Vector2D{T}(x, y)

# Преобразование в кортеж
Base.Tuple(v::Vector2D) = (v.x, v.y)

# 3. Структура для отрезка
struct Segment{T<:Real}
    a::Vector_2D{T}
    b::Vector_2D{T}
end

# 4. Функции для геометрических операций

# Проверка, лежат ли точки по одну сторону от прямой (способ 1 - через векторное произведение)
function same_side_1(p1::Vector_2D{T}, p2::Vector_2D{T}, seg::Segment{T}) where T<:Real
    # Вектор направления прямой
    dir = seg.b - seg.a
    
    # Векторы от точки a отрезка к проверяемым точкам
    v1 = p1 - seg.a
    v2 = p2 - seg.a
    
    # Знаки векторных произведений
    cross1 = cross(dir, v1)
    cross2 = cross(dir, v2)
    
    # Точки лежат по одну сторону, если произведения имеют одинаковый знак
    return cross1 * cross2 >= 0
end

# Проверка, лежат ли точки по одну сторону от прямой (способ 2 - через уравнение прямой)
function same_side_2(p1::Vector_2D{T}, p2::Vector_2D{T}, seg::Segment{T}) where T<:Real
    # Уравнение прямой: (x - xa)(yb - ya) - (y - ya)(xb - xa) = 0
    # Подставляем координаты точек в уравнение прямой
    f1 = (p1.x - seg.a.x) * (seg.b.y - seg.a.y) - (p1.y - seg.a.y) * (seg.b.x - seg.a.x)
    f2 = (p2.x - seg.a.x) * (seg.b.y - seg.a.y) - (p2.y - seg.a.y) * (seg.b.x - seg.a.x)
    
    # Точки лежат по одну сторону, если значения имеют одинаковый знак
    return f1 * f2 >= 0
end

# Проверка выпуклости многоугольника (исправленная сигнатура)
function is_convex(polygon::Vector{Vector_2D{T}}) where T<:Real
    n = length(polygon)
    if n < 3
        return false  # Многоугольник должен иметь хотя бы 3 вершины
    end
    
    # Проверяем знак векторного произведения для всех последовательных ребер
    sign = nothing
    for i in 1:n
        # Текущая, следующая и предыдущая вершины
        a = polygon[i]
        b = polygon[i % n + 1]
        c = polygon[(i + 1) % n + 1]
        
        # Векторы ребер
        v1 = b - a
        v2 = c - b
        
        # Векторное произведение
        cross_val = cross(v1, v2)
        
        # Определяем знак первого ненулевого произведения
        if sign === nothing && !iszero(cross_val)
            sign = cross_val > 0
        end
        
        # Проверяем, что все произведения имеют одинаковый знак
        if !iszero(cross_val) && (cross_val > 0) != sign
            return false
        end
    end
    
    return true
end

# Проверка, лежит ли точка внутри многоугольника (метод трассировки луча)
function point_in_polygon(point::Vector_2D{T}, polygon::Vector{Vector_2D{T}}) where T<:Real
    n = length(polygon)
    inside = false
    
    for i in 1:n
        j = i % n + 1
        p1 = polygon[i]
        p2 = polygon[j]
        
        # Проверяем пересечение горизонтального луча с ребром
        if ((p1.y > point.y) != (p2.y > point.y)) &&
           (point.x < (p2.x - p1.x) * (point.y - p1.y) / (p2.y - p1.y) + p1.x)
            inside = !inside
        end
    end
    
    return inside
end

# Алгоритм Джарвиса (заворачивания подарка) для построения выпуклой оболочки
function jarvis_convex_hull(points::Vector{Vector_2D{T}}) where T<:Real
    n = length(points)
    if n < 3
        return points  # Выпуклая оболочка тривиальна
    end
    
    # Находим самую левую точку
    leftmost = 1
    for i in 2:n
        if points[i].x < points[leftmost].x || 
           (points[i].x == points[leftmost].x && points[i].y < points[leftmost].y)
            leftmost = i
        end
    end
    
    hull = Vector{Vector_2D{T}}()
    current = leftmost
    
    while true
        push!(hull, points[current])
        
        # Ищем следующую точку в оболочке
        next_point = (current % n) + 1
        for i in 1:n
            if i != current
                # Проверяем ориентацию
                cross_val = cross(points[i] - points[current], points[next_point] - points[current])
                if cross_val < 0 || 
                   (iszero(cross_val) && norm(points[i] - points[current]) > norm(points[next_point] - points[current]))
                    next_point = i
                end
            end
        end
        
        current = next_point
        if current == leftmost
            break
        end
    end
    
    return hull
end

# Алгоритм Грехома для построения выпуклой оболочки
function graham_convex_hull(points::Vector{Vector_2D{T}}) where T<:Real
    n = length(points)
    if n < 3
        return points
    end
    
    # Находим точку с минимальной y-координатой (и минимальной x, если совпадение)
    min_idx = 1
    for i in 2:n
        if points[i].y < points[min_idx].y || 
           (points[i].y == points[min_idx].y && points[i].x < points[min_idx].x)
            min_idx = i
        end
    end
    
    # Перемещаем минимальную точку в начало
    points[1], points[min_idx] = points[min_idx], points[1]
    pivot = points[1]
    
    # Сортируем точки по полярному углу относительно pivot
    sorted_points = sort(points[2:end], by=p -> (atan(p.y - pivot.y, p.x - pivot.x), norm(p - pivot)))
    
    # Строим оболочку
    hull = [pivot, sorted_points[1], sorted_points[2]]
    
    for i in 3:length(sorted_points)
        while length(hull) >= 2
            a = hull[end-1]
            b = hull[end]
            c = sorted_points[i]
            
            # Проверяем ориентацию
            if cross(b - a, c - a) <= 0
                pop!(hull)
            else
                break
            end
        end
        push!(hull, sorted_points[i])
    end
    
    return hull
end

# Площадь многоугольника методом трапеций (ориентированная)
function polygon_area_trapezoid(polygon::Vector{Vector_2D{T}}) where T<:Real
    n = length(polygon)
    area = zero(T)
    
    for i in 1:n
        j = i % n + 1
        area += (polygon[i].x + polygon[j].x) * (polygon[j].y - polygon[i].y)
    end
    
    return area / 2
end

# Площадь многоугольника методом треугольников (ориентированная)
function polygon_area_triangles(polygon::Vector{Vector_2D{T}}) where T<:Real
    n = length(polygon)
    area = zero(T)
    
    # Фиксируем первую вершину и разбиваем на треугольники
    for i in 2:n-1
        v1 = polygon[i] - polygon[1]
        v2 = polygon[i+1] - polygon[1]
        area += cross(v1, v2) / 2
    end
    
    return area
end

# Вспомогательные функции
function orientation(a::Vector_2D{T}, b::Vector_2D{T}, c::Vector_2D{T}) where T<:Real
    # Возвращает ориентацию тройки точек (a, b, c)
    # > 0: против часовой стрелки
    # = 0: коллинеарны
    # < 0: по часовой стрелке
    return (b.y - a.y) * (c.x - b.x) - (b.x - a.x) * (c.y - b.y)
end

# Функция для тестирования
function test_geometry_functions()
    println("Тестирование геометрических функций:")
    println("="^50)
    
    # Создаем тестовые точки
    p1 = Vector_2D(0.0, 0.0)
    p2 = Vector_2D(1.0, 0.0)
    p3 = Vector_2D(1.0, 1.0)
    p4 = Vector_2D(0.0, 1.0)
    p5 = Vector_2D(0.5, 0.5)
    p6 = Vector_2D(2.0, 0.5)
    
    # Тестирование базовых операций
    println("1. Базовые операции:")
    println("   p1 + p2 = ", p1 + p2)
    println("   p2 - p1 = ", p2 - p1)
    println("   p1 * 2 = ", p1 * 2)
    println("   dot(p1, p2) = ", dot(p1, p2))
    println("   norm(p3) = ", norm(p3))
    
    # Тестирование проверки сторон
    println("\n2. Проверка расположения точек относительно прямой:")
    seg = Segment(p1, p3)
    println("   Способ 1: ", same_side_1(p2, p4, seg))
    println("   Способ 2: ", same_side_2(p2, p4, seg))
    
    # Тестирование выпуклости (теперь должно работать)
    println("\n3. Проверка выпуклости:")
    square = [p1, p2, p3, p4]
    concave = [p1, p2, Vector_2D(1.5, 0.5), p3, p4]
    println("   Квадрат выпуклый: ", is_convex(square))
    println("   Вогнутый многоугольник выпуклый: ", is_convex(concave))
    
    # Тестирование принадлежности точки
    println("\n4. Принадлежность точки многоугольнику:")
    println("   Точка внутри квадрата: ", point_in_polygon(p5, square))
    println("   Точка вне квадрата: ", point_in_polygon(p6, square))
    
    # Тестирование выпуклых оболочек
    println("\n5. Выпуклые оболочки:")
    random_points = [Vector_2D(rand(), rand()) for _ in 1:10]
    
    jarvis_hull = jarvis_convex_hull(random_points)
    graham_hull = graham_convex_hull(random_points)
    
    println("   Алгоритм Джарвиса: ", length(jarvis_hull), " точек")
    println("   Алгоритм Грехома: ", length(graham_hull), " точек")
    
    # Тестирование площади
    println("\n6. Вычисление площади:")
    area_trap = polygon_area_trapezoid(square)
    area_tri = polygon_area_triangles(square)
    println("   Метод трапеций: ", area_trap)
    println("   Метод треугольников: ", area_tri)
    println("   Погрешность: ", abs(area_trap - area_tri))
    
    # Тестирование Vector2D с FieldVector
    println("\n7. Тестирование Vector2D с FieldVector:")
    v1 = Vector2D(1.0, 2.0)
    v2 = Vector2D(3.0, 4.0)
    println("   v1 + v2 = ", v1 + v2)
    println("   dot(v1, v2) = ", dot(v1, v2))
    println("   norm(v1) = ", norm(v1))
end

# Дополнительные примеры использования
function additional_examples()
    println("\nДополнительные примеры:")
    println("="^50)
    
    # Создание и использование Vector2D с FieldVector
    points = [Vector2D(rand(2)...) for _ in 1:5]
    println("Случайные точки: ", points)
    
    # Матричные операции (благодаря FieldVector)
    M = [1.0 2.0; 3.0 4.0]
    v = Vector2D(1.0, 2.0)
    println("Матричное умножение: ", M * v)
end

# Запуск тестов
if abspath(PROGRAM_FILE) == @__FILE__
    test_geometry_functions()
    additional_examples()
end