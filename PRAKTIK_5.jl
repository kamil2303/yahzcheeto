using BenchmarkTools
using Statistics
using Random

# 1. Реализация функций сортировки с интерфейсом как в стандартной библиотеке

# Сортировка пузырьком
function bubble_sort!(v::AbstractVector)
    n = length(v)
    for i in 1:n-1
        swapped = false
        for j in 1:n-i
            if v[j] > v[j+1]
                v[j], v[j+1] = v[j+1], v[j]
                swapped = true
            end
        end
        !swapped && break
    end
    return v
end

function bubble_sort(v::AbstractVector)
    return bubble_sort!(copy(v))
end

# Сортировка вставками
function insertion_sort!(v::AbstractVector)
    for i in 2:length(v)
        key = v[i]
        j = i - 1
        while j > 0 && v[j] > key
            v[j+1] = v[j]
            j -= 1
        end
        v[j+1] = key
    end
    return v
end

function insertion_sort(v::AbstractVector)
    return insertion_sort!(copy(v))
end

# Функции sortperm (возвращает индексы отсортированного массива)
function bubble_sortperm(v::AbstractVector)
    n = length(v)
    indices = collect(1:n)
    for i in 1:n-1
        for j in 1:n-i
            if v[indices[j]] > v[indices[j+1]]
                indices[j], indices[j+1] = indices[j+1], indices[j]
            end
        end
    end
    return indices
end

function insertion_sortperm(v::AbstractVector)
    n = length(v)
    indices = collect(1:n)
    for i in 2:n
        key_idx = indices[i]
        j = i - 1
        while j > 0 && v[indices[j]] > v[key_idx]
            indices[j+1] = indices[j]
            j -= 1
        end
        indices[j+1] = key_idx
    end
    return indices
end

# Функция issorted (проверка отсортированности)
function is_sorted_custom(v::AbstractVector)
    for i in 1:length(v)-1
        if v[i] > v[i+1]
            return false
        end
    end
    return true
end

# 2. Сортировка расческой
function comb_sort!(v::AbstractVector)
    n = length(v)
    gap = n
    shrink = 1.3
    sorted = false
    
    while !sorted
        gap = floor(Int, gap / shrink)
        if gap <= 1
            gap = 1
            sorted = true
        end
        
        for i in 1:n-gap
            if v[i] > v[i+gap]
                v[i], v[i+gap] = v[i+gap], v[i]
                sorted = false
            end
        end
    end
    return v
end

function comb_sort(v::AbstractVector)
    return comb_sort!(copy(v))
end

# 3. Сортировка Шелла
function shell_sort!(v::AbstractVector)
    n = length(v)
    gap = n ÷ 2
    
    while gap > 0
        for i in gap+1:n
            temp = v[i]
            j = i
            while j > gap && v[j-gap] > temp
                v[j] = v[j-gap]
                j -= gap
            end
            v[j] = temp
        end
        gap ÷= 2
    end
    return v
end

function shell_sort(v::AbstractVector)
    return shell_sort!(copy(v))
end

# 4. Сортировка слиянием
function merge_sort!(v::AbstractVector)
    if length(v) <= 1
        return v
    end
    
    mid = length(v) ÷ 2
    left = merge_sort!(v[1:mid])
    right = merge_sort!(v[mid+1:end])
    
    return merge(left, right)
end

function merge(left::AbstractVector, right::AbstractVector)
    result = similar(left, length(left) + length(right))
    i, j, k = 1, 1, 1
    
    while i <= length(left) && j <= length(right)
        if left[i] <= right[j]
            result[k] = left[i]
            i += 1
        else
            result[k] = right[j]
            j += 1
        end
        k += 1
    end
    
    while i <= length(left)
        result[k] = left[i]
        i += 1
        k += 1
    end
    
    while j <= length(right)
        result[k] = right[j]
        j += 1
        k += 1
    end
    
    return result
end

function merge_sort(v::AbstractVector)
    return merge_sort!(copy(v))
end

# 5. Быстрая сортировка (Хоара)
function quick_sort!(v::AbstractVector, low=1, high=length(v))
    if low < high
        pivot_idx = partition!(v, low, high)
        quick_sort!(v, low, pivot_idx-1)
        quick_sort!(v, pivot_idx+1, high)
    end
    return v
end

function partition!(v::AbstractVector, low, high)
    pivot = v[high]
    i = low - 1
    
    for j in low:high-1
        if v[j] <= pivot
            i += 1
            v[i], v[j] = v[j], v[i]
        end
    end
    v[i+1], v[high] = v[high], v[i+1]
    return i + 1
end

function quick_sort(v::AbstractVector)
    return quick_sort!(copy(v))
end

# 7. Медиана с использованием процедуры Хоара
function quick_select!(v::AbstractVector, k::Int, low=1, high=length(v))
    if low == high
        return v[low]
    end
    
    pivot_idx = partition!(v, low, high)
    
    if k == pivot_idx
        return v[k]
    elseif k < pivot_idx
        return quick_select!(v, k, low, pivot_idx-1)
    else
        return quick_select!(v, k, pivot_idx+1, high)
    end
end

function median_quick(v::AbstractVector)
    n = length(v)
    v_copy = copy(v)
    
    if isodd(n)
        return quick_select!(v_copy, (n+1)÷2)
    else
        mid1 = quick_select!(v_copy, n÷2)
        mid2 = quick_select!(v_copy, n÷2+1)
        return (mid1 + mid2) / 2
    end
end

# 8. Пирамидальная сортировка
function heap_sort!(v::AbstractVector)
    n = length(v)
    
    # Построение max-heap
    for i in n÷2:-1:1
        heapify!(v, n, i)
    end
    
    # Извлечение элементов из кучи
    for i in n:-1:2
        v[1], v[i] = v[i], v[1]
        heapify!(v, i-1, 1)
    end
    
    return v
end

function heapify!(v::AbstractVector, n::Int, i::Int)
    largest = i
    left = 2*i
    right = 2*i + 1
    
    if left <= n && v[left] > v[largest]
        largest = left
    end
    
    if right <= n && v[right] > v[largest]
        largest = right
    end
    
    if largest != i
        v[i], v[largest] = v[largest], v[i]
        heapify!(v, n, largest)
    end
end

function heap_sort(v::AbstractVector)
    return heap_sort!(copy(v))
end

# 9. Сортировка подсчетом (для целых чисел)
function counting_sort(v::AbstractVector{<:Integer}, min_val=minimum(v), max_val=maximum(v))
    range = max_val - min_val + 1
    count = zeros(Int, range)
    
    # Подсчет элементов
    for num in v
        count[num - min_val + 1] += 1
    end
    
    # Накопительная сумма
    for i in 2:range
        count[i] += count[i-1]
    end
    
    # Построение отсортированного массива
    result = similar(v)
    for num in reverse(v)
        idx = num - min_val + 1
        result[count[idx]] = num
        count[idx] -= 1
    end
    
    return result
end

# 10. Быстрый поиск в отсортированном массиве
function binary_search(v::AbstractVector, x)
    low, high = 1, length(v)
    
    while low <= high
        mid = (low + high) ÷ 2
        if v[mid] == x
            return mid
        elseif v[mid] < x
            low = mid + 1
        else
            high = mid - 1
        end
    end
    
    return -1  # Элемент не найден
end

# 11. Сортировка столбцов матрицы по среднеквадратическому отклонению

# Способ 1: Использование sortperm
function sort_columns_method1(A::AbstractMatrix)
    # Вычисление СКО для каждого столбца
    std_devs = [std(A[:, j]) for j in axes(A, 2)]
    
    # Получение индексов отсортированных СКО
    sorted_indices = sortperm(std_devs)
    
    # Перестановка столбцов
    return A[:, sorted_indices]
end

# Способ 2: Сортировка вектора столбцов
function sort_columns_method2(A::AbstractMatrix)
    # Преобразование в вектор столбцов
    columns = [A[:, j] for j in axes(A, 2)]
    
    # Сортировка по СКО
    sort!(columns, by=col -> std(col))
    
    # Обратное преобразование в матрицу
    return hcat(columns...)
end

# Функция для вычисления СКО однопроходным алгоритмом
function std_onepass(v::AbstractVector)
    n = length(v)
    if n < 2
        return 0.0
    end
    
    sum_x = 0.0
    sum_x2 = 0.0
    
    for x in v
        sum_x += x
        sum_x2 += x^2
    end
    
    mean = sum_x / n
    variance = (sum_x2 - n * mean^2) / (n - 1)
    return sqrt(variance)
end

# Функция для тестирования всех алгоритмов
function test_all_algorithms()
    println("Тестирование алгоритмов сортировки:")
    println("="^60)
    
    # Генерация тестовых данных
    Random.seed!(123)
    small_data = rand(100)
    medium_data = rand(1000)
    large_data = rand(10000)
    int_data = rand(1:1000, 1000)
    
    # Тестирование базовых алгоритмов на небольших данных
    println("1. Тестирование базовых алгоритмов (n=100):")
    
    data_copy = copy(small_data)
    @time bubble_sort!(data_copy)
    println("   Пузырьковая сортировка: завершена")
    
    data_copy = copy(small_data)
    @time insertion_sort!(data_copy)
    println("   Сортировка вставками: завершена")
    
    data_copy = copy(small_data)
    @time comb_sort!(data_copy)
    println("   Сортировка расческой: завершена")
    
    # Тестирование эффективности на средних данных
    println("\n2. Тестирование эффективности (n=1000):")
    
    println("   Стандартная сортировка:")
    @time sort(medium_data)
    
    println("   Сортировка расческой:")
    @time comb_sort(medium_data)
    
    println("   Сортировка Шелла:")
    @time shell_sort(medium_data)
    
    println("   Быстрая сортировка:")
    @time quick_sort(medium_data)
    
    # Тестирование на больших данных
    println("\n3. Тестирование на больших данных (n=10000):")
    
    println("   Стандартная сортировка:")
    @time sort(large_data)
    
    println("   Пирамидальная сортировка:")
    @time heap_sort(large_data)
    
    println("   Сортировка слиянием:")
    @time merge_sort(large_data)
    
    # Тестирование сортировки подсчетом
    println("\n4. Сортировка подсчетом (n=1000 целых чисел):")
    @time counting_sort(int_data)
    
    # Тестирование медианы
    println("\n5. Вычисление медианы:")
    test_median = rand(1001)
    println("   Медиана (стандартная): ", median(test_median))
    println("   Медиана (quick select): ", median_quick(test_median))
    
    # Тестирование бинарного поиска
    println("\n6. Бинарный поиск:")
    sorted_data = sort(rand(100))
    target = sorted_data[50]
    println("   Найден элемент $target на позиции: ", binary_search(sorted_data, target))
    
    # Тестирование сортировки столбцов матрицы
    println("\n7. Сортировка столбцов матрицы:")
    A = rand(100, 50)
    
    println("   Способ 1 (sortperm):")
    @time result1 = sort_columns_method1(A)
    
    println("   Способ 2 (sort вектора):")
    @time result2 = sort_columns_method2(A)
    
    # Проверка корректности
    stds1 = [std(result1[:, j]) for j in axes(result1, 2)]
    stds2 = [std(result2[:, j]) for j in axes(result2, 2)]
    
    println("   Результаты корректны: ", is_sorted_custom(stds1) && is_sorted_custom(stds2))
end

# Дополнительная функция для сравнения производительности
function compare_performance()
    println("\nСравнение производительности алгоритмов:")
    println("="^60)
    
    sizes = [100, 500, 1000, 5000]
    algorithms = [
        ("Пузырьковая", bubble_sort!),
        ("Вставками", insertion_sort!),
        ("Расческой", comb_sort!),
        ("Шелла", shell_sort!),
        ("Быстрая", quick_sort!),
        ("Слиянием", merge_sort!),
        ("Пирамидальная", heap_sort!)
    ]
    
    for size in sizes
        println("\nРазмер массива: $size")
        data = rand(size)
        
        for (name, algo) in algorithms
            if size <= 1000 || !(name in ["Пузырьковая", "Вставками"])
                data_copy = copy(data)
                time = @elapsed algo(data_copy)
                println("   $name: ", round(time, digits=6), " сек")
            end
        end
        
        # Стандартная сортировка для сравнения
        data_copy = copy(data)
        time = @elapsed sort!(data_copy)
        println("   Стандартная: ", round(time, digits=6), " сек")
    end
end

# Запуск тестов
if abspath(PROGRAM_FILE) == @__FILE__
    test_all_algorithms()
    compare_performance()
end