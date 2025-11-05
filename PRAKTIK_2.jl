# --- ЗАДАНИЕ 2: РЕАЛИЗАЦИЯ МЕТОДА НЬЮТОНА И РЕШЕНИЕ УРАВНЕНИЯ cos(x) = x ---

# Добавляем необходимые пакеты
using ForwardDiff  # Для автоматического дифференцирования

# Обобщенная функция метода Ньютона
function newthon(r::Function, x; epsilon=1e-8, num_max=10)
    """
    Реализация метода Ньютона для нахождения корня уравнения r(x) = 0.
    
    Параметры:
    r - функция, корень которой ищем (должна возвращать значение и производную)
    x - начальное приближение
    epsilon - точность
    num_max - максимальное количество итераций
    
    Возвращает:
    Приближенное значение корня или nothing, если метод не сошелся
    """
    iteration = 0
    
    while iteration < num_max
        # Вычисляем значение функции и ее производной
        f_val, f_der = r(x)
        
        # Проверяем, не равна ли производная нулю (деление на ноль)
        if abs(f_der) < epsilon
            println("Производная близка к нулю. Метод может не сойтись.")
            return nothing
        end
        
        # Вычисляем следующее приближение по формуле Ньютона
        x_next = x - f_val / f_der
        
        # Проверяем условие сходимости
        if abs(x_next - x) < epsilon
            return x_next
        end
        
        x = x_next
        iteration += 1
    end
    
    println("Метод Ньютона не сошелся за $num_max итераций")
    return nothing
end

# --- СПОСОБ 1: Аналитическое вычисление производной ---

# Функция для уравнения cos(x) = x -> cos(x) - x = 0
function equation_analytic(x)
    """
    Возвращает значение функции f(x) = cos(x) - x и ее производной f'(x) = -sin(x) - 1
    """
    f_val = cos(x) - x
    f_der = -sin(x) - 1
    return (f_val, f_der)
end

# --- СПОСОБ 2: Автоматическое дифференцирование с помощью пакета ForwardDiff ---

function equation_auto_diff(x)
    """
    Использует автоматическое дифференцирование из пакета ForwardDiff
    """
    # Функция, корень которой ищем
    f(x) = cos(x) - x
    
    # Вычисляем значение функции и производной
    f_val = f(x)
    f_der = ForwardDiff.derivative(f, x)
    
    return (f_val, f_der)
end

# --- СПОСОБ 3: Аппроксимация производной конечными разностями ---

function equation_finite_diff(x; h=1e-6)
    """
    Аппроксимирует производную с помощью отношения конечных разностей
    f'(x) ≈ (f(x+h) - f(x)) / h
    """
    f_val = cos(x) - x
    f_der = (cos(x+h) - cos(x)) / h - 1  # Производная от -x равна -1
    return (f_val, f_der)
end

# --- ДЕМОНСТРАЦИЯ РАБОТЫ ---

println("="^50)
println("РЕШЕНИЕ УРАВНЕНИЯ cos(x) = x МЕТОДОМ НЬЮТОНА")
println("="^50)

# Начальное приближение
initial_guess = 0.5

# --- Способ 1: Аналитическая производная ---
println("\n1. АНАЛИТИЧЕСКОЕ ВЫЧИСЛЕНИЕ ПРОИЗВОДНОЙ:")
root_analytic = newthon(equation_analytic, initial_guess)
println("Найденный корень: ", root_analytic)
println("Проверка: cos($root_analytic) - $root_analytic = ", cos(root_analytic) - root_analytic)

# --- Способ 2: Автоматическое дифференцирование ---
println("\n2. АВТОМАТИЧЕСКОЕ ДИФФЕРЕНЦИРОВАНИЕ (ForwardDiff):")
root_auto = newthon(equation_auto_diff, initial_guess)
println("Найденный корень: ", root_auto)
println("Проверка: cos($root_auto) - $root_auto = ", cos(root_auto) - root_auto)

# --- Способ 3: Конечные разности ---
println("\n3. АППРОКСИМАЦИЯ КОНЕЧНЫМИ РАЗНОСТЯМИ:")
root_finite = newthon(equation_finite_diff, initial_guess)
println("Найденный корень: ", root_finite)
println("Проверка: cos($root_finite) - $root_finite = ", cos(root_finite) - root_finite)

# --- Сравнение результатов ---
println("\n" * "="^50)
println("СРАВНЕНИЕ РЕЗУЛЬТАТОВ:")
println("Аналитический метод: ", root_analytic)
println("Авто дифф. (ForwardDiff): ", root_auto)
println("Конечные разности:   ", root_finite)

# Проверка точности решения
exact_solution = 0.7390851332151607  # Известное точное решение
println("\nТочное решение:       ", exact_solution)
println("Погрешность аналитического: ", abs(root_analytic - exact_solution))
println("Погрешность авто дифф.:     ", abs(root_auto - exact_solution))
println("Погрешность конечных разн.: ", abs(root_finite - exact_solution))

# --- ДОПОЛНИТЕЛЬНО: Визуализация сходимости ---
println("\n" * "="^50)
println("ИССЛЕДОВАНИЕ СХОДИМОСТИ:")

function track_newton_convergence(r::Function, x; epsilon=1e-8, num_max=10)
    """
    Версия метода Ньютона, которая отслеживает процесс сходимости
    """
    iterations = []
    x_values = []
    errors = []
    
    iteration = 0
    while iteration < num_max
        f_val, f_der = r(x)
        
        if abs(f_der) < epsilon
            break
        end
        
        x_next = x - f_val / f_der
        error = abs(x_next - x)
        
        push!(iterations, iteration)
        push!(x_values, x)
        push!(errors, error)
        
        if error < epsilon
            break
        end
        
        x = x_next
        iteration += 1
    end
    
    return iterations, x_values, errors
end

# Отслеживаем сходимость для каждого метода
iter_analytic, x_analytic, err_analytic = track_newton_convergence(equation_analytic, initial_guess)
iter_auto, x_auto, err_auto = track_newton_convergence(equation_auto_diff, initial_guess)
iter_finite, x_finite, err_finite = track_newton_convergence(equation_finite_diff, initial_guess)

println("Количество итераций для сходимости:")
println("Аналитический: ", length(iter_analytic))
println("Авто дифф.: ", length(iter_auto))
println("Конечные разности: ", length(iter_finite))

# Выводим процесс сходимости для аналитического метода
println("\nПроцесс сходимости (аналитический метод):")
for i in 1:length(iter_analytic)
    println("Итерация $(iter_analytic[i]): x = $(x_analytic[i]), ошибка = $(err_analytic[i])")
end

sdsd = readline()

