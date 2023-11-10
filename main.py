from math import sin, cos, fabs, log10  # Математичні функції для обчислення.
from tabulate import tabulate  # Створення таблиці.
import matplotlib.pyplot as plt  # Малювання графіку залежності.
import numpy as np  # Для зберігання й обчислення діпазону значення Х і У функції.
from scipy.optimize import fsolve  # Для обрахунку розв'язку системи рівнянь.
from sympy import sin, cos  # Для обрахунку розв'язку системи рівнянь.


# Функція отримання х з першого рівняння системи (phi1).
def getYorPhi1(x=float) -> float:
    return sin(x - 0.6) - 1.6

# Функція отримання y з першого рівняння системи (phi2).
def getXorPhi2(y=float) -> float:
    return (cos(y) + 0.9) / 3.0

# Похідна від функції Phi1.
def getDerivativePhi1(x=float) ->float:
    return cos(x - 0.6)

# Похідна від функції Phi2.
def getDerivativePhi2(y=float) ->float:
    return -sin(y) / 3.0

# Перевірка чи метод підходить для даної системи на заданому проміжку.
def checkStartOrContinue(x=float, y=float):
    return fabs(getDerivativePhi1(x)) <= 1 and fabs(getDerivativePhi2(y)) <= 1

# Метод простої ітерації.
def simpleIteration(x0=float, y0=float, eps=float, maxIteration=int) -> list:
    x, y = x0, y0
    listIterations = []
    messageSuccessStartValues = ""
    if (not checkStartOrContinue(x, y)):
        raise Exception(f"Для заданих значень при цьому методі не виконується обчислення!")
    else:
        messageSuccessStartValues = f"\t*Відбулась перевірка успішна заданих початкових значень ({x0}; {y0}), яка дозволила подальші обчислення."
    for i in range(maxIteration):
        xN, yN = getXorPhi2(y), getYorPhi1(x)
        deltaX, deltaY = fabs(xN - x), fabs(yN - y)
        isTrueDerivativePhi1and2 = checkStartOrContinue(xN, yN)
        isEndIteration = deltaX < eps and deltaY < eps and isTrueDerivativePhi1and2
        if (isTrueDerivativePhi1and2):
            listIterations.append([
                i + 1,
                f"x{i}: {x}",
                f"y{i}: {y}",
                f"x{i + 1}: {xN}",
                f"y{i + 1}: {yN}",
                f"|x{i + 1} - x{i}|: {deltaX}",
                f"|y{i + 1} - y{i}|: {deltaY}",
                f"{getDerivativePhi1(x)} <= 1",
                f"{getDerivativePhi2(y)} <= 1",
                not isEndIteration
            ])
        if isEndIteration and isTrueDerivativePhi1and2:
            return [[xN, yN], i + 1, listIterations, messageSuccessStartValues]
        elif isEndIteration:
            return [[x, y], i + 1, listIterations]
        x, y = xN, yN
    else:
        raise Exception(
            f"Було зроблено максимальну кількість ітерацій: {maxIteration}. Оберіть інші початкові значення (х; у): ({x0}; {y0}).")


# Повернення таблиці кількості ітерацій Методу простої ітерації.
def getTableIterationOfSimpleIteration(x=float, y=float, eps=float, maxIteration=int) -> str:
    dataMethod = simpleIteration(x, y, eps, maxIteration)
    powerEps = int(log10(eps))
    nameTable = f"Таблиця ітерацій при точності: 10 ^ {powerEps}; при початкових (х; у): ({x}; {y})"
    headersTable = [
        "\nІтерація №",
        "\nX(k)",
        "\nY(k)",
        "\nX(k + 1)",
        "\nY(k + 1)",
        "\ndeltaX",
        "\ndeltaY",
        "\n|phi1(X(k + 1))'| <= 1",
        "\n|phi2(Y(k + 1))'| <= 1",
        f"deltaX > eps\ndeltaY > eps\neps: 10 ^ {powerEps} "
    ]
    formatXandY = "{:.11f}" # Бо точність може бути [1; 10] знаків після коми.
    xEps = formatXandY.format(dataMethod[0][0])
    yEps = formatXandY.format(dataMethod[0][1])
    dataMethod[2].append(["=" * 25 for i in range(len(dataMethod[2][0]))])
    dataMethod[2].append([f"Кількість => {dataMethod[1]}", "", "", "", "", f"Точність (eps): 10 ^ {powerEps}", "", "Знайдені корені:", f"x = {xEps}", f"y = {yEps}"])
    strTable = tabulate(dataMethod[2], headersTable, tablefmt="pretty")
    return f"{nameTable}\n{dataMethod[3]}\n{strTable}"


# Графік фукнцій.
def showFunctions(coeficient=float):
    x1 = np.linspace(coeficient * np.pi, coeficient * -np.pi, 400)
    y1 = np.sin(x1 - 0.6) - 1.6
    y2 = x1
    x2 = (np.cos(y2) + 0.9) / 3
    plt.plot(x1, y1, color="yellow", linewidth=2, label="y = sin(x - 0.6) - 1.6")
    plt.plot(x2, y2, color="blue", linewidth=2, label="x = (cos(y) + 0.9) / 3")
    plt.axhline(0, color="red", linewidth=0.5)
    plt.axvline(0, color="red", linewidth=0.5)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend()
    plt.grid(True)
    plt.title(f"ЛР № 5, варіант № 3, Вальчевський П. В., ОІ-11 сп\nГрафіки функцій на межах [-{coeficient} * pi; {coeficient} * pi]")
    plt.show()


# Демонстрація графіку залежності кількості ітерацій та точності.
def showGraphicMethodIterationAndEps(x=float, y=float, maxIteration=int) -> None:
    powEps = [-i for i in range(1, 10 + 1)]
    countIterationInPowEps = [simpleIteration(x, y, 10 ** i, maxIteration)[1] for i in powEps]
    plt.plot(powEps, countIterationInPowEps, label='Залежність')
    plt.xlabel("Степінь значення точності в форматі 10 ^ i")
    plt.ylabel("Кікільсть ітерацій")
    plt.title(f"ЛР № 5, варіант № 3, Вальчевський П. В., ОІ-11 сп\nГрафік залежності кількості ітерацій та точності")
    plt.legend()
    plt.grid(True)
    plt.show()


# Розв'язок системи рівнянь за допомогою бібліотек (для перевірки).
def equations(vars):
    x, y = vars
    eq1 = sin(x - 0.6) - y - 1.6
    eq2 = 3 * x - cos(y) - 0.9
    return [eq1, eq2]


def getInfo():
    strInfo = "ЛР № 5, варіант № 3, Вальчевський П. В., ОІ-11 сп, Метод простої ітерації, система рівнянь:\n"
    strInfo += "\t----\n"
    strInfo += "\t|sin(x - 0.6) - y - 1.6 = 0\n"
    strInfo += "  --|\n"
    strInfo += "\t|3 * x - cos(y) - 0.9 = 0\n"
    strInfo += "\t----\n"
    return strInfo


# Програма виконання.
if __name__ == "__main__":
    x, y, eps, maxIteration = -2.0, -2.0, 10 ** -5, 10 ** 3
    solveSystem = fsolve(equations, [x, y])
    print(getInfo())
    print(getTableIterationOfSimpleIteration(x, y, eps, maxIteration))
    print(f"Розв'язок системи рівнянь знайдений за допомогою бібліотек: ({solveSystem[0]}; {solveSystem[1]})")
    print("Графіки зображені в інших вікнах!")
    showGraphicMethodIterationAndEps(x, y, maxIteration)
    showFunctions(1.5)
