from GaussSeidel import GaussSeidel
from math import pi


def linear(points, toFind):
    if toFind < points[0][0] or toFind > points[len(points)-1][0]:
        print("Error")
        exit()
    p1 = points[0]
    p2 = points[len(points)-1]
    for i in range(len(points)):
        if toFind > points[i][0] > p1[0]:
            p1 = points[i]
        if toFind < points[i][0] < p2[0]:
            p2 = points[i]
    print("*** Linear_Interpolation ***")
    print("Approximate value = ", ((((p1[1] - p2[1]) * toFind) / (p1[0] - p2[0])) + ((p2[1] * p1[0] - p1[1] * p2[0]) / (p1[0] - p2[0]))))
    print("---------------------------------")
    return



def polynomial(points, toFind):
    matrix = [[points[row][0] ** col for col in range(len(points))] for row in range(len(points))]
    vector = [[points[row][1] for _ in range(1)] for row in range(len(points))]
    coefficient = GaussSeidel(matrix, vector)
    print("*** Polynomial_Interpolation ***")
    value = 0.0
    for x in range(len(coefficient)):
        value += coefficient[x][0] * (toFind ** x)
    print("Approximate value = ", value)
    print("---------------------------------")


def Lagrange(points, toFind):
    i = mySum = 0
    for i in range(len(points)):
        xi, yi, myL = points[i][0], points[i][1], 1
        for j in range(len(points)):
            if j != i:
                xj = points[j][0]
                myL *= (toFind - xj) / (xi - xj)
        mySum += myL * yi

    print("*** Lagrange_Interpolation ***")
    print("Approximate value = ", mySum)
    print("---------------------------------")


def Neville(points, toFind):
    if len(points) < 4:
        print("Can't use Neville's algorithm with less than 4 points...")
        exit()
    print("*** Neville's_Algorithm ***")
    print("Approximate value = ", recurssiveNeville(points, 0, len(points) - 1, toFind))
    print("---------------------------------")


def recurssiveNeville(points, m, n, toFind):
    # stop condition
    if m == n:
        return points[m][1]
    xm, xn = points[m][0], points[n][0]
    return (((toFind - xm) * recurssiveNeville(points, m + 1, n, toFind)) - ((toFind - xn) * recurssiveNeville(points, m, n - 1, toFind))) / (xn - xm)


def cubicSpline(points, toFind, derived1, derived2):
    h = [0 for i in range(len(points))]
    lam = [0 for i in range(len(points))]
    u = [0 for i in range(len(points))]
    d = [0 for i in range(len(points))]
    for i in range(len(points) - 1):
        h[i] = (points[i + 1][0] - points[i][0])

    for i in range(1, len(points)):
        lam[i] = h[i] / (h[i - 1] + h[i])
        u[i] = 1 - lam[i]

    for i in range(1, len(points)):
        d[i] = (6 / h[i - 1]) * (derived2- (points[i][1] - points[i - 1][1])/ h[0])

    d[0] = (6 / h[0]) * (((points[1][1]-points[0][1])/h[0]) - derived1)
    u[len(points) - 1] = lam[0] = 1

    matrix = [[0 for col in range(len(points))] for row in range(len(points))]
    vector = [[d[row] for _ in range(1)] for row in range(len(points))]
    matrix[0][0] = matrix[len(points)-1][len(points)-1] = 2
    matrix[0][1], matrix[len(points)-1][len(points)-2] = lam[0], u[len(points)-1]
    for row in range(1, len(matrix)-1):
        for col in range(1, len(matrix)-1):
            if row == col:
                matrix[row][col] = 2
                matrix[row][col - 1] = u[row]
                matrix[row][col+1] = lam[row]
                break

    p1 = 0
    p2 = len(points) - 1
    for i in range(len(points)):
        if toFind > points[i][0] > points[p1][0]:
            p1 = i
        if toFind < points[i][0] < points[p2][0]:
            p2 = i


    M = GaussSeidel(matrix, vector)
    s = ((((points[p2][0] - toFind) ** 3) * M[p1][0] + ((toFind - points[p1][0]) ** 3) * M[p2][0]) / (6 * h[p1])) + (((points[p2][0] - toFind) * points[p1][1] + (toFind - points[p1][0]) * points[p2][1]) / h[p1]) - (((points[p2][0] - toFind) * M[p1][0] + (toFind - points[p1][0]) * M[p2][0]) * h[p1]) / 6
    print("*** Cubic Spline's_Algorithm ***")
    print("Approximate value = ", s)
    print("---------------------------------")


def cubicSplineNatural(points, toFind):
    h = [0 for i in range(len(points))]
    lam = [0 for i in range(len(points))]
    u = [0 for i in range(len(points))]
    d = [0 for i in range(len(points))]
    for i in range(len(points)-1):
        h[i] = (points[i+1][0] - points[i][0])

    for i in range(1, len(points)):
        lam[i] = h[i] / (h[i-1] + h[i])
        u[i] = 1 - lam[i]

    for i in range(1, len(points)-1):
        d[i] = (6/(h[i-1]+h[i]))*((points[i+1][1]-points[i][1])/h[i] - (points[i][1]-points[i-1][1])/h[i-1])

    d[len(points) - 1] = d[0] = u[len(points) - 1] = lam[0] = 0

    matrix = [[0 for col in range(len(points))] for row in range(len(points))]
    vector = [[d[row] for _ in range(1)] for row in range(len(points))]
    for row in range(len(matrix)-1):
        for col in range(len(matrix)-1):
            if row == col:
                matrix[row][col] = 2
                matrix[row][col + 1] = lam[row]
                matrix[row + 1][col] = u[row+1]

    p1 = 0
    p2 = len(points) - 1
    for i in range(len(points)):
        if toFind > points[i][0] > points[p1][0]:
            p1 = i
        if toFind < points[i][0] < points[p2][0]:
            p2 = i


    matrix2 = [[0 for col in range(len(matrix)-2)] for row in range(len(matrix)-2)]
    for row in range(1, len(matrix)-1):
        for col in range(1, len(matrix)-1):
            matrix2[row-1][col-1]=matrix[row][col]
    vector2 = [[d[row] for _ in range(1)] for row in range(1, len(points)-1)]

    coefficient = GaussSeidel(matrix2, vector2)
    M = [[0 for _ in range(1)] for row in range(len(points))]
    for row in range(1, len(points)-1):
        M[row][0] = coefficient[row-1][0]
    s = ((((points[p2][0] - toFind)**3)*M[p1][0] + ((toFind - points[p1][0])**3)*M[p2][0])/(6 * h[p1])) + (((points[p2][0] - toFind)*points[p1][1] + (toFind - points[p1][0])*points[p2][1])/h[p1]) - (((points[p2][0] - toFind)*M[p1][0] + (toFind - points[p1][0])*M[p2][0])*h[p1])/6
    print("*** Cubic Spline's_Algorithm ***")
    print("Approximate value = ", s)
    print("---------------------------------")




points = [[0, 0], [1, 0.8415], [2, 0.9093], [3, 0.1411], [4, -0.7568], [5, -0.9589], [6, -0.2794]]
toFind = 2.5

linear(points, toFind)

points = [[1, 0.8415], [2, 0.9093], [3, 0.1411]]
toFind = 2.5
polynomial(points, toFind)

points = [[1, 1], [2, 0], [4, 1.5]]
toFind = 3
Lagrange(points, toFind)

points = [[1, 0], [1.2, 0.112463], [1.3, 0.167996], [1.4, 0.222709]]
toFind = 1.28
Neville(points, toFind)

points = [[0, 0], [pi/6, 0.5], [pi/4, 0.7072], [pi/2, 1]]
toFind = pi/3
cubicSplineNatural(points , toFind)

points = [[0, 0], [pi/6, 0.5], [pi/4, 0.7072], [pi/2, 1]]
toFind = pi/3
cubicSpline(points , toFind, 1, 0)
