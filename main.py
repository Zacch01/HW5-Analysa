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
