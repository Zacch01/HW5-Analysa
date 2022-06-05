# Solving Linear Equation Using Gauss Seidel Method
from GaussSeidel import GaussSeidel
from InverseMatrix import InverseMatrix
from math import pi

# Global Variable To Store The Machine Precision, (Set the accuracy of the solution)
ACCURACY = 0.00001


def Linear(pointsList, xToFind):
    """
    Method for finding a Point based on the x value

    :param pointsList: List of point represent the points on the graph
    :param xToFind: value on the axis X, that we are searching for
    """
    if len(pointsList) < 2:
        print('Linear approximation must have minimum of  2 points')
        return

    # In case we are using interpolation
    if pointsList[0][0] <= xToFind <= pointsList[len(pointsList) - 1][0]:
        print("*** Linear Interpolation ***")
        foundY = LinearInterpolation(pointsList, xToFind)

    # In case we are using extrapolation
    else:
        print("*** Linear Extrapolation ***")
        foundY = LinearExtrapolation(pointsList, xToFind)

    # The point approximation
    print(f'Point Approximation --> ({xToFind}, {int(foundY * 10 ** 5) / 10 ** 5})')
    print("---------------------------------\n")


def LinearInterpolation(pointsList, xToFind):
    """
    Interpolation for finding a Point based on the x value

    :param pointsList: List of point represent the points on the graph
    :param xToFind: value on the axis X, that we are searching for
    """
    # Loop to find the nearest points to the wanted one
    for i in range(len(pointsList) - 1):

        # In case the needed action is interpolation
        if pointsList[i][0] <= xToFind <= pointsList[i + 1][0]:

            # Return the Y approximation
            return ((xToFind - pointsList[i + 1][0]) * pointsList[i][1] + (pointsList[i][0] - xToFind) * pointsList[i + 1][1]) / (pointsList[i][0] - pointsList[i + 1][0])


def LinearExtrapolation(pointsList, xToFind):
    """
    Extrapolation for finding a Point based on the x value

    :param pointsList: List of point represent the points on the graph
    :param xToFind: value on the axis X, that we are searching for
    """
    # In case the point we search is on the left side
    if xToFind < pointsList[0][0]:
        Index = 0

    # In case the point we search is on the right side
    else:
        Index = len(pointsList) - 2

    # Return the Y approximation
    return pointsList[Index][1] + (xToFind - pointsList[Index][0]) / (pointsList[Index + 1][0] - pointsList[Index][0]) * (pointsList[Index + 1][1] - pointsList[Index][1])


def polynomial(points, toFind):

    if len(points) < 2:
        print('Polynomial approximation must have minimum of  2 points')
        return

    if toFind < points[0][0] or toFind > points[-1][0]:
        print("The wanted point isn't suitable for interpolation")
        return

    matrix = [[points[row][0] ** col for col in range(len(points))] for row in range(len(points))]
    vector = [[points[row][1] for _ in range(1)] for row in range(len(points))]

    print("*** Polynomial Interpolation ***")

    coefficient = GaussSeidel(matrix, vector)
    if coefficient is None:
        coefficient = InverseMatrix(matrix, vector)
        if coefficient is None:
            print('Cant solve the polynomial')
            return

    value = 0.0
    for x in range(len(coefficient)):
        value += coefficient[x] * (toFind ** x)
    print(f'Point Approximation --> ({toFind}, {int(value * 10 ** 5) / 10 ** 5})')
    print("---------------------------------\n")


def Lagrange(points, toFind):

    if len(points) < 2:
        print('Lagrange approximation must have minimum of  2 points')
        return

    if toFind < points[0][0] or toFind > points[-1][0]:
        print("The wanted point isn't suitable for interpolation")
        return

    i = mySum = 0
    for i in range(len(points)):
        xi, yi, myL = points[i][0], points[i][1], 1
        for j in range(len(points)):
            if j != i:
                xj = points[j][0]
                myL *= (toFind - xj) / (xi - xj)
        mySum += myL * yi

    print("*** Lagrange_Interpolation ***")
    print(f'Point Approximation --> ({toFind}, {int(mySum * 10 ** 5) / 10 ** 5})')
    print("---------------------------------\n")


def Neville(pointsList, xToFind):
    """
    Method for finding a Point based on the x value

    :param pointsList: List of point represent the points on the graph
    :param xToFind: value on the axis X, that we are searching for
    """
    if len(points) < 4:
        print('Neville approximation must have minimum of 4 points')
        return

    if toFind < points[0][0] or toFind > points[-1][0]:
        print("The wanted point isn't suitable for interpolation")
        return

    print("*** Neville's Algorithm ***")
    # The point approximation
    print(f'Point Approximation --> ({xToFind}, {int(recursiveNeville(pointsList, xToFind, 0, len(pointsList) - 1) * 10 ** 5) / 10 ** 5})')
    print("---------------------------------\n")


def recursiveNeville(pointsList, xToFind, i, j):
    """
    Recursive method for finding a Point based on the x value

    :param pointsList: List of point represent the points on the graph
    :param xToFind: value on the axis X, that we are searching for
    :param i: Index that represent Xi
    :param j: Index that represent Xj
    :return: The approximation of Y based on the X
    """
    # Stop condition
    if i == j:
        return pointsList[i][1]

    # Saving the calculation of P[i + 1][j]
    if P[i + 1][j] is None:
        P[i + 1][j] = recursiveNeville(pointsList, xToFind, i + 1, j)

    # Saving the calculation of P[i][j - 1]
    if P[i][j - 1] is None:
        P[i][j - 1] = recursiveNeville(pointsList, xToFind, i, j - 1)

    # Create a sub calculating
    return ((xToFind - pointsList[i][0]) * P[i + 1][j] - (xToFind - pointsList[j][0]) * P[i][j - 1]) / (pointsList[j][0] - pointsList[i][0])


def cubicSpline(points, toFind, derived1, derived2):
    if len(points) < 4:
        print('Cubic spline approximation must have minimum of 4 points')
        return

    if toFind < points[0][0] or toFind > points[-1][0]:
        print("The wanted point isn't suitable for interpolation")
        return

    h = [0 for i in range(len(points))]
    lam = [0 for i in range(len(points))]
    u = [0 for i in range(len(points))]
    d = [0 for i in range(len(points))]
    for i in range(len(points) - 1):
        h[i] = (points[i + 1][0] - points[i][0])

    for i in range(1, len(points)):
        lam[i] = h[i] / (h[i - 1] + h[i])
        u[i] = 1 - lam[i]

    for i in range(1, len(h) - 1):
        d[i] = (6 / (h[i - 1] + h[i])) * (((points[i + 1][1] - points[i][1]) / h[i]) - ((points[i][1] - points[i - 1][1]) / h[i - 1]))

    d[-1] = (6 / h[-2]) * (derived2 - ((points[-1][1] - points[-2][1]) / h[0]))

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
    s = ((((points[p2][0] - toFind) ** 3) * M[p1] + ((toFind - points[p1][0]) ** 3) * M[p2]) / (6 * h[p1])) + (((points[p2][0] - toFind) * points[p1][1] + (toFind - points[p1][0]) * points[p2][1]) / h[p1]) - (((points[p2][0] - toFind) * M[p1] + (toFind - points[p1][0]) * M[p2]) * h[p1]) / 6
    print("*** Cubic Spline's_Algorithm ***")
    print(f'Point Approximation --> ({toFind}, {int(s * 10 ** 5) / 10 ** 5})')
    print("---------------------------------\n")


def cubicSplineNatural(points, toFind):
    if len(points) < 4:
        print('Cubic spline approximation must have minimum of 4 points')
        return

    if toFind < points[0][0] or toFind > points[-1][0]:
        print("The wanted point isn't suitable for interpolation")
        return

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
            matrix2[row-1][col-1] = matrix[row][col]
    vector2 = [[d[row] for _ in range(1)] for row in range(1, len(points)-1)]

    coefficient = GaussSeidel(matrix2, vector2)

    M = [[0 for _ in range(1)] for row in range(len(points))]
    for row in range(1, len(points)-1):
        M[row][0] = coefficient[row-1]
    s = ((((points[p2][0] - toFind)**3)*M[p1][0] + ((toFind - points[p1][0])**3)*M[p2][0])/(6 * h[p1])) + (((points[p2][0] - toFind)*points[p1][1] + (toFind - points[p1][0])*points[p2][1])/h[p1]) - (((points[p2][0] - toFind)*M[p1][0] + (toFind - points[p1][0])*M[p2][0])*h[p1])/6
    print("*** Cubic Spline's_Algorithm ***")
    print(f'Point Approximation --> ({toFind}, {int(s * 10 ** 5) / 10 ** 5})')
    print("---------------------------------\n")


"""
points = [[0, 0], [1, 0.8415], [2, 0.9093], [3, 0.1411], [4, -0.7568], [5, -0.9589], [6, -0.2794]]
toFind = 2.5

Linear(points, toFind)

points = [[1, 0.8415], [2, 0.9093], [3, 0.1411]]
toFind = 2.5
polynomial(points, toFind)

points = [[1, 1], [2, 0], [4, 1.5]]
toFind = 3
Lagrange(points, toFind)

points = [[1, 0], [1.2, 0.112463], [1.3, 0.167996], [1.4, 0.222709]]
toFind = 1.28
P = [[None for col in range(len(points))] for row in range(len(points))]
Neville(points, toFind)

points = [[0, 0], [pi/6, 0.5], [pi/4, 0.7072], [pi/2, 1]]
toFind = pi/3
cubicSplineNatural(points, toFind)
"""

points = [[0, 0], [pi/6, 0.5], [pi/4, 0.7072], [pi/2, 1]]
toFind = pi/3
derived1 = 1
derived2 = 0


method = -1
while method > 5 or method < 0:
    method = int(input('0 --> Linear\n1 --> Polynomial\n2 --> Lagrange\n3 --> Neville\n4 --> Natural Cubic Spline\n5 --> Full Cubic Spline\nInput --> '))

    if method > 5 or method < 0:
        print('Invalid input, Try [Zero to Five]')

if method == 0:
    Linear(points, toFind)

elif method == 1:
    polynomial(points, toFind)

elif method == 2:
    Lagrange(points, toFind)

elif method == 3:
    P = [[None for col in range(len(points))] for row in range(len(points))]
    Neville(points, toFind)

elif method == 4:
    cubicSplineNatural(points, toFind)

elif method == 5:
    cubicSpline(points, toFind, derived1, derived2)
