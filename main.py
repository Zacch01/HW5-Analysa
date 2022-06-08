from math import pi

# Global Variable To Store The Machine Precision, (Set the accuracy of the solution)
ACCURACY = 0.00001

# Solving Linear Equation Using Inverse Matrix Method


def InverseMatrix(originMatrix, originVectorB):
    """
    Solving linear equation in the Inverse Matrix method

    :param originMatrix: NxN Matrix
    :param originVectorB: Nx1 Vector
    """
    # Check if the matrix is Quadratic matrix, and check if the vector is in appropriate size
    if len(originMatrix) == len(originMatrix[0]) and len(originVectorB) == len(originMatrix) and len(originVectorB[0]) == 1:

        # In case the matrix has one solution
        if determinantMatrix(originMatrix):

            # Organize the matrix pivots
            originMatrix, originVectorB = organizeMatrix(originMatrix, originVectorB)

            # Getting the inverse matrix of originMatrix
            inverseMatrix = findInverse(originMatrix)

            # Getting the Linear Equation solution
            vectorSolution = finalSolution(originMatrix, originVectorB, multiplyMatrix(inverseMatrix, originVectorB))

            # Saving the Linear Equation final solution
            return list(map(lambda x: x[0], vectorSolution))

        # According message In case there is more or less than one solution
        else:
            return None

    # In case the input Linear Equation isn't meet the demands
    else:
        return None


def organizeMatrix(originMatrix, originVectorB):
    """
    Taking care that the pivot in the every row will be the highest possible, and return the updated Linear Equation

    :param originMatrix: NxN matrix
    :param originVectorB: Nx1 vector
    :return: The updated Linear Equation
    """
    # Loop to get the highest pivots possible
    for i in range(len(originMatrix)):

        # Variable to store the highest value for the pivot
        maxPivot = abs(originMatrix[i][i])

        # Variable to store the new pivot row
        pivotRow = -1

        # Searching the highest potential Pivot for originMatrix[i][i]
        for j in range(i + 1, len(originMatrix)):

            # In case there's a higher pivot (on the Column[i])
            if abs(originMatrix[j][i]) > maxPivot:
                maxPivot = abs(originMatrix[j][i])
                pivotRow = j

        # In case there was a higher pivot, change the matrix so the Pivot will be the maximum
        if maxPivot != abs(originMatrix[i][i]):
            originVectorB[i], originVectorB[pivotRow] = originVectorB[pivotRow], originVectorB[i]
            originMatrix[i], originMatrix[pivotRow] = originMatrix[pivotRow], originMatrix[i]

    # Return the updated Linear Equation
    return originMatrix, originVectorB


def findInverse(matrix):
    """
    Solve the matrix into an Identity matrix, and return the inverse matrix

    :param matrix: NxN matrix
    :return: Inverse NxN matrix
    """
    # Initialize inverseMatrix into an Identity matrix
    inverseMatrix = [[1.0 if row == col else 0.0 for col in range(len(matrix))] for row in range(len(matrix))]

    # Solving matrix into an Identity matrix, and get alongside the Inverse Matrix (Lower part)
    for i in range(len(matrix)):

        # In case the pivot isn't one, we will make sure it will be
        if matrix[i][i] != 1.0:
            inverseMatrix = multiplyMatrix(initElementaryMatrix(len(matrix), i, i, 1 / matrix[i][i]), inverseMatrix)
            matrix = multiplyMatrix(initElementaryMatrix(len(matrix), i, i, 1 / matrix[i][i]), matrix)

        for j in range(i + 1, len(matrix)):
            if matrix[j][i] != 0.0:
                inverseMatrix = multiplyMatrix(initElementaryMatrix(len(matrix), j, i, - matrix[j][i]), inverseMatrix)
                matrix = multiplyMatrix(initElementaryMatrix(len(matrix), j, i, - matrix[j][i]), matrix)

    # Solving matrix into an Identity matrix, and get alongside the Inverse Matrix (Upper part)
    for i in reversed(range(len(matrix))):
        for j in reversed(range(i)):
            if matrix[j][i] != 0:
                inverseMatrix = multiplyMatrix(initElementaryMatrix(len(matrix), j, i, - matrix[j][i]), inverseMatrix)
                matrix = multiplyMatrix(initElementaryMatrix(len(matrix), j, i, - matrix[j][i]), matrix)

    # Return the inverse matrix
    return inverseMatrix


def finalSolution(originMatrix, originVectorB, vectorSolution):
    """
    Getting the Linear equation components, check the accuracy of the solution, if the accuracy isn't precise
    calculate the precise solution and return it

    :param originMatrix: NxN matrix
    :param originVectorB: Nx1 vector
    :param vectorSolution: Nx1 vector semi solution (not surly accurate)
    :return: Nx1 vector, the precise Linear Equation solution
    """
    # Solve r = Ax0 - b (Vector r represent the accuracy of the solution we found)
    vectorR = multiplyMatrix(originMatrix, vectorSolution)
    for i in range(len(vectorR)):
        vectorR[i][0] = vectorR[i][0] - originVectorB[i][0]

    # Update to the correct solution
    for i in range(len(vectorSolution)):
        if abs(vectorSolution[i][0] - round(vectorSolution[i][0])) <= max(1e-09 * max(abs(vectorSolution[i][0]), abs(round(vectorSolution[i][0]))), 0.0):
            vectorSolution[i][0] = round(vectorSolution[i][0])

    # Return the final solution of the Linear Equation
    return vectorSolution


def multiplyMatrix(matrixA, matrixB):
    """
    Multiplying two matrices and return the outcome matrix

    :param matrixA: NxM Matrix
    :param matrixB: NxM Matrix
    :return: NxM matrix
    """
    # Initialize NxM matrix filled with zero's
    matrixC = [[0.0] * len(matrixB[0]) for _ in range(len(matrixA))]

    # Multiply the two matrices and store the outcome in matrixC
    for i in range(len(matrixA)):
        for j in range(len(matrixB[0])):
            for k in range(len(matrixB)):
                matrixC[i][j] = matrixC[i][j] + matrixA[i][k] * matrixB[k][j]

    # Return the outcome matrix
    return matrixC


def initElementaryMatrix(size, row, col, value):
    """
    Initialize elementary matrix, from identity matrix, and a specific value, and return it

    :param size: Matrix size
    :param row: Row index
    :param col: Column index
    :param value: Value parameter
    :return: Return the elementary matrix
    """
    # Initialize the desire elementary matrix
    elementaryMatrix = [[1.0 if row == col else 0.0 for col in range(size)] for row in range(size)]
    elementaryMatrix[row][col] = value

    # Return the elementary matrix
    return elementaryMatrix


def determinantMatrix(matrix):
    """
    Calculate the matrix determinant and return the result

    :param matrix: NxN Matrix
    :return: Matrix determinant
    """
    # Simple case, The matrix size is 2x2
    if len(matrix) == 2:
        value = matrix[0][0] * matrix[1][1] - matrix[1][0] * matrix[0][1]
        return value

    # Initialize our sum variable
    determinantSum = 0

    # Loop to traverse each column of the matrix
    for current_column in range(len(matrix)):
        sign = (-1) ** current_column

        # Calling the function recursively to get determinant value of sub matrix obtained
        determinant_sub = determinantMatrix(
            [row[: current_column] + row[current_column + 1:] for row in (matrix[: 0] + matrix[0 + 1:])])

        # Adding the calculated determinant value of particular column matrix to total the determinantSum
        determinantSum = determinantSum + (sign * matrix[0][current_column] * determinant_sub)

    # Returning the final Sum
    return determinantSum


def GaussSeidel(originMatrix, originVectorB):
    """
    Solving Linear Equation in the Gauss Seidel method

    :param originMatrix: NxN Matrix
    :param originVectorB: Nx1 Vector
    """
    # Check if the matrix is Quadratic matrix, and check if the vector is in appropriate size
    if len(originMatrix) == len(originMatrix[0]) and len(originVectorB) == len(originMatrix) and len(originVectorB[0]) == 1:

        # In case the matrix has one solution
        if determinantMatrix(originMatrix):

            # Organize the matrix pivots
            originMatrix, originVectorB = organizeMatrix(originMatrix, originVectorB)

            # Our lists for the Prev iteration values, and our Current iteration values
            prevIteration = [[0 for _ in range(1)] for _ in range(len(originMatrix))]
            currentIteration = [[0 for _ in range(1)] for _ in range(len(originMatrix))]

            # Loop for finding the solution
            for _ in range(300):

                # Calculate the next guess
                for i in range(len(originMatrix)):
                    rowSum = 0
                    for j in range(len(originMatrix)):
                        if i != j:
                            rowSum = rowSum + originMatrix[i][j] * currentIteration[j][0]
                    currentIteration[i][0] = (originVectorB[i][0] - rowSum) / originMatrix[i][i]

                # In case we found the solution, Stop the program
                if all([False if abs(currentIteration[row][0] - prevIteration[row][0]) > ACCURACY else True for row in range(len(currentIteration))]):
                    break

                # Update the previous solution to be the current solution
                prevIteration = [[currentIteration[row][0] for _ in range(1)] for row in range(len(currentIteration))]

                # According message in case the Matrix is not converge
                if _ == 299:
                    return None
            # Saving the Linear Equation final solution
            return list(map(lambda x: x[0], currentIteration))

        # According message In case there is more or less than one solution
        else:
            return None

    # In case the input Linear Equation isn't meet the demands
    else:
        return None


def isDiagonalDominant(matrix):
    """
    Check if the pivot in every row is bigger than the sum of the whole row (without the pivot),
    If yes return True, else False

    """
    for i in range(len(matrix)):

        # Variable to store, the summation of absolute row [i]
        rowSum = 0
        for j in range(len(matrix)):
            if i != j:
                rowSum = rowSum + abs(matrix[i][j])

        # If the summation of the row is bigger than the pivot, return False (The matrix is not diagonal dominant)
        if rowSum > abs(matrix[i][i]):
            return False

    # The matrix is Diagonal Dominant
    return True


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
    print("Vector coefficients")
    print(coefficient)
    print()
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

    mySum = 0
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
    if len(pointsList) < 4:
        print('Neville approximation must have minimum of 4 points')
        return

    if xToFind < pointsList[0][0] or xToFind > pointsList[-1][0]:
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
        print("Recursive Neville")
        print(P[i + 1][j])
    # Saving the calculation of P[i][j - 1]
    if P[i][j - 1] is None:
        P[i][j - 1] = recursiveNeville(pointsList, xToFind, i, j - 1)
        print("Recursive Neville")
        print(P[i][j - 1])


    # Create a sub calculating
    return ((xToFind - pointsList[i][0]) * P[i + 1][j] - (xToFind - pointsList[j][0]) * P[i][j - 1]) / (pointsList[j][0] - pointsList[i][0])


def cubicSpline(points, toFind, derived1, derived2):
    if len(points) < 4:
        print('Cubic spline approximation must have minimum of 4 points')
        return

    if toFind < points[0][0] or toFind > points[-1][0]:
        print("The wanted point isn't suitable for interpolation")
        return

    h = [0 for _ in range(len(points))]
    lam = [0 for _ in range(len(points))]
    u = [0 for _ in range(len(points))]
    d = [0 for _ in range(len(points))]
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

    matrix = [[0 for _ in range(len(points))] for _ in range(len(points))]
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

    h = [0 for _ in range(len(points))]
    lam = [0 for _ in range(len(points))]
    u = [0 for _ in range(len(points))]
    d = [0 for _ in range(len(points))]
    for i in range(len(points)-1):
        h[i] = (points[i+1][0] - points[i][0])

    for i in range(1, len(points)):
        lam[i] = h[i] / (h[i-1] + h[i])
        u[i] = 1 - lam[i]

    for i in range(1, len(points)-1):
        d[i] = (6/(h[i-1]+h[i]))*((points[i+1][1]-points[i][1])/h[i] - (points[i][1]-points[i-1][1])/h[i-1])

    d[len(points) - 1] = d[0] = u[len(points) - 1] = lam[0] = 0

    matrix = [[0 for _ in range(len(points))] for _ in range(len(points))]
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

    matrix2 = [[0 for _ in range(len(matrix)-2)] for _ in range(len(matrix)-2)]
    for row in range(1, len(matrix)-1):
        for col in range(1, len(matrix)-1):
            matrix2[row-1][col-1] = matrix[row][col]
    vector2 = [[d[row] for _ in range(1)] for row in range(1, len(points)-1)]

    coefficient = GaussSeidel(matrix2, vector2)

    M = [[0 for _ in range(1)] for _ in range(len(points))]
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

table_points = [[1.2, 1.5095], [1.3, 1.6984], [1.4, 1.9043], [1.5, 2.1293], [1.6, 2.3756]]
XtoFind = 1.47
functionDerived1 = 1
functionDerived2 = 0


method = -1
while method > 5 or method < 0:
    method = int(input('0 --> Linear\n1 --> Polynomial\n2 --> Lagrange\n3 --> Neville\n4 --> Natural Cubic Spline\n5 --> Full Cubic Spline\nInput --> '))

    if method > 5 or method < 0:
        print('Invalid input, Try [Zero to Five]')

if method == 0:
    Linear(table_points, XtoFind)

elif method == 1:
    polynomial(table_points, XtoFind)

elif method == 2:
    Lagrange(table_points, XtoFind)

elif method == 3:
    P = [[None for col in range(len(table_points))] for row in range(len(table_points))]
    Neville(table_points, XtoFind)

elif method == 4:
    cubicSplineNatural(table_points, XtoFind)

elif method == 5:
    cubicSpline(table_points, XtoFind, functionDerived1, functionDerived2)
