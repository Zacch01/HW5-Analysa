# Solving Linear Equation Using Gauss Seidel Method


# Global Variable To Store The Machine Precision, (Set the accuracy of the solution)
ACCURACY = 0.00001


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