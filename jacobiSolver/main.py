import sys
import numpy as np
import numpy.linalg as al
import copy

# coefficient = np.array([[12., 3., -5.], [1., 5., 3.], [3., 7., 13.]])
# b = np.array([1., 28., 76.])

coefficient = np.array([[1., -2., 3.], [2., 0., -2.], [-7., 3., 0.]])
b = np.array([10., 11., 3.])

# coefficient = np.array([[4., 2., 1.], [-1., 2., 0.], [2., 1., 4.]])
# b = np.array([11., 3., 16.])

# coefficient = np.array([[8., 4., -1.], [-2., 3., 1.], [2., -1., 6.]])
# b = np.array([11., 4., 7.])

# ?????????coefficient = np.array([[1., -3., 2., 1.], [-2., -6., 1., 4.], [-1., 2., 3., 4.], [0., -1., 1., 1.]])
# ????????b = np.array([-4., 1., 12., 0.])

numOfVar = 3
initialGuess = np.array([0., 0., 0.])
numOfIterations = 100000
absolute_relative_error = 0.05


def check_zero(coefficient):
    if(np.any(np.diag(coefficient) == 0)):
        return True
    return False

def avoid_zero_on_diagonal(coefficient, b, numOfVar):
    if (not check_zero(coefficient)): return

    for i in range(numOfVar):
        if coefficient[i][i] != 0:
            continue
        for j in range(numOfVar):
            if coefficient[j][i] != 0:
                coefficient[[i, j]] = coefficient[[j, i]]
                b[i], b[j] = b[j], b[i]
                break
        if(not check_zero(coefficient)):
            return [coefficient, b]
    if (check_zero(coefficient)):
        print('Cannot be solved')
        sys.exit()




def jacobi(coefficient, b, initialGuess, numOfIterations, absolute_relative_error):
    k = 1
    print(coefficient)
    oldX = initialGuess
    x = copy.deepcopy(oldX)
    if check_zero(coefficient):
        modified = avoid_zero_on_diagonal(coefficient, b, numOfVar)
        coefficient, b = modified[0], modified[1]
        print(coefficient)
        print(b)
    diagonal_matrix = np.diag(np.diag(coefficient))
    inv_diagonal = np.linalg.inv(diagonal_matrix)
    a = (coefficient-diagonal_matrix)

    while k <= numOfIterations:
        x = np.dot(inv_diagonal, b-(np.dot(a, x)))
        # absError = abs(np.divide((x-oldX), x))*100
        if np.any(x == 0):
            x_with_no_zeros = np.where(x == 0, 1e-10, x)
            absError = abs(np.divide((x - oldX), x_with_no_zeros)) * 100
        else:
            absError = abs(np.divide((x - oldX), x)) * 100
        if np.all(absError <= absolute_relative_error):
            print(x)
            break
        oldX = copy.deepcopy(x)
        k = k+1

    if k == (numOfIterations+1):
        print('Cannot be solved')


# Seidel method

def seidel(coefficient, b, initialGuess, numOfIterations, absolute_relative_error, numOfVar):
    k = 1
    oldX = initialGuess
    x = copy.deepcopy(oldX)
    if check_zero(coefficient):
        modified = avoid_zero_on_diagonal(coefficient, b, numOfVar)
        coefficient = modified[0]
        b = modified[1]
        print(coefficient)
        print(b)

    while k<=numOfIterations:
        absError =np.zeros(numOfVar)
        for i in range(numOfVar):
            sum = 0.0
            for j in range(numOfVar):
                if (j!=i):
                    sum = sum + coefficient[i][j]*x[j]
            x[i] = (b[i]-sum)/coefficient[i][i]
            #absError[i] = abs((x[i] - oldX[i]) / x[i]) * 100
            if np.any(x == 0):
                x_with_no_zeros = np.where(x == 0, 1e-10, x)
                absError = abs(np.divide((x - oldX), x_with_no_zeros)) * 100
            else:
                absError = abs(np.divide((x - oldX), x)) * 100
        if np.all(absError <= absolute_relative_error):
            print(x)
            break
        oldX = copy.deepcopy(x)
        k=k+1

    if k == (numOfIterations+1):
        print('Cannot be solved')


jacobi(coefficient, b, initialGuess, numOfIterations, absolute_relative_error)
seidel(coefficient, b, initialGuess, numOfIterations, absolute_relative_error, numOfVar)







