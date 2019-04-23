# https://www.python-course.eu/matrix_arithmetic.php
# https://docs.scipy.org/doc/numpy/user/quickstart.html
import numpy.matlib
import numpy.matmul
import numpy as np

a = np.array([1, 2, 3])
b = np.array([3, 4, 5])
np.dot(a, b)

c = np.matlib.zeros((2, 4))
c

d = np.matlib.full((2, 2), 10)
np.matrix(d)

a = np.matrix(np.identity(3))
b = np.matlib.zeros((3, 4))
# multiplicação de matrizes
np.dot(a, b)
