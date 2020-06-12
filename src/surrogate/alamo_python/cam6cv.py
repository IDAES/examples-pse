def f(X, B):
    import numpy as np
    x1= X[0]
    x2= X[1]
    return  B[0] * x1**2 + B[1] * x2**2 + B[2] * x1**4 + B[3] * x2**4 + B[4] * x1**6 + B[5] * x1*x2 
