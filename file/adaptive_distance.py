import numba 
import numpy as np
import warnings
warnings.filterwarnings("ignore")

@numba.njit(fastmath=True)
def adaptive_euclidean(x, y):
    """Standard euclidean distance.

    ..math::
        D(x, y) = \sqrt{\sum_i (x_i - y_i)^2}
    """
    result = 0.0
    npcs = int((len(x)-1)/2)
    if x[0] != y[0]:
        for i in range(1,(npcs+1)):
            result += (x[i] - y[i]) ** 2
        d = np.sqrt(result)
            
    else:
        for i in range((npcs+1),(2*npcs+1)):
            result += (x[i] - y[i]) ** 2
        d = np.sqrt(result)
        
    return d


@numba.njit(fastmath=True)
def adaptive_euclidean_grad(x, y):
    """Standard euclidean distance and its gradient.

    ..math::
        D(x, y) = \sqrt{\sum_i (x_i - y_i)^2}
        frac{dD(x, y)}{dx} = (x_i - y_i)/D(x,y)
    """
    result = 0.0
    npcs = int((len(x)-1)/2)
    if x[0] != y[0]:
        for i in range(1,(npcs+1)):
            result += (x[i] - y[i]) ** 2
        d = np.sqrt(result)
        grad = [(x[i]-y[i])/ (1e-6 + d) for i in range(1,(npcs+1))]
    else:
        for i in range((npcs+1),(2*npcs+1)):
            result += (x[i] - y[i]) ** 2
        d = np.sqrt(result)
        grad = [(x[i]-y[i])/ (1e-6 + d) for i in range((npcs+1),(2*npcs+1))]
    return d, grad




@numba.njit(fastmath=True)
def adaptive_euclidean2(x, y):
    """Standard euclidean distance.

    ..math::
        D(x, y) = \sqrt{\sum_i (x_i - y_i)^2}
    """
    result = 0.0
    npcs = int((len(x)-2)/3)
    if x[0] != y[0]:
        for i in range(2,(npcs+2)):
            result += (x[i] - y[i]) ** 2
        d = np.sqrt(result)
            
    else:
        if x[1] != y[1] or x[1] == -1:
            for i in range((npcs+2),(2*npcs+2)):
                result += (x[i] - y[i]) ** 2
            d = np.sqrt(result)
        else:
            for i in range((2*npcs+2),(3*npcs+2)):
                result += (x[i] - y[i]) ** 2
            d = np.sqrt(result)
        
    return d


@numba.njit(fastmath=True)
def adaptive_euclidean2_grad(x, y):
    """Standard euclidean distance and its gradient.

    ..math::
        D(x, y) = \sqrt{\sum_i (x_i - y_i)^2}
        frac{dD(x, y)}{dx} = (x_i - y_i)/D(x,y)
    """
    result = 0.0
    npcs = int((len(x)-2)/3)
    if x[0] != y[0]:
        for i in range(2,(npcs+2)):
            result += (x[i] - y[i]) ** 2
        d = np.sqrt(result)
        grad = [(x[i]-y[i])/ (1e-6 + d) for i in range(2,(npcs+2))]
            
    else:
        if x[1] != y[1] or x[1] == -1:
            for i in range((npcs+2),(2*npcs+2)):
                result += (x[i] - y[i]) ** 2
            d = np.sqrt(result)
            grad = [(x[i]-y[i])/ (1e-6 + d) for i in range((npcs+2),(2*npcs+2))]
        else:
            for i in range((2*npcs+2),(3*npcs+2)):
                result += (x[i] - y[i]) ** 2
            d = np.sqrt(result)
            grad = [(x[i]-y[i])/ (1e-6 + d) for i in range((2*npcs+2),(3*npcs+2))]

    return d, grad




