import numpy as np 
from matplotlib import pyplot as plt

def hatfunction(x):
    y = np.copy(x)
    y[np.abs(y)<0.5] = 1
    y[np.abs(y)>0.5] = 0
    y[np.abs(y)==0.5] = 0.5
    return y

def main(x):
    
    pass

if __name__ == '__main__':
    main()
