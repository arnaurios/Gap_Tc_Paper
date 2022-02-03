# coding: utf-8
import numpy as np

def fitfunction(x,k0,k1,k2,k3,d0) :
    conds = [x < k0, (x > k0) & (x < k2), x > k2]

    funcs = [lambda x: 0.,
             lambda x: d0 * np.power(x-k0,2) * np.power(x-k2,2) / ( np.power(x-k0,2) + np.power(k1,2) ) / ( np.power(x-k2,2) + np.power(k3,2) ),
             lambda x: 0.]

    f = np.piecewise(x, conds, funcs)

    return f


def derfitfunction(x,k0,k1,k2,k3,d0) :
    ff = np.zeros([len(x),5])
    conds = [x < k0, (x > k0) & (x < k2), x > k2]

#    funcs = [lambda x: [0.,0.,0.,0.,0.],
#             lambda x:
#             [ 2. * d0 * np.power(x-k0,3) * np.power(x-k2,2) / np.power( ( np.power(x-k0,2) + np.power(k1,2) ),2 ) / ( np.power(x-k2,2) + np.power(k3,2) ) - 2. * d0 * (x-k0) * np.power(x-k2,2) / ( np.power(x-k0,2) + np.power(k1,2) ) / ( np.power(x-k2,2) + np.power(k3,2) ),
#             -2. * d0 * k1 * np.power(x-k0,2) * np.power(x-k2,2) / np.power( ( np.power(x-k0,2) + np.power(k1,2) ),2 ) / ( np.power(x-k2,2) + np.power(k3,2) ),
#               2. * d0 * np.power(x-k0,2) * np.power(x-k2,3) / ( np.power(x-k0,2) + np.power(k1,2) ) / np.power( ( np.power(x-k2,2) + np.power(k3,2) ),2 ) - 2. * d0 * np.power(x-k0,2) * (x-k2) / ( np.power(x-k0,2) + np.power(k1,2) ) / ( np.power(x-k2,2) + np.power(k3,2) ),
#             -2. * d0 * k3 * np.power(x-k0,2) * np.power(x-k2,2) / ( np.power(x-k0,2) + np.power(k1,2) ) / np.power( ( np.power(x-k2,2) + np.power(k3,2) ), 2),
#             np.power(x-k0,2) * np.power(x-k2,2) / ( np.power(x-k0,2) + np.power(k1,2) ) / ( np.power(x-k2,2) + np.power(k3,2) )],
#             lambda x: [0.,0.,0.,0.,0.]]

    func0 = [lambda x: 0.,
             lambda x: 2. * d0 * np.power(x-k0,3) * np.power(x-k2,2) / np.power( ( np.power(x-k0,2) + np.power(k1,2) ),2 ) / ( np.power(x-k2,2) + np.power(k3,2) ) - 2. * d0 * (x-k0) * np.power(x-k2,2) / ( np.power(x-k0,2) + np.power(k1,2) ) / ( np.power(x-k2,2) + np.power(k3,2) ),
             lambda x: 0.]

    func1 = [lambda x: 0.,
             lambda x: -2. * d0 * k1 * np.power(x-k0,2) * np.power(x-k2,2) / np.power( ( np.power(x-k0,2) + np.power(k1,2) ),2 ) / ( np.power(x-k2,2) + np.power(k3,2) ),
             lambda x: 0.]

    func2 = [lambda x: 0.,
             lambda x:  2. * d0 * np.power(x-k0,2) * np.power(x-k2,3) / ( np.power(x-k0,2) + np.power(k1,2) ) / np.power( ( np.power(x-k2,2) + np.power(k3,2) ),2 ) - 2. * d0 * np.power(x-k0,2) * (x-k2) / ( np.power(x-k0,2) + np.power(k1,2) ) / ( np.power(x-k2,2) + np.power(k3,2) ),
             lambda x: 0.]

    func3 = [lambda x: 0.,
             lambda x: -2. * d0 * k3 * np.power(x-k0,2) * np.power(x-k2,2) / ( np.power(x-k0,2) + np.power(k1,2) ) / np.power( ( np.power(x-k2,2) + np.power(k3,2) ), 2),
             lambda x: 0.]

    func4 = [lambda x: 0.,
            lambda x:  np.power(x-k0,2) * np.power(x-k2,2) / ( np.power(x-k0,2) + np.power(k1,2) ) / ( np.power(x-k2,2) + np.power(k3,2) ) ,
            lambda x: 0.]

    ff[:,0] = np.piecewise(x, conds, func0)
    ff[:,1] = np.piecewise(x, conds, func1)
    ff[:,2] = np.piecewise(x, conds, func2)
    ff[:,3] = np.piecewise(x, conds, func3)
    ff[:,4] = np.piecewise(x, conds, func4)

    return ff

def fitmodel(beta,x) :
    k0,k1,k2,k3,d0 = beta
    conds = [x < k0, (x > k0) & (x < k2), x > k2]

    funcs = [lambda x: 0.,
             lambda x: d0 * np.power(x-k0,2) * np.power(x-k2,2) / ( np.power(x-k0,2) + np.power(k1,2) ) / ( np.power(x-k2,2) + np.power(k3,2) ),
             lambda x: 0.]

    #if( x < k0 ) :
    #    f=0
    #elif( x > k0 & x < k2 ) :
    #    f=d0 * np.power(x-k0,2) * np.power(x-k2,2) / ( np.power(x-k0,2) + np.power(k1,2) ) / ( np.power(x-k2,2) + np.power(k4,2) )
    #else :
    #    f=0

    f = np.piecewise(x, conds, funcs)
    return f
