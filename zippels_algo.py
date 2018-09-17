# @author: Ativ Joshi

import random
import math
from sympy import *

prime = 10 ** 9 + 7

"""
pulverizer
"""
def pulverizer(a, b):
    if a == 0:
        return 0, 1

    [x, y] = pulverizer(b % a, a)

    return y - (b / a) * x, x

'''
GCD
'''

def gcd(a,b):

    if a == 0:
        return b
    else:
        return gcd(b%a,a)


'''
 dense interpolation algorithm 'D'. Gives polynomial that fits all the given points.
 polynomial is such that,
 f(p_i) = m_i
'''
def dense_interpolation(p, m):
    x = symbols('x')
    f = m[0] * x ** 0
    # print type(f)
    q = x - p[0]

    for i in range(1, len(p)):
        temp = q.subs(x, p[i])
        qpi_inv = q.subs(x, p[i]) ** (-1)
        f += qpi_inv * q * (m[i] - f.subs(x, p[i]))
        q = (x - p[i]) * q
        # print 'f: ', f, '\nq: ', q, '\nqpi_inv: ', qpi_inv, '\n', 'temp: ', temp, '\n'
    if f == nan:
        return sympify(1)
    else:
        return expand(f)





'''ORACLE function 'F' which gives output of function at given points'''
def oracle(l):
    x, y, z = symbols('x y z')
    # func = 5*(x**4)*(y**2)*(z**3) + (x**3*y*z**2) + x + 8*y*z*x**2 + 3*x**5*y**5*z**5
    func = (y**2*z**2 + y**2*z + 2*x**2*y*z + x*z)*y**2*z**2 + y**2*z + 2*x**2*y*z + x*z
    return func.subs({x: l[0], y: l[1], z: l[2]})

def oracle_gcd(l):
    x, y, z = symbols('x y z')

    d = y**2*z**2 + y**2*z + 2*x**2*y*z + x*z
    f = z**2 + y**2*z + x**2*y*z + x*z + x**2*y**2
    g = y*z + 2*x*z + z + x

    f1 = f*d
    f2 = g*d

    return gcd(f1.subs({x: l[0], y: l[1], z: l[2]}), f2.subs({x: l[0], y: l[1], z: l[2]}))



'''sparse interpolation algorithm 'S' '''
def sparse_interpolation(x_var, a, d):
    S = [sympify(1)]
    p0 = sympify(oracle(a))
    # print p0, ' p0'
    L = []

    '''iterate through each variable'''
    for i in range(0, len(x_var)):
        # print 'i: ', i, '\n'
        r = []
        X = x_var[i]
        P = []

        '''iterate 'd' times corresponding to number of degrees '''
        for j in range(0, d):
            # print '         j=',j,'\n'

            '''pick rj'''
            temp = random.randint(1, 10000)
            # print 'r_j', temp
            r.append(temp)

            '''list for linear equations'''
            L = []
            t = len(S)
            skeletal = sympify(0)
            G = symbols('g0:%d'%(t+1))

            '''generate a skeletal polynomial from S'''
            for k in range(0, t):
                skeletal += G[k] * S[k]
            # print 'skeletal: ', skeletal

            '''iterate 't' times where t is number of monomials in skeletal polynomial'''
            for k in range(0, t):
                # print A, 'A..............................'

                '''
                for i=0, there is no need to solve system of linear equations
                we append the output of oracle [F] in the list 'P'

                otherwise, we need to create a list 'L' of 't' linear equations
                '''
                if i == 0:
                    oracle_out = oracle([r[j]] + a[1:len(a)])
                    # print oracle_out, ' oracle'
                    P.append(oracle_out)
                    sub_list = [(x_var[tem], A[tem]) for tem in range(i)]
                    L.append(skeletal.subs(sub_list) - oracle_out)
                else:
                    A = random.sample(xrange(10000), i)
                    # print 'A: ', A
                    or_lis = A + [r[j]] + a[i+1:len(a)]
                    oracle_out = oracle(or_lis)
                    sub_list = [(x_var[tem], A[tem]) for tem in range(i)]
                    L.append(skeletal.subs(sub_list)-oracle_out)
                    # print 'L: ', L
                    # sol1 = linsolve(L, G[0], G[1], G[2])
                    # print 'sol1: ', sol1
                    # solution = next(iter(sol1))
                    # print 'solution: ', len(solution[0].free_symbols), '\n\n'

                    # while len(solution[0].free_symbols) != 0:
                    #     print 'while'
                    #     A = random.sample(xrange(10000), i)
                    #     or_lis = A + [r[j]] + a[i + 1:len(a)]
                    #     oracle_out = oracle(or_lis)
                    #     sub_list = [(x_var[tem], A[tem]) for tem in range(i)]
                    #     L.append(skeletal.subs(sub_list) - oracle_out)
                    #     solution = next(iter(linsolve(L, G[0], G[1], G[2])))
                    #     print 'solution: ', solution, '\n\n'
            # print 'L: ', L

            '''
            if i != 0, then solve the system of linear equation and produce a polynomial pj
            append the produced polynomial to the list P
            '''
            if i != 0:
                sol1 = linsolve(L, G)
                solution = next(iter(sol1))
                # print 'solution: ', len(solution[0].free_symbols)
                # print solution, '\n'
                temp_p = sympify(0)

                for mon in range(len(S)):
                    temp_p += S[mon]*solution[mon]
                P.append(temp_p)


        '''
        For i = 0,  pass the list p directly to the dense_interpolation algorithm 'D'.
        Store the coefficients to the list 'S' and the resulting polynomial in p0

        If i is not 0, then for each monomial is 'S', pass the corresponding coefficients from
        list 'P' and from p0 to the dense_interpolation algorithm 'D'.
        Merge each monomial to its corresponding result from 'D' and simplify
        '''
        if i == 0:
            # print [p0] + P
            # print [p0] + P, 'P\n', [a[0]] +r, 'r'
            p0 = sympify(dense_interpolation([a[0]] + r, [p0] + P).subs(symbols('x'), x_var[0]))
            S = p0.as_coefficients_dict().keys()
            # print 'p0:', expand(p0), '\n', 'S: ', S, '\n'
            # print 'p0: ', p0, '\nS: ', S, '\n'
            # print p0.as_coefficients_dict(), ' as coeff dict'
        else:
            P = [p0] + P
            r = [a[i]]+r
            p0 = sympify(0)
            # print P, 'P\n', r, 'r'
            for mon in range(len(S)):
                p_i_temp = []
                # m_i_temp = []

                for pol in range(len(P)):
                    p_i_temp += [P[pol].as_coefficients_dict()[S[mon]]]
                # print p_i_temp, 'p_i_\n', r, 'r\n'
                interp_temp =  dense_interpolation(r,p_i_temp).subs(symbols('x'), x_var[i])
                # print 'interp_temp: ',  interp_temp
                p0 += S[mon]*interp_temp
            p0 = expand(p0)
            S = p0.as_coefficients_dict().keys()
            # print 'p0:', expand(p0), '\n', 'S: ', S, '\n'
    return p0



X = symbols('x0:3')
a = [1,2,3]
print sparse_interpolation(X, a, 4)
# print random.sample(xrange(1000),0)
# f_int = dense_interpolation([1, 2, 3, 4, 5, 6], [2, 9, 28, 65, 126, 217])
# print f_int
# coeff = f_int.as_coefficients_dict()
#
# for t in coeff:
#     coeff[t] %= prime
#
# f_new = sympify(0)
#
# for t, c in coeff.iteritems():
#     f_new += c*t
#
# print f_new
