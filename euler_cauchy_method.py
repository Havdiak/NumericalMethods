'''

Euler-Cauchy method, Euler-Cauchy method with refinement
 
One of the methods of solving the Cauchy
problem for ordinary first-order differential equations
on [a,b]

Returns a sequence that coincides with the integral curve solution 
of differential equation 1 order y' = f(x,y)
'''

import math


def euler_cauchy(f, x_0, y_0, a, b, N=10.0):
    h = (b-a)/N
    res_sequence=[(x_0, y_0)]
    N+=1

    i = a+h
    while i<=b:
        res_sequence.append((round(i,7),))
        i+=h

    i = 1
    while i<N:
        x_n = res_sequence[i-1][0]
        y_n = res_sequence[i-1][1]
        
        y_with_hat = y_n + h*y_derivative(x_n, y_n)
        x_with_hat = res_sequence[i][0]

        y_new = y_n + h* ((y_derivative(x_n, y_n) + y_derivative(x_with_hat,y_with_hat))/2.0)
        res_sequence[i] += (round(y_new,6),)

        i+=1
    return res_sequence


def euler_cauchy_iterative_refinement(f, x_0, y_0, a, b, N=15):
    res_sequence=[(x_0, y_0)]
    N+=1
    h = (b-a)/N

    i = a+h
    while i<=b:
        res_sequence.append((round(i,7),))
        i+=h
    
    i = 1
    while i<N:
        x_n = res_sequence[i-1][0]
        y_n = res_sequence[i-1][1]
        
        x_with_hat = res_sequence[i][0]
        y_s_with_hat = y_n + h*y_derivative(x_n, y_n)

        for s in range(4): # 4 is enough, if result is bad make h smaller
            y_s1_with_hat = y_n + (h/2.0) * (y_derivative(x_with_hat, y_s_with_hat) + y_derivative(x_n, y_n))
            if abs(y_s1_with_hat - y_s_with_hat) <= epsilon:
                res_sequence[i] += (round(y_s1_with_hat, 6),)
                break
            else:
                y_s_with_hat = y_s1_with_hat

        y_with_hat = y_n + h*y_derivative(x_n, y_n)
        x_with_hat = res_sequence[i][0]

        y_new = y_n + h* ((y_derivative(x_n, y_n) + y_derivative(x_with_hat,y_with_hat))/2.0)
        res_sequence[i] += (round(y_new,6),)
        i+=1

    return res_sequence


y_derivative = lambda x,y: round(y**2 * math.e**x - 2*y,6)  # Test function
y_solution = lambda x: 1/(math.e**x + math.e**(2*x))  # y* - solution. To compare results

# Euler-Cauchy method
N = 10.0
y_0 = 0.5  #  differential equation initial conditions
x_0 = 0  # differential equation initial conditions

a = 0
b = 2.0
print(f"[a, b] = [{a}, {b}]")
print("\n\t-- res Euler-Cauchy --")

res_euler_cauchy = euler_cauchy(y_derivative, x_0, y_0, a, b, N)
for point in res_euler_cauchy:
    print(f"x = {point[0]}\t\ty = {point[1]}\ty_solution = {round(y_solution(point[0]),6)}")


# Euler-Cauchy  method with refinement. Iterative refinement of each found value

N = 19.0
epsilon = 0.000005

print("\n\t-- res Euler-Cauchy method with refinement --")
res_sequence = euler_cauchy_iterative_refinement(y_derivative, x_0, y_0, a, b, N)
for point in res_sequence:
    print(f"x = {point[0]}\t\ty_my = {point[1]}\t\ty_solution = {round(y_solution(point[0]),6)}")
