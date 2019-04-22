import math
from decimal import Decimal as D

'''solution of the Cauchy problem using Runge-Kutta method
    - solution of the ordinary differential equation 1 and higher orders
    - solution of the system of ordinary differential equations first order
'''


def runge_kutta_method_order_4(f, x_0, y_0, b, N):

    h = (b-x_0)/N
    res_sequence = [(x_0, [y_0, y_0, y_0])]  # contains x points and possible y values

    a = x_0 + h
    while a<=b:

        y_list = []
        x_n = res_sequence[-1][0]

        # 1. y_n+1 = y_n + 1/6 (k1 + 2k2 + 2k3 + k4)
        y_n = res_sequence[-1][1][0]

        # k1 = h*f2(x_n, y_n)
        k1 = h*f(x_n, y_n)
        k2 = h*f(x_n, y_n+(k1/2.0))
        k3 = h*f(x_n, y_n+(k2/2.0))
        k4 = h*f(x_n+h, y_n+k3)
        y_list.append(round(y_n + (k1+2*k2+2*k3+k4)/6.0, 7))

        # 2. y_n+1 = y_n + 1/6(k1+4k3 + k4)
        y_n = res_sequence[-1][1][1]
        k1 = h*f(x_n, y_n)
        k2 = h*f(x_n+(h/4.0), y_n+(k1/4.0))
        k3 = h*f(x_n+(h/2.0), y_n+(k2/2.0))
        k4 = h*f(x_n+h, y_n+k1-2*k2+2*k3)
        y_list.append(round(y_n + (k1+4*k3+k4)/6.0,7))

        # 3. y_n+1 = y_n + 1/8(k1+3k2 + 3k3 + k4)
        y_n = res_sequence[-1][1][2]
        k1 = h*f(x_n, y_n)
        k2 = h*f(x_n + (h/3), y_n + (k1/3))
        k3 = h*f(x_n + (2*h/3), y_n - k1/3 + k2)
        k4 = h*f(x_n + h, y_n+k1-k2+k3)
        y_list.append(round(y_n + (k1+3*k2+3*k3+k4)/8, 7))

        res_sequence.append((a, y_list))
        a = round(a+h, 6)
    return res_sequence

y_derivative = lambda x,y: round(y**2 * math.e**x - 2*y,6)  # Test function
y_solution = lambda x: 1/(math.e**x + math.e**(2*x))  # y* - solution. To compare results

# Runge-Kutta method
N = 20.0
y_0 = 0.5  #  differential equation initial conditions
x_0 = 0  # differential equation initial conditions

a = x_0
b = 2.0
h = (b-a)/N
print(f"[a, b] = [{a}, {b}] h = {h}")
print("\n\t-- Res Runge-Kutta method --")

res = runge_kutta_method_order_4(y_derivative, x_0, y_0, b, N)
for i in res:
    print (f"x={i[0]}\ty1={i[1][0]}\t\ty2={i[1][1]}\t\ty3={i[1][2]}\t\t y*={round(y_solution(float(i[0])),7)} ")
