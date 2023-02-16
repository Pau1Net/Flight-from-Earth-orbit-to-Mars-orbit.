'''
Файл без гомановской(идеальной) траектории. 
Напишите на почту: e49212@voenmeh.ru с домена Военмеха или МГТУ им.Баумана для получения полной версии.
Разработчик: Коцько П.А.(Paul Kotsko)
'''
from numpy import sin, cos, sqrt, sign, arange
from math import pi, asin, acos, sinh, cosh, asinh, acosh, atan
from scipy.optimize import fsolve
import matplotlib.pyplot as plt

mu = 132.712e+9
r0 = 149.6e+6
rk = 227.9e+6
F = 2.20
tp = 270
tp = tp * 24 * 60 * 60

S = sqrt((r0 ** 2) + (rk ** 2) - 2 * r0 * rk * cos(F))
tpar = (((r0 + rk + S) ** (3 / 2)) - ((r0 + rk - S) ** (3 / 2)) * sign(sin(F))) / (6 * sqrt(mu))

Alpha = lambda a: 2 * asinh(sqrt((r0 + rk + S) / (4 * abs(a))))
Betta = lambda a: 2 * asinh(sqrt((r0 + rk - S) / (4 * abs(a))))


def t_el(a):
    global r0, rk, S, tp, mu, F, tm

    d = 2 * asin(sqrt((r0 + rk - S) / (4 * abs(a))))
    eps = 2 * asin(sqrt((r0 + rk + S) / (4 * abs(a))))

    return ((a ** (3 / 2)) * (pi + sign(tm - tp) * (eps - sin(eps) - pi) - sign(sin(F)) * (d - sin(d)))) / sqrt(mu)


def t_hyp(a):
    global mu, F

    alfa = Alpha(a)
    beta = Betta(a)

    return sqrt(((abs(a)) ** 3) / mu) * (sinh(alfa) - alfa - sign(sin(F)) * (sinh(beta) - beta))


def approx(fun_ptr):
    global a1, a2, a3, tp

    while (abs(a2 - a1) > 0.01):
        a3 = (a2 + a1) / 2;

        if (((fun_ptr(a3) - tp) * (fun_ptr(a1) - tp)) >= 0):
            a1 = a3;
        else:
            a2 = a3;


if (tp > tpar):
    print('Эллиптическая орбита:');
    a1 = (r0 + rk + S) / 4;
    if (r0 <= rk):
        a2 = 6 * rk;
    else:
        a2 = 6 * r0;

    dm = 2 * asin(sqrt((r0 + rk - S) / (r0 + rk + S)));
    tm = (((r0 + rk + S) ** (3 / 2)) * (pi - sign(sin(F)) * (dm - sin(dm)))) / (8 * sqrt(mu));

    approx(t_el);

    alfa = Alpha(a3);
    beta = Betta(a3);
    p = (a3 / (S ** 2)) * 4 * r0 * rk * ((sin(F / 2)) ** 2) * (
        sin((alfa + sign(sin(F)) * sign(tm - tp) * beta) / 2)) ** 2;
    e = sqrt(1 - (p / a3));

elif (tp < tpar):
    print('Гиперболическая орбита\n');
    if (r0 <= rk):
        a1 = r0 / 50;
        a2 = rk * 10;
    else:
        a1 = rk / 50;
        a2 = r0 * 10;

    approx(t_hyp);

    alfa = Alpha(a3);
    beta = Betta(a3);
    p = (abs(a3) / (S ** 2)) * 4 * r0 * rk * ((sin(F / 2)) ** 2) * sinh((alfa + sign(sin(F)) * beta) / 2);
    e = sqrt(1 + p / abs(a3));

else:
    print('Параболическая орбита\n');
    p = (sqrt(2) * r0 * rk / (S ** 2)) * ((sin(F / 2)) ** 2) * (sqrt(r0 + rk + S) + sign(sin(F)) * sqrt(r0 + rk - S));
    e = 1;


def func(x):
    return [r0 - x[0] / (1 + x[1] * cos(x[2])),
            rk - x[0] / (1 + x[1] * cos(x[2] + F)),
            x[0] - a3 * (1 - x[1] ** 2)]

print('a = ' + str(a3))
root = fsolve(func, [p, e, 0])

print("p = " + str(root[0]) + ' км')
print("e = " + str(root[1]))
print("Nu = " + str(root[2]))

v0 = sqrt(mu * (2 / r0 + (root[1] ** 2 - 1) / root[0]))
print('v0 = ' + str(v0) + ' км/с');
vk = sqrt(mu * (2 / rk + (root[1] ** 2 - 1) / root[0]))
print('vk = ' + str(vk) + ' км/с');

vz = 30 #Скорость Земли

vr = sqrt(mu / root[0]) * root[1] * sin(root[2])

vt = sqrt(mu / root[0]) * (1 + root[1] * cos(root[2]))

TH = atan(vr / vt)
print('TH = ' + str(TH))

dv = sqrt(v0 ** 2 + vz ** 2 - 2 * v0 * vz * cos(TH))
print('\nОтносительная скорость для перевода КА на межпланетную траекторию:')
print('dv = ' + str(dv) + ' км/с')

plt.plot(0, 0, 'o', markerfacecolor='y', c='yellow', markersize=40, label='Солнце')
plt.plot(r0, 0, 'o', markerfacecolor='g', c = 'green', markersize = 10 , label='нач. положение Земли');
plt.plot(rk * cos(F - (2 * pi / (687 * 86400)) * tp), rk * sin(F - (2 * pi / (687 * 86400)) * tp), 'o',
         markerfacecolor='r', label='нач. положение Марса');
plt.plot(r0 * cos((2 * pi / (365 * 86400)) * tp), r0 * sin((2 * pi / (365 * 86400)) * tp), 'bo',
         label='кон. положение Земли')
plt.plot(rk * cos(F), rk * sin(F), 'o', c = 'red', markersize = 15,  label='кон. положение Марса')

phi = arange(root[2], F + root[2] + 0.01, 0.01)

r = []
for i in range(len(phi)):
    r.append(root[0] / (1 + root[1] * cos(phi[i])));

phi = phi - root[2];
plt.plot(r * cos(phi), r * sin(phi), 'b-', label='перелетная орбита');

phi = arange(0, 2 * pi + 0.01, 0.01)
plt.plot(r0 * cos(phi), r0 * sin(phi), color='gray', linestyle='dashdot');
plt.plot(rk * cos(phi), rk * sin(phi), color='gray', linestyle='dashdot');

plt.grid()
plt.axis('equal')
# plt.legend(loc='lower right', fontsize=6)
# leg = plt.legend(title='Легенда графика', loc='lower right', fontsize=10, markerscale=0.5)
plt.title('Земля-Марс. Гелиоцентрическая система координат')
plt.show()
