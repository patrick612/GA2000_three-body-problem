import numpy as np
import matplotlib.pyplot as plt

def two_body_time(phi, e):
    f_x = 1/np.pi * (2*np.arctan(np.power((1-e)/(1+e), 1/2) * np.tan(phi/2)) -
                     (e * np.power((1 - e**2), 1/2)*np.sin(phi))/(1 + e*np.cos(phi)))
    return f_x

def two_body_time_deriv(phi, e):
    frac = np.power((1-e)/(1+e), 1/2)
    e1 = e * np.power(1-e**2, 1/2)
    deriv = 1 / np.pi * (2/(np.power((frac)*np.tan(phi/2), 2)+1) * 1/2 * frac * np.power(1/np.cos(phi/2), 2)
                 - e1 * np.cos(phi)/(1 + e* np.cos(phi)) - e1 * np.sin(phi)/(np.power(1 + e* np.cos(phi), 2))*e*np.sin(phi)
                 )
    return deriv

y_list = []
for i in np.linspace(0, 2*np.pi, 1000):
    y = two_body_time(i, 0.2)
    y_list.append(y)

plt.plot(np.linspace(0, 2*np.pi, 1000), y_list)
plt.show

y_list = []
for i in np.linspace(0, 2*np.pi, 1000):
    y = two_body_time_deriv(i, 0.2)
    y_list.append(y)

plt.plot(np.linspace(0, 2*np.pi, 1000), y_list)
plt.show

def newton_root(t, e):
    epsilon = 1/np.power(10, 5)
    phi = np.random.rand(1)[0] * np.pi
    while True:
        diff = np.absolute(two_body_time(phi, e) -t)
        if epsilon>= diff:
            break
        phi = phi - (two_body_time(phi, e)-t)/two_body_time_deriv(phi, e)

    return phi

def radius(t, a, e):
    root_phi = newton_root(t, e)
    rad = a * (1 - np.power(e, 2))/(1 + e * np.cos(root_phi))
    return rad
