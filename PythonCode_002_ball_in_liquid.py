"""
Program: The dynamic simulation of Simple harmonic oscillator and simple pendulum 

Author: Yingli Niu
Email: ylniu@bjtu.edu.cn

license: BSD
Please feel free to use and modify this, but keep the above information. Thanks!
"""
#===============================================================================
import numpy as np
import copy
import math
from matplotlib import pyplot as plt
from matplotlib import animation
#===============================================================================
# 1. 弹簧振子：显式欧拉法
# 2. 弹簧振子：隐式欧拉法
# 3. 弹簧振子：梯形欧拉法
# 4. 单　　摆：速度Verlet蛙跳法
itype = 1
#-------------------------------------------------------------------------------
m     = 15      # 质量
k     = 0.5     # 劲度系数
x     = [1, 0]  # 初始位置 x
v     = [0, 0]  # 初始速度 x
w     = (k/m)**0.5
dt    = 0.1     # 时间步长
nt    = 1
g     = 9.8
l     = 50
#-------------------------------------------------------------------------------
pi    = 3.1415926
rho0  = 1000.0
rho   = 7850.0
R     = 0.01
g     = 9.8
eta   = 1.0*10**(-3)
V     = 4/3*pi * R**2
m1    = rho  * V
m0    = rho0 * V
#-------------------------------------------------------------------------------
bb    = (m1*g - m0*g)/(6*pi * eta * R)
aa    = 6*pi*eta * R/m1
#-------------------------------------------------------------------------------
xx    = []
vv    = []
#-------------------------------------------------------------------------------
# Constants
a = dt
b = -k/m * dt
c = 1.0 / (1.0 + (w*dt)**2)
p = 1.0 / (1.0 + (w*dt/2.0)**2)
q =       (1.0 - (w*dt/2.0)**2)
#-------------------------------------------------------------------------------
if (itype==1):
   title="Explicit Euler Method"
   xlable="x"
   ylable="y"
elif (itype==2):
   title="Implicit Euler Method"
   xlable="x"
   ylable="y"
elif (itype==3):
   title="Trapezoid Euler Method"
   xlable="x"
   ylable="y"
elif (itype==4):
   title="Velocity Verlet Leap-frog Method"
   xlable="theta"
   ylable="omega"
#-------------------------------------------------------------------------------
# fig = plt.figure(1)
#-------------------------------------------------------------------------------
# ax1 = fig.add_subplot(211, xlim = (-2,2), ylim = (-2,2))
# line1, = ax1.plot(x[0], x[1], 'og', markersize=10)
# ax1.set_title(title)
#-------------------------------------------------------------------------------
# ax2 = fig.add_subplot(212, xlim = (-2,2), ylim = (-0.5,0.5))
# ax2.set_xlabel(xlable)
# ax2.set_ylabel(ylable)
# line2, = ax2.plot(xx, vv, linestyle="-", label='u(t)', markersize=1)
#===============================================================================
def fun(x):
   return [-g/l * math.sin(x[0]), -g/l * math.sin(x[1])]
#===============================================================================
def dynamics(fun, k, m, x, v, dt, nt, itype):
   for n in range(nt):
      x0 = copy.deepcopy(x)
      v0 = copy.deepcopy(v)
      f0 = fun(x0)
      for i in range(len(x)):
         if (itype==1):
            x[i] =      x0[i] + a * v0[i]
            v[i] =  b * x0[i] +     v0[i]
         elif(itype==2):
            x[i] = (    x0[i] + a * v0[i]) * c
            v[i] = (b * x0[i] +     v0[i]) * c
         elif(itype==3):
            x[i] = (q * x0[i] + a * v0[i]) * p
            v[i] = (b * x0[i] + q * v0[i]) * p
         elif(itype==4):
            x[i] = x0[i] + v0[i] * dt + 0.5 * f0[i] * dt**2
            x1   = copy.deepcopy(x)
            f1   = fun(x1)
            v[i] = v0[i] + 0.5 * (f0[i]+f1[i]) * dt
#===============================================================================
def plot_ball(i, x):
   line1.set_data(x[0], x[1])
   return line1,
#===============================================================================
def plot_line(i, x):
   x0 = copy.deepcopy(x)
   v0 = copy.deepcopy(v)
   xx.append(x0[0])
   vv.append(v0[0])
   line2.set_data(xx, vv)
   return line2,
#===============================================================================
def updateALL(i):
   dynamics(fun, k, m, x, v, dt, nt, itype)
   a = plot_ball(i, x)
   b = plot_line(i, x)
   return a + b
#===============================================================================
# anim  = animation.FuncAnimation(fig, updateALL, frames=1, interval=10, blit=True, repeat=True)
# plt.show()
#===============================================================================
#
v0 = 0.0
dt = 0.1
tt0=[]

vv1=[]
vv2=[]
vvd=[]
vv1.append(v0)
vv2.append(v0)
vd0=0.0
vvd.append(vd0)
tt0.append(0.0)
tmax = 200
nmax = 100000
dt = tmax / nmax
for i in range(nmax):
   t  = i * dt
   v1 = v0 + dt * a * (b - v0)
   v2 = b * (1-math.exp(-a*t))
   tt0.append(t)
   vd0=v2 - v1
   vv2.append(v2)
   vv1.append(v1)
   vvd.append(vd0)
   v0 = v1


#plt.plot(tt0,vv1,label="Euler forward method")
#plt.plot(tt0,vv2,label="Analytical result")
plt.plot(tt0,vvd,label="Error")
plt.legend(loc=1)
plt.xlabel('t')
plt.ylabel('v')

plt.show()