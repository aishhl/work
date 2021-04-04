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
itype = 4
#-------------------------------------------------------------------------------
m     = 15        # 质量
k     = 0.5       # 劲度系数
#x01   = [1, 0]   # 初始位置 x01
#v01   = [0, 0]   # 初始速度 x01
x01   = [1, 1]  # 初始位置 x01
v01   = [-1, 0]   # 初始速度 x01
x02   = [-1, -1]  # 初始位置 x01
v02   = [1, 0]    # 初始速度 x01
v01   = [-0.15, 0]
v02   = [ 0, 0]
w     = (k/m)**0.5
dt    = 0.1       # 时间步长
nt    = 1
g     = 9.8
l     = 50
G0    = 6.67*10**(-11)
mass  = []
mass.append(10**(2))
mass.append(10**(9))
#-------------------------------------------------------------------------------
xx1    = []
vv1    = []
xx2    = []
vv2    = []
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
   xlable="x01"
   ylable="y"
elif (itype==2):
   title="Implicit Euler Method"
   xlable="x01"
   ylable="y"
elif (itype==3):
   title="Trapezoid Euler Method"
   xlable="x01"
   ylable="y"
elif (itype==4):
   title="Velocity Verlet Leap-frog Method"
   xlable="theta"
   ylable="omega"
#-------------------------------------------------------------------------------
fig = plt.figure(1)
#-------------------------------------------------------------------------------
xmax=10
ax  = fig.add_subplot(111, xlim = (-xmax,xmax), ylim = (-xmax,xmax))
b1, = ax.plot(x01[0], x01[1], 'og', markersize=10)
b2, = ax.plot(x02[0], x02[1], 'og', markersize=10)
ax.set_title(title)
#===============================================================================
def fun(x01, x02):
   r=math.sqrt((x01[0]-x02[0])**2 + (x01[1]-x02[1])**2)
   dr=[]
   er=[]
   dr.append(x02[0]-x01[0])
   dr.append(x02[1]-x01[1])
   er.append(dr[0] / r)
   er.append(dr[1] / r)
   fr = G0 * mass[0] * mass[1] / r**2
   return [er[0]*fr, er[1]*fr]
#===============================================================================
def dynamics(fun, k, m, x01, v01, x02, v02, dt, nt, itype):
   for n in range(nt):
      tmp_x01 = copy.deepcopy(x01)
      tmp_v01 = copy.deepcopy(v01)
      tmp_x02 = copy.deepcopy(x02)
      tmp_v02 = copy.deepcopy(v02)
      tmp_f01 = fun(tmp_x01, tmp_x02)
      tmp_f02 = fun(tmp_x02, tmp_x01)
      for i in range(len(x01)):
         if (itype==1):
            x01[i] =      tmp_x01[i] + a * tmp_v01[i]
            v01[i] =  b * tmp_x01[i] +     tmp_v01[i]
         elif(itype==2):
            x01[i] = (    tmp_x01[i] + a * tmp_v01[i]) * c
            v01[i] = (b * tmp_x01[i] +     tmp_v01[i]) * c
         elif(itype==3):
            x01[i] = (q * tmp_x01[i] + a * tmp_v01[i]) * p
            v01[i] = (b * tmp_x01[i] + q * tmp_v01[i]) * p
         elif(itype==4):
            x01[i] = tmp_x01[i] + tmp_v01[i] * dt + 0.5 * tmp_f01[i]/mass[0] * dt**2
            x02[i] = tmp_x02[i] + tmp_v02[i] * dt + 0.5 * tmp_f02[i]/mass[1] * dt**2
            tmp_x01= copy.deepcopy(x01)
            tmp_x02= copy.deepcopy(x02)
            f1     = fun(tmp_x01, tmp_x02)
            f2     = fun(tmp_x02, tmp_x01)
            v01[i] = tmp_v01[i] + 0.5 * (tmp_f01[i]+f1[i])/mass[0] * dt
            v02[i] = tmp_v02[i] + 0.5 * (tmp_f02[i]+f2[i])/mass[1] * dt
#===============================================================================
def plot_ball(i, x01, x02):
   b1.set_data(x01[0], x01[1])
   b2.set_data(x02[0], x02[1])
   return b1, b2,
#===============================================================================
def updateALL(i):
   dynamics(fun, k, m, x01, v01, x02, v02, dt, nt, itype)
   a = plot_ball(i, x01, x02)
   return a
#===============================================================================
anim  = animation.FuncAnimation(fig, updateALL, frames=1, interval=10, blit=True, repeat=True)
plt.show()