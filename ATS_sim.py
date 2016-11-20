#加载软件包

from qutip import *
from pylab import *
import numpy as np
import matplotlib.pyplot as plt
from scipy import io
#退相干常数
t1_01=500e-9
t1_12=500e-9
tphy_11=500e-9
tphy_22=500e-9
tlist=linspace(0,10e-6,31)
psi0=basis(3,0)

#微波驱动项
omega_p=0.5e6*(2*pi)
omega_c=linspace(0.16e6*2*pi,43.8e6*2*pi,31)
delta_p=linspace(-20e6,20e6,31)
#布居数算符
state0=basis(3,0)*basis(3,0).dag()
state1=basis(3,1)*basis(3,1).dag()
state2=basis(3,2)*basis(3,2).dag()
#主方程的塌缩算符
c_ops=[]
c_ops.append(sqrt(1/t1_01)*Qobj([[0,1,0],[0,0,0],[0,0,0]]))
c_ops.append(sqrt(1/t1_12)*Qobj([[0,0,0],[0,0,1],[0,0,0]]))
c_ops.append(sqrt(2/tphy_11)*Qobj([[0,0,0],[0,1,0],[0,0,0]]))
c_ops.append(sqrt(2/tphy_22)*Qobj([[0,0,0],[0,0,0],[0,0,1]]))
#布居数
result0=zeros([31,31])
result1=zeros([31,31])
result2=zeros([31,31])
#开始计算
i=0
for omega_c_temp in omega_c:
    j=0
    for delta_p_temp in delta_p:
        H=Qobj([[0,omega_c_temp/2,0],[omega_c_temp/2,0,omega_p/2],[0,omega_p/2,-delta_p_temp*2*pi]])
        output=mesolve(H,psi0,tlist,c_ops,[state0,state1,state2])
        result0[i][j]=abs(output.expect[0][-1])
        result1[i][j]=abs(output.expect[1][-1])
        result2[i][j]=abs(output.expect[2][-1])
        j=j+1
    i=i+1
data=result1
io.savemat('E:/python3.4/ATS.mat',{'data':data})
ax=plt.subplot()
p=ax.pcolor(delta_p/1e6,omega_c/2e6/pi,result1+result2,edgecolors='none')
p.set_cmap('jet')
ax.set_ylabel(r'$\Omega_c/MHz$')
ax.set_xlabel(r'$\Delta_p/MHz$')
ax.axis('tight')
plt.colorbar(p)
plt.show()
