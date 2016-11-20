#加载软件包

from qutip import *
from pylab import *
import numpy as np
import matplotlib.pyplot as plt
#退相干常数
t1_10=500e-9
t1_01=500e-9
tphy_10=500e-9
tphy_01=500e-9
tlist=linspace(0,10e-6,1001)
psi0=basis(3,0)
#微波驱动项
omega_p=0.1e6*(2*pi)
omega_c=0.5e6*(2*pi)
delta_p=linspace(-5e6,5e6,51)
delta_c=linspace(-5e6,5e6,51)
g=2e6*2*pi
#布居数算符
state0=basis(3,0)*basis(3,0).dag()
state1=basis(3,1)*basis(3,1).dag()
state2=basis(3,2)*basis(3,2).dag()
#主方程的塌缩算符
c_ops=[]
c_ops.append(sqrt(1/t1_10)*Qobj([[0,1,0],[0,0,0],[0,0,0]]))
c_ops.append(sqrt(1/t1_01)*Qobj([[0,0,0],[0,0,1],[0,0,0]]))
c_ops.append(sqrt(2/tphy_10)*Qobj([[0,0,0],[0,1,0],[0,0,0]]))
c_ops.append(sqrt(2/tphy_01)*Qobj([[0,0,0],[0,0,0],[0,0,1]]))
#布居数
result0=zeros([51,51])
result1=zeros([51,51])
result2=zeros([51,51])
#开始计算
i=0
for delta_p_temp in delta_p:
    j=0
    for delta_c_temp in delta_c:
        H=Qobj([[0,omega_p/2,omega_c/2],[omega_p/2,-delta_p_temp*2*pi,g],[omega_c/2,g,-delta_c_temp*2*pi]])
        output=mesolve(H,psi0,tlist,c_ops,[state0,state1,state2])
        result0[i][j]=abs(output.expect[0][-1])
        result1[i][j]=abs(output.expect[1][-1])
        result2[i][j]=abs(output.expect[2][-1])
        j=j+1
    i=i+1
ax=plt.subplot()
p=ax.pcolor(delta_c/1e6,delta_p/1e6,result1+result2,edgecolors='none')
p.set_cmap('jet')
ax.set_xlabel(r'$\Delta_p/MHz$')
ax.set_ylabel(r'$\Delta_c/MHz$')
ax.axis('tight')
plt.colorbar(p)
plt.show()
