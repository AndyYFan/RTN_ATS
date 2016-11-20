from qutip import *
from pylab import *
import numpy as np
import matplotlib.pyplot as plt
from scipy import io
import time
import random

time_paras = 1501
kai = 1e8
def rtnGen(timelist, kai):
    seq = zeros(time_paras+100)

    sign = np.random.randint(0, 2)  # generate random 0,1

    tj = random.expovariate(kai)
    index = 0

    for i in range(time_paras):
        if tj < timelist[i]:
            sign = ~sign + 2
            tj = tj + random.expovariate(kai)
        seq[index] = sign
        index = index + 1
    return seq


# 退相干常数
t1_01 = 500e-9
t1_12 = 500e-9
tphy_11 = 500e-9
tphy_22 = 500e-9
tlist = linspace(0, 10e-6, time_paras)
psi0 = basis(3, 0)

# 微波驱动项
omega_p = 0.1e6 * (2 * pi)
omega_c = 14e6 * 2 * pi
delta_p = linspace(-10e6, 10e6, 101)
# 布居数算符
state0 = basis(3, 0) * basis(3, 0).dag()
state1 = basis(3, 1) * basis(3, 1).dag()
state2 = basis(3, 2) * basis(3, 2).dag()
# 主方程的塌缩算符
c_ops = []
c_ops.append(sqrt(1 / t1_01) * Qobj([[0, 1, 0], [0, 0, 0], [0, 0, 0]]))
c_ops.append(sqrt(1 / t1_12) * Qobj([[0, 0, 0], [0, 0, 1], [0, 0, 0]]))
c_ops.append(sqrt(2 / tphy_11) * Qobj([[0, 0, 0], [0, 1, 0], [0, 0, 0]]))
c_ops.append(sqrt(2 / tphy_22) * Qobj([[0, 0, 0], [0, 0, 0], [0, 0, 1]]))
# 布居数
result0 = zeros([len(delta_p), time_paras])
result1 = zeros([len(delta_p), time_paras])
result2 = zeros([len(delta_p), time_paras])

start_time = time.time()
# 开始计算
# i=0
# for detuning_temp in detuning:
#     delta_c1=detuning_temp
#     delta=-2*delta_c1
j = 0
sequence = rtnGen(tlist, kai)
for delta_p_temp in delta_p:

    H0 = Qobj([[0, omega_p / 2, 0], [omega_p / 2, -delta_p_temp * 2 * pi, 0], [0, 0, -(delta_p_temp) * 2 * pi]])
    H1 = Qobj([[0, 0, 0], [0, 0, omega_c / 2], [0, omega_c / 2, 0]])


    def H1_coeff(tilist, args):
        return sequence[int(tilist*(time_paras-1)*1e5)]


    H = [H0, [H1, H1_coeff]]
    output = mesolve(H, psi0, tlist, c_ops, [state0, state1, state2])
    result0[j] = abs(output.expect[0])
    result1[j] = abs(output.expect[1])
    result2[j] = abs(output.expect[2])
    percent = (j + 1) / len(delta_p) * 100
    end_time = time.time()
    timediff = end_time - start_time
    timeremain = (100 - percent) / percent * timediff
    sys.stdout.write("%.2f" % percent + '% 已进行' + "%.2f" % timediff + 's 剩余' + "%.2f" % timeremain + 's\r')
    sys.stdout.flush()
    j = j + 1
# i=i+1
data = result1 + result2
io.savemat('E:/python3.4/ATS_RTN_test/'+str(kai/1e6)+'M.mat', {'data': data})

