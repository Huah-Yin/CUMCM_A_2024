## 2022 CUMCM -A 数学建模国赛A题—训练

#### 该项目针对学习数学建模 A类赛题的队伍学习参考

## Abstract

 ![](https://img.shields.io/badge/模型-常微分方程组-blue.svg?style=flat)![](https://img.shields.io/badge/解法-RK4 and GA-green.svg?style=flat)

### 问题一

求解得到一组微分方程组(其中$x_1为浮子的位移,x_2为振子的位移$)
$$
\begin{equation}
\left\{
\begin{aligned}
    &(m_1 + m_c) \frac{d^2 x_1}{dt^2} + k_1 x_1 + k_2 (x_1 - x_2) + k_3 \left( \frac{dx_1}{dt} - \frac{dx_2}{dt} \right) + k_4 \frac{dx_2}{dt} + A \cos w = 0 \\
    &m_2 \frac{d^2 x_2}{dt^2} + k_2 (x_2 - x_1) + k_3 \left( \frac{dx_2}{dt} - \frac{dx_1}{dt} \right) = 0
\end{aligned}
\right.

\end{equation}
$$
对方程降阶：
$$
\frac{dX}{dt} =
\begin{bmatrix}
    \frac{dx_1}{dt} \\
    \frac{dx_2}{dt} \\
    \frac{dx_3}{dt} \\
    \frac{dx_4}{dt}
\end{bmatrix}
= f(t, x) =
\begin{cases}
    \begin{array}{c} 
            x_3 \\
            x_4 \\
            \frac{-k_1 x_1 - k_2 (x_1 - x_2) - k_3 \left( \frac{dx_1}{dt} - \frac{dx_2}{dt} \right) - k_4 \frac{dx_1}{dt} + A \cos \omega t}{(m_1 + m_c)} \\
            \frac{k_2 (x_1 - x_2) + k_3 \left( \frac{dx_2}{dt} - \frac{dx_1}{dt} \right)}{m_2}
    \end{array}
\end{cases}
$$


这个方程可以通过$ode45\text{ 自适应步长函数，求解快速方便}$，其中，我们有在网站上找到有队伍通过$Laplace \text{ 变化}$的方式求解，可以求得一组解析解，我们使用$Mathematica$ 来求解，得出的解析解含有复数域，有兴趣的可以尝试

![image-20241112180026970](C:/Users/Mr.yin/AppData/Roaming/Typora/typora-user-images/image-20241112180026970.png)

### 问题二

要使得PTO系统的平均输出功率最大，我们需要得到系统稳定后的状态，以此得到最大的功率.
$$
\overline{P_p} = \frac{1}{T} \int_{t_0}^{t_0 + T} M_{cp} [X'_z(t) - X'_f(t)] \, dt
= \frac{1}{T} \int_{t_0}^{t_0 + T} C_p |X'_z(t) - X'_f(t)|^2 \, dt
$$
阻尼系数和阻尼系数非恒定的状态:
1.恒定状态只有一个参数，可以直接遍历，也可以融合GA算法，但是GA耗费时间较长，遍历更为简洁。

2.系数非恒定状态，一共有两个参数，目标函数为二元函数，将时间离散化后($\;dt=0.08\;$)，目标函数选择$(3)$

积分区间选择$[80,100]$ s,可以再对GA进行一些改进，从而使得算法可以收敛到目标函数最高点:$229\; W$。

### 问题三

纵摇力学分析，太难太难太难，难上加难…………………………………………………!，**but:但是它对最终功率的影响占比较小**，这里不叙述,有兴趣的同学可以自行推导(我们的推导存在问题，果断放弃！)

### 问题四

类比问题二，两个参数的取值，目标函数为二元函数，将时间仍旧离散化，采取$GA+trapz()\text{组合技}$

顺利得出$\;P_{max}=329.5494W$.

请注意：我们队伍当时GA采取小步长$dt=0.005$的 步长，算法会耗费大量的时间，最终我们采取了$\text{大步长+大区间}$的计算方法，算法历经$\;326.5\;s$顺利得到最优解.小步长时间的搜索在直观上会显得准确，但是实际上会耗费大量的时间，使得你每隔$10min$就得查看结果是否得出. 

