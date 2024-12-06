import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import scipy.stats as stats
import math

matplotlib.rcParams.update({'font.size': 14})


col_blue = '#6C8EBF'
col_green = '#82B366'
col_yellow = '#E9CA2C'
col_orange = '#D79B00'
col_red = '#B85450'
col_purple = '#9673A6'
col_grey = '#999999'
col_mid_grey = '#C7C7C7'
col_light_grey = '#F5F5F5'


lw1 = 1.5
lw2 = 1.5

mu = 0
variance = 1
sigma = math.sqrt(variance)
x = np.linspace(mu - 4*sigma, mu + 4*sigma, 100)

fig, ax = plt.subplots(figsize = (10,4), layout='constrained')
ax.set_xticks([-4,-3,-2,-1,0,1,2,3,4], labels = ['μ - 4σ', 'μ - 3σ', 'μ - 2σ', 'μ - 1σ', 'μ', 'μ + 1σ', 'μ + 2σ', 'μ + 3σ', 'μ + 4σ'])

i_label = ['68.27', '95.45', '99.73', '99.99']
i_color = [col_blue, col_orange, col_green, col_red]

for i in range(1,5):
    ax.plot([-i, i],[stats.norm.pdf(-i, mu, sigma), stats.norm.pdf(i, mu, sigma)], label = i_label[i-1]+'%', linewidth = lw1, color = i_color[i-1])
    ax.plot([-i, -i],[stats.norm.pdf(-i, mu, sigma), -0.04], label = i_label[i-1]+'%', linewidth = lw2, linestyle = '--', color = i_color[i-1])
    ax.plot([i, i],[stats.norm.pdf(-i, mu, sigma), -0.04], label = i_label[i-1]+'%', linewidth = lw2, linestyle = '--', color = i_color[i-1])

ax.plot([0, 0],[stats.norm.pdf(0, mu, sigma), -0.04], label = i_label[i-1]+'%', linewidth = lw2, linestyle = '--', color = col_purple)

ax.plot(x, stats.norm.pdf(x, mu, sigma), color = 'black', linewidth = lw2)

plt.text(-0.8,  0.25, '68.27 %', color = col_blue)
plt.text(-0.8,  0.065, '95.45 %', color = col_orange)
plt.text(-0.8,  0.015, '99.73 %', color = col_green)
plt.text(-0.8, -0.022, '99.99 %', color = col_red)

plt.text(-0.3,  0.41, '12.5 ns', color = col_purple)
plt.text(-1.65 ,  0.235, '9.4 ns', color = col_blue)
plt.text( 1.1 ,  0.235, '15.6 ns', color = col_blue)
plt.text(-2.65 ,  0.05, '6.3 ns', color = col_orange)
plt.text( 2.1 ,  0.05, '18.7 ns', color = col_orange)
plt.text(-3.65 ,  0.015, '3.2 ns', color = col_green)
plt.text( 3.1 ,  0.015, '21.8 ns', color = col_green)
plt.text(-4.65 ,  -0.02   , '0.1 ns', color = col_red)
plt.text( 4.1 ,  -0.02  , '24.9 ns', color = col_red)



plt.xlim(-4.8,4.8)
plt.ylim(-0.03, 0.45)
plt.xlabel('x')
plt.ylabel('f(x)')

plt.grid(color=col_mid_grey, linestyle='--', linewidth=1)

# plt.legend()
plt.savefig('annotated_gauss.png', dpi=200)
plt.show()
