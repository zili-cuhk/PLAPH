"""Kinked-MIC geometry figure (main text, Figure 1): normalized penalty tanh(a u^2 + xi |u|)
and its right derivative for u > 0 at a = 90.  Usage: python3 kmic_geometry.py OUT.pdf"""
import sys
import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
a, xis = 90.0, [5, 10, 15, 20, 25.6685, 30, 40]
u = np.linspace(-0.3, 0.3, 1201); up = np.linspace(1e-6, 0.3, 1200)
fig, ax = plt.subplots(1, 2, figsize=(8, 3.2))
cm = plt.get_cmap('viridis')
for i, xi in enumerate(xis):
    col, lw = ('#c0504d', 2.4) if xi == 25.6685 else (cm(i / len(xis)), 1.0)
    ax[0].plot(u, np.tanh(a * u**2 + xi * np.abs(u)), color=col, lw=lw, label=f'$\\xi={xi:g}$')
    ax[1].plot(up, (2 * a * up + xi) / np.cosh(a * up**2 + xi * up)**2, color=col, lw=lw)
ax[0].plot(u, np.tanh(a * u**2), color='gray', ls='--', lw=1, label=r'$\tanh(au^2)$')
ax[1].plot(up, 2 * a * up / np.cosh(a * up**2)**2, color='gray', ls='--', lw=1)
ax[0].set(xlabel='$u$', ylabel=r'$p_K(u)/\lambda$', title='(a) normalized penalty')
ax[1].set(xlabel='$u>0$', ylabel=r"$p_K'(u)/\lambda$", title='(b) right derivative')
ax[0].legend(fontsize=7, ncol=2)
fig.tight_layout(); fig.savefig(sys.argv[1], bbox_inches='tight')
