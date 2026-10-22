"""Plot p_L vs p for one Steane EC round (flagged extraction) vs a bare qubit."""
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

res = json.load(open('../results/steane_ec_results.json'))
A = res['counts']['A']['any']
mc = res['monte_carlo']

INK, INK2, MUTED, GRID = '#0b0b0b', '#52514e', '#8a8984', '#e6e5e1'
STEANE, BARE = '#2a78d6', '#52514e'

plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11,
                     'axes.edgecolor': MUTED, 'axes.labelcolor': INK2,
                     'xtick.color': INK2, 'ytick.color': INK2})
fig, ax = plt.subplots(figsize=(7.5, 5.2), dpi=200)
fig.patch.set_facecolor('#fcfcfb')
ax.set_facecolor('#fcfcfb')

p = np.logspace(-4, -1, 200)
ax.plot(p, p, color=BARE, lw=2, ls=(0, (5, 3)), label='Bare physical qubit  ($p_L = p$)')
ax.plot(p, A * p**2, color=STEANE, lw=2, label=f'Steane EC round, pair counting  ($p_L \\approx {A:.0f}\\,p^2$)')
ps = np.array([m['p'] for m in mc])
pl = np.array([m['p_L'] for m in mc])
er = np.array([m['stderr'] for m in mc])
ax.errorbar(ps, pl, yerr=er, fmt='o', ms=6, color=STEANE, mec='#fcfcfb', mew=1.5,
            ecolor=STEANE, elinewidth=1.2, capsize=0, label='Steane EC round, Monte Carlo', zorder=5)

pth = 1 / A
ax.axvline(pth, color=MUTED, lw=1, ls=':', zorder=0)
ax.annotate('pseudothreshold\n$p^* \\approx %.1f \\times 10^{-3}$' % (pth*1e3), xy=(pth, 2e-5), xytext=(pth * 1.15, 2e-5),
            color=INK2, fontsize=10, va='center')
ax.text(2.2e-4, 3.5e-6, 'Steane better', color=STEANE, fontsize=10)
ax.text(1.5e-2, 1.6e-3, 'bare qubit better', color=INK2, fontsize=10)

ax.set_xscale('log'); ax.set_yscale('log')
ax.set_xlim(1e-4, 1e-1); ax.set_ylim(1e-6, 1)
ax.set_xlabel('Physical error rate per operation, $p$')
ax.set_ylabel('Logical error rate per EC round, $p_L$')
ax.set_title('One Steane [[7,1,3]] EC round (flagged extraction, 9 qubits)',
             color=INK, fontsize=12, loc='left', pad=12)
ax.grid(True, which='major', color=GRID, lw=0.8)
for s in ('top', 'right'):
    ax.spines[s].set_visible(False)
leg = ax.legend(loc='upper left', frameon=False, fontsize=9.5, labelcolor=INK2)
fig.text(0.01, 0.01, 'Circuit-level depolarizing noise, same p for CNOT, reset and measurement; '
         'no idle or shuttling noise yet.', fontsize=8, color=MUTED)
fig.tight_layout(rect=(0, 0.03, 1, 1))
fig.savefig('../figures/steane_ec_threshold.png', facecolor=fig.get_facecolor())
print('saved')
