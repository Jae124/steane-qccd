import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
INK='#1b2433'; MUTED='#6b7280'; ACC='#2a78d6'; FLAG='#d95926'; BG='#fbfbf8'
fig, ax = plt.subplots(figsize=(14, 8.0), dpi=200)
fig.patch.set_facecolor(BG); ax.set_facecolor(BG)
rows = ['q4','q5','q6','q7','a (syndrome)','f (flag)']
ys = {r: -i for i, r in enumerate(rows)}
x0, x1 = 0.6, 10.2
for r in rows:
    c = FLAG if r.startswith('f') else (ACC if r.startswith('a') else INK)
    ax.plot([x0, x1], [ys[r]]*2, color=c, lw=2, zorder=1)
    ax.text(x0-0.15, ys[r], r, ha='right', va='center', fontsize=17, color=c)
def prep(x, r, label, col):
    ax.add_patch(plt.Rectangle((x-0.35, ys[r]-0.28), 0.7, 0.56, fc=BG, ec=col, lw=2, zorder=3))
    ax.text(x, ys[r], label, ha='center', va='center', fontsize=14, color=col, zorder=4)
def cnot(x, c, t, col=INK):
    ax.plot([x, x], [ys[c], ys[t]], color=col, lw=2, zorder=2)
    ax.add_patch(plt.Circle((x, ys[c]), 0.11, color=col, zorder=3))
    ax.add_patch(plt.Circle((x, ys[t]), 0.24, fc=BG, ec=col, lw=2, zorder=3))
    ax.plot([x-0.24, x+0.24], [ys[t]]*2, color=col, lw=2, zorder=4)
    ax.plot([x, x], [ys[t]-0.24, ys[t]+0.24], color=col, lw=2, zorder=4)
prep(1.3, 'a (syndrome)', '|0⟩', ACC); prep(1.3, 'f (flag)', '|+⟩', FLAG)
cnot(2.6, 'q4', 'a (syndrome)')
cnot(3.8, 'f (flag)', 'a (syndrome)', FLAG)
cnot(5.0, 'q5', 'a (syndrome)')
cnot(6.2, 'q6', 'a (syndrome)')
cnot(7.4, 'f (flag)', 'a (syndrome)', FLAG)
cnot(8.6, 'q7', 'a (syndrome)')
prep(9.6, 'a (syndrome)', 'M_Z', ACC); prep(9.6, 'f (flag)', 'M_X', FLAG)
ax.axvspan(4.4, 6.8, ymin=0.05, ymax=0.97, color=FLAG, alpha=0.07, zorder=0)
ax.text(5.6, 0.55, 'hook window: one fault on a here can spread to q6 and q7', ha='center', fontsize=14, color=FLAG)
ax.set_xlim(-0.9, 10.5); ax.set_ylim(-5.6, 0.95); ax.set_aspect('equal'); ax.axis('off')
fig.tight_layout(); fig.savefig('../figures/flagged_circuit.png', facecolor=BG)
