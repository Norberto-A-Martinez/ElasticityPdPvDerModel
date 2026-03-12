#%%
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

load = pd.read_csv("input-data/parametros_load_curvas.csv")
load_prob = pd.read_csv("input-data/parametros_load_probs.csv")
pv = pd.read_csv("input-data/parametros_pv_curvas.csv")
pv_prob = pd.read_csv("input-data/parametros_pv_probs.csv")

month = 1
load = load[load['M'] == month]
load_prob = load_prob[load_prob['M'] == month]
pv = pv[pv['M'] == month]
pv_prob = pv_prob[pv_prob['M'] == month]

fig, ax = plt.subplots(2, 1, figsize=(8, 5))
ax[0].set_title(f'a)', fontsize='small', loc='left')
for s in range(1, 4):
    load_scenario = load[load['D'] == s]
    ax[0].plot(load_scenario['T'], load_scenario['curve_load'], label=f'{s} - $\epsilon$ = {load_prob[load_prob["D"] == s]["prob_D"].values[0]:.2f}')
ax[0].set_ylabel('$p.u.$')
ax[0].set_xticks(range(1, 25))  # Set x-ticks to be from 1 to 24
ax[0].set_xticklabels([])
ax[0].grid(alpha=0.3)
ax[0].legend(title='Scenario:', loc='upper left')

ax[1].set_title(f'b)', fontsize='small', loc='left')
for s in range(1, 4):
    pv_scenario = pv[pv['S'] == s]
    ax[1].plot(pv_scenario['T'], pv_scenario['curve_pv'], label=f'{s} - $\epsilon$ = {pv_prob[pv_prob["S"] == s]["prob_S"].values[0]:.2f}')
ax[1].set_xlabel('Time ($t$)')
ax[1].set_ylabel('$p.u.$')
ax[1].set_xticks(range(1, 25))  # Set x-ticks to be from 1 to 24
ax[1].grid(alpha=0.3)
ax[1].legend(title='Scenario:', loc='upper left')
plt.tight_layout()
plt.savefig('figures/curvas_scenarios.svg') 
plt.show()
