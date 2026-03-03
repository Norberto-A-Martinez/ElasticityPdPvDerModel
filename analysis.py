#%% INITIAL IMPORTS
import pandas as pd
import matplotlib.pyplot as plt

#%% SUBSTATION POWER DELIVERY
sub = pd.read_csv('results/Substation_WithoutAllocation_.csv')

t = sub['T'].unique()
s = sub['S'].unique()
d = sub['D'].unique()

for i in s:
    for j in d:
        sub_sd = sub[(sub['S'] == i) & (sub['D'] == j)]
        plt.figure(figsize=(8,4))
        plt.plot(sub_sd['T'], sub_sd['Pg'], label='Pg')
        plt.plot(sub_sd['T'], sub_sd['Qg'], label='Qg')
        plt.ylabel('Power (kW)')
        plt.xlabel('Time')
        plt.xticks(t,rotation=90)
        plt.title(f'Substation Power Delivery for S{i} D{j}')
        plt.grid(alpha=0.3)
        plt.legend(ncols=2)
        plt.tight_layout()
        plt.savefig(f"figures/se_power_delivery_{i}{j}.svg")
        plt.close()
#%% VOLTAGE
data = pd.read_csv('Voltage_WithoutAllocation_.csv')