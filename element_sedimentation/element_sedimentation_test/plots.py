import matplotlib.pyplot as plt
plt.rcParams["legend.frameon"] = False
plt.rcParams["font.size"] = 14
import mesa_helper as mh

sims = mh.SimulationSeries('LOGS')

x_history = 'star_age'
history_keys = ['radius_Jup', 'Teff', 'e_thermal','total_entropy']

x_profile = 'mass_Jup'
profile_keys = ['y', 'entropy']

fig, axs = plt.subplots(2, 3, figsize=(16, 8))
axs = axs.flatten()

i = 0
for key in history_keys:   
    sims.history_plot(x_history, key, ax=axs[i], set_label=True, set_axes_labels=True)
    axs[i].legend()
    i += 1

for key in profile_keys:


    if key == 'y':
        sims.profile_composition_plot(x_profile, ['x', key], function_y=lambda X, Y: Y/(X+Y), ax=axs[i], set_label=True, set_axes_labels=True)
        axs[i].axhline(0.238-0.005, color='gray', linestyle='--')
        axs[i].axhline(0.238, color='gray', label='Y_Jup')
        axs[i].axhline(0.238+0.005, color='gray', linestyle='--')
    else:
        sims.profile_plot(x_profile, key, ax=axs[i], set_label=True, set_axes_labels=True)

    axs[i].legend()
    i += 1

fig.tight_layout()
# plt.legend()
plt.show()