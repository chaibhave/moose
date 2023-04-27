# EXAMPLE FOR LINE PLOT ON EXODUS FILES


import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
import sys
path_root = Path(__file__).parents[2]
sys.path.append(str('/Users/bhavcv/projects/TIGER'))
print(sys.path)
from MultiExodusReader import MultiExodusReader

# CLOSE EXISTING PLOTS
plt.close('all')

# SET PLOT PARAMS

plt.rcParams.update({'font.family': 'Arial'})
plt.rc('font', family='sans-serif', weight='bold')

# NAME OF EXODUS FILE TO PLOT
filenames = 'kks_example_nested.e'

# OPEN MULTI EXODUS READER
MF = MultiExodusReader(filenames)

# ALL TIME STEPS IN EXODUS FILES
times = MF.global_times

# GETTING CLOSEST TIME STEP TO DESIRED SIMULATION TIME FOR RENDER --> TYPICALLY 200 FRAMES WITH 20 FPS GIVES A GOOD 10 S LONG VIDEO
n_frames = 201
t_max = 8000.0 ##times[-1]
t_frames = np.linspace(0, t_max, n_frames)
idx_frames = [np.where(times-t_frames[i] == min(times -
                       t_frames[i], key=abs))[0][0] for i in range(n_frames)][1:]

# LOOP OVER EACH TIME STEP IN idx_frames
for (i, time_step) in enumerate(idx_frames):
    print("Rendering frame no. ", i+1)
    # GENERATE FIGURE WINDOW
    # fig = plt.figure()
    fig, ax = plt.subplots(figsize=(3.0, 1.85), dpi=500)
    # Adding Twin Axes to plot using dataset_2
    ax2 = ax.twinx()
    # READ c_Ni DATA at FIRST TIMESTEP
    x, y, z, c = MF.get_data_at_time('c_Ni', times[i])
    # For 1D, EXODUS STORES X ELEMENTS TWICE, SO WE SELECT FIRST SET USING x[:,0]
    ax.plot(x[:, 0], c[:], 'r-', linewidth=1.5)

    # READ c_Cr DATA at FIRST TIMESTEP
    x, y, z, c = MF.get_data_at_time('c_Cr', times[i])
    # For 1D, EXODUS STORES X ELEMENTS TWICE, SO WE SELECT FIRST SET USING x[:,0]
    ax.plot(x[:, 0], c[:], 'b-', linewidth=1.5)

    # READ w_Ni DATA at FIRST TIMESTEP
    x, y, z, c = MF.get_data_at_time('w_Ni', times[i])
    # For 1D, EXODUS STORES X ELEMENTS TWICE, SO WE SELECT FIRST SET USING x[:,0]
    ax2.plot(x[:, 0], c[:], 'r--', linewidth=1.0)
    # READ w_Cr DATA at FIRST TIMESTEP
    x, y, z, c = MF.get_data_at_time('w_Cr', times[i])
    # For 1D, EXODUS STORES X ELEMENTS TWICE, SO WE SELECT FIRST SET USING x[:,0]
    ax2.plot(x[:, 0], c[:], 'b--', linewidth=1.0)

    ax.set_yscale('log')
    #FORMATTING PLOT
    xmin = 0.0
    xmax = 51.0
    ymin = 1e-10
    ymax = 1.5

    ax.set_xlim(xmin,xmax)
    ax.set_ylim(ymin,ymax)
    ax2.set_ylim(-1.0,-0.4)

    ax.margins(x=0, y=0)
    ax.set_ylabel('Atomic fraction',fontsize=7,fontweight='bold')
    ax.set_xlabel('X ($\mu$m)',fontsize=7,fontweight='bold')
    ax2.set_ylabel('Chemical potential (eV)',fontsize=7,fontweight='bold')

    #ADD BORDER LINES FOR AXES
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(1.0)
    legend_properties = {'weight': 'bold', 'size': 5}
    lgd = plt.legend([r"$\bf{c_{Ni}}$", r"$\bf{c_{Cr}}$",r"$\bf{\mu_{Ni}}$", r"$\bf{\mu_{Cr}}$"],loc= 'lower left', prop=legend_properties, ncol=2, framealpha=0)

    fig.tight_layout(pad=0.3)

    # SAVEFIG ALLOWS YOU TO SAVE FIGURE IN DESIRED DPI AND TRANSPARENCY
    plt.savefig('temp/two_phase_two_component_' +
                str(i)+'.png', dpi=500)
    # CLOSE FIGURE AFTER YOU ARE DONE WITH IT. OTHERWISE ALL GENERATED FIGURES WILL BE HELD IN MEMORY TILL SCRIPT FINISHES RUNNING
    plt.close(fig)

"""


xticks=[ int(tick) for tick in np.linspace(xmin,xmax, 7)]
ax.set_xticks(xticks)  # arbitrary chosen
ax.set_xticklabels(xticks,fontsize=6)


yticks= np.around(np.linspace(ymin,ymax, 6),1)
ax.set_yticks(yticks)  # arbitrary chosen
ax.set_yticklabels(yticks,fontsize=6)


ax.tick_params(axis='x',direction='in',length=5)
ax.tick_params(axis='y',direction='in',which='major',length=5)




"""
