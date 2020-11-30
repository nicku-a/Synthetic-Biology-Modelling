# Importing useful libraries
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from datetime import datetime
import scipy.stats.mstats as mstats
import os
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as cm
from google.colab import drive
drive.mount('/content/drive')

import seaborn as sns
sns.set_style("dark")
sns.set()

from IPython.display import display, HTML

display(HTML(data="""
<style>
    div#notebook-container    { width: 75%; }
    div#menubar-container     { width: 75%; }
    div#maintoolbar-container { width: 79%; }
</style>
"""))

# Commented out IPython magic to ensure Python compatibility.
# %%javascript
# IPython.OutputArea.auto_scroll_threshold = 9999

from IPython.display import Markdown, display
def printmd(string):
    display(Markdown(string))

base_dir = '/content/drive/My Drive/LabResults/'
os.chdir(base_dir)

mask_data = pd.read_csv('mask.csv')
OD_data = pd.read_csv('OD.csv')
FL_data = pd.read_csv('FL.csv')
#mask_data.head()

t_min = 0
t_max = 300
# after 300 , we see  strange bump on columns 4->12 (set t_max to 500 to see it)


# slicing frames
riboJ_exps = mask_data[(mask_data['Name'].notnull()) & (mask_data['RiboJ'] == 'Y')][['Name', 'Label', 'RiboJ','IPTG', 'Arabinose']].drop_duplicates()
#riboJ_exps.sort_values(by=['Name'], inplace=True)

non_riboJ_exps = mask_data[(mask_data['Name'].notnull()) & (mask_data['RiboJ'] == 'N')][['Name', 'Label', 'RiboJ','IPTG', 'Arabinose']].drop_duplicates()
#non_riboJ_exps.sort_values(by=['Name'], inplace=True)

control_exps = mask_data[(mask_data['Name'].notnull()) & (mask_data['RiboJ'].isnull())][['Name', 'Label', 'RiboJ','IPTG', 'Arabinose']].drop_duplicates()
control_exps.sort_values(by=['Name'], inplace=True)

# creating dicts as it will be easier to iterate over
controls = [{"Name":name,"Label":label,"IPTG":iptg,"Ara":ara}
            for name,label,iptg,ara in zip(control_exps['Name'],control_exps['Label'],control_exps['IPTG'],control_exps['Arabinose'])]
non_riboJ = [{"Name":name,"Label":label,"IPTG":iptg,"Ara":ara}
             for name,label,iptg,ara in zip(non_riboJ_exps['Name'],non_riboJ_exps['Label'],non_riboJ_exps['IPTG'],non_riboJ_exps['Arabinose'])]
riboJ = [{"Name":name,"Label":label,"IPTG":iptg,"Ara":ara}
             for name,label,iptg,ara in zip(riboJ_exps['Name'],riboJ_exps['Label'],riboJ_exps['IPTG'],riboJ_exps['Arabinose'])]


printmd("<span style='font-size:200%;font-weight:900'><br>Plotting the controls</span>")
printmd("<div style='background-color:#4682B4;width:100%;height:20px;'></div>")


for control in controls:

    printmd("<div style='background-color:#4682B4;width:100%;height:10px;'></div>")

    # from mask we retrieve the corresponding wells
    sub_frame = mask_data[mask_data['Name'] == control['Name']]
    wells = [x[0]+str(x[1]) for x in zip(sub_frame['Row'],sub_frame['Col'])]
    OD_frame = OD_data[wells]
    time = OD_data[(OD_data['Time']>=t_min) & (OD_data['Time']<=t_max)]['Time']
    all_ODs = np.array(OD_data[(OD_data['Time']>=t_min) & (OD_data['Time']<=t_max)][wells])
    all_FLs = np.array(FL_data[(FL_data['Time']>=t_min) & (FL_data['Time']<=t_max)][wells])
    geo_mean_OD = mstats.gmean(all_ODs,axis=1)
    geo_mean_FL = mstats.gmean(all_FLs,axis=1)
    fig = plt.figure(figsize=(18,9))
    ax = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)

    for well in wells:
        ax.plot(time, OD_data[(OD_data['Time']>=t_min) & (OD_data['Time']<=t_max)][well],label = well)


    ax.plot(time, geo_mean_OD,lw=8,label = 'geo mean')

    ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    ax.set_title("Optical Density", fontsize=14)
    ax.set_xlabel('time (min)')
    ax.set_ylabel('OD (Arbitrary Units)')

    for well in wells:
        ax2.plot(time, FL_data[(FL_data['Time']>=t_min) & (FL_data['Time']<=t_max)][well],label = well)

    ax2.plot(time, geo_mean_FL,lw=8,label = 'geo mean')
    ax2.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    ax2.set_title("GFP Fluorescence", fontsize=14)
    ax2.set_xlabel('time (min)')
    ax2.set_ylabel('FL (Arbitrary Units)')
    plt.subplots_adjust(left=None, bottom=None, right=1.25, top=None,
                wspace=None, hspace=None)
    printmd("<span style='font-size:180%;font-weight:900;color:#4169E1'><br>{} ({})</span>".format(control['Name'],control['Label']))
    plt.show()
    printmd("<div style='background-color:#4682B4;width:100%;height:10px;'></div>")

printmd("<span style='font-size:200%;font-weight:900'><br>Plotting the non-RiboJ Group</span>")
printmd("<div style='background-color:#4682B4;width:100%;height:20px;'></div>")


for exp in non_riboJ:
    printmd("<div style='background-color:#4682B4;width:100%;height:10px;'></div>")

    # from mask we retrieve the corresponding wells
    sub_frame = mask_data[mask_data['Name'] == exp['Name']]
    wells = [x[0]+str(x[1]) for x in zip(sub_frame['Row'],sub_frame['Col'])]
    #print(wells)
    OD_frame = OD_data[wells]
    time = OD_data[(OD_data['Time']>=t_min) & (OD_data['Time']<=t_max)]['Time']

    fig = plt.figure(figsize=(18,9))
    ax = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)

    for well in wells:
        ax.plot(time, OD_data[(OD_data['Time']>=t_min) & (OD_data['Time']<=t_max)][well],label = well)
    ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    ax.set_title("Optical Density", fontsize=14)
    ax.set_xlabel('time (min)')
    ax.set_ylabel('OD (Arbitrary Units)')

    for well in wells:
        ax2.plot(time, FL_data[(OD_data['Time']>=t_min) & (OD_data['Time']<=t_max)][well],label = well)
    ax2.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    ax2.set_title("GFP Fluorescence", fontsize=14)
    ax2.set_xlabel('time (min)')
    ax2.set_ylabel('FL (Arbitrary Units)')

    printmd("<span style='font-size:180%;font-weight:900;color:#4169E1'><br>{}: no RiboJ // IPTG: {} mM // Ara : {} mM</span>"
                    .format(exp['Name'],exp['IPTG'],exp['Ara']))

    plt.show()
    printmd("<div style='background-color:#4682B4;width:100%;height:10px;'></div>")


printmd("<span style='font-size:200%;font-weight:900'><br>Plotting the RiboJ Group</span>")
printmd("<div style='background-color:#4682B4;width:100%;height:20px;'></div>")


for exp in riboJ:
    printmd("<div style='background-color:#4682B4;width:100%;height:10px;'></div>")

    # from mask we retrieve the corresponding wells
    sub_frame = mask_data[mask_data['Name'] == exp['Name']]
    wells = [x[0]+str(x[1]) for x in zip(sub_frame['Row'],sub_frame['Col'])]
    time = OD_data[(OD_data['Time']>=t_min) & (OD_data['Time']<=t_max)]['Time']

    fig = plt.figure(figsize=(18,9))
    ax = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)

    for well in wells:
        ax.plot(time, OD_data[(OD_data['Time']>=t_min) & (OD_data['Time']<=t_max)][well],label = well)
    ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    ax.set_title("Optical Density", fontsize=14)
    ax.set_xlabel('time (min)')
    ax.set_ylabel('OD (Arbitrary Units)')

    for well in wells:
        ax2.plot(time, FL_data[(OD_data['Time']>=t_min) & (OD_data['Time']<=t_max)][well],label = well)
    ax2.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    ax2.set_title("GFP Fluorescence", fontsize=14)
    ax2.set_xlabel('time (min)')
    ax2.set_ylabel('FL (Arbitrary Units)')

    printmd("<span style='font-size:180%;font-weight:900;color:#4169E1'><br>{}: with RiboJ // IPTG: {} mM // Ara : {} mM</span>"
                    .format(exp['Name'],exp['IPTG'],exp['Ara']))

    plt.show()
    printmd("<div style='background-color:#4682B4;width:100%;height:10px;'></div>")

# We are going to use the medium for OD and FL
# This is far from ideal, but we have too few reliable controls for the FL
sub_frame = mask_data[mask_data['Name'] == 'MED']
wells = [x[0]+str(x[1]) for x in zip(sub_frame['Row'],sub_frame['Col'])]
all_ODs = np.array(OD_data[wells])
all_FLs = np.array(FL_data[wells])
geo_mean_OD = mstats.gmean(all_ODs,axis=1)
geo_mean_FL = mstats.gmean(all_FLs,axis=1)

# Now we can apply the corrections ie subtract these estimates to the previous data (except to the time column)

corrected_OD_data = OD_data.copy(deep=True)
corrected_FL_data = FL_data.copy(deep=True)
well_list = [col for col in OD_data.columns if col != 'Time']

for col in well_list:
    corrected_OD_data[col] = corrected_OD_data[col] - geo_mean_OD
    corrected_FL_data[col] = corrected_FL_data[col] - geo_mean_FL

printmd("<span style='font-size:200%;font-weight:900'><br>Plotting the (corrected) non-RiboJ Group</span>")
printmd("<div style='background-color:#4682B4;width:100%;height:20px;'></div>")


for exp in non_riboJ:
    printmd("<div style='background-color:#4682B4;width:100%;height:10px;'></div>")

    # from mask we retrieve the corresponding wells
    sub_frame = mask_data[mask_data['Name'] == exp['Name']]
    wells = [x[0]+str(x[1]) for x in zip(sub_frame['Row'],sub_frame['Col'])]
    time = corrected_OD_data['Time']

    fig = plt.figure(figsize=(18,9))
    ax = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)

    for well in wells:
        ax.plot(time, corrected_OD_data[well],label = well)
    ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    ax.set_title("Optical Density", fontsize=14)
    ax.set_xlabel('time (min)')
    ax.set_ylabel('OD (Arbitrary Units)')

    for well in wells:
        ax2.plot(time, corrected_FL_data[well],label = well)
    ax2.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    ax2.set_title("GFP Fluorescence", fontsize=14)
    ax2.set_xlabel('time (min)')
    ax2.set_ylabel('FL (Arbitrary Units)')

    printmd("<span style='font-size:180%;font-weight:900;color:#4169E1'><br>{}: no RiboJ // IPTG: {} mM // Ara : {} mM</span>"
                    .format(exp['Name'],exp['IPTG'],exp['Ara']))

    plt.show()
    printmd("<div style='background-color:#4682B4;width:100%;height:10px;'></div>")
printmd("<div style='background-color:#4682B4;width:100%;height:20px;'></div>")

printmd("<span style='font-size:200%;font-weight:900'><br>Plotting the (corrected) RiboJ Group</span>")
printmd("<div style='background-color:#4682B4;width:100%;height:20px;'></div>")


for exp in riboJ:
    printmd("<div style='background-color:#4682B4;width:100%;height:10px;'></div>")

    # from mask we retrieve the corresponding wells
    sub_frame = mask_data[mask_data['Name'] == exp['Name']]
    wells = [x[0]+str(x[1]) for x in zip(sub_frame['Row'],sub_frame['Col'])]
    time = corrected_OD_data['Time']

    fig = plt.figure(figsize=(18,9))
    ax = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)

    for well in wells:
        ax.plot(time, corrected_OD_data[well],label = well)
    ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    ax.set_title("Optical Density", fontsize=14)
    ax.set_xlabel('time (min)')
    ax.set_ylabel('OD (Arbitrary Units)')

    for well in wells:
        ax2.plot(time, corrected_FL_data[well],label = well)
    ax2.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    ax2.set_title("GFP Fluorescence", fontsize=14)
    ax2.set_xlabel('time (min)')
    ax2.set_ylabel('FL (Arbitrary Units)')

    printmd("<span style='font-size:180%;font-weight:900;color:#4169E1'><br>{}: no RiboJ // IPTG: {} mM // Ara : {} mM</span>"
                    .format(exp['Name'],exp['IPTG'],exp['Ara']))

    plt.show()
    printmd("<div style='background-color:#4682B4;width:100%;height:10px;'></div>")
printmd("<div style='background-color:#4682B4;width:100%;height:20px;'></div>")

# So it looks like some time intervals are no good
# 0-> 120 and 300+ for instance
# So let's plot the ratio only for the values in between

display_t_min = 120
display_t_max = 300

estimation_t_min = 150
estimation_t_max = 300

printmd("<span style='font-size:200%;font-weight:900'><br>Analysing the non-RiboJ Group</span>")
printmd("<div style='background-color:#4682B4;width:100%;height:20px;'></div>")


ss_data = np.zeros(12,)

i=0

for exp in non_riboJ:
    printmd("<div style='background-color:#4682B4;width:100%;height:10px;'></div>")

    estimates = []
    sub_frame = mask_data[mask_data['Name'] == exp['Name']]
    wells = [x[0]+str(x[1]) for x in zip(sub_frame['Row'],sub_frame['Col'])]

    time = corrected_OD_data[(corrected_OD_data['Time']>=display_t_min) & (corrected_OD_data['Time']<=display_t_max)]['Time']

    fig = plt.figure(figsize=(9,4), constrained_layout=True)
    ax = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)

    for well in wells:
        ax.plot(time, corrected_OD_data[(corrected_OD_data['Time']>=display_t_min) & (corrected_OD_data['Time']<=display_t_max)][well],label = well)
    ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    ax.set_title("Optical Density\n[IPTG]: {}mM  [ARAB]: {}mM".format(exp['IPTG'], exp['Ara']), fontsize=14)
    ax.set_xlabel('time (min)')
    ax.set_ylabel('OD (Arbitrary Units)')

    for well in wells:
        od = corrected_OD_data[(OD_data['Time']>=display_t_min) & (OD_data['Time']<=display_t_max)][well]
        fl = corrected_FL_data[(OD_data['Time']>=display_t_min) & (OD_data['Time']<=display_t_max)][well]
        short_od = corrected_OD_data[(OD_data['Time']>=estimation_t_min) & (OD_data['Time']<=estimation_t_max)][well]
        short_fl = corrected_FL_data[(OD_data['Time']>=estimation_t_min) & (OD_data['Time']<=estimation_t_max)][well]
        ss = mstats.gmean(short_fl/short_od,0)
        estimates.append(ss)
        ax2.plot(time, fl/od, label = well)
    ax2.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    ax2.set_title("Normalised GFP Fluorescence\n[IPTG]: {}mM  [ARAB]: {}mM".format(exp['IPTG'], exp['Ara']), fontsize=14)
    ax2.set_xlabel('time (min)')
    ax2.set_ylabel('FL (Arbitrary Units)')

    printmd("<span style='font-size:180%;font-weight:900;color:#4169E1'><br>{}: with RiboJ // IPTG: {} mM // Ara : {} mM</span>"
                    .format(exp['Name'],exp['IPTG'],exp['Ara']))

    plt.tight_layout()
    plt.show()

    # Final Results
    print('Estimated Steady States are: ')
    for ss,well in zip(estimates,wells):
        print('Well {}: {}'.format(well,str(ss)))
    print('\nCorresponding for experiment {} to a geometric mean of: {}'.format(exp['Name'],mstats.gmean(estimates,0)))
    printmd("<div style='background-color:#4682B4;width:100%;height:10px;'></div>")


    ss_data[i] = mstats.gmean(estimates)
    i+=1

ss_data

non_riboJ[0]['IPTG']

IPTGmeshgrid = [[0, 0.063, 1],[0, 0.063, 1],[0, 0.063, 1],[0, 0.063, 1]]
ARABmeshgrid = [[0,0,0],[0.0013,0.0013,0.0013],[0.083,0.083,0.083],[5.3,5.3,5.3]]


rj3Dplot = np.zeros([4,3])
rj3Dplot[0,0] = ss_data[0]
rj3Dplot[0,1] = ss_data[1]
rj3Dplot[0,2] = ss_data[2]
rj3Dplot[1,0] = ss_data[3]
rj3Dplot[1,1] = ss_data[4]
rj3Dplot[1,2] = ss_data[5]
rj3Dplot[2,0] = ss_data[6]
rj3Dplot[2,1] = ss_data[7]
rj3Dplot[2,2] = ss_data[8]
rj3Dplot[3,0] = ss_data[9]
rj3Dplot[3,1] = ss_data[10]
rj3Dplot[3,2] = ss_data[11]

print(rj3Dplot)

np.savetxt('Imesh.csv', IPTGmeshgrid, delimiter=',')
np.savetxt('Amesh.csv', ARABmeshgrid, delimiter=',')
np.savetxt('rj3d.csv', rj3Dplot, delimiter=',')

printmd("<span style='font-size:200%;font-weight:900'><br>Analysing the RiboJ Group</span>")
printmd("<div style='background-color:#4682B4;width:100%;height:20px;'></div>")

pss_data = np.zeros(12,)

i=0

for exp in riboJ:
    printmd("<div style='background-color:#4682B4;width:100%;height:10px;'></div>")

    estimates = []
    sub_frame = mask_data[mask_data['Name'] == exp['Name']]
    wells = [x[0]+str(x[1]) for x in zip(sub_frame['Row'],sub_frame['Col'])]

    time = corrected_OD_data[(corrected_OD_data['Time']>=display_t_min) & (corrected_OD_data['Time']<=display_t_max)]['Time']

    fig = plt.figure(figsize=(9,4))
    ax = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)

    for well in wells:
        ax.plot(time, corrected_OD_data[(corrected_OD_data['Time']>=display_t_min) & (corrected_OD_data['Time']<=display_t_max)][well],label = well)
    ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    ax.set_title("Optical Density\n[IPTG]: {}mM  [ARAB]: {}mM".format(exp['IPTG'], exp['Ara']), fontsize=14)
    ax.set_xlabel('time (min)')
    ax.set_ylabel('OD (Arbitrary Units)')

    for well in wells:
        od = corrected_OD_data[(OD_data['Time']>=display_t_min) & (OD_data['Time']<=display_t_max)][well]
        fl = corrected_FL_data[(OD_data['Time']>=display_t_min) & (OD_data['Time']<=display_t_max)][well]
        short_od = corrected_OD_data[(OD_data['Time']>=estimation_t_min) & (OD_data['Time']<=estimation_t_max)][well]
        short_fl = corrected_FL_data[(OD_data['Time']>=estimation_t_min) & (OD_data['Time']<=estimation_t_max)][well]
        ss = mstats.gmean(short_fl/short_od,0)
        estimates.append(ss)
        ax2.plot(time, fl/od, label = well)
    ax2.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    ax2.set_title("Normalised GFP Fluorescence\n[IPTG]: {}mM  [ARAB]: {}mM".format(exp['IPTG'], exp['Ara']), fontsize=14)
    ax2.set_xlabel('time (min)')
    ax2.set_ylabel('FL (Arbitrary Units)')

    printmd("<span style='font-size:180%;font-weight:900;color:#4169E1'><br>{}: with RiboJ // IPTG: {} mM // Ara : {} mM</span>"
                    .format(exp['Name'],exp['IPTG'],exp['Ara']))

    plt.tight_layout()
    plt.show()

    print(estimates)
    print(ss)

    # Final Results
    print('Estimated Steady States are: ')
    for ss,well in zip(estimates,wells):
        print('Well {}: {}'.format(well,str(ss)))
    print('\nCorresponding for experiment {} to a geometric mean of: {}'.format(exp['Name'],mstats.gmean(estimates,0)))
    printmd("<div style='background-color:#4682B4;width:100%;height:10px;'></div>")

    pss_data[i] = mstats.gmean(estimates)
    i+=1

prj3Dplot = np.zeros([4,3])
prj3Dplot[0,0] = pss_data[0]
prj3Dplot[0,1] = pss_data[1]
prj3Dplot[0,2] = pss_data[2]
prj3Dplot[1,0] = pss_data[3]
prj3Dplot[1,1] = pss_data[4]
prj3Dplot[1,2] = pss_data[5]
prj3Dplot[2,0] = pss_data[6]
prj3Dplot[2,1] = pss_data[7]
prj3Dplot[2,2] = pss_data[8]
prj3Dplot[3,0] = pss_data[9]
prj3Dplot[3,1] = pss_data[10]
prj3Dplot[3,2] = pss_data[11]

print(prj3Dplot)

np.savetxt('prj3d.csv', prj3Dplot, delimiter=',')

fig = plt.figure(figsize=(25,8), facecolor='w')
ax = fig.add_subplot(121, projection='3d')
surf = ax.plot_surface(IPTGmeshgrid, ARABmeshgrid, rj3Dplot, cmap=cm.tab20c)
ax.set_xlabel('[IPTG] (mM)')
ax.set_xlim(1, 0)
ax.set_ylabel('[Arabinose] (mM)')
ax.set_zlabel('Normalised Fluorescence (AU)')
ax.set_title('AND Gate without RiboJ')
ax.set_facecolor('w')
fig.colorbar(surf)

ax = fig.add_subplot(122, projection='3d')
surf = ax.plot_surface(IPTGmeshgrid, ARABmeshgrid, prj3Dplot, cmap=cm.tab20c)
ax.set_xlabel('[IPTG] (mM)')
ax.set_xlim(1, 0)
ax.set_ylabel('[Arabinose] (mM)')
ax.set_zlabel('Normalised Fluoorescence (AU)')
ax.set_title('AND Gate with RiboJ')
ax.set_facecolor('w')
fig.colorbar(surf)
plt.show()

print(rj3Dplot)
print(prj3Dplot)
