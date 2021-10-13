import os
from scipy import stats
import matplotlib.pyplot as plt
import numpy as np
import csv
#--------------------------------------------------------------------------

# IN TE VULLEN
filename = "02072021.csv"    # file name met data
time = 36                   # time in minutes between subsequent measurements
theo_time = 6.0067          # theoretical half life of isotope in hours, Tc99m = 6.0067 h

#--------------------------------------------------------------------------

# Initialisation
rows = []
nmeas = []
date = []
activity_all = []
if not os.path.exists('Verwerkte QCs'):
  os.mkdir('Verwerkte QCs')


#--------------------------------------------------------------------------

# reading csv file
with open(filename, 'r') as csvfile:
        csvreader = csv.reader(csvfile,delimiter = ';')

        # extracting each data row one by one
        for row in csvreader:
                rows.append(row)


nval = '1'
rowvalue = 0
for row in rows:
 if not row == []:
   if row[0] == nval:
        nmeas.append(int(nval)*time/60)
        for value in row:
                if len(value) == 15 and value[1] == '/':
                        date.append(value)
                if len(value) > 6 and value[-3:] == 'GBq':
                        nactivity = float(value[:5].replace(',', '.'))*1000
                        activity_all.append(nactivity)
                if len(value) > 6 and value[-3:] == 'MBq':
                        nactivity = float(value[:5].replace(',', '.'))
                        activity_all.append(nactivity)
        nval = str(int(nval)+1)
   rowvalue+=1

# splitting activity in measured and expected
activity_meas  = []
activity_theor = []

for i in range(len(activity_all)):
        if (i % 2 == 0):
              activity_meas.append(activity_all[i])
        else:
              activity_theor.append(activity_all[i])

#--------------------------------------------------------------------------
# Making plot

fig1, ax1  = plt.subplots(1, 1, sharex=True,sharey=True)
x1 = np.linspace(0, nmeas[-1], len(nmeas))
slope, intercept, r_value, p_value, std_err = stats.linregress(x1, np.log(activity_meas))
slope2, intercept2, r_value2, p_value2, std_err2 = stats.linregress(x1, np.log(activity_theor))
pract_time = np.log(2)/-slope

rsqrt = r_value ** 2

# Voor de berekening van de afwijking per punt
delta_x = []
delta_xx = []
x_reflist = []
for i in range(len(x1)):
    x_ref = np.exp(slope * x1[i] + intercept)
    x_reflist.append(x_ref)
    delta_x.append(np.abs((x_ref-activity_meas[i])/(x_ref))*100)
    delta_xx.append(np.abs((activity_theor[i]-activity_meas[i])/(activity_theor[i]))*100)
    max_x = np.argmax(delta_x)

ax1.semilogy(nmeas, activity_meas, 'purple', label=f'Gemeten activiteit met t1/2 = {"{:0.3f}".format(pract_time)} h'
                '\n' r'$\rightarrow$' f' afwijking van {"{:0.2f}".format(np.abs(theo_time-pract_time)/theo_time*100)} % tov {theo_time} h')
#ax1.semilogy(nmeas, activity_theor, 'b.', label='Theoretical activity')
#ax1.semilogy(nmeas, np.exp(slope2*x1+intercept2), 'r--', label='Theoretical activity')
ax1.plot(nmeas, np.exp(slope*x1+intercept),color='orange',linestyle='dashed',label=f'Lineare regressie (op semilog) met r² = {"{:0.3f}".format(rsqrt)}')
ax1.scatter(nmeas[max_x], np.exp(slope*x1[max_x]+intercept), color='r', label=f'Maximal deviation of {"{:0.3f}".format(np.max(delta_x))} %'
            f'\n Extra: Avg. deviation {"{:0.3f}".format(np.average(delta_x))} %')

ax1.set_title('QC Dosis Calibrator')
ax1.set_ylabel('Activiteit [MBq]')
ax1.set_xlabel('Tijd [h]')
ax1.legend()

fig1.savefig(os.path.join('Verwerkte QCs',f'{filename}.pdf'))
