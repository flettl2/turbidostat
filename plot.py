#plot.py
#Creates a simple OD vs time plot in Python
#Creates csv files to import into R for better visualization

#!/usr/bin/env python
#Above line will ensure python is envoked on mac/linux

import matplotlib
matplotlib.use('TkAgg')
import json
import csv #This is code I'm adding to export time, od and dilution rate data to csv files
from pylab import * #This is more of the same


#if we're running in canopy then everything is already imported AND we're in 
# interactive mode for matplotlib (meaning you don't need the show() command).
#try:
#    np
#    no_show = True~\Desktop\School\Turbidostat\Flexostat-interface\Flexostat-interface-master
#except NameError: 
#    from pylab import *
#    no_show = False
    
#Some constants
log_file = "log.dat" #assumes log is named log.dat and in the current directory.

#note: there are more advanced ways to do the stuff below, but it has been 
#      been written this way to make it easier to understand.
data = []
with open(log_file,'r') as f: #open log.dat and assign that to the variable f
    for line in f: #loop over each line in f.
        data.append(json.loads(line))

#go through each data point in the collected data and construct time and 
#OD vectors
t0 = data[0]["timestamp"] #experiment start time
t = []
u = []
od = []
for datum in data:
    #convert time to hours since start
    t.append((datum["timestamp"]-t0)/60.0/60.0) 
    u.append(datum["u"])
    od.append(datum["ods"])

#plotting code
legend_names = []
for chamber in range(1,9):
    legend_names.append("Chamber " + str(chamber))

plot(t,od), #can replace od with u for plotting dilution rate
legend(legend_names,"best")
xlabel("Time hr")
ylabel("OD$_{650}$")
    
#if not no_show:
show()
    
    
#The following is code added by Lucas FLett to export time, od and dilution rate date to csv files
#time
t_file = open("time.csv", "wb")

wr = csv.writer(t_file, dialect='excel')

for item1 in t:
    wr.writerow([item1,])
    
t_file.close()

#od
with open ("od.csv", "wb") as f:
    wtr = csv.writer(f, delimiter=',')
    wtr.writerows(od)


#dilution
with open ("dilution.csv", "wb") as f:
    wtr = csv.writer(f, delimiter=',')
    wtr.writerows(u)
    




