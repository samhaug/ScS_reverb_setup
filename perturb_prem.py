#!/home/samhaug/anaconda2/bin/python

import numpy as np

def write_output(prem_array,title):
    with open(title+'.tvel','a') as f:
        f.write('prem_5km.tvel P\n')
        f.write('prem_5km.tvel S\n')
        np.savetxt(f,prem_array)


perturb = np.round(np.linspace(-0.03,0.03,num=12),2)
print len(perturb)

for ii in perturb:
    prem = np.genfromtxt('./models/prem.tvel',skip_header=2)
    prem[0:40,:] *= 1+ii
    write_output(prem,'prem_'+str(ii))




