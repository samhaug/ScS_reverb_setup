#!/home/samhaug/anaconda2/bin/python

import numpy as np

def write_output(prem_array,title):
    with open(title+'.tvel','a') as f:
        f.write('prem_5km.tvel P\n')
        f.write('prem_5km.tvel S\n')
        np.savetxt(f,prem_array,fmt='%.3f')


perturb = np.round(np.linspace(-0.3,0.3,num=20),3)
print len(perturb)

for ii in perturb:
    prem = np.genfromtxt('./models/prem.tvel',skip_header=2)
    prem[2:40,1::] *= 1+ii
    write_output(prem,'prem_'+str(ii))




