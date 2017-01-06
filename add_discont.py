#!/home/samhaug/anaconda2/python

from sys import argv
import numpy as np
from scipy.signal import blackman
from scipy.interpolate import interp1d
from matplotlib import pyplot as plt

'''
argv[1]: path to mineos model txt file. Must be in tabular format.
'''

def make_discont(rad,amp):
    '''
    rad: location of discontinuity
    amp: amplitude of discontinuity
    '''
    r = np.linspace(3480000,6371000,250)
    kink = blackman(15)*amp
    kink[0:7]*=-1.
    #kink +=1
    kink_ind = np.argmin(np.abs(r-rad))
    out = np.hstack((np.zeros(kink_ind-7),kink,np.zeros(len(r)-kink_ind-len(kink)/2-1)))
    return out,r

def read_model(filename):
    a = open(filename)
    f = a.read()
    meta = f.strip().split('\n')[2].split()
    knots,nic,noc = meta[0],meta[1],meta[2]
    print noc
    a.close()
    model = np.genfromtxt(filename,skip_header=3)
    cores = model[0:int(noc),:]
    mantle = model[int(noc)::,:]
    return meta,cores,mantle

def write_output(out_array,discont,amp,**kwargs):
    title= kwargs.get('title','test.txt')
    with open(title, 'w') as f:
        f.write('Radius: {} m, Amplitude: {} %vs\n'.format(str(discont),
                str(amp)))
        f.write('1    1.0000   1\n')
        f.write('{}   {}   {}\n'.format(str(out_array.shape[0]),33,66))
        np.savetxt(f, out_array, delimiter="   ",fmt='%.1f')

out,r = make_discont(4700000,0.01)
meta,cores,mantle = read_model(argv[1])

f = interp1d(mantle[:,0],mantle[:,1::],axis=0)
int_mant = np.hstack((np.transpose([r]),f(r)))
int_mant[:,3]*=(1+out)
int_mant[:,7]*=(1+out)
plt.plot(int_mant[:,0],int_mant[:,3])
plt.show()

output_prem = np.vstack((cores,int_mant))
#output_prem = np.around(output_prem,3)

write_output(output_prem,2700000,0.01)




