#!/home/samhaug/anaconda2/python

from sys import argv
import numpy as np
from scipy.signal import blackman
from scipy.interpolate import interp1d
from matplotlib import pyplot as plt
from scipy.interpolate import UnivariateSpline

'''
argv[1]: path to mineos model txt file. Must be in tabular format.

This script writes an axisem MESHER .bm file. This is to make a mesh file that
has a 1D discontinuity in it.
'''

def make_discont(rad,amp):
    '''
    rad: location of discontinuity
    amp: amplitude of discontinuity
    '''
    r = np.linspace(3480000,6371000,20)
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
    title= kwargs.get('title','{}_km_{}.bm'.format(discont,amp))
    with open(title, 'w') as f:
        f.write('ANELASTIC T\n')
        f.write('ANISOTROPIC F\n')
        f.write('UNITS m\n')
        f.write('COLUMNS radius rho vpv vsv qka qmu vph vsh eta\n')
        np.savetxt(f, out_array, delimiter="   ",fmt='%.3f')

discont = 4700000
amp = 0.000
#out,r = make_discont(discont,amp)
meta,cores,mantle = read_model(argv[1])

f = interp1d(mantle[:,0],mantle[:,1::],axis=0)
int_mant = np.hstack((np.transpose([r]),f(r)))

'''
int_mant[:,3]*=(1+out)
int_mant[:,7]*=(1+out)
'''

rad = mantle[:,0]
props = mantle[:,1::]
long_int = np.linspace(int_mant[:,0].min(),int_mant[:,0].max(),num=100)
new_mantle = np.zeros((len(long_int),mantle.shape[1]))
new_mantle[:,0] = long_int
for ii in range(1,mantle.shape[1]):
    if ii == 4:
       f = interp1d(mantle[:,0],mantle[:,4])
       new_mantle[:,4] = f(long_int)
       continue
    if ii == 5:
       f = interp1d(mantle[:,0],mantle[:,5])
       new_mantle[:,5] = f(long_int)
       continue
    if ii == 8:
       f = interp1d(mantle[:,0],mantle[:,8])
       new_mantle[:,8] = f(long_int)
       continue
    s = UnivariateSpline(int_mant[:,0],int_mant[:,ii],s=10000)
    new_mantle[:,ii] = s(long_int)
    #plt.plot(long_int,new_mantle[:,ii])
    #s = UnivariateSpline(int_mant[:,0],int_mant[:,ii],s=100)
    #plt.plot(long_int,s(long_int))
    #plt.plot(int_mant[:,0],int_mant[:,ii])

    #plt.show()


f_m = interp1d(mantle[:,0],mantle[:,1::],axis=0)
r = np.linspace(mantle[:,0].min(),mantle[:,0].max(),num=new_mantle.shape[0])
mantle = np.hstack((np.transpose([r]),f_m(r)))

for ii in range(1,new_mantle.shape[1]):
    diff = mantle[:,ii].sum()-new_mantle[:,ii].sum()
    new_mantle[:,ii] += diff/new_mantle.shape[0]

output_prem = np.vstack((cores,new_mantle))
prem = np.vstack((cores,mantle))
write_output(output_prem,0,0,title='smooth_prem.bm')





