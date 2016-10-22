from pyma_v1 import *

a = np.array([1,2,3,4,5,6,7,8,9])
b = np.array([2,3,5,6,3,8,6,5,4])
axis = np.linspace(-2,len(a)-2,len(a)+1)
# a2d = np.repeat(a,len(a)).reshape((len(a),len(a)))
a2d, b2d = np.meshgrid(a,b,indexing='ij')
a2d = a2d*b2d

plt.figure()

plt.subplot(3,2,1)
plt.title('original')
plt.pcolormesh(axis,axis,a2d,vmin=a2d.min(), vmax=a2d.max())
print "original:"
print a2d

plt.subplot(3,2,2)
plt.title('EitoEf')
plt.pcolormesh(axis,axis,EitoEf(a2d,axis),vmin=a2d.min(), vmax=a2d.max())
print "EitoEf:"
print EitoEf(a2d,axis)

plt.subplot(3,2,3)
plt.title('EftoEi(EitoEf)')
plt.pcolormesh(axis,axis,EftoEi(EitoEf(a2d,axis),axis),vmin=a2d.min(), vmax=a2d.max())
print "EftoEi(EitoEf):"
print EftoEi(EitoEf(a2d,axis),axis)

plt.subplot(3,2,4)
plt.title('EitoEi(EitoEf(EftoEi))')
plt.pcolormesh(axis,axis,EitoEf(EftoEi(EitoEf(a2d,axis),axis),axis),vmin=a2d.min(), vmax=a2d.max())
print "EitoEi(EitoEf(EftoEi)):"
print EitoEf(EftoEi(EitoEf(a2d,axis),axis),axis)

plt.subplot(3,2,5)
plt.title('..and again..')
plt.pcolormesh(axis,axis,EftoEi(EitoEf(EftoEi(EitoEf(a2d,axis),axis),axis),axis),vmin=a2d.min(), vmax=a2d.max())

plt.subplot(3,2,6)
plt.title('..and again..')
plt.pcolormesh(axis,axis,EitoEf(EftoEi(EitoEf(EftoEi(EitoEf(a2d,axis),axis),axis),axis),axis),vmin=a2d.min(), vmax=a2d.max())


plt.show()

