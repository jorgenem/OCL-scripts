from pyma_v1 import *

a = np.array([1,5,8,3,5,7,9,2,10])
axis = np.linspace(-2,len(a)-2,len(a)+1)
a2d = np.repeat(a,len(a)).reshape((len(a),len(a)))

plt.figure()

plt.subplot(2,2,1)
plt.title('original')
plt.pcolormesh(axis,axis,a2d)
print "original:"
print a2d

plt.subplot(2,2,2)
plt.title('EitoEf')
plt.pcolormesh(axis,axis,EitoEf(a2d,axis))
print "EitoEf:"
print EitoEf(a2d,axis)

plt.subplot(2,2,3)
plt.title('EftoEi(EitoEf)')
plt.pcolormesh(axis,axis,EftoEi(EitoEf(a2d,axis),axis))
print "EftoEi(EitoEf):"
print EftoEi(EitoEf(a2d,axis),axis)

plt.subplot(2,2,4)
plt.title('EitoEi(EitoEf(EftoEi))')
plt.pcolormesh(axis,axis,EitoEf(EftoEi(EitoEf(a2d,axis),axis),axis))
print "EitoEi(EitoEf(EftoEi)):"
print EftoEi(EitoEf(a2d,axis),axis)

plt.show()

