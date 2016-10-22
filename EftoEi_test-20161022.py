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
plt.title('EftoEi')
plt.pcolormesh(axis,axis,EftoEi(a2d,axis))
print "EftoEi:"
print EftoEi(a2d,axis)

plt.subplot(2,2,3)
plt.title('EitoEf(EftoEi)')
plt.pcolormesh(axis,axis,EitoEf(EftoEi(a2d,axis),axis))
print "EitoEf:"
print EitoEf(EftoEi(a2d,axis),axis)

plt.subplot(2,2,4)
plt.title('EftoEi(EitoEf(EftoEi))')
plt.pcolormesh(axis,axis,EftoEi(EitoEf(EftoEi(a2d,axis),axis),axis))

plt.show()

