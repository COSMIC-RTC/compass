import numpy as np
import chakra as ch



a=np.ones((6),dtype=np.float64)*np.pi
b=np.ones((2,3),dtype=np.float32)*np.pi
c=ch.chakra_context()


#reshape mandatory: data must be a np.ndarray(ndim=3)
print "constructor double, data"
sh=[a.shape[0],1,1]
oa=ch.chakra_obj_Double(c,data=np.reshape(a,sh),ndim=1)

#reshape mandatory: data must be a np.ndarray(ndim=3)
print "constructor float, data"
sh=[b.shape[0],b.shape[1],1]
ob=ch.chakra_obj_Float(c,data=np.reshape(b,sh),ndim=2)

print "constructor copy double"
o_copy=ch.chakra_obj_Double(obj=oa)
print "constructor copy float"
o_copy2=ch.chakra_obj_Float(ctxt=c,obj=ob)


b_out=ob.device2host(2)
cp_out=o_copy.device2host(1)
cp2_out=o_copy2.device2host(2)

##reshape mandatory
a2=np.ones((6),dtype=np.float64)*(np.pi+1)
sh=o_copy.get_Dims()
sh=np.append(sh[1:],[1,1])
o_copy.host2device(np.reshape(a2,sh[0:3]))

a2=np.ones((6),dtype=np.float32)*(np.pi+1)
sh=o_copy2.get_Dims()
sh=np.append(sh[1:],[1,1])
o_copy2.host2device(np.reshape(a2,sh[0:3]))


a_out=oa.device2host(1)
print "oa:"
print a_out
o_copy.copyInto(oa)
a_out=oa.device2host(1)
print " copy into oa:"
print a_out

b_out=ob.device2host(2)
print "ob:"
print b_out
ob.copyFrom(o_copy2)
b_out=oa.device2host(1)
print " copy from cp2:"
print b_out

