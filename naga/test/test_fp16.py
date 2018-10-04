import naga
import numpy as np

c = naga.naga_context.get_instance()

m = 5000
n = 24000

A = np.random.random((m, n))
x = np.random.random((n, 1))
y = np.zeros((m, 1))

mat = naga.make_naga_obj_half(c, A)
matF = naga.naga_obj_float(c, A)
X = naga.make_naga_obj_half(c, x)
Y = naga.make_naga_obj_half(c, y)
XF = naga.naga_obj_float(c, x)
YF = naga.naga_obj_float(c, y)
