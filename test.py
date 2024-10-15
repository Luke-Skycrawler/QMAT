from pyqmat import qmat
import numpy as np
a = qmat("build/bug.off", "build/bug.ma")
a.simplify_slab(35)
a.export_ply("output/bug")
a.export_ma("output/bug")
h = a.export_hausdorff_distance()
print(h)