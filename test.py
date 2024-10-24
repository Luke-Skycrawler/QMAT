from pyqmat import qmat, PointAdder
import numpy as np


record_start = 100
record_end = 20 
record_step = 10

model = "hand"

off_name = "data/spider.off"
complete_ma = "data/spider.ma"
fine_ma = "data/spider_v100.ma"
coarse_ma = "data/spider_v25.ma"

def test_reentrant():
    a = qmat(off_name, complete_ma)    
    a.simplify_slab(100)
    a.export_ply("output/spider_to100")

    a.init_mergelist()
    a.init_collapse_queue()
    a.simplify(75)

    a.export_ply("output/spider_100to25")
    M = a.export_mergelist()
    # print(M)

    return M

def test_add_sphere():
    a = qmat(off_name, fine_ma)
    # a.simplify_slab(100)
    a.adjust_storage()

    b = qmat(a, coarse_ma)
    sphere = np.array([-0.2176294, 0.06370435, 0.09288197, 0.04448237])

    diag = a.get_diagonal()
    adder = PointAdder(diag, b, a)
    L = test_reentrant()
    adder.set_mergelist(L)
    adder.add_new_node(*sphere)
    adder.export_ply("output/spider_v26")

if __name__ == "__main__":

    # test_reentrant()
    test_add_sphere()
