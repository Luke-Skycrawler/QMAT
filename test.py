from pyqmat import qmat
import numpy as np


record_start = 100
record_end = 20 
record_step = 10

model = "hand"

off_name = "data/spider.off"
complete_ma = "data/spider.ma"
fine_ma = "data/spider_v100.ma"
coarse_ma = "data/spider_v20.ma"

def test_reentrant():
    a = qmat(off_name, complete_ma)    
    a.simplify_slab(100)
    a.export_ply("output/spider_to100")

    a.init_mergelist()
    a.init_collapse_queue()
    a.simplify(75)

    a.export_ply("output/spider_100to25")
    M = a.export_mergelist()
    print(M)


def test_nqmat():
    a = nqmat(f"data/{model}.off", f"data/{model}.ma")

    frame = 0
    for i in range(record_start, record_end, -record_step): 
        # a.simplify(i)
        if frame == 0:
            a.simplify(i)
        else :
            a.delete_n(record_step)
        # a.clean_up()
        a.export_ply(f"output/{model}_{frame}")
        frame += 1 

def test_qmat():
    a = qmat(f"data/{model}.off", f"data/{model}.ma")
    frame = 0
    for i in range(record_start, record_end, -record_step): 
        a.simplify_slab(i)
        a.export_ply(f"output/{model}_{frame}")
        frame += 1

def test_qmat_simple():
    a = qmat(f"data/{model}.off", f"data/{model}.ma")
    a.simplify_slab(35)
    a.export_ply(f"output/{model}_0")

if __name__ == "__main__":

    test_reentrant()
    # test_nqmat()

    # test_qmat()
    # test_qmat_simple()
# a.simplify_slab(35)
# a.clean_up()
# a.export_ply(f"output/{model}")
# a.export_ma("output/bug")
# h = a.hausdorff()
# print(h)