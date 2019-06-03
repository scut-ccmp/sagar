import math
import os
import numpy

from sagar.crystal.structure import Cell
from sagar.io.vasp import  *
from sagar.crystal.utils import non_dup_hnfs, _is_hnf_dup, _hnfs
from itertools import combinations


def match(pcell_A, pcell_B, n_max, distance, err_l, err_s, err_theta, n_min=1, MS=0, translation_vectors=numpy.array([[0, 0]])):
    s_A =  abs(numpy.linalg.det(numpy.array(pcell_A._lattice[0:2, 0:2])))
    s_B =  abs(numpy.linalg.det(numpy.array(pcell_B._lattice[0:2, 0:2])))
    NM = match_S(s_A, s_B, err_s, n_min, n_max)######得到匹配的面积比NM
    print(NM)
    #######去除非Z轴的对称性：
    pcell_A_X = add_atom_bottom(pcell_A)
    if MS==0:
        pcell_B_X = pcell_B
    else:
        pcell_B_X = add_atom_bottom(pcell_B)
    ###存储匹配对应的AB的胞矩阵：
    A_b = []
    B_a = []
    err = [] ##err_s err_l1 err-l2 err_theta
    for nm in NM:
        ######得到所有固定面积超胞:
        slatt_A = non_dup_hnfs(pcell_A_X, nm[0], symprec=1e-2, comprec=1e-2, dimension=2)
        slatt_B = non_dup_hnfs(pcell_B_X, nm[1], symprec=1e-2, comprec=1e-2,  dimension=2)
        ######超胞变换为唯一形状:
        slatt_A_d = []
        slatt_B_d = []
        for l_A in slatt_A:
            lat = delaunay_reduce_2D(pcell_A._lattice, l_A)
            slatt_A_d.append(lat)
        for l_B in slatt_B:
            lat = delaunay_reduce_2D(pcell_B._lattice, l_B)
            slatt_B_d.append(lat)
        for l_A_d in slatt_A_d:
            ######比较边长夹角得到匹配的超胞
            dis_A = numpy.linalg.norm(l_A_d[0:2,0:2], axis=1)
            theta_A = (l_A_d[0,0]*l_A_d[1,0]+l_A_d[0,1]*l_A_d[1,1])/(dis_A[0]*dis_A[1])
            for l_B_d in slatt_B_d:
                dis_B = numpy.linalg.norm(l_B_d[0:2,0:2], axis=1)
                theta_B = (l_B_d[0,0]*l_B_d[1,0]+l_B_d[0,1]*l_B_d[1,1])/(dis_B[0]*dis_B[1])
                error_s = abs(nm[0]*s_A-nm[1]*s_B)/(nm[0]*s_A)
                error_l1 = abs(dis_A[0] - dis_B[0])/dis_A[0]
                error_l2 = abs(dis_A[1] - dis_B[1])/dis_A[1]
                error_theta = abs(theta_A - theta_B)
                if (error_l1 < err_l and error_l2 < err_l and error_theta < err_theta):
                    A_b.append(numpy.round(l_A_d*numpy.matrix(pcell_A._lattice).I).astype(int))
                    B_a.append(numpy.round(l_B_d*numpy.matrix(pcell_B._lattice).I).astype(int))
                    err.append([error_s, error_l1, error_l2, error_theta])
    with open('error',"w") as f:
        f.writelines('error_area\terror_l1\terror_l2\terror_theta\n')
        for ii in range(len(err)):
            f.writelines('\t'.join([str(i) for i in err[ii]])+'\n')
    # numpy.savetxt('latticematching_files/output/error',err)
    # numpy.savetxt('error',err)
    n_poscar = 0
    print(A_b, B_a)
    for i in range(0, len(A_b)):
        for translation_vector in translation_vectors:
            ######产生匹配好的POSCAR
            scell_A = Cell.extend(pcell_A, A_b[i])
            scell_B = Cell.extend(pcell_B, B_a[i])
            latt_layer = scell_A._lattice
            translation_vector = numpy.matrix(translation_vector)*\
                numpy.matrix(pcell_A._lattice[0:2,0:2])*numpy.linalg.inv(latt_layer[0:2,0:2])
            positions_B = scell_B._positions\
                +numpy.repeat([[translation_vector[0,0], translation_vector[0,1],
                    numpy.amax(scell_A._positions[:,2])-numpy.amin(scell_B._positions[:,2])\
                +distance/latt_layer[2,2]]], len(scell_B._atoms), axis=0)###A上面扣上B
            positions_layer = numpy.vstack((scell_A._positions,positions_B))
            atoms_layer = list(scell_A._atoms) + list(scell_B._atoms)
            if numpy.linalg.det(scell_A._lattice)<0:
                ###保证右手系
                latt_layer[[0,1], :] = latt_layer[[1,0], :]
                positions_layer[:, [0,1]] = positions_layer[:, [1,0]]
            cell_layer =  Cell(latt_layer, positions_layer, atoms_layer)
            # fname = 'latticematching_files/output/POSCAR'+str(n_poscar+1)
            fname = 'POSCAR'+str(n_poscar+1)
            n_poscar += 1
            write_vasp(cell_layer, filename=fname, suffix='.vasp', long_format=True)


def add_atom_bottom(pcell):
    ###通过在底部加一层和上一层位置相同的特殊原子来去除非Z轴对称性
    X_pos = numpy.concatenate((pcell._positions, numpy.array([[0, 0, 0 ]])),axis=0)
    X_pos[-1][2] = X_pos[-1][2]+0.5
    X_atoms = numpy.append(pcell._atoms,numpy.array([666]))
    return Cell(pcell._lattice, X_pos, X_atoms)


def match_S(s_A, s_B, err_s, n_min, n_max):
    ###寻找面积匹配的公倍数
    nms = []
    for n in range(n_min, n_max+1):
        for m in range(math.floor(s_A/s_B), math.floor(n_max*s_A/s_B)+2):
            if abs(n*s_A-m*s_B) < err_s*n*s_A:
                nms.append([n, m])
    NM = numpy.array(nms)
    return NM


def delaunay_reduce_2D(platt, hnfm):
    ###变换基矢到饱满,算法来自1957<<Remarks on the Delaunay Reduction>>
    latt = numpy.matrix(hnfm)*numpy.matrix(platt)
    s = abs(numpy.linalg.det(numpy.array(latt)))
    latt_2D = latt[0:2, 0:2]
    a = numpy.row_stack((latt_2D, -(latt_2D[0]+latt_2D[1])))
    com = list(combinations([0, 1, 2], 2))
    for t in range(1,1000):
        k = 0
        for ij in com:
            hij = a[ij[0]]*a[ij[1]].T
            if hij>0.00001:
                a[0] = -a[ij[0]]
                a[1] = a[ij[1]]
                a[2] = -(a[0]+a[1])
                break
            k += 1
        if k == 3:
            ###三个夹角都是非锐角时,结束
            break
    dis = numpy.linalg.norm(a[0:3,0:2], axis=1)
    a_dis = numpy.column_stack((a, dis))
    b = numpy.matrix(a_dis[numpy.lexsort(a_dis.T)])####选取最短两条
    latt[0:2, 0:2] = b[0:2, 0:2]
    if (t < 999 and abs(abs(numpy.linalg.det(numpy.array(latt)))-s) < 0.1):
        return latt
    else:
        raise ValueError("This cell can't transform to a full cell.")


def del_file(path):
    ###删除文件
    ls = os.listdir(path)
    for i in ls:
        c_path = os.path.join(path, i)
        if os.path.isdir(c_path):
            del_file(c_path)
        else:
            os.remove(c_path)

def main_match(poscar_a,poscar_b,n_max=16,distance=2.2,err_l=0.08,err_s=0.08,
            err_theta=0.08,n_min=1,MS=1,translation_vectors=numpy.array([[0, 0]])):
    pcell_A = read_vasp(poscar_a)
    pcell_B = read_vasp(poscar_b)
    # del_file('./latticematching_files/output')
    match(pcell_A,pcell_B,n_max=n_max,distance=distance,err_l=err_l,err_s=err_s,
        err_theta=err_theta,n_min=n_min,MS=MS,translation_vectors=translation_vectors)
#match(pcell_A, pcell_B, n_max=16, distance=2.2, err_l = 0.08, err_s = 0.08, err_theta = 0.08)
##distance:吸附层和衬底间距离(挨)
##MS=0：单层薄膜匹配（可翻转）MS=1：块块界面匹配（不支持翻转）
##err_theta 是余弦值的差（不是相对值，因为可能是0）
#translation_vectors : 错位向量
#[0, 0],[0.5, 0],[0, 0.5],[0.5, 0.5],[0.33333,0.33333],[0.66666,0.66666]

if __name__ == '__main__':
    main_match('/home/hecc/Desktop/input/Ag111.vasp','/home/hecc/Desktop/input/silicene.vasp')
