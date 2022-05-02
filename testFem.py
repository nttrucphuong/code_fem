from fem import fem_calculate 
from fpdf import FPDF
from tabulate import tabulate
import os


def writeFile(x,y, stress=True):
    if stress:
        so_chia=9
    else:
        so_chia=11
    K1_stress, K2_stress, K3_stress, K_tt, Se, = fem_calculate(x,y, stress)
    K1_stress=tabulate(K1_stress)
    K2_stress=tabulate(K2_stress)
    K3_stress=tabulate(K3_stress)
    K_tt=tabulate(K_tt)
    Se=tabulate(Se)
    f.write(f'X={x} Y ={y}')
    f.write(f'\n \n')
    f.write(f'K1= (10^{so_chia}) *  \n')
    f.write(f'\n')
    f.write(K1_stress)
    f.write(f'\n \n')
    f.write(f'K2= (10^{so_chia}) *   \n')
    f.write(K2_stress)
    f.write(f'\n \n')
    f.write(f'K3= (10^{so_chia}) *   \n')
    f.write(K3_stress)
    f.write(f'\n \n')
    f.write(f'K_tt= (10^{so_chia}) *   \n')
    f.write(K_tt)
    f.write(f'\n \n')
    f.write(f'Se= (10^10) *  \n')
    f.write(Se)
    f.write(f'\n \n')



X=[1,2,3,4,5]
Y=[1,2,3,4,5,6,7,8,9,10]


file_name='plane_stress.txt'
if os.path.exists(file_name):
    os.remove(file_name)
global f
f = open(file_name, "a")
for x in X:
    for y in Y:
        writeFile(x,y, stress=True)
        f.write('\n \n \n')
        f.write('============================================ \n')
        

f.close()

file_name='plane_strain.txt'
if os.path.exists(file_name):
    os.remove(file_name)
f = open(file_name, "a")
for x in X:
    for y in Y:
        writeFile(x,y, stress=False)
        f.write('\n \n \n')
        f.write('============================================ \n')
        

f.close()
