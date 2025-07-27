

""" PACOTES """
import numpy as np

""" INPUTS """
N1 = int(180*1) #numero total de atomos de c1
N2 = int(180*1) #numero total de atomos de c1b
N3 = int(180*1) #numero total de atomos de c2
N4 = int(180*1) #numero total de atomos de c3
N5 = int(180*1) #numero total de atomos de c2b
N6 = int(180*1) #numero total de atomos de c3b
N7 = int(180*1) #numero total de atomos de n1
N8 = int(180*1) #numero total de atomos de n2
N9 = int(180*1) #numero total de atomos de n1b
N10 = int(180*1) #numero total de atomos de n2b
Nt =  int(180*10) #numero total de atomos
nframes = 1

""" INICIANDO MATRIZES E VETORES """
rv1 = np.zeros((nframes, N1, 3), dtype=np.float_)
rv2 = np.zeros((nframes, N2, 3), dtype=np.float_)
rv3 = np.zeros((nframes, N3, 3), dtype=np.float_)
rv4 = np.zeros((nframes, N4, 3), dtype=np.float_)
rv5 = np.zeros((nframes, N5, 3), dtype=np.float_)
rv6 = np.zeros((nframes, N6, 3), dtype=np.float_)
rv7 = np.zeros((nframes, N7, 3), dtype=np.float_)
rv8 = np.zeros((nframes, N8, 3), dtype=np.float_)
rv9 = np.zeros((nframes, N9, 3), dtype=np.float_)
rv10 = np.zeros((nframes, N10, 3), dtype=np.float_)

""" LEITURA DO CENTRO DE MASSA DAS MOLECULAS """
f = open("TCNE.pdb", "r")

for frame in range(1):
    i1 = 0
    i2 = 0
    i3 = 0
    i4 = 0
    i5 = 0
    i6 = 0
    i7 = 0
    i8 = 0
    i9 = 0
    i10 = 0
    for n in range(1):
        f.readline()
    for n in range(Nt):
        a1, a2, ida, a4, a5, a6, rx, ry, rz, b1, b2 = f.readline().split() 
        if (int(ida) == 1): 
            rv1[frame, i1, :] = float(rx), float(ry), float(rz)
            i1 += 1
        elif (int(ida) == 2): 
            rv2[frame, i2, :] = float(rx), float(ry), float(rz)
            i2 += 1
        elif (int(ida) == 3):
            rv3[frame, i3, :] = float(rx), float(ry), float(rz)
            i3 += 1           
        if (int(ida) == 4): 
            rv4[frame, i4, :] = float(rx), float(ry), float(rz)
            i4 += 1
        elif (int(ida) == 5): 
            rv5[frame, i5, :] = float(rx), float(ry), float(rz)
            i5 += 1
        elif (int(ida) == 6):
            rv6[frame, i6, :] = float(rx), float(ry), float(rz)
            i6 += 1   
        elif (int(ida) == 7):
            rv7[frame, i7, :] = float(rx), float(ry), float(rz)
            i7 += 1   
        elif (int(ida) == 8):
            rv8[frame, i8, :] = float(rx), float(ry), float(rz)
            i8 += 1   
        elif (int(ida) == 9):
            rv9[frame, i9, :] = float(rx), float(ry), float(rz)
            i9 += 1   
        elif (int(ida) == 10):
            rv10[frame, i10, :] = float(rx), float(ry), float(rz)
            i10 += 1   
f.close()

f = open('positions_c1.txt', 'w')
for n in range(N1):
    f.write('create_atoms \t %u \t single \t %f \t %f \t %f \t remap yes \t units box\n' % (11, rv1[0,n,0], rv1[0,n,1], rv1[0,n,2]))
f.close()
f = open('positions_c1b.txt', 'w')
for n in range(N2):
    f.write('create_atoms \t %u \t single \t %f \t %f \t %f \t remap yes \t units box\n' % (12, rv2[0,n,0], rv2[0,n,1], rv2[0,n,2]))
f.close()
f = open('positions_c2.txt', 'w')
for n in range(N3):
    f.write('create_atoms \t %u \t single \t %f \t %f \t %f \t remap yes \t units box\n' % (13, rv3[0,n,0], rv3[0,n,1], rv3[0,n,2]))
f.close()
f = open('positions_c3.txt', 'w')
for n in range(N4):
    f.write('create_atoms \t %u \t single \t %f \t %f \t %f \t remap yes \t units box\n' % (14, rv4[0,n,0], rv4[0,n,1], rv4[0,n,2]))
f.close()
f = open('positions_c2b.txt', 'w')
for n in range(N5):
    f.write('create_atoms \t %u \t single \t %f \t %f \t %f \t remap yes \t units box\n' % (15, rv5[0,n,0], rv5[0,n,1], rv5[0,n,2]))
f.close()
f = open('positions_c3b.txt', 'w')
for n in range(N6):
    f.write('create_atoms \t %u \t single \t %f \t %f \t %f \t remap yes \t units box\n' % (16, rv6[0,n,0], rv6[0,n,1], rv6[0,n,2]))
f.close()
f = open('positions_n1.txt', 'w')
for n in range(N7):
    f.write('create_atoms \t %u \t single \t %f \t %f \t %f \t remap yes \t units box\n' % (17, rv7[0,n,0], rv7[0,n,1], rv7[0,n,2]))
f.close()
f = open('positions_n2.txt', 'w')
for n in range(N8):
    f.write('create_atoms \t %u \t single \t %f \t %f \t %f \t remap yes \t units box\n' % (18, rv8[0,n,0], rv8[0,n,1], rv8[0,n,2]))
f.close()
f = open('positions_n1b.txt', 'w')
for n in range(N9):
    f.write('create_atoms \t %u \t single \t %f \t %f \t %f \t remap yes \t units box\n' % (19, rv9[0,n,0], rv9[0,n,1], rv9[0,n,2]))
f.close()
f = open('positions_n2b.txt', 'w')
for n in range(N10):
    f.write('create_atoms \t %u \t single \t %f \t %f \t %f \t remap yes \t units box\n' % (20, rv10[0,n,0], rv10[0,n,1], rv10[0,n,2]))
f.close()
"""
1 c1
2 c1b
3 c2
4 c3
5 c2b
6 c3b
7 n1
8 n2
9 n1b
10 n2b
"""



