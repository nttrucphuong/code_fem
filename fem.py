import numpy as np
def MT_cung_pt_tamgiac(xi, xj, xk, yi, yj, yk):

    xij = xi - xj
    xjk = xj - xk
    xik = xi - xk
    yij = yi - yj
    yjk = yj - yk
    yik = yi - yk


    k11 = yjk**2+ lamda*xjk**2
    k12 = -C2*xjk*yjk - lamda*yjk*xjk
    k13 = -yik*yjk - lamda*xjk*xik
    k14 = C2*xik*yjk + lamda*yik*xjk
    k15 = yjk*yij + lamda*xjk*xij
    k16 = -C2*yjk*xij - lamda*xjk*yij

    k21 = k12
    k22 = xjk**2 + lamda*yjk**2
    k23 = C2*xjk*yik + lamda*xik*yjk
    k24 = -xjk*xik - lamda*yjk*yik
    k25 = -C2*xjk*yij - lamda*yjk*xij
    k26 = xij*xjk + lamda*yjk*yij

    k31 = k13
    k32 = k23
    k33 = yik**2 + lamda*xik**2
    k34 = -C2*xik*yik - lamda*xik*yik
    k35 = -yik*yij - lamda*xik*xij
    k36 = C2*xij*yik + lamda*xik*yij

    k41 = k14
    k42 = k24
    k43 = k34
    k44 = xik**2 + lamda*yik**2
    k45 = C2*xik*yij + lamda*yik*xij
    k46 = -xik*xij - lamda*yik*yij

    k51 = k15
    k52 = k25
    k53 = k35
    k54 = k45
    k55 = yij**2 + lamda*xij**2
    k56 = -C2*xij*yij - lamda*xij*yij

    k61 = k16
    k62 = k26
    k63 = k36
    k64 = k46
    k65 = k56
    k66 = xij**2 + lamda*yij**2

    K_pt = np.array([[k11, k12, k13, k14, k15, k16],
                          [k21, k22, k23, k24, k25, k26],
                          [k31, k32, k33, k34, k35, k36],
                          [k41, k42, k43, k44, k45, k46],
                          [k51, k52, k53, k54, k55, k56],
                          [k61, k62, k63, k64, k65, k66]])
    return np.round_((C1*h/(4*Ae))*K_pt/(so_chia),decimals = 4)

def MT_cung_pt_tugiac(a,b,lamda,C1,C2):

    # 2a: tong chieu dai
    # 2b: tong chieu rong

    k11 = 4*(b**2+ lamda*a**2)/3
    k12 = a*b*(lamda + C2)
    k13 = 2*(lamda*a**2 - 2*b**2)/3
    k14 = (C2 - lamda)*a*b
    k15 = -2*(b**2 + lamda*a**2)/3
    k16 = -a*b*(lamda + C2)
    k17 = 2*(b**2 - 2*lamda*a**2)/3
    k18 = a*b*(lamda - C2)

    k21 = k12
    k22 = 4*(a**2 + lamda*b**2)/3
    k23 = k18
    k24 = 2*(a**2 - 2*lamda*b**2)/3
    k25 = k16
    k26 = -2*(a**2 + lamda*b**2)/3
    k27 = k14
    k28 = 2*(lamda*b**2 - 2*a**2)/3

    k31 = k13
    k32 = k23
    k33 = k11
    k34 = k16
    k35 = k17
    k36 = k14
    k37 = k15
    k38 = k12

    k41 = k14
    k42 = k24
    k43 = k34
    k44 = k22
    k45 = k18
    k46 = k28
    k47 = k12
    k48 = k26

    k51 = k15
    k52 = k25
    k53 = k35
    k54 = k45
    k55 = k11
    k56 = k12
    k57 = k13
    k58 = k14

    k61 = k16
    k62 = k26
    k63 = k36
    k64 = k46
    k65 = k56
    k66 = k22
    k67 = k18
    k68 = k24

    k71 = k17
    k72 = k27
    k73 = k37
    k74 = k47
    k75 = k57
    k76 = k67
    k77 = k11
    k78 = k16

    k81 = k18
    k82 = k28
    k83 = k38
    k84 = k48
    k85 = k58
    k86 = k68
    k87 = k78
    k88 = k22


    K_pt = np.array([[k11, k12, k13, k14, k15, k16, k17, k18],
                      [k21, k22, k23, k24, k25, k26, k27, k28],
                      [k31, k32, k33, k34, k35, k36, k37, k38],
                      [k41, k42, k43, k44, k45, k46, k47, k48],
                      [k51, k52, k53, k54, k55, k56, k57, k58],
                      [k61, k62, k63, k64, k65, k66, k67, k68],
                      [k71, k72, k73, k74, k75, k76, k77, k78],
                      [k81, k82, k83, k84, k85, k86, k87, k88]])
    return np.round_((C1*h/(4*a*b))*K_pt/(so_chia),decimals = 4)

def MT_Se(a,b,xo,yo,lamda,C1,C2):
    x1_ = x2 - xo
    x2_ = x3 - xo
    x3_ = x5 - xo
    x4_ = x6 - xo

    y1_ = y2 - yo
    y2_ = y3 - yo
    y3_ = y5 - yo
    y4_ = y6 - yo
    Se = np.array([[-y4_,       -C2*x2_,    y4_,       C2*x1_,    -y1_,       -C2*x1_,    y1_,       C2*x2_],
                     [-C2*y4_,    -x2_,       C2*y4_,    x1_,       -C2*y1_,    -x1_,       C2*y1_,    x2_],
                     [-lamda*x2_, -lamda*y4_, lamda*x1_, lamda*y4_, -lamda*x1_, -lamda*y1_, lamda*x2_, lamda*y1_]])
    return np.round_((C1/(4*a*b))*Se/(10**10),decimals = 4)


def fem_calculate(x,y, stress=True):
        
    global nuy, E, Ae, lamda, C1, C2, h,x1, x2, x3, x4, x5, x6, y1, y2, y3, y4, y5, y6, so_chia
    if stress:
        so_chia=10**9
    else:
        so_chia=10**11
    nuy = 0.3
    E = 2e11

    # x = 1
    # y =1
    x1 = 0
    x2 = x1 + x
    x3 = x2 + y
    x4 = x3 + x
    x5 = x3
    x6 = x2

    y1 = y2 = y3 = y4 =0
    y5 =y1 + x
    y6 = y5

    Ae = 0.5*x**2

    if stress:
        lamda = (1-nuy)/2
        C1 = E/(1-nuy**2)
        C2 = nuy
        h = 0.01
    else:
        lamda = (1-2*nuy)/(2*(1-nuy))
        C1 = ((1-nuy)*E)/((1+nuy)*(1-2*nuy))
        C2 = nuy/(1-nuy)
        h=1


    K1_stress = MT_cung_pt_tamgiac(x1,x2,x6,y1,y2,y6)
    K3_stress = MT_cung_pt_tamgiac(x3,x4,x5,y3,y4,y5)


    a = x/2
    b = y/2
    K2_stress = MT_cung_pt_tugiac(a,b,lamda,C1,C2)

    K_tt = np.zeros((12,12))

    Bool_1 = [0,1,2,3,10,11]
    n =0
    m =0
    # print(K1_stress[6,7])
    for i in Bool_1:
        for j in Bool_1:
            K_tt[i,j] = K_tt[i,j] + K1_stress[n,m]
            m += 1
            # print(K_stress)
        m = 0
        n += 1

    Bool_2 = [2,3,4,5,8,9,10,11]

    m = 0
    n = 0

    for i in Bool_2:
        for j in Bool_2:
            K_tt[i,j] = K_tt[i,j] + K2_stress[n,m]
            m += 1
            # print(K_stress)
        m = 0
        n += 1

    Bool_3 = [4,5,6,7,8,9]

    m = 0
    n = 0

    for i in Bool_3:
        for j in Bool_3:
            K_tt[i,j] = K_tt[i,j] + K3_stress[n,m]
            m += 1
            # print(K_stress)
        m = 0
        n += 1

    Se = MT_Se(x/2,y/2,x+y/2,x/2,lamda,C1,C2)

    # print('10e9*','\n', K1_stress)
    # print('10e9*','\n',K2_stress)
    # print('10e9*','\n',K3_stress)
    # print('10e9*','\n',K_tt) 
    return K1_stress, K2_stress, K3_stress, K_tt, Se
