import numpy as np
import matplotlib.pyplot as plt
import math

# THE PROGRAM IS RESTRICTED TO SOLVE ONLY CIRCULAR CROSS SECTIONS

#..........FUNCTION TO CALCULATE TRANSFORMED LOCAL STIFFNESS MATRIX..........#
def Kg(L,thta,E,R):
  n = 6
  K_e = np.zeros((n,n))
  T_e = np.zeros((n,n))
  r = np.radians(thta)
  l = np.cos(r)
  mu = np.sin(r)
  A = (math.pi)*R**2
  I = ((math.pi)*R**4)/4

  T_e[0,0] = T_e[1,1] = T_e[3,3] = T_e[4,4] = l
  T_e[2,2] = T_e[5,5] = 1
  T_e[0,1] = T_e[3,4] = mu
  T_e[1,0] = T_e[4,3] = -mu

  t1 = E*A/L
  t2 = E*I/(L**3)
  t3 = E*I/(L**2)
  t4 = E*I/L

  K_e[0,0] = K_e[3,3] = t1
  K_e[0,3] = K_e[3,0] = -t1
  K_e[1,1] = K_e[4,4] = 12*t2
  K_e[1,4] = K_e[4,1] = -12*t2
  K_e[1,2] = K_e[2,1] = K_e[1,5] = K_e[5,1] = 6*t3
  K_e[2,4] = K_e[4,2] = K_e[4,5] = K_e[5,4] = -6*t3
  K_e[2,5] = K_e[5,2] = 2*t4
  K_e[5,5] = K_e[2,2] = 4*t4

  T_e_T = np.transpose(T_e)
  K_g = np.dot(T_e_T,np.dot(K_e,T_e))
  return K_g

#..........FUNCTION TO CALCULATE GLOBAL STIFFNESS MATRIX..........#
def Km(x,y,K_g,num_nodes):
  m = (3*x) - 2 - 1
  n = (3*y) - 2 - 1
  K = np.zeros(((num_nodes*3),(num_nodes*3)))

  K[m,m]   = K_g[0,0] 
  K[m,m+1] = K[m+1,m] = K_g[0,1]
  K[m,m+2] = K[m+2,m] = K_g[0,2]
  K[m,n]   = K[n,m]   = K_g[0,3]
  K[m,n+1] = K[n+1,m] = K_g[0,4]
  K[m,n+2] = K[n+2,m] = K_g[0,5]

  K[m+1,m+1] = K_g[1,1]
  K[m+1,m+2] = K[m+2,m+1] = K_g[1,2]
  K[m+1,n]   = K[n,m+1]   = K_g[1,3]
  K[m+1,n+1] = K[n+1,m+1] = K_g[1,4]
  K[m+1,n+2] = K[n+2,m+1] = K_g[1,5]

  K[m+2,m+2] = K_g[2,2]
  K[m+2,n]   = K[n,m+2]   = K_g[2,3]
  K[m+2,n+1] = K[n+1,m+2] = K_g[2,4]
  K[m+2,n+2] = K[n+2,m+2] = K_g[2,5]

  K[n,n]   = K_g[3,3]
  K[n,n+1] = K[n+1,n] = K_g[3,4]
  K[n,n+2] = K[n+2,n] = K_g[3,5]

  K[n+1,n+1] = K_g[4,4]
  K[n+1,n+2] = K[n+2,n+1] = K_g[4,5]

  K[n+2,n+2] = K_g[5,5]
  return K

#..........Function to Extract Input Data from Input File..........#
def extract_data(a):

  #.....Extraction of Nodal Data and Storing them into a NODE Matrix .....#

  p = 0                           
  num_nodes = 0
  num_ele = 0
  for i in range (p,len(a)):
    if (a[i])[0] == '*ELEMENT':
      break
    num_nodes = i
  #print(num_nodes)
  NODE = np.zeros((num_nodes,4))
  for i in range (p,num_nodes):
    for j in range (0,4):
      ((a[i+1])[j]) = ((a[i+1])[j]).replace(',','')
      NODE[i][j] = float(((a[i+1])[j]))
  np.set_printoptions(suppress=True)
  NODE = NODE[np.argsort(NODE[:,0],kind = 'mergesort')]
  #print("Node Matrix = ")
  #print(NODE)

  #.....Extraction of Elemental Data and Storing them into a ELE Matrix .....#

  p = p+num_nodes+1
  for i in range (p,len(a)):
    if (a[i])[0] == '*MATERIAL':
      break
    num_ele = i - (p)
  #print(num_ele)
  ELE = np.zeros((num_ele,3))
  for i in range (0,num_ele):
    for j in range (0,3):
      ((a[p+1+i])[j]) = ((a[p+1+i])[j]).replace(',','')
      ELE[i][j] = float(((a[p+1+i])[j]))
  ELE = ELE.astype(int)
  ELE = ELE[np.argsort(ELE[:,0],kind = 'mergesort')]
  #print("Element Matrix = ")
  #print(ELE)

  #.....Extraction of Material Data and Storing them into a MAT Matrix .....#

  p = p + num_ele + 1
  num_mat = 0
  for i in range (p,len(a)):
    if (a[i])[0] == '*SECTION':
      break
    num_mat = i - (p)
  #print(num_mat)
  MAT = np.zeros((num_mat,3))
  for i in range (0,num_mat):
    for j in range (0,3):
      ((a[p+1+i])[j]) = ((a[p+1+i])[j]).replace(',','')
      MAT[i][j] = float(((a[p+1+i])[j]))
  #print("Material Matrix = ")
  np.set_printoptions(suppress=True)
  MAT = MAT[np.argsort(MAT[:,0],kind = 'mergesort')]
  #print(MAT)

  #.....Extraction of Section Propeties and Storing them into a SEC_P Matrix .....#

  p = p + num_mat + 1
  SEC_P = np.zeros((num_ele,3))
  for i in range (0,num_ele):
    for j in range (0,3):
      ((a[p+1+i])[j]) = ((a[p+1+i])[j]).replace(',','')
      SEC_P[i][j] = float(((a[p+1+i])[j]))
  SEC_P = SEC_P[np.argsort(SEC_P[:,0],kind = 'mergesort')]
  #print("Section Properties Matrix = ")
  #print(SEC_P)

  #.....Extraction of Boundary Conditions and Storing them into a BC Matrix .....#

  p = num_ele+1+p
  num_bc = 0
  for i in range (p,len(a)):
    if (a[i])[0] == '*FORCE':
      break
    num_bc = i - (p)
  #print(num_bc)
  BC = np.zeros((num_bc,4))
  for i in range (0,num_bc):
    for j in range (0,4):
      ((a[p+1+i])[j]) = ((a[p+1+i])[j]).replace(',','')
      BC[i][j] = float(((a[p+1+i])[j]))
  BC = BC.astype(int)
  BC = BC[np.argsort(BC[:,0],kind = 'mergesort')]
  #print("Boundary Condition Matrix = ")
  #print(BC)

  #.....Extraction of Forces and Storing them into a F Matrix .....#

  p = num_bc+1+p
  num_f = 0
  for i in range (p,len(a)):
    if (a[i])[0] == '*MOMENT':
      break
    num_f = i - (p)
  #print(num_f)
  F = np.zeros((num_f, 4))
  for i in range (0,num_f):
    for j in range (0,4):
      ((a[p+1+i])[j]) = ((a[p+1+i])[j]).replace(',','')
      F[i][j] = float(((a[p+1+i])[j]))
  F = F.astype(int)
  F = F[np.argsort(F[:,0],kind = 'mergesort')]
  #print("Force Matrix = ")
  #print(F)

  #.....Extraction of Moments and Storing them into a M Matrix .....#

  p = num_f + 1 + p
  num_m = 0
  for i in range (p,len(a)):
    num_m = i - (p)
  #print(num_m)
  M = np.zeros((num_m, 4))
  for i in range (0,num_m):
    for j in range (0,4):
      ((a[p+1+i])[j]) = ((a[p+1+i])[j]).replace(',','')
      M[i][j] = float(((a[p+1+i])[j]))
  M = M.astype(int)
  M = M[np.argsort(M[:,0],kind = 'mergesort')]
  #print("Moment Matrix = ")
  #print(M)
  return NODE, num_nodes, ELE, num_ele, MAT, num_mat, SEC_P, BC, num_bc, F, num_f, M, num_m

#..........FUNCTION TO CALCULATE LENGTH and THETA MATRIX..........#
def len_thta(NODE, ELE, num_ele, num_nodes):
  LEN = np.zeros((num_ele,2))
  THTA = np.zeros((num_ele,2))
  for i in range (0,num_ele):
    g = ELE[i,1] 
    h = ELE[i,2]
    for j in range (0,num_nodes):
        if NODE[j,0] == g:
          x1 = NODE[j,1]
          y1 = NODE[j,2]
          for k in range (0, num_nodes):
            if NODE[k,0] == h:
              x2 = NODE[k,1]
              y2 = NODE[k,2]
              LEN[i,1] = ((y2-y1)**2 + (x2-x1)**2)**(0.5)
              LEN[i,0] = THTA[i,0] = ELE[i,0]
              THTA[i,1] = np.degrees(math.atan2((y2-y1),(x2-x1)))
  print("LENGTH MATRIX = ")
  print(LEN)
  print("THETA MATRIX = ")
  print(THTA)
  return LEN, THTA

#..........FUNCTION TO CALCULATE GLOBAL STIFFNESS MATRIX..........#
def KGlo(NODE,ELE,MAT,SEC_P,LEN,THTA,num_ele,num_nodes,num_mat):
  K_glo = np.zeros(((num_nodes*3),(num_nodes*3)))

  for i in range (0,num_ele):
    for j in range (0,num_mat):
        if SEC_P[i,1] == MAT[j,0]:
            mt = MAT[j,1]
            K_eg = Kg(LEN[i,1], THTA[i,1], mt, SEC_P[i,2])
            K_ebig = Km(ELE[i,1],ELE[i,2],K_eg,num_nodes)
            K_glo = K_glo + K_ebig

  print("GLOBAL STIFFNESS MATRIX = ")
  print(K_glo)
  return K_glo

#.......FUNCTION TO CALCULATE DEGREE OF FREEDOM MATRIX.......#
def solveDof(b_c, F, M, K_glo, num_nodes):
    A = []
    B = []
    C = []
    m = (3 * num_nodes)
    u_refi = [0] * m
    i = 1
    for z in range(1, num_nodes+1):
        if i > len(b_c):
            i = i-1

        if (b_c[i - 1][0]) == z:
            for j in range(1, 4):
                if (b_c[i - 1, j]) != 0:
                    q = (3 * z - 3) + (j - 1)
                    A.append(q)
                    u_refi[q] = b_c[i - 1, j]
            i = i+1
        else:
            t = (3 * z - 3)
            A.append(t)
            A.append(t + 1)
            A.append(t + 2)
    i = 1
    s = 1
    for z in range(1, num_nodes+1):
        h = 1
        if i > len(F):
            i = i - 1
        if (F[i - 1][0]) == z:

            C.append(3 * z - 3)
            C.append(3 * z - 2)
            B.append(F[i - 1][1])
            B.append(F[i - 1][2])
            i = i+1
            h = 0
        else:
            for t in range(1, len(b_c) + 1):
                if (b_c[t - 1][0]) == z:
                    h = 0

                    for j in range(1, 3):
                        if (b_c[t - 1, j]) != 0:
                            q = (3 * z - 3) + (j - 1)
                            C.append(q)
                            B.append(0)

            for t in range(1, len(M) + 1):
                if (M[t - 1, 0]) == z:
                    if (M[t - 1, 3]) != 0:
                        C.append(3 * z - 3)
                        C.append(3 * z - 2)
                        B.append(0)
                        B.append(0)
                        h = 0

        if (s > len(M)):
            s = s - 1

        if (M[s - 1][0]) == z:

            C.append(3 * z - 1)
            B.append(M[s - 1][3])
            s = s+1
        else:
            for t in range(1, len(F) + 1):
                if (F[t - 1][0]) == z:
                    q = (3 * z - 1)
                    C.append(q)
                    B.append(0)
        if (h == 1):
         C.append(3 * z - 3)
         C.append(3 * z - 2)
         C.append(3 * z - 1)
         B.append(0)
         B.append(0)
         B.append(0)

    k = np.zeros((len(C), len(A)))
    k_L = np.zeros((len(C), len(A)))
    for i in range(0, len(C)):
        for j in range(0, len(A)):
            k[i][j] = K_glo[C[i]][A[j]]

    k_L = np.linalg.inv(k)
    D = np.dot(k_L, B)
    for i in range(len(D)):
        u_refi[A[i]] = D[i]

    return u_refi

#.......FUNCTION TO CALCULATE INTERNAL FORCE MATRIX.......#
def solveFglobal(u_refi, K_glo):
    F_intf = np.dot(K_glo, u_refi)
    return F_intf

#.......FUNCTION TO CALCULATE NODAL STRESS and STRAIN MATRIX.......#
def strain_stress(u_refi,F_intf,NODE,ELE,MAT,SEC_P,LEN,THTA,num_ele,num_nodes,num_mat):
  u = u_refi
  Q = F_intf
  e_elem = np.zeros((num_ele,5))
  e_node = np.zeros((num_nodes,2))
  S_node = np.zeros((num_nodes,2))
  #.......SUB FUNCTION TO CALCULATE NODAL STRESS and STRAIN MATRIX.......#
  def ex_bar(L,thta,node1,node2):
    r = np.radians(thta)
    l = np.cos(r)
    m = np.sin(r)
    T0=np.array([-l,-m,l,m])
    U0=np.array([[u[node1*3-3,0]],[u[node1*3-2,0]],[u[node2*3-3,0]],[u[node2*3-2,0]]])
    e1=1/L*(np.dot(T0,U0))
    return e1
  def ex_beam(x,L,thta,E,I,R,node1,node2):
    r = np.radians(thta)
    l = np.cos(r)
    m = np.sin(r)
    T1=np.array([m,-l,1,-m,l,1])
    U1=np.array([[Q[node1*3-3,0]*x],[Q[node1*3-2,0]*x],[Q[node1*3-1,0]],[Q[node2*3-3,0]],[Q[node2*3-2,0]*(L-x)],[Q[node2*3-1,0]*(L-x)]])
    e2=abs(R/(I*E)*(np.dot(T1,U1)))
    return e2
  for i in range (0,num_ele):  
    for j in range (0,num_mat):
      if SEC_P[i,1] == MAT[j,0]:
        E = MAT[j,1]
        R=SEC_P[i,2]
        I=math.pi*R**4/4
        A=math.pi*R**2
        e_elem[i,0] = ELE[i,0]
        e_elem[i,1] = ex_bar(LEN[i,1], THTA[i,1], ELE[i,1], ELE[i,2])
        e_elem[i,2] = ex_beam(0,LEN[i,1], THTA[i,1], E, I, R, ELE[i,1], ELE[i,2])
        e_elem[i,3] = ex_bar(LEN[i,1], THTA[i,1], ELE[i,1], ELE[i,2])
        e_elem[i,4] = ex_beam(LEN[i,1],LEN[i,1], THTA[i,1], E, I, R, ELE[i,1], ELE[i,2])
  for i in range (0,num_nodes):  
    n=0
    node=NODE[i,0]
    ebar=0
    ebeam=0
    Sbar=0
    Sbeam=0
    for j in range (0,num_ele):
      for h in range (0,num_mat):
        if SEC_P[j,1] == MAT[h,0]:
          if (ELE[j,1]==node):
            ebar=ebar+e_elem[j,1]
            ebeam=ebeam+e_elem[j,2]
            Sbar=Sbar+MAT[h,1]*e_elem[j,1]
            Sbeam=Sbeam+MAT[h,1]*e_elem[j,2]
            n=n+1
          elif (ELE[j,2]==node):
            ebar=ebar+e_elem[j,3]
            ebeam=ebeam+e_elem[j,4]
            Sbar=Sbar+MAT[h,1]*e_elem[j,3]
            Sbeam=Sbeam+MAT[h,1]*e_elem[j,4]
            n=n+1
    e_node[i,0]=node
    S_node[i,0]=node
    if ebar>0 :
     e_node[i,1]=(ebar+abs(ebeam))/n
    else:
     e_node[i,1]=(ebar-abs(ebeam))/n
    if Sbar>0 :
      S_node[i,1]=(Sbar+abs(Sbeam))/n
    else:
      S_node[i,1]=(Sbar-abs(Sbeam))/n
  print("STRAIN MATRIX =")
  print(e_node)
  print("STRESS MATRIX =")
  print(S_node)
  return e_node,S_node

#............FUNCTION TO PLOT UNDEFORMED and DEFORMED STRUCTURE...........#
def plot(NODE,num_nodes,ELE,num_ele,u_refi):
  Xud = np.zeros(num_nodes)
  Yud = np.zeros(num_nodes)
  Xd = np.zeros(num_nodes)
  Yd = np.zeros(num_nodes)
  x_ud = np.zeros(num_nodes)
  y_ud = np.zeros(num_nodes)
  for i in range(0,num_nodes):
    Xud[i] = NODE[i,1]
    Yud[i] = NODE[i,2]

  sf = 1000       # SCALING FACTOR TO SEE THE DEFORMATION QUALITATIVELY, IT CAN BE CHANGED AS PER REQUIREMENTS
  
  for i in range(0,num_nodes):
    Xd[i] = NODE[i,1] + sf*(u_refi[(3*i),0])
    Yd[i] = NODE[i,2] + sf*(u_refi[(3*i)+1,0])  

  for i in range (0,num_ele):
    x_ud = [Xud[ELE[i,1]-1],Xud[ELE[i,2]-1]]
    y_ud = [Yud[ELE[i,1]-1],Yud[ELE[i,2]-1]]
    plt.plot(x_ud,y_ud,color = 'blue')

  for i in range (0,num_ele):
    x_d = [Xd[ELE[i,1]-1],Xd[ELE[i,2]-1]]
    y_d = [Yd[ELE[i,1]-1],Yd[ELE[i,2]-1]]
    plt.plot(x_d,y_d,color = 'red')

  plt.scatter(Xud,Yud,color = 'black',label = 'UnDeformed')
  plt.scatter(Xd,Yd,color = 'yellow', label = 'Deformed')
  plt.xlabel('X')
  plt.ylabel('Y')
  plt.legend()
  plt.grid()
  plt.show()
  return

#..................MAIN CODE STARTS.................#

print("1D-FRAME ELEMENT PROGRAM: TEAM A --- FEM (ME 615) ")
print("THE PROGRAM IS RESTRICTED TO SOLVE ONLY CIRCULAR CROSS SECTIONS")
name = input("Enter the Input File Path with its Name and Extension or Just the Name if it is in the same Directory: ")
file = open(name, 'r')
line = file.readlines()
a = []
for lin in line:
	fields = lin.split()
	a.append(fields)          # GETTING ALL THE DATA IN SINGLE DYNAMIC ARRAY 'a'
file.close()   

#.....CALLING THE extract_data FUNCTION TO EXTRACT INPUT DATA.....#
NODE, num_nodes, ELE, num_ele, MAT, num_mat, SEC_P, BC, num_bc, F, num_f, M, num_m = extract_data(a)  

#.....CALLING THE len_thta FUNCTION TO GET VALUES OF LENGTH and THETA OF EACH ELEMENT.....#
LEN, THTA = len_thta(NODE, ELE, num_ele, num_nodes)

#.....CALLING THE KGlo FUNCTION TO GET GLOBAL STIFFNESS MATRIX.....#
K_glo = KGlo(NODE,ELE,MAT,SEC_P,LEN,THTA,num_ele,num_nodes,num_mat)

#.....CALLING THE solveDof & solveFglobal FUNCTION TO GET DEGREE OF FREEDOM and INTERNAL FORCE  MATRIX.....#
u_refi = solveDof(BC, F, M, K_glo, num_nodes)
u_refi = np.array(u_refi)
print("GLOBAL DEGREE OF FREEDOM MATRIX = ")
tolerance = 0
u_refi[np.isclose(u_refi, 0, atol=tolerance)] = 0
u_refi = u_refi.reshape((len(u_refi),1))
print(u_refi)

F_intf = solveFglobal(u_refi, K_glo)
print("GLOBAL INTERNAL FORCE MATRIX = ")
tolerance = 0
F_intf[np.isclose(F_intf, 0, atol=tolerance)] = 0
F_intf = F_intf.reshape((len(F_intf),1))
print(F_intf)

#.....CALLING THE strain_stress FUNCTION TO GET NODAL STRESS and STRAIN MATRIX.....#
strain,stress = strain_stress(u_refi,F_intf,NODE,ELE,MAT,SEC_P,LEN,THTA,num_ele,num_nodes,num_mat)

#.....CALLING THE plot FUNCTION TO PLOT DEFORMED and UNDEFORMED STRUCTURE.....#
plot(NODE,num_nodes,ELE,num_ele,u_refi)

# WRITING THE OUPUT FILE
file1 = open('output.txt', 'w')
u_str=''
e_str=''
S_str=''
for i in range (0,num_nodes*3):
  u_str=u_str+str(i+1)+','+' '+str(u_refi[i,0])+'\n'
for i in range (0,num_nodes):
  e_str=e_str+str(int(strain[i,0]))+','+' '+str(strain[i,1])+'\n'
  S_str=S_str+str(int(stress[i,0]))+','+' '+str(stress[i,1])+'\n'
L = ["*DISPLACEMENT \n",u_str, "\n*STRAIN \n",e_str, "\n*STRESS \n",S_str]
file1.writelines(L)
file1.close