import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri

#Função que cria uma malha irregular a das coordenadas dos cantos da malha
def irregular_mesh(x,y,perturbation=0.25):
  """This function returns three arrays: two arrays with shape (ny, nx) with the x and y coordinates of the points of an irregular mesh;
  xy array with shape (ny*nx,2) containing all the (x,y) points
  """
  dx = np.abs(x[1]-x[0])
  dy = np.abs(y[1]-y[0])
  x_mesh = np.ndarray(shape=(len(y),len(x)),dtype=float)
  y_mesh = np.ndarray(shape=(len(y),len(x)),dtype=float)
  x_mesh[0:len(y),:]=x

  for i in range(len(x)):
    y_mesh[:,i]=y
  cond = ((x_mesh>x.min()) & (x_mesh<x.max()) & (y_mesh>y.min()) & (y_mesh<y.max()))
  np.random.seed(1)
  x_mesh[cond] += dx*perturbation*((np.random.rand(len(y),len(x))-0.5)*2)[cond]
  y_mesh[cond] += dy*perturbation*((np.random.rand(len(y),len(x))-0.5)*2)[cond]
  
  return x_mesh, y_mesh


topo = np.loadtxt("../ETOPO10.txt")

topo = np.reshape(topo,(1801,3601))

latmin = -15.0
latmax = 15.0

longmin = -85.0
longmax = -40.0


ilongmin = int((longmin+180)*10)
ilongmax = int((longmax+180)*10)

ilatmin = int((90-latmin)*10)
ilatmax = int((90-latmax)*10)

topo = topo[ilatmax:ilatmin+1,ilongmin:ilongmax+1]


plt.figure(figsize=(10,10))
plt.imshow(topo)
plt.savefig("topo_Amazon.png",dpi=300)
plt.close()


#Resolução 0.1 x 0.1 graus
#nx = (longmax-longmin)*10+1
#x = np.linspace(longmin,longmax,nx)
#ny = (latmax-latmin)*10+1
#y = np.linspace(latmin,latmax,ny)
#topo = topo[::-1,:]

#Resolução 0.2 x 0.2 graus
#nx = (longmax-longmin)*5+1
#x = np.linspace(longmin,longmax,nx)
#ny = (latmax-latmin)*5+1
#y = np.linspace(latmin,latmax,ny)
#y -= np.min(y)
#x -= np.min(x)
#topo = topo[::-2,::2]

#Resolução 0.5 x 0.5 graus
nx = int((longmax-longmin)*2+1)
x = np.linspace(longmin,longmax,nx)
ny = int((latmax-latmin)*2+1)
y = np.linspace(latmin,latmax,ny)
y -= np.min(y) - 10.0
x -= np.min(x) - 10.0
topo = topo[::-5,::5]


print(np.size(x),np.size(y),np.shape(topo))
print(int(nx*ny),int((nx-1)*(ny-1)*2))

#X,Y = np.meshgrid(x,y)
X,Y = irregular_mesh(x,y)


plt.figure(figsize=(10,10))
plt.contourf(x,y,topo,100)
plt.axis("equal")
plt.savefig("topo_Amazon_mesh.png",dpi=300)
plt.close()


X = np.reshape(X,np.size(X))
Y = np.reshape(Y,np.size(Y))

topo = np.reshape(topo,np.size(topo))

moho = X*0 + 35000.0
lithology = X*0 + 800000.0
Te = X*0 + 20000.0


#Pontos - graus em m
np.savetxt(
    "pontos.txt",
    np.c_[X*100000,Y*100000],
    fmt="%.2f"
)


#########
#Malha triangular
#Topo, moho, lito, Te
print(np.shape(topo), np.shape(moho), np.shape(lithology), np.shape(Te))
np.savetxt("topo_moho_lito_Te.txt",np.c_[topo,moho,lithology,Te],fmt="%.2f")

# Criar triangulação
triang = tri.Triangulation(X, Y)

# Obter os índices dos triângulos
indices = triang.triangles   # array (ntri, 3)

# Número de triângulos
ntri = indices.shape[0]

# Exportar com cabeçalho = número de triângulos
np.savetxt(
    "malha_externa.txt",
    indices,
    fmt="%d",
    header=str(ntri),
    comments=''   # evita o "#" no começo
)

# Plotar apenas a malha
plt.figure(figsize=(20,10))
plt.triplot(triang, color="k")
#plt.plot(X, Y, "o", color="red")
plt.gca().set_aspect("equal")  # mesma proporção
plt.savefig("Triagular_mesh.png",dpi=500)
plt.close()


#Falhas
#open("falhas.txt", "w").close()
with open("falhas.txt", "w") as f:
    f.write("0")



#Parametros gerais

params = f"""
maxy {np.max(Y)*100000}
miny {np.min(Y)*100000}
ny {ny}
nx {nx}

nx_flexural {(ny-1)//2+1}
ny_flexural {(nx-1)//2+1}

axis_stream 500000.0

Te_rigida 70000.0
Te_offshore 15000.0

vR 2.0
time_ofchangevR 20000000.0
vR2 2.0

vRandes 2.0
time_ofchangevRandes 20000000.0
vR2andes 2.0

Kf 0.08
Km 1000.0

ls 200.0
lb 800000.0
lb2 4000000.0

uplift_scale 1.0
time_ofchangeu 1000000.0
uplift_scale2 1.0

tempo_max 40000000.0
dt 200

"""

# Create the parameter file
with open("param_OrogSedFlex_1.1.txt", "w") as f:
    for line in params.split("\n"):
        line = line.strip()
        if len(line):
            f.write(" ".join(line.split()) + "\n")


##############
##############
##############
##############
# Create the parameter file for dynamic topography
params = f"""
1
1.0
40.0
40.0
0.0
"""

# Create the parameter file
with open("param_topo_din.txt", "w") as f:
    for line in params.split("\n"):
        line = line.strip()
        if len(line):
            f.write(" ".join(line.split()) + "\n")


#dynamic topography (dt1 dt2)
np.savetxt("topo_din_map.txt",np.c_[topo*0,topo*0],fmt="%.2f")



##############
##############
##############
##############
# Create the parameter file for uplift
params = f"""
1
1.0 0.0 40000000.0
"""

# Create the parameter file
with open("param_uplift.txt", "w") as f:
    for line in params.split("\n"):
        line = line.strip()
        if len(line):
            f.write(" ".join(line.split()) + "\n")

#uplift map
np.savetxt("uplift_map.txt",np.c_[topo*0],fmt="%.2f")
