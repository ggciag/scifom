import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri

#Função que cria uma malha irregular a das coordenadas dos cantos da malha
def irregular_mesh(x,y,perturbation=0.2):
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

def nearest(Xn,Yn,Map,x,y):
  i = 0
  mapa = x*0
  for xi,yi in zip(x,y):
    dist = (Xn-xi)**2 + (Yn-yi)**2
    j = np.argmin(dist)
    mapa[i] = Map[j]
    i+=1
  return mapa


topo = np.loadtxt("../ETOPO10.txt")

topo = np.reshape(topo,(1801,3601))

latmin = -60.0
latmax = 15.0

longmin = -85.0
longmax = -30.0


ilongmin = int((longmin+180)*10)
ilongmax = int((longmax+180)*10)

ilatmin = int((90-latmin)*10)
ilatmax = int((90-latmax)*10)

topo = topo[ilatmax:ilatmin+1,ilongmin:ilongmax+1]


plt.figure(figsize=(10,10))
plt.imshow(topo)
plt.savefig("topo_Amazon.png",dpi=300)
plt.close()

resolutions = [0.1,0.2,0.5]
resol = resolutions[2]


#Resolução 0.1 x 0.1 graus
if (resol==resolutions[0]):
    nx = int(longmax-longmin)*10+1
    x = np.linspace(longmin,longmax,nx)
    ny = int(latmax-latmin)*10+1
    y = np.linspace(latmin,latmax,ny)
    topo = topo[::-1,:]

#Resolução 0.2 x 0.2 graus
elif (resol==resolutions[1]):
    nx = int(longmax-longmin)*5+1
    x = np.linspace(longmin,longmax,nx)
    ny = int(latmax-latmin)*5+1
    y = np.linspace(latmin,latmax,ny)
    topo = topo[::-2,::2]

#Resolução 0.5 x 0.5 graus
elif (resol==resolutions[2]):
    nx = int((longmax-longmin)*2+1)
    x = np.linspace(longmin,longmax,nx)
    ny = int((latmax-latmin)*2+1)
    y = np.linspace(latmin,latmax,ny)
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

lat = np.copy(Y)
lon = np.copy(X)

X = (X-np.min(x)+10.0)
Y = (Y-np.min(y)+10.0)

topo = np.reshape(topo,np.size(topo))
topo_copy = np.copy(topo)

##### Remove Andes

cond =               (lon<-63.0) & (lat<-15.0) & (topo>0.0)
topo[cond] *= 0.001

cond = (lon>-81.0) & (lon<-66.0) & (lat>=-15.0) & (lat<0.0) & (topo>0.0)
topo[cond] *= 0.001

cond = (lon>-81.0) & (lon<-69.0) & (lat>=-15.0) & (topo>0.0)
topo[cond] *= 0.001

uplift_rate = (topo_copy - topo)
uplift_rate[topo_copy<800.0] = 0.0
uplift_rate *=0.18/1.0E6



moho = X*0 + 35000.0
lithology = X*0 + 800000.0

#####################
#####################
######### Te map

Te = X*0

x_Te_all,y_Te_all,Te_all = np.loadtxt("../te_global.xyz",unpack=True)
cond = np.isnan(Te_all)
Te_all[cond] = 10
Te = nearest(x_Te_all,y_Te_all,Te_all,lon,lat)
Te *=1000.0


plt.figure()
plt.scatter(lon,lat,c=Te,cmap="turbo")
plt.colorbar()
plt.savefig("Te_map.png")
plt.close()


# translation in x and y (temporary) !!! see (X-np.min(x)+10.0)*100000,(Y-np.min(y)+10)*100000
#Pontos - graus em m
np.savetxt(
    "pontos.txt",
    np.c_[X*100000,Y*100000],#np.c_[X*100000,Y*100000],
    fmt="%.2f"
)

print(np.min(X),np.max(X),np.min(Y),np.max(Y))

#########
#Malha triangular
#Topo, moho, lito, Te
print(np.shape(topo), np.shape(moho), np.shape(lithology), np.shape(Te))
np.savetxt("topo_moho_lito_Te.txt",np.c_[topo,moho,lithology,Te],fmt="%.2f")

# Criar triangulação
triang = tri.Triangulation(X, Y)

# Compute the area of each triangle
xt = X[triang.triangles]
yt = Y[triang.triangles]

# Shoelace formula for triangle area (absolute value)
areas = 0.5 * np.abs(
    (xt[:,0]*(yt[:,1]-yt[:,2]) +
     xt[:,1]*(yt[:,2]-yt[:,0]) +
     xt[:,2]*(yt[:,0]-yt[:,1]))
)

print(np.min(areas),np.max(areas))

plt.figure()
plt.hist(areas)
plt.savefig("Tri_areas.png")
plt.close()


# Mask degenerate (zero or tiny area) triangles
tol = 1e-12
mask = areas < tol
triang.set_mask(mask)

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

#Plota topo
cores=[(0.0, 0, 1.0),
       (0.09803921568627451, 0.5137254901960786, 0.1568627450980392),
       (0.15, 0.596078431372549, 0.19607843137254902),
       (0.18, 0.7, 0.3),
       (0.2, 0.8, 0.4),
       (0.5, 0.85, 0.47),
       (0.7333333333333333, 0.8941176470588236, 0.5725490196078431),
       (0.83, 0.87, 0.65),
       (1.0, 0.8627450980392157, 0.7254901960784313),
       (0.9529411764705882, 0.792156862745098, 0.5372549019607843),
       (0.9019607843137255, 0.7215686274509804, 0.34509803921568627),
       (0.8509803921568627, 0.6509803921568628, 0.15294117647058825),
       (0.6588235294117647, 0.6039215686274509, 0.12156862745098039),
       (0.6431372549019608, 0.5647058823529412, 0.09803921568627451),
       (0.6352941176470588, 0.5254901960784314, 0.07450980392156863),
       (0.6235294117647059, 0.4823529411764706, 0.050980392156862744),
       (0.617, 0.46, 0.038),
       (0.611764705882353, 0.44313725490196076, 0.027450980392156862),
       (0.6, 0.4, 0.0),
       (0.617, 0.37, 0.18),
       (0.6352941176470588, 0.34901960784313724, 0.34901960784313724),
       (0.6980392156862745, 0.4627450980392157, 0.4627450980392157),
       (0.705, 0.51, 0.51),
       (0.7176470588235294, 0.5764705882352941, 0.5764705882352941),
       (0.7607843137254902, 0.6901960784313725, 0.6901960784313725),
       (0.8, 0.8, 0.8),
       (0.8980392156862745, 0.8980392156862745, 0.8980392156862745),
       (0.9490196078431372, 0.9490196078431372, 0.9490196078431372),
       (1.0, 1.0, 1.0)]

levels=[-10000,      0,    200,  400,    800,   1200,   1600,
   1800,   2000,   2200,   2400,   2600,   2800,   3100,   3400,   3700,
   4000,   4300,   4600,   4900,   5200,   5500,   5800,   8000]

plt.figure(figsize=(20,10))
plt.axis("equal")
plt.tricontourf(triang, topo,levels=levels,colors=cores)
plt.colorbar()
plt.savefig("Topo.png")
plt.close()

plt.figure(figsize=(20,10))
plt.axis("equal")
plt.tricontourf(triang, Te,levels=30)
plt.colorbar()
plt.savefig("Te.png")
plt.close()

plt.figure(figsize=(20,10))
plt.axis("equal")
plt.tricontourf(triang, lithology,levels=30)
plt.colorbar()
plt.savefig("Lithology.png")
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

nx_flexural {(ny-1)//4+1}
ny_flexural {(nx-1)//4+1}

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
n_sub_dt 1000

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
np.savetxt("uplift_map.txt",np.c_[uplift_rate],fmt="%lg")

plt.figure(figsize=(20,10))
plt.axis("equal")
plt.tricontourf(triang, uplift_rate*1.0E6,levels=30)
plt.colorbar()
plt.savefig("Uplift.png")
plt.close()


Prec = X*0 + 1.0
#cond = (lat>-10) & (lat<5) & (lon>-80) & (lon<-72)
#Prec[cond] = 3.0

for i in range(41):
    np.savetxt("Prec_%d.txt"%(i),np.c_[Prec, Prec, Prec],fmt="%.2f")

plt.figure(figsize=(20,10))
plt.axis("equal")
plt.tricontourf(triang, Prec,levels=30)
plt.colorbar()
plt.savefig("Prec.png")
plt.close()
