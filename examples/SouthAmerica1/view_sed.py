import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from matplotlib.colors import from_levels_and_colors
import sys

x,y = np.loadtxt("m_xy.txt",unpack=True,skiprows=1)
Tri = np.loadtxt("m_Tri.txt",unpack=True,skiprows=1)
Tri = np.transpose(Tri)

def make_rivers(a,b):
    n = np.size(a)
    f = a*np.nan
    c = np.append(a,b)
    c = np.append(c,f)
    c = np.reshape(c,(3,n))
    c = np.transpose(c)

    c = np.reshape(c,(3*n))
    return(c)


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

print(np.shape(levels))
print(np.shape(cores))


triang = mtri.Triangulation(x, y, Tri)

Qr_limits = [1000,5.0E10,5.0E11]
width_limits = [0.3,0.8,2]

argv = sys.argv


if len(argv)>1:
    time_start = float(sys.argv[1])
    time_stop = float(sys.argv[2])
    time_step = float(sys.argv[3])
else:
    time_start = 2.0
    time_stop = 9.0
    time_step = 1.0

for t in np.arange(time_start,time_stop,time_step):
    print(t)
    plt.figure(figsize=(30,30))
    plt.axis("equal")
    topo,bed,Qr,direc,lake,prec = np.loadtxt("m_Topo_%.3f.txt"%(t),unpack=True)
    direc = direc.astype("int")
    index = np.arange(np.size(direc)).astype("int")
    

    dist = (x[index]-x[direc])**2 + (y[index]-y[direc])**2

    for i in range(3):
        cond = (Qr>Qr_limits[i]) & (dist<200.0E3**2)

        id1_a = index[cond]
        id2_a = direc[cond]

        llo = make_rivers(x[id1_a],x[id2_a])
        lla = make_rivers(y[id1_a],y[id2_a])
        plt.plot(llo,lla,"b",linewidth=width_limits[i])


    cond = lake==0

    #plt.tricontourf(triang, topo,levels=100, vmin=0,vmax=4)
    plt.tricontourf(triang, topo-bed,levels=100,vmin=0.0,vmax=500.0)
    plt.colorbar()

    plt.savefig("Sed_%07.3f.png"%(t))

    plt.close()
