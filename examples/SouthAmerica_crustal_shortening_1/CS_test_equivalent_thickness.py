import numpy as np
import matplotlib.pyplot as plt

A = np.loadtxt("crustal_shortening.txt",skiprows=2)

time = A[::-1,0]


short = A[:,1:20:2]

lat_vec = np.array([38,36,33.5,30.5,25,20.5,18,15,10,5])
lat_vec *= -1.0


print(np.shape(short))
print(np.shape(lat_vec))

lat_interp = np.arange(-40,5,0.1)

cumulat = lat_interp * 0.0

hf_all = []
short_all = []
cumulat_all = []
hcr_all = []

plt.figure()

hc0 = 38.0
li = 350.0 + 0.0*lat_interp
tmax = 50
cont = 0
dt = 1.0

hf = hc0 + 0.0*lat_interp

for t in range(51):
    selec_short = short[tmax-t,:]
    short_interp = np.interp(lat_interp,lat_vec, selec_short)

    cumulat += short_interp

    dl = np.copy(short_interp)


    hcr = hc0*dl/li/dt
    hf += hcr*dt

    short_all = np.append(short_all,short_interp)
    cumulat_all = np.append(cumulat_all,cumulat)
    hf_all = np.append(hf_all,hf)
    hcr_all = np.append(hcr_all,hcr)

    plt.plot(lat_interp,short_interp)
    cont+=1

plt.savefig("short_interp.png")
plt.close()

linhas = np.size(short_all)//cont

short_all = np.reshape(short_all,(cont,linhas))
plt.contourf(lat_interp, time, short_all,levels=100)
plt.ylim(50,0)
plt.colorbar()
plt.savefig("short_all.png")
plt.close()

cumulat_all = np.reshape(cumulat_all,(cont,linhas))
plt.contourf(lat_interp, time, cumulat_all,levels=100)
plt.ylim(50,0)
plt.colorbar()
plt.savefig("cumulat_all.png")
plt.close()

hf_all = np.reshape(hf_all,(cont,linhas))
plt.contourf(lat_interp, time, hf_all,levels=100)
plt.ylim(50,0)
plt.colorbar()
plt.savefig("hf_all.png")
plt.close()

hcr_all = np.reshape(hcr_all,(cont,linhas))
plt.contourf(lat_interp, time, hcr_all,levels=100)
plt.ylim(50,0)
plt.colorbar()
plt.savefig("hcr_all.png")
plt.close()



"""
hc = 38.0
dl = 288.0
li = 350.0
dt = 50.0
hcr = hc*dl/(li-dl)/dt
hf = hc*li/(li-dl)
print("hcr = %.3f km/Myr, hf = %.3f km"%(hcr,hf))



###################################################


hc = 38.0
dl = 20.0
li = 100.0
dt = 5.0
hcr = hc*dl/(li-dl)/dt
hf = hc*li/(li-dl)
print("hcr = %.3f km/Myr, hf = %.3f km"%(hcr,hf))


hc = 38.0
dl = 4.0
li = 100.0
dt = 1.0
hcr = hc*dl/(li-dl)/dt
hf = hc*li/(li-dl)
print("hcr = %.3f km/Myr, hf = %.3f km"%(hcr,hf))

hc = hf
dl = 4.0
li -= dl
dt = 1.0
hcr = hc*dl/(li-dl)/dt
hf = hc*li/(li-dl)
print("hcr = %.3f km/Myr, hf = %.3f km"%(hcr,hf))

hc = hf
dl = 4.0
li -= dl
dt = 1.0
hcr = hc*dl/(li-dl)/dt
hf = hc*li/(li-dl)
print("hcr = %.3f km/Myr, hf = %.3f km"%(hcr,hf))

hc = hf
dl = 4.0
li -= dl
dt = 1.0
hcr = hc*dl/(li-dl)/dt
hf = hc*li/(li-dl)
print("hcr = %.3f km/Myr, hf = %.3f km"%(hcr,hf))

hc = hf
dl = 4.0
li -= dl
dt = 1.0
hcr = hc*dl/(li-dl)/dt
hf = hc*li/(li-dl)
print("hcr = %.3f km/Myr, hf = %.3f km"%(hcr,hf))





#print(t)

"""

