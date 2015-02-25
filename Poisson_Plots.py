#!/usr/bin/env python

#Author: Craig Lage, NYU; 
#Date: 26-Jan-15

#This program plots the Poisson equation solutions from the C++ Poisson solver
import matplotlib
matplotlib.use("PDF")
from pylab import *
import sys, time, h5py

#****************SUBROUTINES*****************
class Array3dHDF5(object):
    def __init__(self, dir, filebase, LogEField):
        phifile = dir+'/'+filebase+'_phi'
        rhofile = dir+'/'+filebase+'_rho'
        hdfphi = h5py.File(phifile,'r')
        Dimension = hdfphi[hdfphi.items()[0][0]].attrs[u'Dimension']
        self.nx=Dimension[0]
        self.ny=Dimension[1]
        self.nz=Dimension[2]
        
        Lower_Left = hdfphi[hdfphi.items()[0][0]].attrs[u'Lower_Left']
        self.xmin=Lower_Left[0]
        self.ymin=Lower_Left[1]
        self.zmin=Lower_Left[2]

        Upper_Right = hdfphi[hdfphi.items()[0][0]].attrs[u'Upper_Right']
        self.xmax=Upper_Right[0]
        self.ymax=Upper_Right[1]
        self.zmax=Upper_Right[2]
        
        self.dx=(self.xmax-self.xmin)/self.nx
        self.dy=(self.ymax-self.ymin)/self.ny
        self.dz=(self.zmax-self.zmin)/self.nz
        self.volume = self.dx * self.dy * self.dz
        
        self.x=linspace(self.xmin+self.dx/2,self.xmax-self.dx/2,self.nx)
        self.y=linspace(self.ymin+self.dy/2,self.ymax-self.dy/2,self.ny)
        self.z=linspace(self.zmin+self.dz/2,self.zmax-self.dz/2,self.nz)

        self.phi=array(hdfphi[hdfphi.items()[0][0]])
        hdfrho = h5py.File(rhofile,'r')
        self.rho=array(hdfrho[hdfrho.items()[0][0]])
        if LogEField == 1:
            Exfile = dir+'/'+filebase+'_Ex'
            Eyfile = dir+'/'+filebase+'_Ey'
            Ezfile = dir+'/'+filebase+'_Ez'
            hdfEx = h5py.File(Exfile,'r')
            self.Ex=array(hdfEx[hdfEx.items()[0][0]])
            hdfEy = h5py.File(Eyfile,'r')
            self.Ey=array(hdfEy[hdfEy.items()[0][0]])
            hdfEz = h5py.File(Ezfile,'r')
            self.Ez=array(hdfEz[hdfEz.items()[0][0]])

def ReadConfigFile(filename):
    # This reads the config file for the necessary settings
    # and returns a dictionary with the values
    file = open(filename,'r')
    lines=file.readlines()
    file.close()

    ParamNames = ['outputfiledir', 'outputfilebase', 'ScaleFactor', 'GridsPerPixel', 'PixelBoundaryLowerLeft', 'PixelBoundaryUpperRight', 'PixelBoundaryStepSize', 'LogPixels', 'LogPixelPaths', 'LogEField', 'EdgePlot', 'PlotEField']
    ConfigData = {}
    
    for ParamName in ParamNames:
        try:
            for line in lines:
                    ThisLine=line.strip().split()
                    ThisLineLength=len(ThisLine)
                    if ThisLineLength < 3:
                        continue
                    if list(ThisLine[0])[0]=='#' or ThisLine[0]=='\n':
                        continue
                    if ThisLine[0] !=  ParamName:
                        continue
                    else:
                        try:
                            ThisLine.remove(ThisLine[0])
                            for counter,item in enumerate(ThisLine):
                                if list(item)[0] == '#':
                                    del ThisLine[counter:] # Strip the rest of the line as a comment
                                    continue
                                if item == '=':
                                    ThisLine.remove(item)
                                    continue
                            if len(ThisLine) == 0:
                                continue
                            elif len(ThisLine) == 1:
                                ThisParam = ThisLine[0]
                                try: ConfigData[ParamName] = int(ThisParam)
                                except ValueError:
                                    try: ConfigData[ParamName] = float(ThisParam)
                                    except ValueError:
                                        ConfigData[ParamName] = ThisParam
                            else:
                                ThisParam = []
                                for item in ThisLine:
                                    try: ThisParam.append(int(item))
                                    except ValueError:
                                        try: ThisParam.append(float(item))
                                        except ValueError:
                                            ThisParam.append(item)
                                ConfigData[ParamName] = ThisParam

                        except:
                            continue
        except:
            print "Error reading .cfg file"

    return ConfigData


#****************MAIN PROGRAM*****************

# First, read the .cfg file

configfile = sys.argv[1]
ConfigData = ReadConfigFile(configfile)
outputfilebase = ConfigData["outputfilebase"]
outputfiledir = ConfigData["outputfiledir"]

dat = Array3dHDF5(outputfiledir, outputfilebase, ConfigData["LogEField"]) 
# This holds all of the data
ScaleFactor = ConfigData["ScaleFactor"]
GridsPerPixel = ConfigData["GridsPerPixel"]
nxx = dat.nx - 1
nyy = dat.ny - 1
nzz = dat.nz - 1

# A couple of things to customize the plots
PlotEField = bool(ConfigData["PlotEField"])
EdgePlot = bool(ConfigData["EdgePlot"])

if EdgePlot:
    NumPixelsPlotted = 8 # Special for edge plots
    nxcenter = nxx/2 + NumPixelsPlotted/2 * ScaleFactor + ScaleFactor * 6 # Special for edge plots
    nycenter = nyy/2 + NumPixelsPlotted/2 * ScaleFactor - ScaleFactor * 2 # Special for edge plots
else:
    NumPixelsPlotted = 4
    nxcenter = nxx/2 + NumPixelsPlotted/2 * ScaleFactor
    nycenter = nyy/2 + NumPixelsPlotted/2 * ScaleFactor

nxmin = nxcenter - NumPixelsPlotted/2 * ScaleFactor * GridsPerPixel
nxmax = nxcenter + NumPixelsPlotted/2 * ScaleFactor * GridsPerPixel
nymin = nycenter - NumPixelsPlotted/2 * ScaleFactor * GridsPerPixel
nymax = nycenter + NumPixelsPlotted/2 * ScaleFactor * GridsPerPixel
nzmin = 0
nzmax = 8 * ScaleFactor


rcParams['contour.negative_linestyle'] = 'solid'
rcParams.update({'font.size': 6})

print "Making array edge potential plots\n"
figure()
suptitle("Array Edge Potentials. Grid = %d*%d*%d."%(nxx,nyy,nzz),fontsize = 18)
subplot(2,2,1)
title("Front Edge")
ylim(-20.0, 10.0)
plot(dat.x[:],dat.phi[:,0,0])
plot(dat.x[:],dat.phi[:,0,1])
plot(dat.x[:],dat.phi[:,0,2])
plot(dat.x[:],dat.phi[:,0,3])
plot(dat.x[:],dat.phi[:,0,10])

subplot(2,2,2)
title("Back Edge")
ylim(-20.0, 10.0)
plot(dat.x[:],dat.phi[:,dat.ny-1,0])
plot(dat.x[:],dat.phi[:,dat.ny-1,1])
plot(dat.x[:],dat.phi[:,dat.ny-1,2])
plot(dat.x[:],dat.phi[:,dat.ny-1,3])
plot(dat.x[:],dat.phi[:,dat.ny-1,10])

subplot(2,2,3)
title("Left Edge")
ylim(-20.0, 10.0)
plot(dat.y[:],dat.phi[0,:,0])
plot(dat.y[:],dat.phi[0,:,1])
plot(dat.y[:],dat.phi[0,:,2])
plot(dat.y[:],dat.phi[0,:,3])
plot(dat.y[:],dat.phi[0,:,10])

subplot(2,2,4)
title("Right Edge")
ylim(-20.0, 10.0)
plot(dat.y[:],dat.phi[dat.nx-1,:,0])
plot(dat.y[:],dat.phi[dat.nx-1,:,1])
plot(dat.y[:],dat.phi[dat.nx-1,:,2])
plot(dat.y[:],dat.phi[dat.nx-1,:,3])
plot(dat.y[:],dat.phi[dat.nx-1,:,10])

savefig("plots/"+outputfilebase+"_Edge_Potentials.pdf")

print "Making 1D potential plots\n"
figure()

suptitle("1D Potential Slices. Grid = %d*%d*%d."%(nxx,nyy,nzz),fontsize = 18)
plotcounter = 1
subplots_adjust(hspace=0.3, wspace=0.3)

nxcenter2 = nxcenter
nycenter2 = nycenter
subplot(2,3,1)
title("Phi-Collect Gate, x = %.2f, y = %.2f"%(((dat.x[nxcenter2]+dat.x[nxcenter2-1])/2.0),((dat.y[nycenter2]+dat.y[nycenter2-1])/2.0)))
plot(dat.z[0:20],(dat.phi[nxcenter2,nycenter2,0:20]+dat.phi[nxcenter2-1,nycenter2,0:20]+dat.phi[nxcenter2,nycenter2-1,0:20]+dat.phi[nxcenter2-1,nycenter2-1,0:20])/4.0)

xlabel("Z-Dimension (microns)")
ylim(-20.0, 10.0)

subplot(2,3,4)
title("Rho-Collect Gate, x = %.2f, y = %.2f"%(((dat.x[nxcenter2]+dat.x[nxcenter2-1])/2.0),((dat.y[nycenter2]+dat.y[nycenter2-1])/2.0)))
plot(dat.z[0:20],(dat.rho[nxcenter2,nycenter2,0:20]+dat.rho[nxcenter2-1,nycenter2,0:20]+dat.rho[nxcenter2,nycenter2-1,0:20]+dat.rho[nxcenter2-1,nycenter2-1,0:20])/4.0)
xlabel("Z-Dimension (microns)")
ylim(-60.0, 20.0)

nxcenter2 = nxcenter
nycenter2 = nycenter + GridsPerPixel * ScaleFactor / 2
subplot(2,3,2)
title("Phi-ChanStop, x = %.2f, y = %.2f"%(((dat.x[nxcenter2]+dat.x[nxcenter2-1])/2.0),((dat.y[nycenter2]+dat.y[nycenter2-1])/2.0)))
plot(dat.z[0:20],(dat.phi[nxcenter2,nycenter2,0:20]+dat.phi[nxcenter2-1,nycenter2,0:20]+dat.phi[nxcenter2,nycenter2-1,0:20]+dat.phi[nxcenter2-1,nycenter2-1,0:20])/4.0)
xlabel("Z-Dimension (microns)")
ylim(-20.0, 10.0)

subplot(2,3,5)
title("Rho-ChanStop, x = %.2f, y = %.2f"%(((dat.x[nxcenter2]+dat.x[nxcenter2-1])/2.0),((dat.y[nycenter2]+dat.y[nycenter2-1])/2.0)))
plot(dat.z[0:20],(dat.rho[nxcenter2,nycenter2,0:20]+dat.rho[nxcenter2-1,nycenter2,0:20]+dat.rho[nxcenter2,nycenter2-1,0:20]+dat.rho[nxcenter2-1,nycenter2-1,0:20])/4.0)
xlabel("Z-Dimension (microns)")
ylim(-60.0, 20.0)

nxcenter2 = nxcenter + GridsPerPixel * ScaleFactor / 2
nycenter2 = nycenter
subplot(2,3,3)
title("Phi-Barrier Gate, x = %.2f, y = %.2f"%(((dat.x[nxcenter2]+dat.x[nxcenter2-1])/2.0),((dat.y[nycenter2]+dat.y[nycenter2-1])/2.0)))
plot(dat.z[0:20],(dat.phi[nxcenter2,nycenter2,0:20]+dat.phi[nxcenter2-1,nycenter2,0:20]+dat.phi[nxcenter2,nycenter2-1,0:20]+dat.phi[nxcenter2-1,nycenter2-1,0:20])/4.0)
xlabel("Z-Dimension (microns)")
ylim(-20.0, 10.0)

subplot(2,3,6)
title("Rho-Barrier Gate, x = %.2f, y = %.2f"%(((dat.x[nxcenter2]+dat.x[nxcenter2-1])/2.0),((dat.y[nycenter2]+dat.y[nycenter2-1])/2.0)))
plot(dat.z[0:20],(dat.rho[nxcenter2,nycenter2,0:20]+dat.rho[nxcenter2-1,nycenter2,0:20]+dat.rho[nxcenter2,nycenter2-1,0:20]+dat.rho[nxcenter2-1,nycenter2-1,0:20])/4.0)
xlabel("Z-Dimension (microns)")
ylim(-60.0, 20.0)
savefig("plots/"+outputfilebase+"_1D_Potentials.pdf")

print "Making summary plots\n"
figure()

suptitle("CCD Charge Collection. Grid = %d*%d*%d."%(nxx,nyy,nzz),fontsize = 18)
plotcounter = 1
subplots_adjust(hspace=0.3, wspace=0.3)
[yy,xx] = meshgrid(dat.y[nymin:nymax],dat.x[nxmin:nxmax])

slicez = 0
subplot(2,3,1, aspect = 1)
title("Phi, z = 0.0")
levels = linspace(-10.0, 10.0, 51)
contourf(xx,yy,dat.phi[nxmin:nxmax,nymin:nymax,slicez],levels)
xlabel("X-Dimension (microns)")
ylabel("Y-Dimension (microns)")
colorbar()

slicez = 1 * ScaleFactor
subplot(2,3,2, aspect = 1)
title("Rho, z = %.2f"%dat.z[slicez])
contourf(xx,yy,(dat.rho[nxmin:nxmax,nymin:nymax,slicez]+dat.rho[nxmin:nxmax,nymin:nymax,slicez+1])/2.0)
xlabel("X-Dimension (microns)")
ylabel("Y-Dimension (microns)")
colorbar()

[yy,xx] = meshgrid(dat.y[nymin:nymax],dat.x[nxmin:nxmax])
slicez = 1 * ScaleFactor
subplot(2,3,3, aspect = 1)
title("Phi, z = %.2f"%dat.z[slicez])
levels = linspace(-30.0, 10.0, 81)
contourf(xx,yy,dat.phi[nxmin:nxmax,nymin:nymax,slicez],levels)
xlabel("X-Dimension (microns)")
ylabel("Y-Dimension (microns)")
plot([dat.x[nxmin+1],dat.x[nxmax-1]],[((dat.y[nycenter]+dat.y[nycenter-1])/2.0),((dat.y[nycenter]+dat.y[nycenter-1])/2.0)],ls = "-", color="k")
plot([((dat.x[nxcenter]+dat.x[nxcenter-1])/2.0),((dat.x[nxcenter]+dat.x[nxcenter-1])/2.0)],[dat.y[nymin+1],dat.y[nymax-1]],ls = "-", color="k")
colorbar()

slicez = 4 * ScaleFactor
subplot(2,3,4, aspect = 1)
title("Phi, z = %.2f"%dat.z[slicez])
contourf(xx,yy,dat.phi[nxmin:nxmax,nymin:nymax,slicez],levels)
xlabel("X-Dimension (microns)")
ylabel("Y-Dimension (microns)")
colorbar()

subplot(2,3,5)
title("Phi and (-)E in Gate Region. y = %.2f"%((dat.y[nycenter]+dat.y[nycenter-1])/2.0))
xlabel("X-Dimension (microns)")
ylabel("Z-Dimension (microns)")
nzmin = 0
[zz,xx] = meshgrid(dat.z[nzmin:nzmax],dat.x[nxmin:nxmax])
contourf(xx,zz,(dat.phi[nxmin:nxmax,nycenter,nzmin:nzmax]+dat.phi[nxmin:nxmax,nycenter-1,nzmin:nzmax])/2.0,levels)
colorbar()
if ConfigData["LogEField"] == 1 and PlotEField:
    nzmin = 1
    [zz,xx] = meshgrid(dat.z[nzmin:nzmax],dat.x[nxmin:nxmax])
    quiver(xx, zz, (dat.Ex[nxmin:nxmax,nycenter,nzmin:nzmax]+dat.Ex[nxmin:nxmax,nycenter-1,nzmin:nzmax])/2.0, (dat.Ez[nxmin:nxmax,nycenter,nzmin:nzmax]+dat.Ez[nxmin:nxmax,nycenter-1,nzmin:nzmax])/2.0, color='b', scale = 150.0)#, scale_units="width")
nzmin = 0
[zz,xx] = meshgrid(dat.z[nzmin:nzmax],dat.x[nxmin:nxmax])
ylim(zz[0,0], zz[-1,-1])
xlim(xx[0,0], xx[-1,-1])


subplot(2,3,6)
title("Phi and (-)E in Gate Region. x = %.2f"%((dat.x[nxcenter]+dat.x[nxcenter-1])/2.0))
xlabel("Y-Dimension (microns)")
ylabel("Z-Dimension (microns)")
nzmin = 0
[zz,yy] = meshgrid(dat.z[nzmin:nzmax],dat.y[nymin:nymax])
contourf(yy,zz,(dat.phi[nxcenter,nymin:nymax,nzmin:nzmax]+dat.phi[nxcenter-1,nymin:nymax,nzmin:nzmax])/2.0,levels)
colorbar()
if ConfigData["LogEField"] == 1 and PlotEField:
    nzmin = 1
    [zz,yy] = meshgrid(dat.z[nzmin:nzmax],dat.y[nymin:nymax])
    quiver(yy, zz, (dat.Ey[nxcenter,nymin:nymax,nzmin:nzmax]+dat.Ey[nxcenter-1,nymin:nymax,nzmin:nzmax])/2.0, (dat.Ez[nxcenter,nymin:nymax,nzmin:nzmax]+dat.Ez[nxcenter-1,nymin:nymax,nzmin:nzmax])/2.0, color='b', scale = 150.0)#, scale_units="width")
nzmin = 0
[zz,yy] = meshgrid(dat.z[nzmin:nzmax],dat.y[nymin:nymax])
ylim(zz[0,0], zz[-1,-1])
xlim(yy[0,0], yy[-1,-1])

savefig("plots/"+outputfilebase+"_Summary.pdf")

if ConfigData["LogPixels"] == 1:
    # Next, plots of the pixel boundaries
    print "Making pixel plots\n"
    figure()
    rcParams['contour.negative_linestyle'] = 'solid'
    #rcParams.update({'font.size': 18})

    suptitle("CCD Pixel Plots. Grid = %d*%d*%d."%(nxx,nyy,nzz),fontsize = 24)
    plotcounter = 1
    subplots_adjust(hspace=0.3, wspace=0.1)

    filename = outputfiledir+"/"+outputfilebase+"_Pts"
    file = open(filename,"r")
    lines = file.readlines()
    file.close()

    redsx=[]
    redsy=[]
    blacksx=[]
    blacksy=[]
    plottedxin = -1000.0
    plottedyin = -1000.0
    for line in lines:
        if line.split()[0] == "xin":
            continue

        values = line.split()
        zout = float(values[5])    
        if zout > dat.z[2]:
            continue # Skip these in case of LogTracePaths = 1
        xin = float(values[0])
        yin = float(values[1])
        xout = float(values[3])
        yout = float(values[4])
        pixxout = int(xout/10.0)
        pixyout = int(yout/10.0)    

        if (xin > plottedxin - .001) and (xin < plottedxin + .001) and \
           (yin > plottedyin - .001) and (yin < plottedyin + .001):
            continue # Skips it if it is already plotted
        if (pixxout + pixyout) % 2 == 0:
            redsx.append(xin)
            redsy.append(yin)
        else:
            blacksx.append(xin)
            blacksy.append(yin)

        plottedxin = xin
        plottedyin = yin

    subplot(1,1,1,aspect=1)
    title("Pixel Boundaries",fontsize = 12)
    scatter(redsx,redsy,s=5,color="red")
    scatter(blacksx,blacksy,s=5,color="black")

    if EdgePlot:
        for linex in linspace(120.0,200.0,9):
            plot((linex,linex),(20.0,70.0),linewidth=1.0, color='blue')

    xlabel("X(microns)",fontsize = 18)
    ylabel("Y(microns)",fontsize = 18)
    xlim(ConfigData["PixelBoundaryLowerLeft"][0], ConfigData["PixelBoundaryUpperRight"][0])
    ylim(ConfigData["PixelBoundaryLowerLeft"][1], ConfigData["PixelBoundaryUpperRight"][1])


    savefig("plots/"+outputfilebase+"_Pixels.pdf")

if ConfigData["LogPixelPaths"] == 1 and ConfigData["LogPixels"] == 1:
    # Last, plots of the electron paths
    print "Making array electron path plots\n"

    yline = (ConfigData["PixelBoundaryLowerLeft"][1] + ConfigData["PixelBoundaryUpperRight"][1] + ConfigData["PixelBoundaryStepSize"][1]) / 2.0
    xline = (ConfigData["PixelBoundaryLowerLeft"][0] + ConfigData["PixelBoundaryUpperRight"][0] + ConfigData["PixelBoundaryStepSize"][0]) / 2.0

    figure()
    vertical_zoom = 4
    suptitle("Electron Path Plot - Vertical Zoom = %d"%vertical_zoom, fontsize = 24)
    subplots_adjust(wspace=0.2)
    subplot(1,2,1,aspect=vertical_zoom)
    oldxin = 1000000.0
    for line in lines:
        if line.split()[0] == "xin":
            continue
        values = line.split()
        xin = float(values[0])
        yin = float(values[1])
        if (yin < yline - .001) or (yin > yline + .001):
            continue

        pixxin = int(xin/10.0)
        if pixxin % 2 == 0:
            color = "red"
        else:
            color = "black"
        if abs(xin - oldxin) > .001:
            if oldxin < 100000.0:
                plot(xpaths, zxpaths, color = color, linewidth = 0.1)
            oldxin = xin
            xpaths=[]
            zxpaths=[]
            
        xout = float(values[3])
        yout = float(values[4])
        zout = float(values[5])
        xpaths.append(xout)
        zxpaths.append(zout)
    ylabel("Z(microns")
    xlabel("X (microns)")
    ylim(0.0,10.0)
    xlim(ConfigData["PixelBoundaryLowerLeft"][0], ConfigData["PixelBoundaryUpperRight"][0])

    subplot(1,2,2,aspect=vertical_zoom)
    oldyin = 1000000.0
    for line in lines:
        if line.split()[0] == "xin":
            continue
        values = line.split()
        xin = float(values[0])
        yin = float(values[1])
        if (xin < xline - .001) or (xin > xline + .001):
            continue

        pixyin = int(yin/10.0)
        if pixyin % 2 == 0:
            color = "red"
        else:
            color = "black"
        if abs(yin - oldyin) > .001:
            if oldyin < 100000.0:
                plot(ypaths, zypaths, color = color, linewidth = 0.1)
            oldyin = yin
            ypaths=[]
            zypaths=[]
            
        xout = float(values[3])
        yout = float(values[4])
        zout = float(values[5])
        ypaths.append(yout)
        zypaths.append(zout)
    ylabel("Z(microns")
    xlabel("Y (microns)")
    ylim(0.0,10.0)
    xlim(ConfigData["PixelBoundaryLowerLeft"][1], ConfigData["PixelBoundaryUpperRight"][1])
    savefig("plots/"+outputfilebase+"_Paths.pdf")
