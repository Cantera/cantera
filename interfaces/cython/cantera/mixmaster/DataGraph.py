import sys

if sys.version_info[0] == 3:
    from tkinter import *
else:
    from Tkinter import *

import math
from numpy import *

def plotLimits(ypts, f=0.0, ndiv=5, logscale=0):
    """Return plot limits that"""
    if logscale:
        threshold = 1.0e-19
    else:
        threshold = -1.0e20
    ymax = -1.e20
    ymin = 1.e20
    for y in ypts:
        if y > ymax: ymax = y
        if y < ymin and y > threshold: ymin = y

    dy = abs(ymax - ymin)

    if logscale:
        ymin = math.floor(math.log10(ymin))
        ymax = math.floor(math.log10(ymax))+1
        fctr = 1.0

##         if dy < 0.2*ymin:
##             ymin = ymin*.9
##             ymax = ymax*1.1
##             dy = abs(ymax - ymin)
##         else:
    else:
        ymin = ymin - f*dy
        ymax = ymax + f*dy
        dy = abs(ymax - ymin)

        try:
            p10 = math.floor(math.log10(0.1*dy))
            fctr = math.pow(10.0, p10)
        except:
            return (ymin -1.0, ymax + 1.0, 1.0)
        mm = [2.0, 2.5, 2.0]
        i = 0
        while dy/fctr > ndiv:
            fctr = mm[i % 3]*fctr
            i = i + 1
        ymin = fctr*math.floor(ymin/fctr)
        ymax = fctr*(math.floor(ymax/fctr+0.999))

    return (ymin, ymax, fctr)


class DataGraph(Frame):
    def __init__(self,master,
                 data, ix=0, iy=0,
                 title='',
                 label = ('x-axis','y-axis'),
                 logscale = (0,0),
                 pixelX=500,
                 pixelY=500):
        self.logscale = logscale
        self.data = data
        self.ix = ix
        self.iy = iy
        self.minX, self.maxX, self.dx = plotLimits(data[ix,:],
                                                   logscale=self.logscale[0])
        self.minY, self.maxY, self.dy = plotLimits(data[iy,:],
                                                   logscale=self.logscale[1])

        Frame.__init__(self,master, relief=RIDGE, bd=2)
        self.title = Label(self,text=' ')
        self.title.grid(row=0,column=1,sticky=W+E)
        self.graph_w, self.graph_h = pixelX - 120, pixelY - 70
        self.origin = (100, 20)
        self.canvas = Canvas(self,
                             width=pixelX,
                             height=pixelY,
                             relief=SUNKEN,bd=1)
        id = self.canvas.create_rectangle(self.origin[0],self.origin[1],
                                          pixelX-20,pixelY-50)
        self.canvas.grid(row=1,column=1,rowspan=2,sticky=N+S+E+W)
        self.last_points=[]
        self.ticks(self.minX, self.maxX, self.dx,
                   self.minY, self.maxY, self.dy, 10)
        self.screendata()
        self.draw()
        self.canvas.create_text(self.origin[0] + self.graph_w/2,
                                self.origin[1] + self.graph_h + 30,
                                text=label[0],anchor=N)
        self.canvas.create_text(self.origin[0] - 50,
                                self.origin[1] + self.graph_h/2,
                                text=label[1],anchor=E)


    def writeValue(self, y):
        yval = '%15.4f' % (y)
        self.title.config(text = yval)

    def delete(self, ids):
        for id in ids:
            self.canvas.delete(id)

    def screendata(self):
        self.xdata = array(self.data[self.ix,:])
        self.ydata = array(self.data[self.iy,:])
        npts = len(self.ydata)
        if self.logscale[0] > 0:
            self.xdata = log10(self.xdata)
        if self.logscale[1] > 0:
            self.ydata = log10(self.ydata)
        f = float(self.graph_w)/(self.maxX-self.minX)
        self.xdata = (self.xdata - self.minX)*f + self.origin[0]
        f = float(self.graph_h)/(self.maxY-self.minY)
        self.ydata = (self.maxY - self.ydata)*f + self.origin[1]

    def toscreen(self,x,y):
        if self.logscale[0] > 0:
            x = log10(x)
        if self.logscale[1] > 0:
            y = log10(y)
        f = float(self.graph_w)/(self.maxX-self.minX)
        xx = (x - self.minX)*f + self.origin[0]
        f = float(self.graph_h)/(self.maxY-self.minY)
        yy = (self.maxY - y)*f + self.origin[1]
        return (xx, yy)

    def move(self, id, newpos, oldpos):
        dxpt = (newpos[0] - oldpos[0])/(self.maxX-self.minX)*self.graph_w
        dypt = -(newpos[1] - oldpos[1])/(self.maxY-self.minY)*self.graph_h
        self.canvas.move(id, dxpt, dypt)
        self.writeValue(newpos[1])

    def plot(self,n,color='black'):
        xpt, ypt = self.toscreen(self.data[self.ix,n],
                                 self.data[self.iy,n])
        #xpt = (x-self.minX)/(self.maxX-self.minX)*float(self.graph_w) + self.origin[0]
        #ypt = (self.maxY-y)/(self.maxY-self.minY)*float(self.graph_h) + self.origin[1]
        id_ycross = self.canvas.create_line(xpt,self.graph_h+self.origin[1],xpt,self.origin[1],fill = 'gray')
        id_xcross = self.canvas.create_line(self.origin[0],ypt,self.graph_w+self.origin[0],ypt,fill = 'gray')
        id = self.canvas.create_oval(xpt-2,ypt-2,xpt+2,ypt+2,fill=color)
        #self.writeValue(y)
        s = '(%g, %g)' % (self.data[self.ix,n],self.data[self.iy,n])
        if n > 0 and self.data[self.iy,n] > self.data[self.iy,n-1]:
            idt = self.canvas.create_text(xpt+5,ypt+5,text=s,anchor=NW)
        else:
            idt = self.canvas.create_text(xpt+5,ypt-5,text=s,anchor=SW)

        return [id,id_xcross,id_ycross, idt]

    def draw(self,color='red'):
        npts = len(self.xdata)
        for n in range(1,npts):
            self.canvas.create_line(self.xdata[n-1],self.ydata[n-1],
                            self.xdata[n],self.ydata[n],fill=color)

    def addLabel(self, y, orient=0):
        if orient==0:
            xpt, ypt = self.toscreen(y, 1.0)
            ypt = self.origin[1] + self.graph_h + 5
            self.canvas.create_text(xpt,ypt,text=y,anchor=N)
        else:
            xpt, ypt = self.toscreen(self.minX, y)
            xpt = self.origin[0] - 5
            self.canvas.create_text(xpt,ypt,text=y,anchor=E)
    def addLegend(self,text,color=None):
        m=Message(self,text=text,width=self.graph_w-10)
        m.pack(side=BOTTOM)
        if color:
            m.config(fg=color)

    def pauseWhenFinished(self):
        self.wait_window()

    def minorTicks(self, x0, x1, y, n, size, orient=0):
        xtick = x0
        dx = (x1 - x0)/float(n)
        if orient == 0:
            while xtick <= x1:
                xx, yy = self.toscreen(xtick, y)
                self.canvas.create_line(xx,yy,
                                        xx,yy-size)
                xtick += dx
        else:
            while xtick <= x1:
                xx, yy = self.toscreen(y, xtick)
                self.canvas.create_line(xx,yy,
                                        xx+size,yy)
                xtick += dx


    def ticks(self, xmin, xmax, dx, ymin, ymax, dy, size):

        if self.logscale[0]:
            xmin = math.pow(10.0,xmin)
            xmax = math.pow(10.0,xmax)
        if self.logscale[1]:
            ymin = math.pow(10.0,ymin)
            ymax = math.pow(10.0,ymax)

        n = 5
        ytick = ymin
        while ytick <= ymax:
            xx, yy = self.toscreen(xmin, ytick)
            self.canvas.create_line(xx, yy, xx + size,yy)
            self.addLabel(ytick,1)
            xx, yy = self.toscreen(xmax, ytick)
            self.canvas.create_line(xx, yy, xx - size,yy)
            ytick0 = ytick
            if self.logscale[1]:
                ytick *= 10.0
                n = 10
            else: ytick = ytick + dy
            if ytick <= ymax:
                self.minorTicks(ytick0, ytick, xmin, n, 5, 1)
                self.minorTicks(ytick0, ytick, xmax, n, -5, 1)

        n = 5
        xtick = xmin
        while xtick <= xmax:
            xx, yy = self.toscreen(xtick, ymin)
            self.canvas.create_line(xx, yy, xx, yy - size)
            self.addLabel(xtick,0)
            xx, yy = self.toscreen(xtick, ymax)
            self.canvas.create_line(xx, yy, xx, yy + size)
            if self.logscale[0]:
                xtick *= 10.0
                n = 10
            else: xtick = xtick + dx
            if xtick <= xmax:
                self.minorTicks(xtick - dx, xtick, ymin, n, 5, 0)
                self.minorTicks(xtick - dx, xtick, ymax, n, -5, 0)
