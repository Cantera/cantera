
from Tkinter import *
import math

class Graph(Frame):
    def __init__(self,master,title,minX,maxX,minY,maxY,pixelX=250,pixelY=250):
        Frame.__init__(self,master, relief=RIDGE, bd=2)
        #       self.pack()
        self.title = Label(self,text=' ')
        self.title.grid(row=0,column=1,sticky=W+E)
        self.graph_w, self.graph_h = pixelX, pixelY
        self.maxX, self.maxY = maxX, maxY #float(math.floor(maxX + 1)), \
                               #float(math.floor(maxY + 1))
        self.minX, self.minY = minX, minY # float(math.floor(minX)), float(math.floor(minY))
        self.canvas = Canvas(self,
                             width=self.graph_w,
                             height=self.graph_h,
                             relief=SUNKEN,bd=1)
        ymintext = "%8.1f" % (self.minY)
        ymaxtext = "%8.1f" % (self.maxY)
        self.ml=Label(self, text=ymintext)
        self.mr=Label(self, text=ymaxtext)
        self.ml.grid(row=2,column=0,sticky=S+E)
        self.mr.grid(row=1,column=0,sticky=N+E)
        self.canvas.grid(row=1,column=1,rowspan=2,sticky=N+S+E+W)
        self.last_points=[]


    def writeValue(self, y):
        yval = '%15.4f' % (y)
        self.title.config(text = yval)

    def delete(self, ids):
        for id in ids:
            self.canvas.delete(id)

    def move(self, id, newpos, oldpos):
        dxpt = (newpos[0] - oldpos[0])/(self.maxX-self.minX)*self.graph_w
        dypt = -(newpos[1] - oldpos[1])/(self.maxY-self.minY)*self.graph_h
        self.canvas.move(id, dxpt, dypt)
        self.writeValue(newpos[1])

    def plot(self,x,y,color='black'):
        xpt = (x-self.minX)/(self.maxX-self.minX)*float(self.graph_w) + 1.5
        ypt = (self.maxY-y)/(self.maxY-self.minY)*float(self.graph_h) - 1.5
        id_ycross = self.canvas.create_line(xpt,self.graph_h,xpt,0,fill = 'gray')
        id_xcross = self.canvas.create_line(0,ypt,self.graph_w,ypt,fill = 'gray')
        id = self.canvas.create_oval(xpt-2,ypt-2,xpt+2,ypt+2,fill=color)
        self.writeValue(y)
        return [id,id_xcross,id_ycross]

    def reset(self,minX,maxX,minY,maxY):
        self.maxX, self.maxY = maxX, maxY
        self.minX, self.minY = minX, minY
        self.canvas.destroy()
        self.canvas = Canvas(self,
                             width=self.graph_w,
                             height=self.graph_h,
                             relief=SUNKEN,bd=1)
        self.canvas.create_text(4,2,text=self.maxY,anchor=NW)
        self.canvas.create_text(4,self.graph_h,text=self.minY,anchor=SW)
        self.ml["text"] = `minX`
        self.mr["text"] = `maxX`
        self.canvas.pack()
        self.last_points = []

    def join(self,point_list):
        i = 0
        for pt in point_list:
            x, y, color = pt
            if self.last_points == []:
                last_x, last_y, last_color = pt
            else:
                last_x, last_y, last_color = self.last_points[i]
            i = i + 1
            xpt = (x - self.minX)/(float(self.maxX - self.minX)/self.graph_w) + 1.5
            ypt = (self.maxY-y)/(float(self.maxY - self.minY)/self.graph_h) - 1.5
            last_xpt = (last_x - self.minX)/(float(self.maxX - self.minX)/self.graph_w) + 1.5
            last_ypt = (self.maxY-last_y)/(float(self.maxY - self.minY)/self.graph_h) - 1.5
            self.canvas.create_line(last_xpt,last_ypt,
                            xpt,ypt,fill=color)
        self.last_points = point_list
        self.canvas.update()
        return


    def addLegend(self,text,color=None):
        m=Message(self,text=text,width=self.graph_w-10)
        m.pack(side=BOTTOM)
        if color:
            m.config(fg=color)

    def pauseWhenFinished(self):
        self.wait_window()


if __name__=='__main__':
    root= Tk()
    g = Graph(root,'graph1',0,10,0.01,120)
    h = Graph(root,'graph2',0,15,0,20000)
    g.pack(side=LEFT)
    h.pack(side=RIGHT)

    #root.protocol("WM_DELETE_WINDOW", root.destroy())
    j = Graph(root,'Graph',0,1000,0,2000)
    j.pack()

    j.plot(0, 0, color='red')
    j.last_points = [ (0, 0, 'red') ]
    for i in range(100):
        j.join( [ ( (i*10),(i*10+500), 'red' ) ] )


    g.addLegend('An example of the GraphFrame')
    h.addLegend('This is where the legend goes')
    for i in range(100):
        if root:
            x,y = float(i)/10, i
            g.plot(x,y,color='red')
            h.plot(i,i**2)#(0,0)
            #h.join([(i,i**2,'black')])
        else:
            break

    #print("finished")
    g.pauseWhenFinished()
    h.pauseWhenFinished()
    print g
