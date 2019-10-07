# in ipython notebook, enable inline plotting with:


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from matplotlib.font_manager import FontProperties
font = FontProperties()
font.set_weight('bold')
font.set_size('medium')
# create some circles
#def draw_circle(x,y):
#    if ((x>-1.5 and x<1.5) and (y>-1.5 and y<1.5)):
#        if (x>0 and y>0):
#            circle1 = plt.Circle((x,y), 0.3, color='r', alpha=.4)
#        if (x>0 and y<0):
#            circle1 = plt.Circle((x,y), 0.3, color='b', alpha=.4)
#        if (x<0 and y<0):
#            circle1 = plt.Circle((x,y), 0.3, color='g', alpha=.4)
#        if (x<0 and y>0):
#            circle1 = plt.Circle((x,y), 0.3, color='y', alpha=.4)
#        ax = plt.gca()
#        plt.axhline(y=0, color='r', linestyle='-')
#        plt.axvline(x=0, color='r', linestyle='-')
#        ax.add_artist(circle1)
#        ax.set_xlim(-2.5, 2.5); ax.set_ylim(-2.5, 2.5)
#        ax.set_aspect('equal')
#
#for i in range (80):
##    x = np.random.uniform(0,2)
##    y = np.random.uniform(0,2)
#    x = np.random.normal(0,2)
#    y = np.random.normal(0,2)
#    draw_circle(x,y)
#
#plt.plot()
#plt.show()





NUM = 50
w =1
h =2

ells = []
fig, ax = plt.subplots(subplot_kw={'aspect': 'equal'})   
plt.axhline(y=0, color='r', linestyle='-')
plt.axvline(x=0, color='r', linestyle='-')  
ax.set_xlim(-3, 3); ax.set_ylim(-3, 3)
       
for i in range(15):
    x = np.random.uniform(-2,2)
    y = np.random.uniform(-2,2)
    if (x>0 and y>0):
        ells = Ellipse(xy=(x,y),
            width=w, height=h,
            angle=np.random.rand() * 360, gid = 'ID')
               
        ax.add_artist(ells)
        ax.annotate(i,xy=(x,y),fontproperties=font)
        ells.set_clip_box(ax.bbox)
        ells.set_alpha(0.5)
        ells.set_facecolor('b')
    if (x>0 and y<0):
        ells = Ellipse(xy=(x,y),
            width=w, height=h,
            angle=np.random.rand() * 360)
               
        ax.add_artist(ells)
        ax.annotate(i,xy=(x,y),fontproperties=font)
        ells.set_clip_box(ax.bbox)
        ells.set_alpha(0.5)
        ells.set_facecolor('r')
    if (x<0 and y>0):
        ells = Ellipse(xy=(x,y),
            width=w, height=h,
            angle=np.random.rand() * 360)
               
        ax.add_artist(ells)
        ax.annotate(i,xy=(x,y),fontproperties=font)
        ells.set_clip_box(ax.bbox)
        ells.set_alpha(0.5)
        ells.set_facecolor('y')
    if (x<0 and y<0):
        ells = Ellipse(xy=(x,y),
            width=w, height=h,
            angle=np.random.rand() * 360)
               
        ax.add_artist(ells)
        ax.annotate(i,xy=(x,y),fontproperties=font)
        ells.set_clip_box(ax.bbox)
        ells.set_alpha(0.5)
        ells.set_facecolor('g')
plt.show()
