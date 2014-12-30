# -*- coding: utf-8 -*-

from mpl_toolkits.mplot3d import proj3d
import matplotlib.patches as mpatches
from constant import *

def make_3d_plot_with_speed(ax, position, time, steps, n, plot_type='SPEED'):
    ax.set_xlim3d([0, n[0]])
    ax.set_ylim3d([0, n[1]])
    ax.set_zlim3d([0, n[2]])
    xStepGrid, yStepGrid, zStepGrid = steps
    gridDimension = (int(n[0]/xStepGrid), int(n[1]/yStepGrid), int(n[2]/zStepGrid))
    n = position.shape
    labels = [[[False for i in range(gridDimension[0])] for j in range(gridDimension[1])] for k in range(gridDimension[2])]
    flags = labels
    for num in range(n[1]):
        xBig, yBig, zBig = position[time][num][0], position[time][num][1], position[time][num][2]
        # print('{}: {} {} {}'.format(num, xBig, yBig, zBig))
        # print('{} {} {}'.format(xStepGrid, yStepGrid, zStepGrid))
        # print('{} {} {}'.format(xBig/xStepGrid, yBig/yStepGrid, zBig/zStepGrid))
        i, j, k = int(xBig/xStepGrid), int(yBig/yStepGrid), int(zBig/zStepGrid)
        x, y, z = i*xStepGrid, j*yStepGrid, k*zStepGrid
        ax.scatter(xBig, yBig, zBig, color='red')
        if (j == k and not flags[0][j][k]):
            flags[0][j][k] = True
            x2, y2, _ = proj3d.proj_transform(xBig, yBig, zBig, ax.get_proj())
            str_label = ''
            if plot_type=='SPEED':
                str_label = '{:5.7f},{:5.7f},{:5.7f}'.format(position[time][num][5]*Mnuc, position[time][num][6]*Mnuc, position[time][num][7]*Mnuc)
            if plot_type=='COORD':
                str_label = '{:5.7f},{:5.7f},{:5.7f}'.format(position[time][num][0], position[time][num][1], position[time][num][2])
            labels[i][j][k] = ax.annotate (str_label, xy=(x2, y2), xytext=(-15, 15),
                textcoords='offset points', ha='right', va='bottom',
                bbox=dict(boxstyle='round,pad=0.5', fc='yellow', alpha=0.5),
                arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0'))
    return ax

def make_2d_plot_with_tracks(ax, position, time, steps, n):
    ax.set_xlim([0, n[0]])
    ax.set_ylim([0, n[1]])
    xStepGrid, yStepGrid, zStepGrid = steps
    gridDimension = (int(n[0]/xStepGrid), int(n[1]/yStepGrid), int(n[2]/zStepGrid))
    n = position.shape
    labels = [[[False for i in range(gridDimension[0])] for j in range(gridDimension[1])] for k in range(gridDimension[2])]
    for num in range(n[1]):
        xBig, yBig, zBig = position[time][num][0], position[time][num][1], position[time][num][2]
        i, j, k = int(xBig/xStepGrid), int(yBig/yStepGrid), int(zBig/zStepGrid)
        x, y, z = i*xStepGrid, j*yStepGrid, k*zStepGrid
        ax.scatter(xBig, yBig, color='red')
        ax.plot([position[t][num][0] for t in range(time)], [position[t][num][1] for t in range(time)])
    return ax

