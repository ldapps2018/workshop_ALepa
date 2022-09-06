#!/usr/bin/env python
# coding: utf-8
# %%

# %%


import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib import patches, lines
from ipywidgets import interact, widgets



ti = widgets.FloatSlider(min=0, max=2.1624, step=0.02, value=0.55, description='Zeit [$s$]', continuous_update=False)

def wurf_plot(ti):
    def get_angle_plot(line1, line2, offset = 1, color = None, origin = [0,0], len_x_axis = 2.5, len_y_axis = 1):
        l1xy = line1.get_xydata()
        # Angle between line1 and x-axis
        slope1 = (l1xy[1][1] - l1xy[0][1]) / float(l1xy[1][0] - l1xy[0][0])
        angle1 = abs(math.degrees(math.atan(slope1))) # Taking only the positive angle

        l2xy = line2.get_xydata()

        # Angle between line2 and x-axis
        slope2 = (l2xy[1][1] - l2xy[0][1]) / float(l2xy[1][0] - l2xy[0][0])
        angle2 = abs(math.degrees(math.atan(slope2)))

        theta1 = min(angle1, angle2)
        theta2 = max(angle1, angle2)

        angle = theta2 - theta1

        if color is None:
            color = line1.get_color() # Uses the color of line 1 if color parameter is not passed.

        return patches.Arc(origin, len_x_axis*offset, len_y_axis*offset, 0, theta1, theta2, color=color)
    alpha =  45
    m = 2        # kg
    v0 = 15  # in m/s
    g = 9.81
    tc = 2*v0*math.sin(math.radians(alpha))/g
    t = np.linspace(0, ti, 100)
    vx = v0*math.cos(math.radians(alpha))
    vy = v0*math.sin(math.radians(alpha)) - g*t

    x = (v0*math.cos(math.radians(alpha)))*t
    y = (v0*math.sin(math.radians(alpha)))*t - 0.5*g*t**2
    mt = 1/5
    fig, ax = plt.subplots(figsize=(10,8))
    plt.plot(x,y, zorder = 1)
    plt.ylabel('y in $m$')
    plt.xlabel('x in $m$')
    xf = x[len(x)-1]
    yf = y[len(y)-1]
    vyf = vy[len(vy)-1]
    ax.scatter(xf, yf, s=50, color="r", zorder = 2)
    vxp = ax.quiver(xf, yf, mt*vx, 0, angles='xy', scale_units='xy', scale=1, color='green')
    vyp = ax.quiver(xf, yf, 0, mt*vyf, angles='xy', scale_units='xy', scale=1, color='green')
    vr = ax.quiver(xf, yf, mt*vx, mt*vyf, angles='xy', scale_units='xy', scale=1, color='blue')
    dia1 = ax.plot([xf, xf+mt*vx],[yf+vyf*mt,yf+vyf*mt], color='green', linestyle='dashed')
    dia2 = ax.plot([xf+mt*vx,xf+mt*vx],[yf,yf+vyf*mt], color='green', linestyle='dashed')
    vx_t = plt.text(xf+mt*vx+0.5,yf, r'$\vec v_{x}$', fontsize=12, color='green')
    vy_t = plt.text(xf-1.25, yf+vyf*mt+0.15, r'$\vec v_{y}$', fontsize=12, color='green')
    Vr_t = plt.text(xf+mt*vx+0.5, yf+vyf*mt, r'$\vec v$', fontsize=12, color='blue')
    mt_t = plt.text(-4.5, 7.5, "Maßstab 1 $m$: 5 $m/s$", fontsize=12)
    alpha_t = plt.text(1.25, .25, r'$\alpha$', fontsize=12)
    ground = np.array([[-5,-3],[-5,0],[30,0],[30,-3]])
    ground_patch = plt.Polygon(ground, color =  'saddlebrown', alpha = 0.5)
    plt.gca().add_patch(ground_patch)
    line_1 = lines.Line2D([0,1], [0,1], linewidth=0.1, linestyle = "-", color="blue")
    line_2 = lines.Line2D([0,.1], [0,0], linewidth=0.1, linestyle = "-", color="saddlebrown")
    angle_plot = get_angle_plot(line_1, line_2, 1)
    ax.add_patch(angle_plot)
    ax.grid()
    ax.set_xlim([-5, 30])
    ax.set_ylim([-3, 8])
    ax.minorticks_on()
    return plt.show()

