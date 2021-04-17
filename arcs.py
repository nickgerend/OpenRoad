# Written by: Nick Gerend, @dataoutsider
# Viz: "Open Road", enjoy!

import pandas as pd
import os
import math
import scipy.optimize as optimize
import numpy as np
import matplotlib.pyplot as plt

class point:
    def __init__(self, index, item, x, y, path, sub_item = ''): 
        self.index = index
        self.item = item
        self.x = x
        self.y = y
        self.path = path
        self.sub_item = sub_item
    def to_dict(self):
        return {
            'index' : self.index,
            'item' : self.item,
            'x' : self.x,
            'y' : self.y,
            'path' : self.path,
            'sub_item' : self.sub_item }

#region Functions

def angle_func(x, a, c):
    return 2*a*np.sin(x/2.*(np.pi/180.))-c*(x*np.pi/180.)

def circle_angle_radius(arc_length, chord_length):
    t = optimize.bisect(angle_func, 0.01, 359.9, args=(arc_length, chord_length), maxiter=1000)
    return t, arc_length/(t*np.pi/180.)

def circle_origin(x1, y1, x2, y2, r, flip=False):
    xa = 0.5*(x2-x1)
    ya = 0.5*(y2-y1)
    xo = x1+xa
    yo = y1+ya
    a = np.sqrt(xa**2.+ya**2.)
    b = np.sqrt(r**2.-a**2.)
    sign = 1.
    if flip:
        sign = -1.
    x_o = xo+sign*(b*ya/a)
    y_o = yo-sign*(b*xa/a)
    return x_o, y_o

def angle_from_origin(xo, yo, x, y):
    return math.atan2(y-yo, x-xo)*180./math.pi

def DistBtwTwoPnts(x1, y1, x2, y2):
    return math.sqrt((x2-x1)**2+(y2-y1)**2)

def arc(x1, y1, x2, y2, arc_factor, points, flip=False, radius=None, xo_=None, yo_=None, a1=-90, a2=270):
    r = 0
    xo = 0
    yo = 0
    if radius is not None:
        r = radius
    else:
        chord_length = DistBtwTwoPnts(x1, y1, x2, y2)
        arc_length = chord_length*arc_factor
        a,r = circle_angle_radius(arc_length, chord_length)
    if xo_ is not None and xo_ is not None:
        xo = xo_
        yo = yo_
    else:
        xo,yo = circle_origin(x1, y1, x2, y2, r, flip)
    angles = np.linspace(a1, a2, num=points)
    xs_arc = r*np.sin(angles*np.pi/180.)+xo
    ys_arc = r*np.cos(angles*np.pi/180.)+yo

    return list(zip(xs_arc,ys_arc)),r,xo,yo

def line(p1, p2):
    A = (p1[1] - p2[1])
    B = (p2[0] - p1[0])
    C = (p1[0]*p2[1] - p2[0]*p1[1])
    return A, B, -C

def intersection(x1_1, y1_1, x2_1, y2_1, x1_2, y1_2, x2_2, y2_2):
    L1 = line([x1_1,y1_1], [x2_1,y2_1])
    L2 = line([x1_2,y1_2], [x2_2,y2_2])
    D  = L1[0] * L2[1] - L1[1] * L2[0]
    Dx = L1[2] * L2[1] - L1[1] * L2[2]
    Dy = L1[0] * L2[2] - L1[2] * L2[0]
    if D != 0:
        x = Dx / D
        y = Dy / D
        return x,y
    else:
        return False

def line_from_point_angle(x, y, angle):
    m = math.tan((angle-90)*math.pi/180.)
    b = y-m*x
    return m, b

def rescale(x, xmin, xmax, newmin, newmax):
    rescaled = (newmax-newmin)*((x-xmin)/(xmax-xmin))+newmin
    return rescaled

#endregion

#region Load Data
df = pd.read_csv(os.path.dirname(__file__) + '/Interstates_Clean2.csv')
#endregion

#region scale
min_y_x1 = df['S_W_Lat'].min()
max_y_x1 = df['S_W_Lat'].max()
min_x_y1 = df['S_W_Long'].min()
max_x_y1 = df['S_W_Long'].max()

min_y_x2 = df['N_E_Lat'].min()
max_y_x2 = df['N_E_Lat'].max()
min_x_y2 = df['N_E_Long'].min()
max_x_y2 = df['N_E_Long'].max()

buffer = 10

min_x = min(min_x_y1,min_x_y2) - buffer
max_x = min(max_x_y1,max_x_y2) + buffer
min_y = min(min_y_x1,min_y_x2) - buffer
max_y = min(max_y_x1,max_y_x2) + buffer

df['S_W_Lat_norm'] =  [rescale(x, min_y, max_y, 0, 1) for x in df['S_W_Lat']]  
df['S_W_Long_norm'] =  [rescale(x, min_x, max_x, 0, 1) for x in df['S_W_Long']]
df['N_E_Lat_norm'] =  [rescale(x, min_y, max_y, 0, 1) for x in df['N_E_Lat']]
df['N_E_Long_norm'] =  [rescale(x, min_x, max_x, 0, 1) for x in df['N_E_Long']]

d_buffer = 100
min_d = df['Length_Miles'].min()-d_buffer
max_d = df['Length_Miles'].max()+d_buffer
df['width'] =  [rescale(x, min_d, max_d, 0, 1)/5. for x in df['Length_Miles']]
#endregion

#region Algorithm
H = 'S_W_Lat_norm'
V = 'S_W_Long_norm'
box_buffer = .2
points = 200
ratio_scale = 10. #5.
list_xy = []
jx = 0
boarder_ratio = 0.25
line_length_ratio = 10./12.
line_sep_ratio = 30./12.
line_width_ratio = (4./12.)/12.
line_points = 10
for j, row in df.iterrows():
    
    #region Direction
    direction = row['Direction']
    start = 0.
    if direction == 'H':
        start = row[H]
    else:
        start = row[V]
    #endregion

    #region Endpoint
    angle = row['Heading']
    if direction == 'V':
        angle += 90
    m_x,b_x = line_from_point_angle(0, start, angle)
    end = m_x*1. + b_x
    x_i = 1
    y_i = end
    #endregion

    #region Midpoint
    f_m = m_x*0.5 + b_x
    f = True
    y_l = 1
    if f_m < 0.5:
        f = False
        y_l = 0
    #endregion

    #region Adjust for Intersection
    if end > 1 or end < 0:
        x_i,y_i = intersection(0, y_l, 1, y_l, 0, start, 1, end)
    #endregion

    #region Outlines
    ratio = row['arc_ratio']
    width = row['width']
    if ratio > 1:
        ratio = 1.
    row_xy,r,xo,yo = arc(0, start, x_i, y_i, 1+ratio/ratio_scale, points, flip=f)

    #total width
    r1 = r-width
    r2 = r+width
    row_xy = arc(0, start, 1, end, None, points, flip=f, radius=r1, xo_=xo, yo_=yo)[0]
    row_xy2 = arc(0, start, 1, end, None, points, flip=f, radius=r2, xo_=xo, yo_=yo)[0]

    #boarders
    b1_1 = r1
    b1_2 = r1-width*boarder_ratio
    b1_1_row_xy = arc(0, start, 1, end, None, points, flip=f, radius=b1_1, xo_=xo, yo_=yo)[0]
    b1_2_row_xy2 = arc(0, start, 1, end, None, points, flip=f, radius=b1_2, xo_=xo, yo_=yo)[0]
    b2_1 = r2
    b2_2 = r2+width*boarder_ratio
    b2_1_row_xy = arc(0, start, 1, end, None, points, flip=f, radius=b2_1, xo_=xo, yo_=yo)[0]
    b2_2_row_xy2 = arc(0, start, 1, end, None, points, flip=f, radius=b2_2, xo_=xo, yo_=yo)[0]

    #markers
    lanes = row['States'].count(",") + 1
    lane_width = width*2/lanes
    list_lanes = []
    lane_radius = r1 + lane_width
    line_length_scaled = lane_width*2*line_length_ratio
    line_sep_scaled = lane_width*2*line_sep_ratio
    line_fraction = line_length_scaled/(line_length_scaled+line_sep_scaled)
    line_width_scaled = lane_width*line_width_ratio/2.
    line_itr = 0
    for i in range(lanes-1):
        a_start = 66 #-90
        circumference = 2*math.pi*lane_radius
        lines = circumference/(line_length_scaled+line_sep_scaled)
        for k in range(int(lines)):
            angles_1 = np.linspace(a_start, a_start+(360./lines)*line_fraction, num=line_points)
            angles_2 = np.linspace(a_start+(360./lines)*line_fraction, a_start, num=line_points)

            x_line_1 = (lane_radius-line_width_scaled)*np.sin(angles_1*np.pi/180.)+xo
            y_line_1 = (lane_radius-line_width_scaled)*np.cos(angles_1*np.pi/180.)+yo
            line_path_1 = np.linspace(1, line_points, num=line_points)
            line_1 = list(zip(x_line_1, y_line_1, line_path_1))

            x_line_2 = (lane_radius+line_width_scaled)*np.sin(angles_2*np.pi/180.)+xo
            y_line_2 = (lane_radius+line_width_scaled)*np.cos(angles_2*np.pi/180.)+yo
            line_path_2 = np.linspace(line_points+1, line_points*2, num=line_points)
            line_2 = list(zip(x_line_2, y_line_2, line_path_2))

            a_start += 360./lines
            itr = np.linspace(line_itr, line_itr, num=line_points*2)
            line_12 = line_1 + line_2
            
            list_lanes += list(zip(line_12, itr))
            line_itr += 1
        lane_radius += lane_width
    #endregion

    #region Output

    #region total width
    for i in range(len(row_xy)):
        if direction == 'H':
            list_xy.append(point(jx, row['i_Segment'], row_xy[i][0], row_xy[i][1], i, 'road'))
        else:
            list_xy.append(point(jx, row['i_Segment'], row_xy[i][1], row_xy[i][0], i, 'road'))
        jx += 1
    for i in range(len(row_xy2)):
        if direction == 'H':
            list_xy.append(point(jx, row['i_Segment'], row_xy2[i][0], row_xy2[i][1], (len(row_xy)+len(row_xy2))-i, 'road'))
        else:
            list_xy.append(point(jx, row['i_Segment'], row_xy2[i][1], row_xy2[i][0], (len(row_xy)+len(row_xy2))-i, 'road'))
        jx += 1

    for i in range(len(row_xy)):
        if direction == 'H':
            list_xy.append(point(jx, row['i_Segment'], row_xy[i][0], row_xy[i][1], i, 'a_road'))
        else:
            list_xy.append(point(jx, row['i_Segment'], row_xy[i][1], row_xy[i][0], i, 'a_road'))
        jx += 1
    for i in range(len(row_xy2)):
        if direction == 'H':
            list_xy.append(point(jx, row['i_Segment'], row_xy2[i][0], row_xy2[i][1], (len(row_xy)+len(row_xy2))-i, 'a_road'))
        else:
            list_xy.append(point(jx, row['i_Segment'], row_xy2[i][1], row_xy2[i][0], (len(row_xy)+len(row_xy2))-i, 'a_road'))
        jx += 1
    #endregion

    #region boarders
    for i in range(len(b1_1_row_xy)):
        if direction == 'H':
            list_xy.append(point(jx, row['i_Segment'], b1_1_row_xy[i][0], b1_1_row_xy[i][1], i, 'b1'))
        else:
            list_xy.append(point(jx, row['i_Segment'], b1_1_row_xy[i][1], b1_1_row_xy[i][0], i, 'b1'))
        jx += 1
    for i in range(len(b1_2_row_xy2)):
        if direction == 'H':
            list_xy.append(point(jx, row['i_Segment'], b1_2_row_xy2[i][0], b1_2_row_xy2[i][1], (len(b1_1_row_xy)+len(b1_2_row_xy2))-i, 'b1'))
        else:
            list_xy.append(point(jx, row['i_Segment'], b1_2_row_xy2[i][1], b1_2_row_xy2[i][0], (len(b1_1_row_xy)+len(b1_2_row_xy2))-i, 'b1'))
        jx += 1
    for i in range(len(b2_1_row_xy)):
        if direction == 'H':
            list_xy.append(point(jx, row['i_Segment'], b2_1_row_xy[i][0], b2_1_row_xy[i][1], i, 'b2'))
        else:
            list_xy.append(point(jx, row['i_Segment'], b2_1_row_xy[i][1], b2_1_row_xy[i][0], i, 'b2'))
        jx += 1
    for i in range(len(b2_2_row_xy2)):
        if direction == 'H':
            list_xy.append(point(jx, row['i_Segment'], b2_2_row_xy2[i][0], b2_2_row_xy2[i][1], (len(b2_1_row_xy)+len(b2_2_row_xy2))-i, 'b2'))
        else:
            list_xy.append(point(jx, row['i_Segment'], b2_2_row_xy2[i][1], b2_2_row_xy2[i][0], (len(b2_1_row_xy)+len(b2_2_row_xy2))-i, 'b2'))
        jx += 1
    
    for i in range(len(b1_1_row_xy)):
        if direction == 'H':
            list_xy.append(point(jx, row['i_Segment'], b1_1_row_xy[i][0], b1_1_row_xy[i][1], i, 'a_b1'))
        else:
            list_xy.append(point(jx, row['i_Segment'], b1_1_row_xy[i][1], b1_1_row_xy[i][0], i, 'a_b1'))
        jx += 1
    for i in range(len(b1_2_row_xy2)):
        if direction == 'H':
            list_xy.append(point(jx, row['i_Segment'], b1_2_row_xy2[i][0], b1_2_row_xy2[i][1], (len(b1_1_row_xy)+len(b1_2_row_xy2))-i, 'a_b1'))
        else:
            list_xy.append(point(jx, row['i_Segment'], b1_2_row_xy2[i][1], b1_2_row_xy2[i][0], (len(b1_1_row_xy)+len(b1_2_row_xy2))-i, 'a_b1'))
        jx += 1
    for i in range(len(b2_1_row_xy)):
        if direction == 'H':
            list_xy.append(point(jx, row['i_Segment'], b2_1_row_xy[i][0], b2_1_row_xy[i][1], i, 'a_b2'))
        else:
            list_xy.append(point(jx, row['i_Segment'], b2_1_row_xy[i][1], b2_1_row_xy[i][0], i, 'a_b2'))
        jx += 1
    for i in range(len(b2_2_row_xy2)):
        if direction == 'H':
            list_xy.append(point(jx, row['i_Segment'], b2_2_row_xy2[i][0], b2_2_row_xy2[i][1], (len(b2_1_row_xy)+len(b2_2_row_xy2))-i, 'a_b2'))
        else:
            list_xy.append(point(jx, row['i_Segment'], b2_2_row_xy2[i][1], b2_2_row_xy2[i][0], (len(b2_1_row_xy)+len(b2_2_row_xy2))-i, 'a_b2'))
        jx += 1
    #endregion

    #region markers
    for i in range(len(list_lanes)):
        if direction == 'H':
            list_xy.append(point(jx, row['i_Segment'], list_lanes[i][0][0], list_lanes[i][0][1], list_lanes[i][0][2], 'marker_'+str(list_lanes[i][1])))
        else: 
            list_xy.append(point(jx, row['i_Segment'], list_lanes[i][0][1], list_lanes[i][0][0], list_lanes[i][0][2], 'marker_'+str(list_lanes[i][1])))
        jx += 1
    #endregion

    #region reference points
    if direction == 'H':
        list_xy.append(point(jx, row['i_Segment'], 0., start, 1, 'ref_start'))
        jx += 1
        list_xy.append(point(jx, row['i_Segment'], x_i, y_i, 2, 'ref_end'))
        jx += 1
    else:
        list_xy.append(point(jx, row['i_Segment'], start, 0., 1, 'ref_start'))
        jx += 1
        list_xy.append(point(jx, row['i_Segment'], y_i, x_i, 2, 'ref_end'))
        jx += 1
    #endregion

    #endregion

#endregion

df_out = pd.DataFrame.from_records([s.to_dict() for s in list_xy])
df_out.to_csv(os.path.dirname(__file__) + '/arcs.csv', encoding='utf-8', index=False)

print('finished')