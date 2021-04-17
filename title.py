# Written by: Nick Gerend, @dataoutsider
# Viz: "Open Road", enjoy!

import numpy as np 
import cv2
import os

#region prepare data
# load image 
img2 = cv2.imread(os.path.dirname(__file__) + '/title.jpg', cv2.IMREAD_COLOR) 
   
# duplicate to gray scale
img = cv2.imread(os.path.dirname(__file__) + '/title.jpg', cv2.IMREAD_GRAYSCALE) 

_,threshold = cv2.threshold(img, 100, 255, cv2.THRESH_BINARY) #20,50
#threshold = cv2.adaptiveThreshold(img,50,cv2.ADAPTIVE_THRESH_GAUSSIAN_C,cv2.THRESH_BINARY,151,9)
   
# detect shapes 
contours,_ = cv2.findContours(threshold, cv2.RETR_TREE, cv2.CHAIN_APPROX_SIMPLE)
#endregion
   
#region algorithm
# find polygons within regions
poly = []
index = 0
print(len(contours))
for cnt in contours : 
    area = cv2.contourArea(cnt) 
   
    # Shortlisting the regions based on there area. 
    if area >= 10:  
        polygon = cv2.approxPolyDP(cnt, 0.0000001 * cv2.arcLength(cnt, True), True) #.009
   
        # Checking if the no. of sides of the selected region is 7. 
        if(len(polygon) >= 5):  
            
            cv2.drawContours(img2, [polygon], 0, (0, 0, 255), 5)
            poly.append(polygon)
    index += 1
   
# # Showing the image along with outlined arrow. 
# #img2 = cv2.resize(img2, (900, 900))
# cv2.imshow('image2', img2)  
   
# # Exiting the window if 'q' is pressed on the keyboard. 
# if cv2.waitKey(0) & 0xFF == ord('q'):  
#     cv2.destroyAllWindows()
#endregion

#region output
path = 1
it = 1
import csv
with open(os.path.dirname(__file__) + '/title.csv', 'w',) as csvfile:
    writer = csv.writer(csvfile, lineterminator = '\n')
    writer.writerow(['x', 'y', 'path', 'item'])
    for polygon in poly:
        for j in polygon:
            writer.writerow([j[0][0],j[0][1],path,it])
            path += 1
        path = 1
        it += 1
#endregion

print('finish')