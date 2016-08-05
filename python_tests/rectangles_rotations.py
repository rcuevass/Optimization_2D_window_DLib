import numpy as np
import cv2
import math as mt



def f(x,y,xc,yc,r,mu):
	
	###To set a white rectangle
    val = 1.0
    val = val/( 1. + mt.exp(-mu*(x - xc + r))  )
    val = val/( 1. + mt.exp(-mu*(xc - x + r))  )
    val = val/( 1. + mt.exp(-mu*( y - yc + r))  )
    val = val/( 1. + mt.exp(-mu*( yc - y + r))  )
    return val




def g(x,y,xc,yc,r,mu):

	# To generate a diamond

    val = 1.0
    
    # Right-upper part of diamond
    val = val/( 1. + mt.exp(-mu*( (x - xc) + r - (y - yc)  ))  )
    
    # Right-lower part of diamond
    val = val/( 1. + mt.exp(-mu*( - (x - xc) + r - (y - yc) ))  )
    
    # Left-upper pert of diamond
    val = val/( 1. + mt.exp(-mu*( (x - xc) + (y - yc) + r ))  )
    
    # Left-lower pert of diamond
    val = val/( 1. + mt.exp(-mu*( - (x - xc) + (y - yc) + r ))  )
    return val




def J3(x,y,xc,yc,rx,ry,mu,theta):


	### THIS IS IT!
	### THIS IS THE FUNCTION TO WORK WITH!!!
	### KEEP IN MIND YOU ARE ROTATING THE REFERENCE FRAME !!!

    
    
    sVal = mt.sin(mt.pi*(90-theta)/180)
    cVal = mt.cos(mt.pi*(90-theta)/180)
    
    val = 1.

     # Positive (counter-clockwise) rotation on x-xc and y - yc!!!
    xNew = (x - xc)*cVal - (y - yc)*sVal
    yNew = (x - xc)*sVal + (y - yc)*cVal
    

    
    # Right-lower part of diamond
    val = val/( 1. + mt.exp(-mu*( - yNew + ry ))  )
    
    # Left-upper pert of diamond
    val = val/( 1. + mt.exp(-mu*( yNew + ry ))  )
    
    
    # Right-upper part of diamond
    #val = val/( 1. + mt.exp(-mu*mt.sqrt(2)*( (x - xc)/mt.sqrt(2) - (y - yc)/mt.sqrt(2) + r   ))  )
    
    val = val/( 1. + mt.exp(-mu*(  xNew + rx   ))  )
    
    
    # Left-lower pert of diamond
    val = val/( 1. + mt.exp(-mu*( -xNew + rx ))  )
    
    return val


img = np.zeros((512,512), np.uint8)
cv2.imwrite('black_canvas.jpg',img)



# White rectangle
xc = 256
yc = 256
mu = 4.0
r = 100
for j in range(512):
    for i in range(512):
        #print i , j , f(i,j,xc,yc,r,mu)
        img[i][j] = 255*f(i,j,xc,yc,r,mu)

cv2.imwrite('square_01.jpg',img)





img = np.zeros((512,512), np.uint8)
cv2.imwrite('black_canvas.jpg',img)



xc = 256
yc = 256
mu = 1.0
r = 100
for j in range(512):
    for i in range(512):
        #print i , j , f(i,j,xc,yc,r,mu)
        img[i][j] = 255*g(i,j,xc,yc,r,mu)


cv2.imwrite('diamond_01.jpg',img)






img = np.zeros((512,512), np.uint8)
cv2.imwrite('black_canvas.jpg',img)



theta = 30;
xc = 256
yc = 256
mu = 1.0
rx = 100
ry = 50
for j in range(512):
    for i in range(512):
        #print i , j , f(i,j,xc,yc,r,mu)
        img[i][j] = 255*J3(i,j,xc,yc,rx,ry,mu,theta)


cv2.imwrite('rotated_rectangle_01.jpg',img)

