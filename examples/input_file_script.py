import numpy as np
import matplotlib.pyplot as plt

x_rough= np.loadtxt("x.out")
y_rough= np.loadtxt("y.out")
len_arrays= 2*len(x_rough) + 4   #(4 corners + 4 points on the edge of the fault zone )
height_of_domain= 40     #(km)


 # working on the node part of dynosol 2d

x = np.zeros(len_arrays)
y= np.zeros(len_arrays)

x[0]= 0.0              											    		# top left corner
x[1]= x_rough[len(x_rough)-1]													# top right corner

x[2]=  x_rough[0]															# upper center left corner
x[3:3+len(x_rough)-2]=  x_rough[1: (len(x_rough)-1 )]   								# upper level of body of rough fault
x[3+len(x_rough)-2]=  x_rough[len(x_rough)-1]									# upper center right corner	

x[3+len(x_rough)-2+ 1]=	x_rough[0]												# lower center left corner
x[3+len(x_rough)-2+ 2 : 3+ 2* (len(x_rough)-2) +2 ]=  x_rough[1:len(x_rough)-1]    		# lower level of body of rough fault
x[3+ 2* (len(x_rough)-2) +2] = x_rough[len(x_rough)-1]    					 		# lower center right corner

x[3+ 2* (len(x_rough)-2) +3] =  0.0     											# lower left corner
x[3+ 2* (len(x_rough)-2) +4] =  x_rough[len(x_rough)-1]    							# lower right corner



y[0]= height_of_domain               										# top left corner
y[1]= height_of_domain														# top right corner

y[2]= 1.0 + y_rough[0] 														# upper center left corner
y[3:3+len(x_rough)-2]= 1.0 + y_rough[1: (len(x_rough)-1 )]  				 		# upper level of body of rough fault
y[3+len(x_rough)-2]=  1.0 + y_rough[len(x_rough)-1]								# upper center right corner	

y[3+len(x_rough)-2+ 1]=	y_rough[0] -1.0  													# lower center left corner
y[3+len(x_rough)-2+ 2 : 3+ 2* (len(x_rough)-2) +2 ]=  y_rough[1:len(x_rough)-1] -1.0  					# lower level of body of rough fault
y[3+ 2* (len(x_rough)-2) +2] =  y_rough[len(x_rough)-1] -1.0    									# lower center right corner

y[3+ 2* (len(x_rough)-2) +3] =  0.0    															 # lower left corner
y[3+ 2* (len(x_rough)-2) +4] =  0.0     															# lower right corner

# np.savetxt('x_values.txt', x)
# np.savetxt('y_values.txt', y)
x= x*1000.0                          # convert to km
y=y*1000.0				 
 

# working on the element part of dynosol 2d

element= np.zeros(len_arrays+2)
pj0=	np.zeros(len_arrays+2)
pj1=	np.zeros(len_arrays+2)
boundary_flag=  np.zeros(len_arrays+2)




# left boundary from top to bottom ----------------------
pj0[0]= 0 
pj1[0]=2
boundary_flag[0]= 1

pj0[1]= 2 
pj1[1]=3+len(x_rough)-2+ 1
boundary_flag[1]= 1

pj0[2]=3+len(x_rough)-2+ 1
pj1[2]=3+ 2* (len(x_rough)-2) +3
boundary_flag[2]= 1

# Bottom boundary ----------------------

pj0[3]=  3+ 2* (len(x_rough)-2) +3
pj1[3]= 3+ 2* (len(x_rough)-2) +4
boundary_flag[3]= 16

# Right boundary from bottom to top----------------------
pj0[4]=  3+ 2* (len(x_rough)-2) +4
pj1[4]= 3+ 2* (len(x_rough)-2) +2
boundary_flag[4]= 2

pj0[5]=  3+ 2* (len(x_rough)-2) +2
pj1[5]=  3+len(x_rough)-2
boundary_flag[5]= 2

pj0[6]=  3+len(x_rough)-2
pj1[6]=  1
boundary_flag[6]= 2

# Top boundary ----------------------
pj0[7]=  1
pj1[7]=  0
boundary_flag[7]= 32

# Not a boundary-- Body segments -------------------



for ii in range(len(x_rough)-1):
	next_value= 8+ ii
	pj0[8+ ii]=  2+ii
	pj1[8+ii] = 2+ii +1  
	boundary_flag[8+ii] = 0

for ii in range(len(x_rough)-1):	
	pj0[next_value+1+ii] = 3+len(x_rough)-2+ 1 + ii
	pj1[next_value+1+ii] = 3+len(x_rough)-2+ 1 + ii +1
	boundary_flag[next_value+1+ii] = 0 


#Write output
f = open('coupling_input.poly','w')

#need to write this part
# npoints ndims 0 0
#  13      2     0 0
f.write('# input file for dynosol\n')
f.write('# \n')
f.write("{} {} {} {}\n".format( '#npoints' , 'ndims', '0', '0' ) )
f.write("{} {} {} {}\n".format( len_arrays , 2, 0, 0 ) )
f.write("{} {} {} \n".format( '#i' , 'xi', 'yi' ) )
for i in range(len_arrays):

	f.write("{} {:E} {:E}\n".format(i , x[i], y[i] ) ) 

## nsegments 1
#  16        1
f.write('# segments\n')
f.write("{} {}\n".format( '#nsegments' , '1' ) )
f.write("{:d} {:d}\n".format( len_arrays+2 , 1 ) )
f.write("{} {} {} {}\n".format( '#i' , 'pj0', 'pj1', 'boundary_flag' ) )


for i in range(len_arrays+2):
	pj_0= int ( pj0[i] )
	pj_1 =int ( pj1[i] )
	boundary= int( boundary_flag[i] )
	f.write("{} {:d} {:d} {:d}\n".format(i , pj_0, pj_1, boundary ) ) 
#f.write('# author='+author+'\n')

f.write('# #### holes, must be 0 ####\n')
f.write("{:d}\n".format( 0 ) )

f.write('#### regions ####\n')
f.write('# nregions\n')
f.write( "{}\n".format(3) )
f.write("{} {} {} {}\n".format( '#k' , 'xk', 'yk', 'mattype', 'size' ) )

# Working on the regions part --------------------------------------


zones_1_x = 40000.0 
zones_1_y=  10000.0
element_size_1= 400000.0
element_type_1= 0

zones_2_x = 0.0
#zones_2_x = x_rough[1] *1000 +500
zones_2_y=  20100.0
#zones_2_y = y_rough[1]*1000 +500.0
element_size_2 = 10000.0
element_type_2= 1

zones_3_x = 40000.0
zones_3_y=  30000.0
element_size_3= 400000.0
element_type_3= 0

f.write("{} {:E} {:E} {:d} {:E}\n".format(0 , zones_1_x, zones_1_y, element_type_1, element_size_1 ) )
f.write("{} {:E} {:E} {:d} {:E}\n".format(1 , zones_2_x, zones_2_y, element_type_2, element_size_2 ) )
f.write("{} {:E} {:E} {:d} {:E}\n".format(2 , zones_3_x, zones_3_y ,element_type_3, element_size_3) )

f.close()






