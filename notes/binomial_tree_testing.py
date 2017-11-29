import math

num_procs = 16 

# IT WORKS
for i in range(math.ceil(math.log2(num_procs))):
    print("*****[STEP {}]*****".format(i))

    for j in range(num_procs):
        #root
        if j == 0:
            print("{} receive from {}".format(j,j+2**i))
        
        #other
        elif j % (2**(i+1)) == 0:
            print("{} receive from {}".format(j, j+2**i))
        else:
            if j % (2**i) == 0:
                print("{} send to {}".format(j, j-2**i)) 


