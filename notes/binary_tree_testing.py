from math import ceil,log2

num_procs = 8 
start = 0

# use this algorithm to tell each process what level they are on

for i in range(ceil(log2(num_procs)) + 1 + 1):
    bound = start + 2**i
    print("start = {}, bound = {}".format(start, bound))
    for j in range(start, bound):
        if j <= num_procs:
            print("node: {}, level: {}".format(j, i)) 
        else:
            break

    start = bound

