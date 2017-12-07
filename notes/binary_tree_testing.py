from math import ceil,log2

num_procs = 8
root = 1 
start = 0
my_flat_circle = [(k + root) % num_procs for k in range(num_procs)]
print(my_flat_circle)

'''
# use this algorithm to tell each process what level they are on
print("printing each processes's level....")
for i in range(ceil(log2(num_procs)) + 1):
    bound = start + 2**i
    print("start = {}, bound = {}".format(start, bound))
    for j in range(start, bound):
        if j < num_procs: 
            print("node: {}, level: {}".format(my_flat_circle[j], i)) 
        else:
            break

    start = bound
'''

'''
print("")

root = 1

print("printing ring")
for k in range(num_procs):
    print((k + root) % num_procs)
'''

'''
print("testing number of levels")
for j in range(12):
    print("num_procs = {}, num_levels = {}".format(j, ceil(log2(j + 1) - 1)))
'''
