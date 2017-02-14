#Author: Matt Lineback
#CSCI 373 Parallel Programming HW 3

import matplotlib.pyplot as plt
import numpy as np

dat = np.loadtxt('timing_data.txt', skiprows=1)

proc4 = [0,1,2,3,4,5,6,7,8]
proc3 = [9,10,11,12,13,14,15,16,17]
proc2 = [18,19,20,21,22,23,24,25,26]
proc1 = [27,28,29,30,31,32,33,34,35]

plt.plot(dat[proc4, 1], dat[proc4, 2], label = "4 procs")
plt.plot(dat[proc3, 1], dat[proc3, 2], label = "3 procs")
plt.plot(dat[proc2, 1], dat[proc2, 2], label = "2 procs")
plt.plot(dat[proc1, 1], dat[proc1, 2], label = "1 procs")
plt.suptitle("Timing Exercise Showing All Four Processors")
plt.ylabel("Size of Problem")
plt.xlabel("Seconds")
plt.legend(loc = 'lower right')
plt.show()
#plt.savefig("timing_all_procs")

plt.plot(dat[proc4, 1], dat[proc4, 2], marker ='o', linestyle='--', label = "4 procs")
plt.suptitle("Timings for Different Problem Sizes")
plt.ylabel("Size of Problem")
plt.xlabel("Seconds")
plt.legend(loc = 'lower right')
plt.show()
#plt.savefig("timing_diff_n.png")

Cols = [0,9,18,27]
plt.bar(dat[Cols, 0], dat[Cols, 1], .5, align = 'center')
plt.suptitle("Timings for each Number of Processors while n = 1000000000")
plt.ylabel("Seconds")
plt.xlabel("Number of Processors")
#plt.show()
plt.savefig("timing_num_procs.png")
