import matplotlib.pyplot as plt
import sys

pa = []
ho = []

input_file = sys.argv[1]
num_rw = int(sys.argv[2])

first_line = True

with open(input_file) as f:
    for line in f:
        if first_line is True:
            first_line = False
            continue
        items = line.split()
        if num_rw != int(items[4]):
            continue
        if int(items[1]) == 0 and int(items[2]) == 0:
            pa.append(float(items[5]))
        if int(items[1]) == 100 and int(items[2]) == 100:
            ho.append(float(items[5]))


plt.plot(pa, ho, "o", [0, 1], [0, 1], "-")

plt.xlabel("Personalized")
plt.ylabel('High-order')
plt.title(input_file + ', rw = ' + str(num_rw))
plt.axis([0, 1.0, 0, 1.0])
plt.grid(True)

plt.show()
