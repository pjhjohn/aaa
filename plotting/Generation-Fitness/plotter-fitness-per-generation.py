
textlen = 0
with open('output') as text :
	data = text.read().split('\n')

best, average = [], []
for i in range(len(data)) :
	row = data[i].split()
	best.append(float(row[0]))
	average.append(float(row[1]))

axis = range(1, len(data)+1)
import matplotlib.pyplot as plt
plt.plot(axis, best)
plt.plot(axis, average)
plt.axis([0, len(data), 800, 1000])
plt.show()