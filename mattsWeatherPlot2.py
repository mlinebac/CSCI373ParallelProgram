from matplotlib import pyplot
from pylab import genfromtxt

data = genfromtxt('weatherData.txt')

diff = (data[:, 5])-(data[:, 3])
pyplot.plot(diff, label="Difference")
pyplot.suptitle("Difference between Average High and Observed High")
pyplot.ylabel("Difference in Temperature (F)")
pyplot.xlabel("Day of the Year")
pyplot.legend()
pyplot.savefig("DiffPlot.pdf")
