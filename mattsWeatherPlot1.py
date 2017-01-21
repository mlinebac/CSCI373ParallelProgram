from matplotlib import pyplot
from pylab import genfromtxt

data = genfromtxt('weatherData.txt')

pyplot.plot(data[:, 2], 'p', label="ObservedLo")
pyplot.plot(data[:, 3], 'p', label="ObservedHi")
pyplot.plot(data[:, 4], label="normalLo")
pyplot.plot(data[:, 5], label="normalHi")
pyplot.suptitle("Asheville Weather 2015")
pyplot.ylabel("Temperature (F)")
pyplot.xlabel("Day of the Year")
pyplot.legend()
pyplot.savefig("MyPlot.pdf")
