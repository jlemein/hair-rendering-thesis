import numpy as np

data = []
for theta in range(0, 3):
    row = []
    for phi in range(10, 13):
        row.append("theta: " + str(theta) + "phi: " + str(phi))
    data.append(row)

npdata = np.array(data)
print("data: ", data)
print("np data: ", npdata)
print("axis should be phi: ", npdata[:,0])
print("axis should be theta: ", npdata[0,:])