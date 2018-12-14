from possiblefinal import find_nearest
import numpy as np

bins = 10
vacradius = 1.1
vesradius = 1.2
increment = (vesradius - vacradius)/bins

absorb_radius = np.zeros((bins, 2), dtype = float)
for i in range(bins):
    absorb_radius[i,0] = vacradius + (i+1)*increment

result = find_nearest(absorb_radius[:,0], 1.147)
print(absorb_radius)
print(type(result))
print(result)
