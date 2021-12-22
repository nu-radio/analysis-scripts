import numpy as np
import numpy
import gzip
import itertools  

def get_string_heads(file, **kwargs):
    with gzip.GzipFile(file) as f:
        geo = np.loadtxt(f, dtype=np.dtype(
            [('string', int), ('om', int)] + [(c, float) for c in 'xyz']))
        pos = geo[geo['om'] == 1]
        return pos[list('xy')]


def get_sectors():
    sectors = dict()
    sectors["Dark Sector"] = numpy.array([
        [ -6745.5344331313745 , -6666.666666666667 ],
        [  360.24866210820073 , -704.0174330560367 ],
        [   929.3377603199988 ,  732.0277663199995 ],
        [ -1230.6941085521246 ,  6666.666666666667 ],
        [ -1230.6941085521246 ,  6666.666666666667 ],
        [ -11166.666666666666 ,  6666.666666666667 ],
        [            -12500.0 ,  5333.333333333334 ],
        [            -12500.0 , -5333.333333333334 ],
        [ -11166.666666666666 , -6666.666666666667 ],
        [ -6745.5344331313745 , -6666.666666666667 ],
        ])

    sectors["Downwind-West Sector"] = numpy.array([
        [ -2002.6843965330781 ,  -6666.666666666667 ],
        [  -566.7523103221163 , -3043.2210840457365 ],
        [  360.24866210820073 ,  -704.0174330560367 ],
        [ -6745.5344331313745 ,  -6666.666666666667 ],
        [ -6745.5344331313745 ,  -6666.666666666667 ],
        [ -2002.6843965330781 ,  -6666.666666666667 ],
        ])

    sectors["Downwind-East Sector"] = numpy.array([
        [   318.6018363789806 ,  -6666.666666666667 ],
        [   567.9378335400033 , -3986.2232073738505 ],
        [  -566.7523103221163 , -3043.2210840457365 ],
        [ -2002.6843965330781 ,  -6666.666666666667 ],
        [ -2002.6843965330781 ,  -6666.666666666667 ],
        [   318.6018363789806 ,  -6666.666666666667 ],
        ])

    sectors["Clean Air Sector"] = numpy.array([
        [ -1230.6941085521246 ,   6666.666666666667 ],
        [   929.3377603199988 ,   732.0277663199995 ],
        [  1323.2565559199993 ,  -350.2552507199998 ],
        [  2039.2957245914913 ,  -599.4482567216783 ],
        [   7499.999999999999 , -2499.8600956438077 ],
        [   7499.999999999999 , -2499.8600956438077 ],
        [   7499.999999999999 ,   5333.333333333334 ],
        [   6166.666666666666 ,   6666.666666666667 ],
        [ -1230.6941085521246 ,   6666.666666666667 ],
        ])

    sectors["Quiet Sector"] = numpy.array([
        [  7499.999999999999 , -2499.8600956438077 ],
        [ 2039.2957245914913 ,  -599.4482567216783 ],
        [  825.1526280298203 , -1221.0801017904032 ],
        [  567.9378335400033 , -3986.2232073738505 ],
        [  318.6018363789806 ,  -6666.666666666667 ],
        [  318.6018363789806 ,  -6666.666666666667 ],
        [  6166.666666666666 ,  -6666.666666666667 ],
        [  7499.999999999999 ,  -5333.333333333334 ],
        [  7499.999999999999 , -2499.8600956438077 ],
        ])
    return sectors

sectorColors = {
    "Dark Sector" : ('g', 0.1), 
    "Downwind-West Sector" : ('b', 0.1), 
    "Downwind-East Sector" : ('y', 0.1), 
    "Clean Air Sector" : ('m', 0.1), 
    "Quiet Sector" : ('c', 0.1), 
    "Old South Pole Station Hazard Area" : ('r', 0.5),
    "IceCube Perimeter" : ('k', 0.1), 
}

def ray(points):
    x0 = points[0]
    dx = np.diff(points, axis=0)[0, :]
    return lambda x: x0[1] + dx[1]/dx[0]*(x-x0[0])

def half_fantasy_radio_geometry(sector, spacing=1.5e3, nstations=200):
    edge = sector[1:3]
    top = ray(sector[2:4])
    bottom = ray(sector[0:2])

    # basis vectors along the skiway edge of the dark sector
    dx = np.diff(edge, axis=0)[0, :]
    dx /= np.hypot(*dx)
    perpdx = -np.asarray([dx[1], -dx[0]])
    x0 = edge[0]
    print(x0)
    print(dx)
    print(spacing)

    locations = []

    stations = nstations*50
    for row in itertools.count():
        corner = x0 + row*perpdx*spacing
        for col in itertools.count():
            pos = corner + col*dx*spacing
            if pos[1] > top(pos[0]):
                break
            locations.append(pos)
            if len(locations) >= stations:
                break
        for col in itertools.count(1):
            pos = corner - col*dx*spacing
            if pos[1] < bottom(pos[0]):
                break
            locations.append(pos)
            if len(locations) >= stations:
                break
        if len(locations) >= stations:
            break

    locs = np.vstack(locations)

    order = np.hypot(*locs.T).argsort()
    return locs[order][:nstations]


def rotate_point(point_x, point_y, origin_x, origin_y, degrees):
    radians = np.deg2rad(degrees)
    x,y = point_x, point_y
    offset_x, offset_y = origin_x, origin_y
    adjusted_x = (x - offset_x)
    adjusted_y = (y - offset_y)
    cos_rad = np.cos(radians)
    sin_rad = np.sin(radians)
    qx = offset_x + cos_rad * adjusted_x + sin_rad * adjusted_y
    qy = offset_y + -sin_rad * adjusted_x + cos_rad * adjusted_y
    return qx, qy


def trim_to_fit(in_x, in_y, min_rad=2e4, max_rad=1.75e4):
    sectors = get_sectors()
    top = ray(sectors['Dark Sector'][2:4])
    bottom = ray(sectors['Dark Sector'][0:2])
    new_x = []
    new_y = []
    for x, y in zip(in_x, in_y):
        if y > top(x):
            continue
        if y < bottom(x):
            continue
        if np.sqrt(x**2 + y**2) > max_rad:
            continue
        if np.sqrt(x**2 + y**2) < min_rad:
            continue
        new_x.append(x)
        new_y.append(y)
    return np.asarray(new_x), np.asarray(new_y)

def softtrim_to_fit(in_x, in_y, min_rad=2e4, max_rad=1.75e4, b0=np.array([0,0]), b1=np.array([0,0])):
    sectors = get_sectors()
    top = ray(sectors['Dark Sector'][2:4])
    bottom = ray(sectors['Dark Sector'][0:2])
    new_x = []
    new_y = []
    b0mag = np.sqrt((b0*b0).sum())
    b1mag = np.sqrt((b1*b1).sum())
    maxGap = max(b0mag, b1mag) + 1e-8
    for x, y in zip(in_x, in_y):
        if y > top(x):
            continue
        if y < bottom(x):
            continue
        if np.sqrt(x**2 + y**2) > max_rad + maxGap:
            continue
        if np.sqrt(x**2 + y**2) < min_rad:
            continue
        new_x.append(x)
        new_y.append(y)
    return np.asarray(new_x), np.asarray(new_y)
    
def trim_to_fit_edge(in_x, in_y, min_rad=2e4, max_rad=1.75e4, edge=0.):
    sectors = get_sectors()
    
    tds = np.array([[   929.3377603199988 - edge,  732.0277663199995 ],
                    [ -1230.6941085521246 - edge,  6666.666666666667 ]]).T
    inner_top_edge = ray(tds)
    
    bds = np.array([
        [ -6745.5344331313745 - edge , -6666.666666666667 ],
        [  360.24866210820073 - edge, -704.0174330560367 ]])
    inner_bottom_edge = ray(bds)

    top = ray(sectors['Dark Sector'][2:4])
    bottom = ray(sectors['Dark Sector'][0:2])
    
    new_x = []
    new_y = []
    for x, y in zip(in_x, in_y):
        if y > top(x+edge):
            continue
        if y < bottom(x+edge):
            continue
    
        if np.sqrt(x**2 + y**2) > max_rad:
            continue
        if np.sqrt(x**2 + y**2) < min_rad:
            continue
        new_x.append(x)
        new_y.append(y)
    return np.asarray(new_x), np.asarray(new_y)
    
def trim_to_fit_fill_edge(in_x, in_y, min_rad=2e4, max_rad=1.75e4, edge=0.):
    sectors = get_sectors()
    
    tds = np.array([[   929.3377603199988 - edge,  732.0277663199995 ],
                    [ -1230.6941085521246 - edge,  6666.666666666667 ]])
    inner_top_edge = ray(tds)
    
    bds = np.array([
        [ -6745.5344331313745 - edge , -6666.666666666667 ],
        [  360.24866210820073 - edge, -704.0174330560367 ]])
    inner_bottom_edge = ray(bds)

    top = ray(sectors['Dark Sector'][2:4])
    bottom = ray(sectors['Dark Sector'][0:2])
    
    new_x = []
    new_y = []
    for x, y in zip(in_x, in_y):
        if y > top(x):
            continue
        if y < bottom(x):
            continue
    
        if np.sqrt(x**2 + y**2) > max_rad:
            continue
        if np.sqrt(x**2 + y**2) < min_rad:
            # if you're in the band
            if  ( y < 0 and y > inner_bottom_edge(x) ) or (y >= 0 and y < inner_top_edge(x)):
                continue
        
        new_x.append(x)
        new_y.append(y)
    return np.asarray(new_x), np.asarray(new_y)

def softtrim_to_radius(in_x, in_y, min_rad=2e4, max_rad=1.75e4, b0=np.array([0,0]), b1=np.array([0,0])):
    sectors = get_sectors()
    top = ray(sectors['Dark Sector'][2:4])
    maxY = sectors['Dark Sector'][2:4,1].max()
    new_x = []
    new_y = []
    b0mag = np.sqrt((b0*b0).sum())
    b1mag = np.sqrt((b1*b1).sum())
    maxGap = max(b0mag, b1mag) + 1e-6
    for x, y in zip(in_x, in_y):
        if y > top(x) and y < maxY:
            continue
        if np.sqrt(x**2 + y**2) > max_rad + maxGap:
            continue
        if np.sqrt(x**2 + y**2) < min_rad:
            continue
        new_x.append(x)
        new_y.append(y)
    return np.asarray(new_x), np.asarray(new_y)

def get_shallow_ring(b0, b1, xdeep, ydeep, xshallow, yshallow, offset=np.array([0,0])):

    xring = []
    yring = []

    deep_pts = np.array([xdeep, ydeep]).T
    shallow_pts = np.array([xshallow, yshallow]).T
    b0mag = np.sqrt((b0*b0).sum())
    b1mag = np.sqrt((b1*b1).sum())
    minGap = min(b0mag, b1mag) - 1e-6
    maxGap = max(b0mag, b1mag) + 1e-6

    for point in deep_pts:

        for i in range(-5, 5):
            for j in range(-5, 5):

                candidate = point + i*b0 + j*b1 + offset
                rcandidate = np.linalg.norm(candidate - point)
                if(rcandidate > maxGap):
                    continue
                rshallow = np.linalg.norm(candidate - shallow_pts, axis=1)
                rdeep = np.linalg.norm(candidate - deep_pts, axis=1)
                ring_pts = np.array([xring, yring]).T
                rring = np.linalg.norm(candidate - ring_pts, axis=1)
                if( (not np.any(rshallow < minGap/2)) and (not np.any(rdeep < 1e-6)) 
                        and np.any(rdeep < maxGap) and (not np.any(rring < minGap)) ):
                    xring.append(candidate[0])
                    yring.append(candidate[1]) 

    return np.asarray(xring), np.asarray(yring)

# adjust stations so that none would require a ring station in clean air sector
def adjust_stations(b0, b1, x, y, offset=np.array([0,0]), b0ring=None, b1ring=None):

    sectors = get_sectors()
    top = ray(sectors['Dark Sector'][2:4])
    maxY = sectors['Dark Sector'][2:4,1].max()

    if(b0ring is None):
        b0ring = b0
    if(b1ring is None):
        b1ring = b1

    b0ringmag = np.sqrt((b0ring*b0ring).sum())
    b1ringmag = np.sqrt((b1ring*b1ring).sum())
    minGap = min(b0ringmag, b1ringmag) - 1e-6
    maxGap = max(b0ringmag, b1ringmag) + 1e-6
    
    pts = np.array([x,y]).T
    
    for idx, point in enumerate(pts):
        nBadRing = 0
        for i in range(-5, 5):
            for j in range(-5, 5):
                
                # check where possible ring stations would be
                candidate = point + i*b0ring + j*b1ring + offset
                rcandidate = np.linalg.norm(candidate - point)
                if(rcandidate > maxGap):
                    continue
                xx, yy = candidate
                if yy > top(xx) and yy < maxY:
                    nBadRing += 1

        xcandidate = []
        ycandidate = []
        if nBadRing > 0:
            for i in range(-10, 10): 
                for j in range(-10, 10):
                    # find new candidate locations
                    candidate = point + i*b0 + j*b1
                    xx, yy = candidate
                    if yy < maxY:
                        continue
                    #if yy > top(xx) :
                    #    continue
                    r = np.linalg.norm(candidate - pts, axis=1)
                    if np.any(r < minGap):
                        continue
                    xcandidate.append(xx)
                    ycandidate.append(yy)
            candpts = np.array([xcandidate, ycandidate]).T
            boundary = top(xcandidate)
            #r = np.linalg.norm(point - candpts, axis=1)
            r = (ycandidate - boundary)
            imin = r.argmin()
            #pts[idx] = candpts[imin] # moves station to closest possible candidate location
            pts[idx] = candpts[imin] # moves station to candidate location to left of boundary (or as close as possible) --- this could be made better

    x = pts[:,0]
    y = pts[:,1]
    return x, y


def get_ARA():
    # A1: [-2346.96, -345.6432]
    # A2: [-3344.2656, -2078.1264]
    # A3: [-4344.6192, -345.948]
    # A4: [-3388.7664, 1437.4368]
    # A5: [-4342.7904, -3804.5136]
    locations_x = np.asarray([-2346.96, -3344.2656, -4344.6192, -3388.7664, -4342.7904])
    locations_y = np.asarray([-345.6432, -2078.1264, -345.948, 1437.4368, -3804.5136])
    return locations_x, locations_y


