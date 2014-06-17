def sparea(lon, lat, R=1, units='deg'):
    """Calculate the area of a spherical polygon.

    It is assumed that the polygon contains no poles! The method is adapted
    from "Some algorithms for polygons on a sphere" by Chamberlain and
    Duquette (http://trs-new.jpl.nasa.gov/dspace/handle/2014/40409)

    Parameters
    ----------
    lon, lat : array-like
        Longitude and latitude of each vertex in the polygon. The first
        vertex may be listed either only at the beginning of the sequence,
        or both at the beginning and at the end of the sequence.
    R : float, optional
        Radius of the sphere (default is 1).
    units : {'deg', 'rad'}, optional
        Sets the unit system for `lon` and `lat`.

    Returns
    -------
    float
        Total area of the spherical polygon.

    """
    lon, lat = np.array(lon), np.array(lat)

    # Ensure the polygon is closed
    if not lon[-1]==lon[0] or not lat[-1]==lat[0]:
        lon, lat = np.append(lon, lon[0]), np.append(lat, lat[0])

    if units == 'deg':
        lon, lat = lon * np.pi / 180, lat * np.pi / 180

    # Great circle segments between vertices
    l = gcdist(lon[:-1], lat[:-1], lon[1:], lat[1:], deg=False)

    # Semiperimeter of each spherical triangle
    s = 0.5 * (l + np.pi + lat[:-1] + lat[1:])

    # Spherical excess of each spherical triangle from L'Huilier's theorem.
    # Note that none of the terms should be negative (not 100% sure about
    # that); assume that any negative values are within machine precision
    # of 0.
    term1 = (s - (np.pi / 2 + lat[:-1])) /2
    term1[term1<0] = 0
    term2 = (s - (np.pi / 2 + lat[1:])) /2
    term2[term2<0] = 0
    result = np.tan(s/2) * np.tan((s-l)/2) * np.tan(term1) * np.tan(term2)
    E = 4 * np.arctan(np.sqrt(result))

    # Let A<0 for lon[i]<lon[i+1], A>0 for lon[i+1]<lon[i] assuming ccw
    # traversal
    sign = 2*(lon[1:] < lon[:-1]) - 1

    # Total area
    A = np.sum(sign * E) * R**2
    if units == 'deg':
        A = A * (180 / np.pi)**2
    if A < 0:  # Fix the sign in case the vertices are not listed ccw
        A = -A
    return A
