import numpy as np


def parse_ctype(ctype_str):
    """Get the coordinate system and projection for a given value of the
    CTYPEi keyword.

    """
    elements = ctype_str.split('-')
    coordsys, projection = elements[0], elements[-1]
    return coordsys, projection


def pix2proj(x, y, hdr):
    """Convert pixel coordinates into projection plane coordinates
    according to Eq. 3 in Greisen & Calabretta (2002).

    Parameters
    ----------
    x, y : float or array
        x and y pixel coordinates.
    hdr : astropy.io.fits.Header or dictionary
        A FITS header. Required keywords:

        - CRPIX1, CRPIX2
        - CD1_1, CD1_2, CD2_1, CD2_2

        Optional keywords:

        - CUNIT1, CUNIT2: Used for converting the CD matrix into deg/pix.
          If omitted, the CD matrix is assumed to already have units
          deg/pix so that no conversion is necessary.

    Returns
    -------
    float or array
        x and y projection plane coordinates in degrees.

    """
    cunit1, cunit2 = hdr.get('CUNIT1'), hdr.get('CUNIT2')
    if cunit1 is not None:
        pass  # Always assume degrees for now
    if cunit2 is not None:
        pass  # Always assume degrees for now

    # (G&C02 eq.3/pg.1063/pdf.3)
    xp = hdr['CD1_1']*(x-hdr['CRPIX1']) + hdr['CD1_2']*(y-hdr['CRPIX2'])
    yp = hdr['CD2_1']*(x-hdr['CRPIX1']) + hdr['CD2_2']*(y-hdr['CRPIX2'])

    return xp, yp


def proj2pix(xp, yp, hdr):
    """Convert projection plane coordinates into pixel coordinates
    according to Eq. 3 in Greisen & Calabretta (2002).

    Parameters
    ----------
    xp, yp : float or array
        x and y projection plane coordinates in degrees.
    hdr : astropy.io.fits.Header or dictionary
        A FITS header. Required keywords:

        - CRPIX1, CRPIX2
        - CD1_1, CD1_2, CD2_1, CD2_2

        Optional keywords:

        - CUNIT1, CUNIT2: Used for converting the CD matrix into deg/pix.
          If omitted, the CD matrix is assumed to already have units
          deg/pix so that no conversion is necessary.

    Returns
    -------
    float or array
        x and y pixel coordinates.

    """
    cunit1, cunit2 = hdr.get('CUNIT1'), hdr.get('CUNIT2')
    if cunit1 is not None:
        pass  # Always assume degrees for now
    if cunit2 is not None:
        pass  # Always assume degrees for now

    # (G&C02 eq.3/pg.1063/pdf.3)
    c = hdr['CD1_1']*hdr['CD2_2'] - hdr['CD1_2']*hdr['CD2_1']
    x = (hdr['CD2_2']*xp - hdr['CD1_2']*yp) / c + hdr['CRPIX1']
    y = (-hdr['CD2_1']*xp + hdr['CD1_1']*yp) / c + hdr['CRPIX2']

    return x, y


def proj2natsph(xp, yp, hdr):
    """Convert projection plane coordinates into native spherical
    coordinates for a given projection according to Calabretta & Greisen
    (2002).

    Parameters
    ----------
    xp, yp : float or array
        x and y projection plane coordinates in degrees.
    hdr : astropy.io.fits.Header or dictionary
        A FITS header. Required keywords:

        - CTYPE1, CTYPE2

    Returns
    -------
    float or array
        Native spherical longitude (phi) and latitude (theta) in degrees
        for the projection specified by CTYPE1 and CTYPE2.

    """
    # Degrees to radians
    xp = xp * np.pi/180
    yp = yp * np.pi/180

    projection1 = parse_ctype(hdr['CTYPE1'])[1]
    projection2 = parse_ctype(hdr['CTYPE2'])[1]
    projection = projection1  # Always assume projection1 and projection2 are the same?

    if projection == 'TAN':
        phi = np.arctan2(xp, -yp)  # (C&G02 eq.14/pg.1085/pdf.9)
        r = np.sqrt(xp**2 + yp**2)  # (C&G02 eq.15/pg.1085/pdf.9)
        theta = np.arctan(1/r)  # (C&G02 eq.55/pg.1088/pdf.12)
    else:
        # Raise an error?
        phi, theta = None, None

    # Radians to degrees
    phi *= 180/np.pi
    theta *= 180/np.pi

    return phi, theta


def natsph2proj(phi, theta, hdr):
    """Convert native spherical coordinates into projection plane
    coordinates for a given projection according to Calabretta & Greisen
    (2002).

    Parameters
    ----------
    phi, theta : float or array
        Native spherical longitude (phi) and latitude (theta) in degrees
        for the projection specified by CTYPE1 and CTYPE2.
    hdr : astropy.io.fits.Header or dictionary
        A FITS header. Required keywords:

        - CTYPE1, CTYPE2

    Returns
    -------
    float or array
        x and y projection plane coordinates in degrees.

    """
    # Degrees to radians
    phi = phi * np.pi/180
    theta = theta * np.pi/180

    projection1 = parse_ctype(hdr['CTYPE1'])[1]
    projection2 = parse_ctype(hdr['CTYPE2'])[1]
    projection = projection1  # Always assume projection1 and projection2 are the same?

    if projection == 'TAN':
        r = 1/np.tan(theta)  # (C&G02 eq.54/pg.1088/pdf.12)
        xp = r*np.sin(phi)  # (C&G02 eq.12/pg.1085/pdf.9)
        yp = -r*np.cos(phi)  # (C&G02 eq.13/pg.1085/pdf.9)
    else:
        # Raise and error?
        xp, yp = None, None

    # Radians to degrees
    xp *= 180/np.pi
    yp *= 180/np.pi

    return xp, yp


def natsph2celsph(phi, theta, hdr):
    """Convert native spherical coordinates for a given projection into the
    given celestial spherical coordinates according to Calabretta & Greisen
    (2002).

    Parameters
    ----------
    phi, theta : float or array
        Native spherical longitude (phi) and latitude (theta) in degrees
        for the projection specified by CTYPE1 and CTYPE2.
    hdr : astropy.io.fits.Header or dictionary
        A FITS header. Required keywords:

        - CTYPE1, CTYPE2
        - CRVAL1, CRVAL2

        Optional keywords:

        - CUNIT1, CUNIT2: Used for converting CRVAL1 and CRVAL2 into
          degrees and converting celestial longitude and latitude into the
          proper units. If omitted, the everything is assumed to be in
          degrees.

    Returns
    -------
    float or array
        Celestial spherical longitude (lon) and latitude (lat) for the
        celestial coordinate system specified by CTYPE1 and CTYPE2.

    """
    cunit1, cunit2 = hdr.get('CUNIT1'), hdr.get('CUNIT2')
    if cunit1 is not None:
        pass  # Always assume CRVAL1 and longitude are in degrees for now
    if cunit2 is not None:
        pass  # Always assume CRVAL2 and latitude are in degrees for now

    # Degrees to radians
    phi = phi * np.pi/180
    theta = theta * np.pi/180

    coordsys1, projection1 = parse_ctype(hdr['CTYPE1'])
    coordsys2, projection2 = parse_ctype(hdr['CTYPE2'])
    projection = projection1  # Always assume projection1 and projection2 are the same?

    if projection == 'TAN':
        lon_p = hdr['CRVAL1'] * np.pi/180
        lat_p = hdr['CRVAL2'] * np.pi/180
        phi_p = 180.0 * np.pi/180
    else:
        # Raise an error?
        lon_p, lat_p, phi_p = None, None, None

    # (C&G02 eq.2/pg.1079/pdf.3)
    dphi = phi - phi_p
    lon = lon_p + np.arctan2(-np.cos(theta) * np.sin(dphi),
                             np.sin(theta) * np.cos(lat_p) -
                             np.cos(theta) * np.sin(lat_p) * np.cos(dphi))
    lat = np.arcsin(np.sin(theta) * np.sin(lat_p) +
                    np.cos(theta) * np.cos(lat_p) * np.cos(dphi))

    # Radians to degrees
    lon *= 180/np.pi
    lat *= 180/np.pi

    return (lon, lat)


def celsph2natsph(lon, lat, hdr):
    """Convert the given celestial spherical coordinates into native
    spherical coordinates for a given projection according to Calabretta &
    Greisen (2002).

    Parameters
    ----------
    lon, lat : float or array
        Celestial spherical longitude (lon) and latitude (lat) for the
        celestial coordinate system specified by CTYPE1 and CTYPE2.
    hdr : astropy.io.fits.Header or dictionary
        A FITS header. Required keywords:

        - CTYPE1, CTYPE2
        - CRVAL1, CRVAL2

        Optional keywords:

        - CUNIT1, CUNIT2: Used for converting CRVAL1, CRVAL2, and celestial
          longitude and latitude into degrees. If omitted, the everything
          is assumed to be in degrees.

    Returns
    -------
    float or array
        Native spherical longitude (phi) and latitude (theta) in degrees
        for the projection specified by CTYPE1 and CTYPE2.

    """
    cunit1, cunit2 = hdr.get('CUNIT1'), hdr.get('CUNIT2')
    if cunit1 is not None:
        pass  # Always assume CRVAL1 and longitude are in degrees for now
    if cunit2 is not None:
        pass  # Always assume CRVAL2 and latitude are in degrees for now

    # Degrees to radians
    lon = lon * np.pi/180
    lat = lat * np.pi/180

    coordsys1, projection1 = parse_ctype(hdr['CTYPE1'])
    coordsys2, projection2 = parse_ctype(hdr['CTYPE2'])
    projection = projection1  # Always assume projection1 and projection2 are the same?

    if projection == 'TAN':
        lon_p = hdr['CRVAL1'] * np.pi/180
        lat_p = hdr['CRVAL2'] * np.pi/180
        phi_p = 180.0 * np.pi/180
    else:
        # Raise an error?
        lon_p, lat_p, phi_p = None, None, None

    # (C&G02 eq.5/pg.1080/pdf.4)
    dlon = lon - lon_p
    phi = phi_p + np.arctan2(-np.cos(lat) * np.sin(dlon),
                             np.sin(lat) * np.cos(lat_p) -
                             np.cos(lat) * np.sin(lat_p) * np.cos(dlon))
    theta = np.arcsin(np.sin(lat) * np.sin(lat_p) +
                      np.cos(lat) * np.cos(lat_p) * np.cos(dlon))

    # Radians to degrees
    phi *= 180/np.pi
    theta *= 180/np.pi

    return phi, theta


def pix2world(x, y, hdr):
    """Convert pixel coordinates into celestial coordinates accoring to
    Calabretta & Greisen (2002).

    .. note:: This transformation is already supported by the `astropy.wcs`
       module.

    Parameters
    ----------
    x, y : float or array
        x and y pixel coordinates.
    hdr : astropy.io.fits.Header or dictionary
        A FITS header. Required keywords:

        - CTYPE1, CTYPE2
        - CRPIX1, CRPIX2
        - CRVAL1, CRVAL2
        - CD1_1, CD1_2, CD2_1, CD2_2

        Optional keywords:

        - CUNIT1, CUNIT2: Used for converting the CD matrix into deg/pix,
          CRVAL1 and CRVAL2 into degrees, and celestial longitude and
          latitude into the proper units. If omitted, everything is assumed
          to be in degrees (deg/pix for the CD matrix).

    Returns
    -------
    float or array
        Celestial spherical longitude (lon) and latitude (lat) for the
        celestial coordinate system specified by CTYPE1 and CTYPE2.

    """
    xp, yp = pix2proj(x, y, hdr)
    phi, theta = proj2natsph(xp, yp, hdr)
    lon, lat = natsph2celsph(phi, theta, hdr)
    return lon, lat


def world2pix(lon, lat, hdr):
    """Convert celestial coordinates into pixel coordinates accoring to
    Calabretta & Greisen (2002).

    .. note:: This transformation is already supported by the `astropy.wcs`
       module.

    Parameters
    ----------
    lon, lat : float or array
        Celestial spherical longitude (lon) and latitude (lat) for the
        celestial coordinate system specified by CTYPE1 and CTYPE2.
    hdr : astropy.io.fits.Header or dictionary
        A FITS header. Required keywords:

        - CTYPE1, CTYPE2
        - CRPIX1, CRPIX2
        - CRVAL1, CRVAL2
        - CD1_1, CD1_2, CD2_1, CD2_2

        Optional keywords:

        - CUNIT1, CUNIT2: Used for converting the CD matrix into deg/pix
          and CRVAL1, CRVAL2, and celestial longitude and latitude into
          degrees. If omitted, everything is assumed to be in degrees
          (deg/pix for the CD matrix).

    Returns
    -------
    float or array
        x and y pixel coordinates.

    """
    phi, theta = celsph2natsph(lon, lat, hdr)
    xp, yp = natsph2proj(phi, theta, hdr)
    x, y = proj2pix(xp, yp, hdr)
    return x, y
