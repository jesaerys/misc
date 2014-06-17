"""
astropy.wcs already handles these transformations, but I spent some time
figuring out how they work, so here they are! They are mostly limited to
gnomonic (TAN) projections with the celestial coordinates in deg.

"""
def pix2proj(x, y, hdr):
    """Convert pixel coordinates into projection plane coordinates (deg)
    according to Greisen & Calabretta (2002).

    Required header keywords:
    CRPIX1, CRPIX2
    CD1_1, CD1_2, CD2_1, CD2_2

    CDi_j are assumed to have units deg/pix.

    """
    # (G&C02 eq.3/pg.1063/pdf.3)
    u = hdr['cd1_1']*(x-hdr['crpix1']) + hdr['cd1_2']*(y-hdr['crpix2'])
    v = hdr['cd2_1']*(x-hdr['crpix1']) + hdr['cd2_2']*(y-hdr['crpix2'])

    return u, v


def proj2pix(u, v, hdr):
    """Convert projection plane coordinates (deg) into pixel coordinates
    according to Greisen & Calabretta (2002).

    Required header keywords:
    CRPIX1, CRPIX2
    CD1_1, CD1_2, CD2_1, CD2_2

    CDi_j are assumed to have units deg/pix.

    """
    # (G&C02 eq.3/pg.1063/pdf.3)
    c = hdr['cd1_1']*hdr['cd2_2'] - hdr['cd1_2']*hdr['cd2_1']
    x = (hdr['cd2_2']*u - hdr['cd1_2']*v) / c + hdr['crpix1']
    y = (-hdr['cd2_1']*u + hdr['cd1_1']*v) / c + hdr['crpix2']

    return x, y


def proj2natsph(u, v):
    """Convert projection plane coordinates (deg) into native spherical
    coordinates (deg) assuming a gnomonic (TAN) projection according to
    Calabretta & Greisen  (2002).

    CTYPE1 and CTYPE2 are assumed to end with 'TAN'.

    """
    u = u * np.pi/180
    v = v * np.pi/180
    p = np.arctan2(u, -v)  # (C&G02 eq.14/pg.1085/pdf.9)
    r = np.sqrt(u**2 + v**2)  # (C&G02 eq.15/pg.1085/pdf.9)
    t = np.arctan(1/r)  # (C&G02 eq.55/pg.1088/pdf.12)
    p *= 180/np.pi
    t *= 180/np.pi

    return p, t


def natsph2proj(p, t):
    """Convert native spherical coordinates (deg) into projection plane
    coordinates (deg) assuming a gnomonic (TAN) projection according to
    Calabretta & Greisen  (2002).

    CTYPE1 and CTYPE2 are assumed to end with 'TAN'.

    """
    p = p * np.pi/180
    t = t * np.pi/180
    r = 1/np.tan(t)  # (C&G02 eq.54/pg.1088/pdf.12)
    u = r*np.sin(p)  # (C&G02 eq.12/pg.1085/pdf.9)
    v = -r*np.cos(p)  # (C&G02 eq.13/pg.1085/pdf.9)
    u *= 180/np.pi
    v *= 180/np.pi

    return u, v


def natsph2celsph(p, t, hdr):
    """Convert native spherical coordinates (deg) into celestial spherical
    coordinates (deg) assuming a gnomonic (TAN) projection according to
    Calabretta & Greisen (2002). The celestial coordinate system is
    determined by CTYPE1 and CTYPE2. For example, if CTYPE1 and CTYPE2 are
    'RA---TAN' and 'DEC--TAN', then a and d correspond to RA and dec (as do
    CRVAL1 and CRVAL2). The units of CRVAL1 and CRVAL2 are assumed to be
    deg.

    Required header keywords:
    CRVAL1, CRVAL2

    CTYPE1 and CTYPE2 are assumed to end with 'TAN', and CRVAL1 and CRVAL2
    are assumed to have units deg.

    """
    # TAN-specific settings
    ap, dp = hdr['crval1'], hdr['crval2']
    pp = 180.0

    # (C&G02 eq.2/pg.1079/pdf.3)
    p = p * np.pi/180
    t = t * np.pi/180
    ap *= np.pi/180
    dp *= np.pi/180
    pp *= np.pi/180
    a = ap + np.arctan2(-np.cos(t)*np.sin(p-pp),
                        np.sin(t)*np.cos(dp)-np.cos(t)*np.sin(dp)*np.cos(p-pp))
    d = np.arcsin(np.sin(t)*np.sin(dp) + np.cos(t)*np.cos(dp)*np.cos(p-pp))
    a *= 180/np.pi
    d *= 180/np.pi

    return (a, d)


def celsph2natsph(a, d, hdr):
    """Convert celestial spherical coordinates (deg) into native spherical
    coordinates (deg) assuming a gnomonic (TAN) projection according to
    Calabretta & Greisen (2002). The celestial coordinate system is
    determined by CTYPE1 and CTYPE2. For example, if CTYPE1 and CTYPE2 are
    'RA---TAN' and 'DEC--TAN', then a and d correspond to RA and dec (as do
    CRVAL1 and CRVAL2). The units of CRVAL1 and CRVAL2 are assumed to be
    deg.

    Required header keywords:
    CRVAL1, CRVAL2

    CTYPE1 and CTYPE2 are assumed to end with 'TAN', and CRVAL1 and CRVAL2
    are assumed to have units deg.

    """
    # TAN-specific settings
    ap, dp = hdr['crval1'], hdr['crval2']
    pp = 180.0

    # (C&G02 eq.5/pg.1080/pdf.4)
    a = a * np.pi/180
    d = d * np.pi/180
    ap *= np.pi/180
    dp *= np.pi/180
    pp *= np.pi/180
    p = pp + np.arctan2(-np.cos(d)*np.sin(a-ap),
                        np.sin(d)*np.cos(dp)-np.cos(d)*np.sin(dp)*np.cos(a-ap))
    t = np.arcsin(np.sin(d)*np.sin(dp) + np.cos(d)*np.cos(dp)*np.cos(a-ap))
    p *= 180/np.pi
    t *= 180/np.pi

    return p, t


def pix2world(x, y, hdr):
    """Convert pixel coordinates into RA and dec assuming a gnomonic (TAN)
    projection accoring to Calabretta & Greisen (2002).

    Required keywords in ``hdr``:
    CRPIX1, CRPIX2
    CRVAL1, CRVAL2
    CD1_1, CD1_2, CD2_1, CD2_2

    CTYPE1 and CTYPE2 are assumed to be 'RA---TAN' and 'DEC--TAN'. CDi_j
    are assumed to have units deg/pix. CRVAL1 and CRVAL2 are assumed to
    have units deg.

    """
    u, v = pix2proj(x, y, hdr)
    p, t = proj2natsph(u, v)
    a, d = natsph2celsph(p, t, hdr)
    return a, d


def world2pix(a, d, hdr):
    """Convert RA and dec into pixel coordinates assuming a gnomonic (TAN)
    projection accoring to Calabretta & Greisen (2002).

    Required keywords in ``hdr``:
    CRPIX1, CRPIX2
    CRVAL1, CRVAL2
    CD1_1, CD1_2, CD2_1, CD2_2

    CTYPE1 and CTYPE2 are assumed to be 'RA---TAN' and 'DEC--TAN'. CDi_j
    are assumed to have units deg/pix. CRVAL1 and CRVAL2 are assumed to
    have units deg.

    """
    p, t = celsph2natsph(a, d, hdr)
    u, v = natsph2proj(p, t)
    x, y = proj2pix(u, v, hdr)
    return x, y
