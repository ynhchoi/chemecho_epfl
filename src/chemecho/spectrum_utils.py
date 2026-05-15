# Priority order for yunits when multiple spectra are available.
# extract_spectrum_data in get_spectrum.py iterates this list in order.
YUNITS_PRIORITY = [
    'TRANSMITTANCE',               # 1st: direct use, no conversion
    'ABSORBANCE',                  # 2nd: T = 10^(-A) * 100 (Beer-Lambert)
    'MOLAR_ABSORPTIVITY',          # 3rd: same log scale as absorbance (approximation)
]


def _transmittance_fractional(data_y: list) -> list:
    """Convert fractional transmittance (0–1) to % (0–100).

    Args:
        data_y: list of transmittance values in 0–1 scale
    Returns:
        list of transmittance values in %
    """
    return [y * 100 for y in data_y]


def _transmittance_percent(data_y: list) -> list:
    """Return transmittance already in % as-is.

    Args:
        data_y: list of transmittance values in % scale
    Returns:
        list of transmittance values in % (unchanged)
    """
    return list(data_y)


def _absorbance_to_transmittance(data_y: list) -> list:
    """Convert absorbance to transmittance (%).

    Formula: T = 10^(-A) * 100  (Beer-Lambert law)

    Args:
        data_y: list of absorbance values
    Returns:
        list of transmittance values in %
    """
    return [10 ** (-a) * 100 for a in data_y]


def _molar_absorptivity_to_transmittance(data_y: list) -> list:
    """Approximate conversion of molar absorptivity to transmittance (%).

    NIST quantitative IR DB stores values as (micromol/mol)-1m-1 (base 10),
    which is a log-base-10 absorption coefficient. Without concentration and
    path length, an exact conversion is impossible, but peak positions and
    relative intensities are preserved by treating values as absorbance-like.

    Formula (approximation): T ≈ 10^(-α) * 100

    Args:
        data_y: list of molar absorptivity values
    Returns:
        list of approximate transmittance values in %
    """
    return [10 ** (-a) * 100 for a in data_y]


def _detect_transmittance_scale(data_y: list) -> str:
    """Detect whether transmittance values are fractional (0–1) or percent (0–100).

    NIST reports both as yunits = 'TRANSMITTANCE', so we distinguish by value range.

    Args:
        data_y: list of transmittance values
    Returns:
        'fractional' if max <= 1.5, 'percent' otherwise
    """
    return "fractional" if max(data_y) <= 1.5 else "percent"


def classify_yunits(yunits: str) -> str:
    """Map a raw NIST yunits string to one of the keys in YUNITS_PRIORITY.

    Args:
        yunits: raw yunits string from JCAMP-DX header (any case)
    Returns:
        one of 'TRANSMITTANCE', 'ABSORBANCE', 'MOLAR_ABSORPTIVITY',
        or 'UNKNOWN' if not recognised
    """
    yunits_upper = yunits.upper() if yunits else ''

    if 'TRANSMITTANCE' in yunits_upper:
        return 'TRANSMITTANCE'
    if 'ABSORBANCE' in yunits_upper:
        return 'ABSORBANCE'
    if 'MICROMOL' in yunits_upper or 'MOL' in yunits_upper:
        return 'MOLAR_ABSORPTIVITY'
    return 'UNKNOWN'


def to_transmittance(data_y: list, yunits: str) -> list:
    """Convert any supported NIST IR y-axis format to transmittance (%).

    Supported formats (in priority order):
        1. TRANSMITTANCE  → direct use (% or fractional auto-detected)
        2. ABSORBANCE     → T = 10^(-A) * 100
        3. MOLAR_ABSORPTIVITY → T ≈ 10^(-α) * 100  (approximation)

    No normalisation applied — absolute values are preserved.

    Args:
        data_y: list of raw y-axis values from NIST JCAMP-DX
        yunits: y-axis unit string from JCAMP-DX header (case-insensitive)
    Returns:
        list of transmittance values in %
    Raises:
        ValueError: if yunits is not a recognised format

    Example:
        >>> to_transmittance([1.0, 0.5, 0.0], 'ABSORBANCE')
        [10.0, 31.62..., 100.0]
    """
    kind = classify_yunits(yunits)

    if kind == 'TRANSMITTANCE':
        scale = _detect_transmittance_scale(data_y)
        if scale == "fractional":
            return _transmittance_fractional(data_y)
        else:
            return _transmittance_percent(data_y)

    if kind == 'ABSORBANCE':
        return _absorbance_to_transmittance(data_y)

    if kind == 'MOLAR_ABSORPTIVITY':
        return _molar_absorptivity_to_transmittance(data_y)

    raise ValueError(
        f"Unrecognised y-unit format: '{yunits}'. "
        f"Supported: TRANSMITTANCE, ABSORBANCE, (micromol/mol)-1m-1 (base 10)."
    )


def to_wavenumber(data_x: list, xunits: str) -> list:
    """Convert any NIST IR x-axis format to wavenumbers (cm⁻¹).

    Handles two formats:
        1. Already in cm⁻¹  (1/CM, cm-1)
        2. Micrometers (μm) → wavenumber: ν = 10000 / λ

    Args:
        data_x: list of raw x-axis values from NIST JCAMP-DX
        xunits: x-axis unit string from JCAMP-DX header (case-insensitive)
    Returns:
        list of wavenumbers in cm⁻¹
    Raises:
        ValueError: if xunits is not a recognised format

    Example:
        >>> to_wavenumber([10.0, 5.0, 2.5], 'MICROMETERS')
        [1000.0, 2000.0, 4000.0]
    """
    xunits_upper = xunits.upper() if xunits else ''

    if 'MICRON' in xunits_upper or 'MICROMETER' in xunits_upper:
        return [10000 / x for x in data_x]

    if '1/CM' in xunits_upper or 'CM' in xunits_upper:
        return list(data_x)

    raise ValueError(
        f"Unrecognised x-unit format: '{xunits}'. "
        "Supported: 1/CM, MICROMETERS."
    )
