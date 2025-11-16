import time
import vtk
import math
import numpy as np
import pyvista as pv
from pyvista import examples
import random
import logging
from scipy.ndimage import gaussian_filter
from scipy.optimize import newton
from datetime import datetime, timezone, timedelta
from scipy.spatial.transform import Rotation as R


pv.global_theme.show_scalar_bar = False
pv.global_theme.volume_mapper = 'smart'  


class NoScalarVolumeMapper(vtk.vtkSmartVolumeMapper):
    def SetScalarVisibility(self, visibility):
        super().SetScalarVisibility(0)  

pv.global_theme.allow_empty_mesh = True


# Logging & warnings
logging.basicConfig(level=logging.ERROR)
import warnings
warnings.filterwarnings("ignore", message=".*pickpoint.*")
random.seed(42)
pl = pv.Plotter(window_size=[1920, 1080])
pv.global_theme.show_scalar_bar = False
pv.global_theme.title = ""                  
pv.global_theme.axes.show = False




SIM_START_UTC = datetime.now(timezone.utc)
SIM_TIME = SIM_START_UTC.timestamp()         
SIM_SECONDS_PER_REAL_SECOND = 1

# ---------------------------- SYSTEM CONSTANTS ----------------------------
GLOBAL_SCALE_MULTIPLIER = 1  # Reduced for better visualization
DISTANCE_KM_PER_UNIT = 1.0e5
SIZE_KM_PER_UNIT = DISTANCE_KM_PER_UNIT * 0.02
KM_TO_SCENE = 1.0 / DISTANCE_KM_PER_UNIT
SIZEKM_TO_SCENE = 1.0 / SIZE_KM_PER_UNIT
J2000 = 2451545.0
MIN_ZOOM_DISTANCE = 5

# Moon counts
TOTAL_MOONS = 95  # Jupiter
INNER_MOON_COUNT = 15
OUTER_MOON_COUNT = TOTAL_MOONS - INNER_MOON_COUNT
NEPTUNE_MOON_COUNT = 14
SATURN_MOON_COUNT = 83  # Reduced to known moons for performance
URANUS_MOON_COUNT = 27
MARS_MOON_COUNT = 2
EARTH_MOON_COUNT = 1

AU_KM = 149597870.7
AU_SCALE = AU_KM * KM_TO_SCENE
DISTANCE_SCALE = 149.6e6 / AU_SCALE

SPEED_SCALE = DISTANCE_SCALE ** 1.5 * 1e-80

# Radii (scaled relative to Jupiter)
REAL_SUN_RADIUS_KM = 696340.0
REAL_EARTH_RADIUS_KM = 6371.0
REAL_MARS_RADIUS_KM = 3390.0
REAL_JUPITER_RADIUS_KM = 71492.0
REAL_SATURN_RADIUS_KM = 58232.0
REAL_URANUS_RADIUS_KM = 25559.0
REAL_NEPTUNE_RADIUS_KM = 24622.0
REAL_MERCURY_RADIUS_KM = 2440.0
REAL_VENUS_RADIUS_KM = 6052.0
REAL_PLUTO_RADIUS_KM = 1188.3
REAL_ERIS_RADIUS_KM = 1163.0
REAL_HAUMEA_RADIUS_KM = 816.0  
REAL_MAKEMAKE_RADIUS_KM = 715.0
REAL_CERES_RADIUS_KM = 473.0

EARTH_AU = 1.0
MERCURY_AU = 0.387
VENUS_AU = 0.723
MARS_AU = 1.524
JUPITER_AU = 5.2
SATURN_AU = 9.58
URANUS_AU = 19.22
NEPTUNE_AU = 30.07
PLUTO_AU = 39.48
ERIS_AU = 67.78
HAUMEA_AU = 43.13
MAKEMAKE_AU = 45.79
CERES_AU = 2.77

SUN_RADIUS = REAL_SUN_RADIUS_KM * SIZEKM_TO_SCENE
EARTH_RADIUS = REAL_EARTH_RADIUS_KM * SIZEKM_TO_SCENE
MARS_RADIUS = REAL_MARS_RADIUS_KM * SIZEKM_TO_SCENE
JUPITER_RADIUS = REAL_JUPITER_RADIUS_KM * SIZEKM_TO_SCENE
SATURN_RADIUS = REAL_SATURN_RADIUS_KM * SIZEKM_TO_SCENE
URANUS_RADIUS = REAL_URANUS_RADIUS_KM * SIZEKM_TO_SCENE
NEPTUNE_RADIUS = REAL_NEPTUNE_RADIUS_KM * SIZEKM_TO_SCENE
MERCURY_RADIUS = REAL_MERCURY_RADIUS_KM * SIZEKM_TO_SCENE
VENUS_RADIUS = REAL_VENUS_RADIUS_KM * SIZEKM_TO_SCENE
PLUTO_RADIUS = REAL_PLUTO_RADIUS_KM * SIZEKM_TO_SCENE
ERIS_RADIUS = REAL_ERIS_RADIUS_KM * SIZEKM_TO_SCENE
HAUMEA_RADIUS = REAL_HAUMEA_RADIUS_KM * SIZEKM_TO_SCENE
MAKEMAKE_RADIUS = REAL_MAKEMAKE_RADIUS_KM * SIZEKM_TO_SCENE
CERES_RADIUS = REAL_CERES_RADIUS_KM * SIZEKM_TO_SCENE


# Orbital distances (approximate, but positions now computed dynamically)
EARTH_ORBIT_RADIUS = EARTH_AU * AU_SCALE
MERCURY_ORBIT_RADIUS = MERCURY_AU * AU_SCALE
VENUS_ORBIT_RADIUS = VENUS_AU * AU_SCALE
MARS_ORBIT_RADIUS = MARS_AU * AU_SCALE
JUPITER_ORBIT_RADIUS = JUPITER_AU * AU_SCALE
SATURN_ORBIT_RADIUS = SATURN_AU * AU_SCALE
URANUS_ORBIT_RADIUS = URANUS_AU * AU_SCALE
NEPTUNE_ORBIT_RADIUS = NEPTUNE_AU * AU_SCALE
PLUTO_ORBIT_RADIUS = PLUTO_AU * AU_SCALE
ERIS_ORBIT_RADIUS = ERIS_AU * AU_SCALE
HAUMEA_ORBIT_RADIUS = HAUMEA_AU * AU_SCALE
MAKEMAKE_ORBIT_RADIUS = MAKEMAKE_AU * AU_SCALE
CERES_ORBIT_RADIUS = CERES_AU * AU_SCALE


# Gravitational parameters (km^3/s^2)
GM_SUN = 1.327e11
GM_JUPITER = 1.266e8
GM_SATURN = 3.793e7
GM_NEPTUNE = 6.836e6
GM_URANUS = 5.794e6
GM_MARS = 4.282e4
GM_EARTH = 3.986e5
GM_VENUS = 3.249e5
GM_MERCURY = 2.2032e4
GM_PLUTO = 8.71e6
GM_ERIS = 8.31e6
GM_HAUMEA = 4.00e6  
GM_MAKEMAKE = 3.19e6  
GM_CERES = 6.26e4


GM_PLANET = {
    "earth": GM_EARTH,
    "mars": 4.282837e4,
    "jupiter": 1.2668653e8,
    "saturn": 3.7931187e7,
    "uranus": 5.793939e6,
    "neptune": 6.836529e6,
    "pluto": GM_PLUTO,
    "eris": GM_ERIS,
    "haumea": GM_HAUMEA,
    "makemake": GM_MAKEMAKE,
    "ceres": GM_CERES,
}



FPS = 60

# ---------------------------- ROTATION PERIODS (real-world) ----------------------------
# All values in **seconds** – taken from NASA fact sheets
ROTATION_PERIOD = {
    "sun":       25.0 * 86400,          # differential, ~25 d at equator
    "mercury":   58.646 * 86400,
    "venus":    -243.025 * 86400,       # retrograde → negative
    "earth":      0.99726968 * 86400,  # sidereal day
    "mars":       1.025957 * 86400,
    "jupiter":    9.925 * 3600,
    "saturn":    10.656 * 3600,
    "uranus":   -17.24 * 3600,         # retrograde
    "neptune":   16.11 * 3600,
    "pluto": -153.29 * 3600,  # Retrograde
    "eris": 25.9 * 3600,      # Approximate
    "haumea": 3.915 * 3600,   # Fast rotator
    "makemake": 22.48 * 3600, # Approximate
    "ceres": 9.074 * 3600,
    # Moons (tidally locked → same as orbital period)
    "moon":      27.321661 * 86400,
    "phobos":     0.31891 * 86400,
    "deimos":     1.26244 * 86400,
    "charon": 6.387 * 86400,  # Tidally locked to Pluto
    "styx": 3.24 * 86400,     # Approximate orbital
    "nix": 25.9 * 86400,      # Approximate
    "kerberos": 32.2 * 86400, # Approximate
    "hydra": 38.2 * 86400,
    # Add any other moon you want a spin for (most are locked)
}
# Convert to **degrees per simulated second**
ROTATION_DEG_PER_SIMSEC = {k: 360.0 / p if p != 0 else 0 for k, p in ROTATION_PERIOD.items()}

# ---------------------------- AXIAL TILTS (real-world, degrees from NASA/JPL) ----------------------------
AXIAL_TILTS = {
    "sun":      7.25,     # to ecliptic
    "mercury":  0.03,
    "venus":    177.4,    # retrograde, effective 2.6° inverted
    "earth":    23.44,
    "mars":     25.19,
    "jupiter":  3.13,
    "saturn":   26.73,
    "uranus":   97.77,    # nearly sideways, retrograde
    "neptune":  28.32,
    "pluto": 122.5,  # Extreme tilt
    "eris": 0.0,     # Unknown, assume equatorial
    "haumea": 0.0,   # Unknown, assume equatorial
    "makemake": 0.0, # Unknown, assume equatorial
    "ceres": 0.0,    # Minor tilt

    # Moons (relative to parent planet's equator; most tidally locked ~0°)
    "moon":     1.54,     # Earth's Moon
    "phobos":   0.0,
    "deimos":   0.0,
    # Add major Jupiter moons (all ~0°)
    "io":       0.0,
    "europa":   0.0,
    "ganymede": 0.0,
    "callisto": 0.0,
    # Add others as 0° (simplification for irregulars)
}


SECONDS_PER_FRAME = SIM_SECONDS_PER_REAL_SECOND / FPS


# Saturn ring constants (km)
REAL_SATURN_MEAN_RADIUS_KM = 58232.0
RING_C_INNER, RING_C_OUTER = 74658, 92000
RING_B_INNER, RING_B_OUTER = 92000, 117580
RING_CASSINI_INNER, RING_CASSINI_OUTER = 117580, 122170
RING_A_INNER, RING_A_OUTER = 122170, 136775

# Thresholds
POINT_LOD_THRESHOLD = JUPITER_RADIUS * 0.001
POINT_DISPLAY_SIZE = 6
JUPITER_TINY_THRESHOLD = JUPITER_RADIUS * 0.01
SATURN_TINY_THRESHOLD = SATURN_RADIUS * 0.01
URANUS_TINY_THRESHOLD = URANUS_RADIUS * 0.01
NEPTUNE_TINY_THRESHOLD = NEPTUNE_RADIUS * 0.01
EARTH_TINY_THRESHOLD = EARTH_RADIUS * 0.01
VENUS_TINY_THRESHOLD = VENUS_RADIUS * 0.01

TINY_MOON_SPIN_SCALE_MATCH = True

def hex_to_rgb(hex_str):
    """Convert hex color string to RGB tuple (floats 0-1)."""
    hex_str = hex_str.lstrip('#')
    return [int(hex_str[i:i+2], 16) / 255.0 for i in (0, 2, 4)]

SUBLIMATION_DISTANCE = 3 * AU_KM
# Comet data (real JPL Horizons osculating elements for Nov 11, 2025 epoch JD 2460990.5; a in km)
COMET_ELEMENTS = {
    "halley": {  # No tail
        "a_km": 17.84 * AU_KM,
        "e": 0.9671,
        "i_deg": 162.26,
        "Omega_deg": 59.40,
        "omega_deg": 112.05,
        "M0_deg": 345.5,
        "Tp_jd": 2461545.0,
        "epoch_jd": 2460990.5,
        "radius_km": 5.5,
        "color": "lightblue"
    },
    "halebopp": {  # No tail
        "a_km": 177.43 * AU_KM,
        "e": 0.99498,
        "i_deg": 89.29,
        "Omega_deg": 282.73,
        "omega_deg": 130.41,
        "M0_deg": 3.9,
        "Tp_jd": 2450539.64,
        "epoch_jd": 2460990.5,
        "radius_km": 25.0,
        "color": "cyan"
    },
    "enscke": {  # No tail
        "a_km": 2.215 * AU_KM,
        "e": 0.84833,
        "i_deg": 11.76,
        "Omega_deg": 334.6,
        "omega_deg": 186.5,
        "M0_deg": 180.2,
        "Tp_jd": 2460250.0,
        "epoch_jd": 2460990.5,
        "radius_km": 2.0,
        "color": "yellow"
    },
    "lovejoy": {  # No tail
        "a_km": 393.0 * AU_KM,
        "e": 0.998,
        "i_deg": 80.3,
        "Omega_deg": 12.4,
        "omega_deg": 90.3,
        "M0_deg": 0.27,
        "Tp_jd": 2457058.0,
        "epoch_jd": 2460990.5,
        "radius_km": 10.0,
        "color": "green"
    },
}



MOON_SCALE_FACTORS = {
    "earth": [
        {
            "name": "Moon",
            "a_km": 384400,
            "e": 0.0554,
            "i_deg": 5.16,
            "Omega_deg": 125.08,
            "omega_deg": 318.15,
            "M0_deg": 135.27,
            "epoch_jd": 2458849.5,
            "radius_km": 1737.4  # 0.2723 * 6371 km (Earth’s radius)
        }
    ],
    "mars": [
        {
            "name": "Phobos",
            "a_km": 9375,
            "e": 0.015,
            "i_deg": 1.1,
            "Omega_deg": 169.2,
            "omega_deg": 216.3,
            "M0_deg": 189.7,
            "epoch_jd": 2458849.5,
            "radius_km": 11.2  # 0.0033 * 3390 km (Mars’ radius)
        },
        {
            "name": "Deimos",
            "a_km": 23457,
            "e": 0.000,
            "i_deg": 1.8,
            "Omega_deg": 54.3,
            "omega_deg": 0.0,
            "M0_deg": 205.0,
            "epoch_jd": 2458849.5,
            "radius_km": 6.4  # 0.0019 * 3390 km
        }
    ],
    "jupiter": [
        {
            "name": "Metis",
            "a_km": 128000,
            "e": 0.0002,
            "i_deg": 0.060,
            "Omega_deg": 0,
            "omega_deg": 0,
            "M0_deg": 0,
            "epoch_jd": 2458849.5,
            "radius_km": 21.0  # 0.0003 * 69911 km
        },
        {
            "name": "Adrastea",
            "a_km": 129000,
            "e": 0.0015,
            "i_deg": 0.030,
            "Omega_deg": 0,
            "omega_deg": 0,
            "M0_deg": 0,
            "epoch_jd": 2458849.5,
            "radius_km": 7.0  # 0.0001 * 69911 km
        },
        {
            "name": "Amalthea",
            "a_km": 181400,
            "e": 0.0032,
            "i_deg": 0.374,
            "Omega_deg": 0,
            "omega_deg": 0,
            "M0_deg": 0,
            "epoch_jd": 2458849.5,
            "radius_km": 83.9  # 0.0012 * 69911 km
        },
        {
            "name": "Thebe",
            "a_km": 221900,
            "e": 0.0175,
            "i_deg": 1.076,
            "Omega_deg": 0,
            "omega_deg": 0,
            "M0_deg": 0,
            "epoch_jd": 2458849.5,
            "radius_km": 48.9  # 0.0007 * 69911 km
        },
        {
            "name": "Io",
            "a_km": 421800,
            "e": 0.0041,
            "i_deg": 0.050,
            "Omega_deg": 43.977,  # Using example value
            "omega_deg": 84.129,
            "M0_deg": 171.016,
            "epoch_jd": 2458849.5,
            "radius_km": 1821.6  # 0.0255 * 69911 km (matches known radius)
        },
        {
            "name": "Europa",
            "a_km": 671100,
            "e": 0.0090,
            "i_deg": 0.470,
            "Omega_deg": 219.106,  # Using example value
            "omega_deg": 88.970,
            "M0_deg": 317.021,
            "epoch_jd": 2458849.5,
            "radius_km": 1560.8  # 0.0218 * 69911 km (matches known radius)
        },
        {
            "name": "Ganymede",
            "a_km": 1070400,
            "e": 0.0013,
            "i_deg": 0.200,
            "Omega_deg": 0,
            "omega_deg": 0,
            "M0_deg": 0,
            "epoch_jd": 2458849.5,
            "radius_km": 2631.2  # 0.0368 * 69911 km
        },
        {
            "name": "Callisto",
            "a_km": 1882700,
            "e": 0.0074,
            "i_deg": 0.192,
            "Omega_deg": 0,
            "omega_deg": 0,
            "M0_deg": 0,
            "epoch_jd": 2458849.5,
            "radius_km": 2410.3  # 0.0337 * 69911 km
        },
        {
            "name": "Themisto",
            "a_km": 7397000,
            "e": 0.257,
            "i_deg": 44.3,
            "Omega_deg": 0,
            "omega_deg": 0,
            "M0_deg": 0,
            "epoch_jd": 2458849.5,
            "radius_km": 7.0  # 0.0001 * 69911 km
        },
        {
            "name": "Leda",
            "a_km": 11145200,
            "e": 0.162,
            "i_deg": 28.2,
            "Omega_deg": 0,
            "omega_deg": 0,
            "M0_deg": 0,
            "epoch_jd": 2458849.5,
            "radius_km": 14.0  # 0.0002 * 69911 km
        },
        {
            "name": "Ersa",
            "a_km": 11399400,
            "e": 0.117,
            "i_deg": 29.0,
            "Omega_deg": 0,
            "omega_deg": 0,
            "M0_deg": 0,
            "epoch_jd": 2458849.5,
            "radius_km": 1.4  # 0.00002 * 69911 km
        },
        {
            "name": "S/2018 J 2",
            "a_km": 11419700,
            "e": 0.152,
            "i_deg": 28.3,
            "Omega_deg": 0,
            "omega_deg": 0,
            "M0_deg": 0,
            "epoch_jd": 2458849.5,
            "radius_km": 1.4  # 0.00002 * 69911 km
        },
        {
            "name": "Himalia",
            "a_km": 11439000,
            "e": 0.160,
            "i_deg": 28.4,
            "Omega_deg": 0,
            "omega_deg": 0,
            "M0_deg": 0,
            "epoch_jd": 2458849.5,
            "radius_km": 69.9  # 0.0010 * 69911 km
        },
        {
            "name": "Pandia",
            "a_km": 11479600,
            "e": 0.178,
            "i_deg": 28.9,
            "Omega_deg": 0,
            "omega_deg": 0,
            "M0_deg": 0,
            "epoch_jd": 2458849.5,
            "radius_km": 1.4  # 0.00002 * 69911 km
        },
        {
            "name": "Lysithea",
            "a_km": 11699100,
            "e": 0.117,
            "i_deg": 27.7,
            "Omega_deg": 0,
            "omega_deg": 0,
            "M0_deg": 0,
            "epoch_jd": 2458849.5,
            "radius_km": 21.0  # 0.0003 * 69911 km
        },
        {
            "name": "Elara",
            "a_km": 11710700,
            "e": 0.212,
            "i_deg": 27.8,
            "Omega_deg": 0,
            "omega_deg": 0,
            "M0_deg": 0,
            "epoch_jd": 2458849.5,
            "radius_km": 41.9  # 0.0006 * 69911 km
        },
        {
            "name": "S/2011 J 3",
            "a_km": 11716800,
            "e": 0.192,
            "i_deg": 27.6,
            "Omega_deg": 0,
            "omega_deg": 0,
            "M0_deg": 0,
            "epoch_jd": 2458849.5,
            "radius_km": 1.4  # 0.00002 * 69911 km
        },
        {
            "name": "Dia",
            "a_km": 12257900,
            "e": 0.232,
            "i_deg": 29.1,
            "Omega_deg": 0,
            "omega_deg": 0,
            "M0_deg": 0,
            "epoch_jd": 2458849.5,
            "radius_km": 2.1  # 0.00003 * 69911 km
        },
        {
            "name": "S/2018 J 4",
            "a_km": 16328500,
            "e": 0.177,
            "i_deg": 50.2,
            "Omega_deg": 0,
            "omega_deg": 0,
            "M0_deg": 0,
            "epoch_jd": 2458849.5,
            "radius_km": 0.7  # 0.00001 * 69911 km
        },
        {
            "name": "Carpo",
            "a_km": 17039500,
            "e": 0.415,
            "i_deg": 53.3,
            "Omega_deg": 0,
            "omega_deg": 0,
            "M0_deg": 0,
            "epoch_jd": 2458849.5,
            "radius_km": 1.4  # 0.00002 * 69911 km
        },
        {
            "name": "Valetudo",
            "a_km": 18690100,
            "e": 0.217,
            "i_deg": 34.5,
            "Omega_deg": 0,
            "omega_deg": 0,
            "M0_deg": 0,
            "epoch_jd": 2458849.5,
            "radius_km": 0.5  # 0.000007 * 69911 km
        },
        {
            "name": "Euporie",
            "a_km": 19261900,
            "e": 0.148,
            "i_deg": 145.5,
            "Omega_deg": 0,
            "omega_deg": 0,
            "M0_deg": 0,
            "epoch_jd": 2458849.5,
            "radius_km": 0.7  # 0.00001 * 69911 km
        },
        {
            "name": "S/2003 J 18",
            "a_km": 20332800,
            "e": 0.102,
            "i_deg": 145.7,
            "Omega_deg": 0,
            "omega_deg": 0,
            "M0_deg": 0,
            "epoch_jd": 2458849.5,
            "radius_km": 0.7  # 0.00001 * 69911 km
        },
        {
            "name": "Eupheme",
            "a_km": 20763400,
            "e": 0.234,
            "i_deg": 147.9,
            "Omega_deg": 0,
            "omega_deg": 0,
            "M0_deg": 0,
            "epoch_jd": 2458849.5,
            "radius_km": 0.7  # 0.00001 * 69911 km
        },
        {
            "name": "S/2021 J 3",
            "a_km": 20776600,
            "e": 0.239,
            "i_deg": 147.9,
            "Omega_deg": 0,
            "omega_deg": 0,
            "M0_deg": 0,
            "epoch_jd": 2458849.5,
            "radius_km": 0.7  # 0.00001 * 69911 km
        },
        {
            "name": "S/2010 J 2",
            "a_km": 20786900,
            "e": 0.244,
            "i_deg": 148.0,
            "Omega_deg": 0,
            "omega_deg": 0,
            "M0_deg": 0,
            "epoch_jd": 2458849.5,
            "radius_km": 0.5  # 0.000007 * 69911 km
        },
        {
            "name": "S/2016 J 1",
            "a_km": 20796700,
            "e": 0.245,
            "i_deg": 145.1,
            "Omega_deg": 0,
            "omega_deg": 0,
            "M0_deg": 0,
            "epoch_jd": 2458849.5,
            "radius_km": 0.5  # 0.000007 * 69911 km
        },
        {
            "name": "Mneme",
            "a_km": 20815800,
            "e": 0.240,
            "i_deg": 147.8,
            "Omega_deg": 0,
            "omega_deg": 0,
            "M0_deg": 0,
            "epoch_jd": 2458849.5,
            "radius_km": 0.7  # 0.00001 * 69911 km
        },
        {
            "name": "Euanthe",
            "a_km": 20822900,
            "e": 0.243,
            "i_deg": 148.1,
            "Omega_deg": 0,
            "omega_deg": 0,
            "M0_deg": 0,
            "epoch_jd": 2458849.5,
            "radius_km": 1.4  # 0.00002 * 69911 km
        },
        {
            "name": "S/2003 J 16",
            "a_km": 20877500,
            "e": 0.238,
            "i_deg": 147.8,
            "Omega_deg": 0,
            "omega_deg": 0,
            "M0_deg": 0,
            "epoch_jd": 2458849.5,
            "radius_km": 0.7  # 0.00001 * 69911 km
        }
    ],
    "saturn": [
        {
            "name": "S/2009 S 1",
            "a_km": 116900,
            "e": 0.000,
            "i_deg": 0.0,
            "Omega_deg": 0,
            "omega_deg": 0,
            "M0_deg": 0,
            "epoch_jd": 2458849.5,
            "radius_km": 151.4  # 0.0026 * 58232 km
        },
        {
            "name": "Pan",
            "a_km": 133600,
            "e": 0.000,
            "i_deg": 0.0,
            "Omega_deg": 0,
            "omega_deg": 0,
            "M0_deg": 0,
            "epoch_jd": 2458849.5,
            "radius_km": 11.6  # 0.0002 * 58232 km
        },
        {
            "name": "Daphnis",
            "a_km": 136500,
            "e": 0.000,
            "i_deg": 0.0,
            "Omega_deg": 0,
            "omega_deg": 0,
            "M0_deg": 0,
            "epoch_jd": 2458849.5,
            "radius_km": 4.1  # 0.00007 * 58232 km
        },
        {
            "name": "Atlas",
            "a_km": 137700,
            "e": 0.001,
            "i_deg": 0.0,
            "Omega_deg": 0,
            "omega_deg": 0,
            "M0_deg": 0,
            "epoch_jd": 2458849.5,
            "radius_km": 17.5  # 0.0003 * 58232 km
        },
        {
            "name": "Prometheus",
            "a_km": 139400,
            "e": 0.002,
            "i_deg": 0.0,
            "Omega_deg": 0,
            "omega_deg": 0,
            "M0_deg": 0,
            "epoch_jd": 2458849.5,
            "radius_km": 40.8  # 0.0007 * 58232 km
        },
        {
            "name": "Pandora",
            "a_km": 141700,
            "e": 0.004,
            "i_deg": 0.0,
            "Omega_deg": 0,
            "omega_deg": 0,
            "M0_deg": 0,
            "epoch_jd": 2458849.5,
            "radius_km": 40.8  # 0.0007 * 58232 km
        },
        {
            "name": "Epimetheus",
            "a_km": 151400,
            "e": 0.020,
            "i_deg": 0.3,
            "Omega_deg": 0,
            "omega_deg": 0,
            "M0_deg": 0,
            "epoch_jd": 2458849.5,
            "radius_km": 58.2  # 0.0010 * 58232 km
        },
        {
            "name": "Janus",
            "a_km": 151500,
            "e": 0.007,
            "i_deg": 0.2,
            "Omega_deg": 0,
            "omega_deg": 0,
            "M0_deg": 0,
            "epoch_jd": 2458849.5,
            "radius_km": 87.3  # 0.0015 * 58232 km
        },
        {
            "name": "Aegaeon",
            "a_km": 167500,
            "e": 0.000,
            "i_deg": 0.0,
            "Omega_deg": 0,
            "omega_deg": 0,
            "M0_deg": 0,
            "epoch_jd": 2458849.5,
            "radius_km": 0.3  # 0.000006 * 58232 km
        },
        {
            "name": "Mimas",
            "a_km": 186000,
            "e": 0.020,
            "i_deg": 1.6,
            "Omega_deg": 0,
            "omega_deg": 0,
            "M0_deg": 0,
            "epoch_jd": 2458849.5,
            "radius_km": 198.0  # 0.0034 * 58232 km
        },
        {
            "name": "Methone",
            "a_km": 194700,
            "e": 0.002,
            "i_deg": 0.0,
            "Omega_deg": 0,
            "omega_deg": 0,
            "M0_deg": 0,
            "epoch_jd": 2458849.5,
            "radius_km": 1.2  # 0.00002 * 58232 km
        },
        {
            "name": "Anthe",
            "a_km": 198100,
            "e": 0.002,
            "i_deg": 0.0,
            "Omega_deg": 0,
            "omega_deg": 0,
            "M0_deg": 0,
            "epoch_jd": 2458849.5,
            "radius_km": 0.9  # 0.000015 * 58232 km
        },
        {
            "name": "Pallene",
            "a_km": 212300,
            "e": 0.004,
            "i_deg": 0.2,
            "Omega_deg": 0,
            "omega_deg": 0,
            "M0_deg": 0,
            "epoch_jd": 2458849.5,
            "radius_km": 2.3  # 0.00004 * 58232 km
        },
        {
            "name": "Enceladus",
            "a_km": 238400,
            "e": 0.005,
            "i_deg": 0.0,
            "Omega_deg": 0,
            "omega_deg": 0,
            "M0_deg": 0,
            "epoch_jd": 2458849.5,
            "radius_km": 252.1  # 0.0043 * 58232 km
        },
        {
            "name": "Tethys",
            "a_km": 295000,
            "e": 0.001,
            "i_deg": 1.1,
            "Omega_deg": 0,
            "omega_deg": 0,
            "M0_deg": 0,
            "epoch_jd": 2458849.5,
            "radius_km": 531.0  # 0.0091 * 58232 km
        },
        {
            "name": "Telesto",
            "a_km": 295000,
            "e": 0.001,
            "i_deg": 1.2,
            "Omega_deg": 0,
            "omega_deg": 0,
            "M0_deg": 0,
            "epoch_jd": 2458849.5,
            "radius_km": 11.6  # 0.0002 * 58232 km
        },
        {
            "name": "Calypso",
            "a_km": 295000,
            "e": 0.001,
            "i_deg": 1.5,
            "Omega_deg": 0,
            "omega_deg": 0,
            "M0_deg": 0,
            "epoch_jd": 2458849.5,
            "radius_km": 11.6  # 0.0002 * 58232 km
        },
        {
            "name": "Helene",
            "a_km": 377600,
            "e": 0.007,
            "i_deg": 0.2,
            "Omega_deg": 0,
            "omega_deg": 0,
            "M0_deg": 0,
            "epoch_jd": 2458849.5,
            "radius_km": 17.5  # 0.0003 * 58232 km
        },
        {
            "name": "Polydeuces",
            "a_km": 377600,
            "e": 0.019,
            "i_deg": 0.2,
            "Omega_deg": 0,
            "omega_deg": 0,
            "M0_deg": 0,
            "epoch_jd": 2458849.5,
            "radius_km": 1.2  # 0.00002 * 58232 km
        },
        {
            "name": "Dione",
            "a_km": 377700,
            "e": 0.002,
            "i_deg": 0.0,
            "Omega_deg": 0,
            "omega_deg": 0,
            "M0_deg": 0,
            "epoch_jd": 2458849.5,
            "radius_km": 561.7  # 0.0096 * 58232 km
        },
        {
            "name": "Rhea",
            "a_km": 527200,
            "e": 0.001,
            "i_deg": 0.3,
            "Omega_deg": 0,
            "omega_deg": 0,
            "M0_deg": 0,
            "epoch_jd": 2458849.5,
            "radius_km": 764.5  # 0.0131 * 58232 km
        },
        {
            "name": "Titan",
            "a_km": 1221900,
            "e": 0.029,
            "i_deg": 0.3,
            "Omega_deg": 28.06,  # Using example value
            "omega_deg": 186.59,
            "M0_deg": 230.0,
            "epoch_jd": 2458849.5,
            "radius_km": 2575.0  # 0.0442 * 58232 km (matches known radius)
        },
        {
            "name": "Hyperion",
            "a_km": 1481500,
            "e": 0.105,
            "i_deg": 0.6,
            "Omega_deg": 0,
            "omega_deg": 0,
            "M0_deg": 0,
            "epoch_jd": 2458849.5,
            "radius_km": 134.0  # 0.0023 * 58232 km
        },
        {
            "name": "Iapetus",
            "a_km": 3561700,
            "e": 0.028,
            "i_deg": 7.6,
            "Omega_deg": 0,
            "omega_deg": 0,
            "M0_deg": 0,
            "epoch_jd": 2458849.5,
            "radius_km": 734.5  # 0.0126 * 58232 km
        },
        {
            "name": "S/2023 S 1",
            "a_km": 11205400,
            "e": 0.386,
            "i_deg": 48.8,
            "Omega_deg": 0,
            "omega_deg": 0,
            "M0_deg": 0,
            "epoch_jd": 2458849.5,
            "radius_km": 1503.0  # 0.0258 * 58232 km
        },
        {
            "name": "S/2019 S 1",
            "a_km": 11245400,
            "e": 0.384,
            "i_deg": 49.5,
            "Omega_deg": 0,
            "omega_deg": 0,
            "M0_deg": 0,
            "epoch_jd": 2458849.5,
            "radius_km": 2498.8  # 0.0429 * 58232 km
        }
    ],
    "uranus": [
        {
            "name": "Cordelia",
            "a_km": 49800,
            "e": 0.000,
            "i_deg": 0.2,
            "Omega_deg": 0,
            "omega_deg": 0,
            "M0_deg": 0,
            "epoch_jd": 2458849.5,
            "radius_km": 20.3  # 0.0008 * 25362 km
        },
        {
            "name": "Ophelia",
            "a_km": 53800,
            "e": 0.011,
            "i_deg": 0.1,
            "Omega_deg": 0,
            "omega_deg": 0,
            "M0_deg": 0,
            "epoch_jd": 2458849.5,
            "radius_km": 20.3  # 0.0008 * 25362 km
        },
        {
            "name": "S/2025 U 1",
            "a_km": 57800,
            "e": 0.039,
            "i_deg": 4.0,
            "Omega_deg": 0,
            "omega_deg": 0,
            "M0_deg": 0,
            "epoch_jd": 2458849.5,
            "radius_km": 10.1  # 0.0004 * 25362 km
        },
        {
            "name": "Bianca",
            "a_km": 59200,
            "e": 0.001,
            "i_deg": 0.1,
            "Omega_deg": 0,
            "omega_deg": 0,
            "M0_deg": 0,
            "epoch_jd": 2458849.5,
            "radius_km": 25.4  # 0.0010 * 25362 km
        },
        {
            "name": "Cressida",
            "a_km": 61800,
            "e": 0.000,
            "i_deg": 0.1,
            "Omega_deg": 0,
            "omega_deg": 0,
            "M0_deg": 0,
            "epoch_jd": 2458849.5,
            "radius_km": 40.6  # 0.0016 * 25362 km
        },
        {
            "name": "Desdemona",
            "a_km": 62700,
            "e": 0.000,
            "i_deg": 0.1,
            "Omega_deg": 0,
            "omega_deg": 0,
            "M0_deg": 0,
            "epoch_jd": 2458849.5,
            "radius_km": 35.5  # Approximate, based on known size ~35 km
        },
        {
            "name": "Juliet",
            "a_km": 64400,
            "e": 0.001,
            "i_deg": 0.1,
            "Omega_deg": 0,
            "omega_deg": 0,
            "M0_deg": 0,
            "epoch_jd": 2458849.5,
            "radius_km": 46.7  # Approximate, based on known size ~47 km
        },
        {
            "name": "Portia",
            "a_km": 66100,
            "e": 0.000,
            "i_deg": 0.1,
            "Omega_deg": 0,
            "omega_deg": 0,
            "M0_deg": 0,
            "epoch_jd": 2458849.5,
            "radius_km": 67.6  # Approximate, based on known size ~68 km
        },
        {
            "name": "Rosalind",
            "a_km": 69900,
            "e": 0.000,
            "i_deg": 0.1,
            "Omega_deg": 0,
            "omega_deg": 0,
            "M0_deg": 0,
            "epoch_jd": 2458849.5,
            "radius_km": 36.0  # Approximate, based on known size ~36 km
        },
        {
            "name": "Belinda",
            "a_km": 75300,
            "e": 0.000,
            "i_deg": 0.1,
            "Omega_deg": 0,
            "omega_deg": 0,
            "M0_deg": 0,
            "epoch_jd": 2458849.5,
            "radius_km": 45.3  # Approximate, based on known size ~45 km
        },
        {
            "name": "Puck",
            "a_km": 86000,
            "e": 0.000,
            "i_deg": 0.3,
            "Omega_deg": 0,
            "omega_deg": 0,
            "M0_deg": 0,
            "epoch_jd": 2458849.5,
            "radius_km": 81.0  # Approximate, based on known size ~81 km
        },
        {
            "name": "Miranda",
            "a_km": 129900,
            "e": 0.0013,
            "i_deg": 4.3,
            "Omega_deg": 0,
            "omega_deg": 0,
            "M0_deg": 0,
            "epoch_jd": 2458849.5,
            "radius_km": 235.8  # Approximate, based on known size ~236 km
        },
        {
            "name": "Ariel",
            "a_km": 190900,
            "e": 0.0012,
            "i_deg": 0.0,
            "Omega_deg": 0,
            "omega_deg": 0,
            "M0_deg": 0,
            "epoch_jd": 2458849.5,
            "radius_km": 578.9  # Approximate, based on known size ~579 km
        },
        {
            "name": "Umbriel",
            "a_km": 266000,
            "e": 0.0039,
            "i_deg": 0.0,
            "Omega_deg": 0,
            "omega_deg": 0,
            "M0_deg": 0,
            "epoch_jd": 2458849.5,
            "radius_km": 584.7  # Approximate, based on known size ~585 km
        },
        {
            "name": "Titania",
            "a_km": 436300,
            "e": 0.0011,
            "i_deg": 0.0,
            "Omega_deg": 0,
            "omega_deg": 0,
            "M0_deg": 0,
            "epoch_jd": 2458849.5,
            "radius_km": 788.9  # Approximate, based on known size ~789 km
        },
        {
            "name": "Oberon",
            "a_km": 583500,
            "e": 0.0014,
            "i_deg": 0.0,
            "Omega_deg": 0,
            "omega_deg": 0,
            "M0_deg": 0,
            "epoch_jd": 2458849.5,
            "radius_km": 761.4  # Approximate, based on known size ~761 km
        },
        {
            "name": "Francisco",
            "a_km": 4282900,
            "e": 0.145,
            "i_deg": 147.3,
            "Omega_deg": 0,
            "omega_deg": 0,
            "M0_deg": 0,
            "epoch_jd": 2458849.5,
            "radius_km": 11.0  # Approximate, based on known size ~11 km
        },
        {
            "name": "Caliban",
            "a_km": 7231000,
            "e": 0.181,
            "i_deg": 141.7,
            "Omega_deg": 0,
            "omega_deg": 0,
            "M0_deg": 0,
            "epoch_jd": 2458849.5,
            "radius_km": 36.0  # Approximate, based on known size ~36 km
        },
        {
            "name": "Stephano",
            "a_km": 8004000,
            "e": 0.229,
            "i_deg": 143.8,
            "Omega_deg": 0,
            "omega_deg": 0,
            "M0_deg": 0,
            "epoch_jd": 2458849.5,
            "radius_km": 16.0  # Approximate, based on known size ~16 km
        },
        {
            "name": "Trinculo",
            "a_km": 8504000,
            "e": 0.219,
            "i_deg": 166.3,
            "Omega_deg": 0,
            "omega_deg": 0,
            "M0_deg": 0,
            "epoch_jd": 2458849.5,
            "radius_km": 9.0  # Approximate, based on known size ~9 km
        },
        {
            "name": "Sycorax",
            "a_km": 12179000,
            "e": 0.522,
            "i_deg": 159.4,
            "Omega_deg": 0,
            "omega_deg": 0,
            "M0_deg": 0,
            "epoch_jd": 2458849.5,
            "radius_km": 75.0  # Approximate, based on known size ~75 km
        },
        {
            "name": "Margaret",
            "a_km": 14345000,
            "e": 0.677,
            "i_deg": 57.4,
            "Omega_deg": 0,
            "omega_deg": 0,
            "M0_deg": 0,
            "epoch_jd": 2458849.5,
            "radius_km": 10.0  # Approximate, based on known size ~10 km
        },
        {
            "name": "Prospero",
            "a_km": 16268000,
            "e": 0.445,
            "i_deg": 151.8,
            "Omega_deg": 0,
            "omega_deg": 0,
            "M0_deg": 0,
            "epoch_jd": 2458849.5,
            "radius_km": 25.0  # Approximate, based on known size ~25 km
        },
        {
            "name": "Setebos",
            "a_km": 17420000,
            "e": 0.591,
            "i_deg": 158.2,
            "Omega_deg": 0,
            "omega_deg": 0,
            "M0_deg": 0,
            "epoch_jd": 2458849.5,
            "radius_km": 24.0  # Approximate, based on known size ~24 km
        },
        {
            "name": "Ferdinand",
            "a_km": 20430000,
            "e": 0.368,
            "i_deg": 169.8,
            "Omega_deg": 0,
            "omega_deg": 0,
            "M0_deg": 0,
            "epoch_jd": 2458849.5,
            "radius_km": 10.0  # Approximate, based on known size ~10 km
        }
    ],
    "neptune": [
        {
            "name": "Naiad",
            "a_km": 48224,
            "e": 0.0047,
            "i_deg": 4.691,
            "Omega_deg": 0,
            "omega_deg": 0,
            "M0_deg": 0,
            "epoch_jd": 2458849.5,
            "radius_km": 29.5  # 0.0012 * 24622 km
        },
        {
            "name": "Thalassa",
            "a_km": 50074,
            "e": 0.0018,
            "i_deg": 0.2,
            "Omega_deg": 0,
            "omega_deg": 0,
            "M0_deg": 0,
            "epoch_jd": 2458849.5,
            "radius_km": 41.0  # Approximate, based on known size ~41 km
        },
        {
            "name": "Despina",
            "a_km": 52526,
            "e": 0.0002,
            "i_deg": 0.2,
            "Omega_deg": 0,
            "omega_deg": 0,
            "M0_deg": 0,
            "epoch_jd": 2458849.5,
            "radius_km": 75.0  # Approximate, based on known size ~75 km
        },
        {
            "name": "Galatea",
            "a_km": 61953,
            "e": 0.0001,
            "i_deg": 0.1,
            "Omega_deg": 0,
            "omega_deg": 0,
            "M0_deg": 0,
            "epoch_jd": 2458849.5,
            "radius_km": 87.0  # Approximate, based on known size ~87 km
        },
        {
            "name": "Larissa",
            "a_km": 73548,
            "e": 0.0014,
            "i_deg": 0.2,
            "Omega_deg": 0,
            "omega_deg": 0,
            "M0_deg": 0,
            "epoch_jd": 2458849.5,
            "radius_km": 97.0  # Approximate, based on known size ~97 km
        },
        {
            "name": "Hippocamp",
            "a_km": 105283,
            "e": 0.0009,
            "i_deg": 0.0,
            "Omega_deg": 0,
            "omega_deg": 0,
            "M0_deg": 0,
            "epoch_jd": 2458849.5,
            "radius_km": 17.4  # Approximate, based on known size ~17 km
        },
        {
            "name": "Proteus",
            "a_km": 117646,
            "e": 0.0005,
            "i_deg": 0.1,
            "Omega_deg": 0,
            "omega_deg": 0,
            "M0_deg": 0,
            "epoch_jd": 2458849.5,
            "radius_km": 210.0  # Approximate, based on known size ~210 km
        },
        {
            "name": "Triton",
            "a_km": 354759,
            "e": 0.0000,
            "i_deg": 156.8,
            "Omega_deg": 0,
            "omega_deg": 0,
            "M0_deg": 0,
            "epoch_jd": 2458849.5,
            "radius_km": 1353.4  # Approximate, based on known size ~1353 km
        },
        {
            "name": "Nereid",
            "a_km": 5513400,
            "e": 0.7507,
            "i_deg": 7.1,
            "Omega_deg": 0,
            "omega_deg": 0,
            "M0_deg": 0,
            "epoch_jd": 2458849.5,
            "radius_km": 178.5  # Approximate, based on known size ~179 km
        },
        {
            "name": "Halimede",
            "a_km": 16611000,
            "e": 0.264,
            "i_deg": 134.1,
            "Omega_deg": 0,
            "omega_deg": 0,
            "M0_deg": 0,
            "epoch_jd": 2458849.5,
            "radius_km": 31.0  # Approximate, based on known size ~31 km
        },
        {
            "name": "Sao",
            "a_km": 22228000,
            "e": 0.293,
            "i_deg": 49.9,
            "Omega_deg": 0,
            "omega_deg": 0,
            "M0_deg": 0,
            "epoch_jd": 2458849.5,
            "radius_km": 22.0  # Approximate, based on known size ~22 km
        },
        {
            "name": "Laomedeia",
            "a_km": 23567000,
            "e": 0.424,
            "i_deg": 34.7,
            "Omega_deg": 0,
            "omega_deg": 0,
            "M0_deg": 0,
            "epoch_jd": 2458849.5,
            "radius_km": 21.0  # Approximate, based on known size ~21 km
        },
        {
            "name": "Psamathe",
            "a_km": 46695000,
            "e": 0.461,
            "i_deg": 137.7,
            "Omega_deg": 0,
            "omega_deg": 0,
            "M0_deg": 0,
            "epoch_jd": 2458849.5,
            "radius_km": 20.0  # Approximate, based on known size ~20 km
        },
        {
            "name": "Neso",
            "a_km": 49245000,
            "e": 0.424,
            "i_deg": 136.5,
            "Omega_deg": 0,
            "omega_deg": 0,
            "M0_deg": 0,
            "epoch_jd": 2458849.5,
            "radius_km": 30.0  # Approximate, based on known size ~30 km
        }
    ],
    "pluto": [
        {
            "name": "Charon",
            "a_km": 19571,
            "e": 0.0002,
            "i_deg": 0.0,
            "Omega_deg": 0.0,
            "omega_deg": 0.0,
            "M0_deg": 0.0,
            "epoch_jd": 2458849.5,
            "radius_km": 606.0  # Spherical
        },
        {
            "name": "Styx",
            "a_km": 42656,
            "e": 0.000,
            "i_deg": 0.9,
            "Omega_deg": 0.0,
            "omega_deg": 0.0,
            "M0_deg": 0.0,
            "epoch_jd": 2458849.5,
            "radius_km": 6.0  # Irregular, small
        },
        {
            "name": "Nix",
            "a_km": 48694,
            "e": 0.000,
            "i_deg": 0.1,
            "Omega_deg": 0.0,
            "omega_deg": 0.0,
            "M0_deg": 0.0,
            "epoch_jd": 2458849.5,
            "radius_km": 19.0  # Irregular
        },
        {
            "name": "Kerberos",
            "a_km": 57783,
            "e": 0.000,
            "i_deg": 0.4,
            "Omega_deg": 0.0,
            "omega_deg": 0.0,
            "M0_deg": 0.0,
            "epoch_jd": 2458849.5,
            "radius_km": 9.0  # Irregular, small
        },
        {
            "name": "Hydra",
            "a_km": 64738,
            "e": 0.000,
            "i_deg": 0.3,
            "Omega_deg": 0.0,
            "omega_deg": 0.0,
            "M0_deg": 0.0,
            "epoch_jd": 2458849.5,
            "radius_km": 24.0  # Irregular
        }
    ],
    "eris": [
        {
            "name": "Dysnomia",
            "a_km": 37350,
            "e": 0.000,
            "i_deg": 0.0,
            "Omega_deg": 0,
            "omega_deg": 0,
            "M0_deg": 0,
            "epoch_jd": 2458849.5,
            "radius_km": 35.0  # Approximate
        }
    ],
    "haumea": [
        {
            "name": "Hi'iaka",
            "a_km": 49100,
            "e": 0.07,
            "i_deg": 0.0,
            "Omega_deg": 0,
            "omega_deg": 0,
            "M0_deg": 0,
            "epoch_jd": 2458849.5,
            "radius_km": 80.0  # Approximate
        }
    ],
    "makemake": [
        {
            "name": "MK2",
            "a_km": 21000,  # Approximate
            "e": 0.0,
            "i_deg": 0.0,
            "Omega_deg": 0,
            "omega_deg": 0,
            "M0_deg": 0,
            "epoch_jd": 2458849.5,
            "radius_km": 175.0  # Approximate
        }
    ],
}


# Object type classifier for zoom limits
OBJECT_MIN_ZOOM_MULTIPLIER = {
    # Planets: surface skim
    "mercury": 1.1, "venus": 1.1, "earth": 1.1, "mars": 1.1, "haumea": 1.1, "makemake": 1.1, "ceres": 1.1,
    "jupiter": 1.1, "saturn": 1.1, "uranus": 1.1, "neptune": 1.1, "pluto": 1.1, "eris": 1.1,
    "sun": 2.0,  
    "default": 5.0
}

planet_names = ["mercury", "venus", "earth", "mars", "jupiter", "saturn", "uranus", "neptune", "pluto", "eris", "haumea", "makemake", "ceres"] 

def get_object_min_dist(actor_name, radius_scene=0):
    """
    Compute min zoom dist based on actor name and size.
    For moons: Use radius_km from data if available.
    """
    if not isinstance(actor_name, str):
        actor_name = str(getattr(actor_name, 'name', 'default'))
    
    # Extract base name (e.g., "jupiter_io" -> "jupiter", moon="io")
    parts = actor_name.lower().split('_')
    base = parts[0] if parts else "default"
    is_moon = len(parts) > 1 and base in planet_names
    
    if not is_moon:
        # Planet/Sun: multiplier * radius
        radius = {
            "mercury": MERCURY_RADIUS, "venus": VENUS_RADIUS, "earth": EARTH_RADIUS,
            "mars": MARS_RADIUS, "jupiter": JUPITER_RADIUS, "saturn": SATURN_RADIUS,
            "uranus": URANUS_RADIUS, "neptune": NEPTUNE_RADIUS, "sun": SUN_RADIUS
        }.get(base, 1.0)
        return OBJECT_MIN_ZOOM_MULTIPLIER.get(base, 5.0) * radius
    
    # Moon: dynamic based on size category
    moon_name = '_'.join(parts[1:])  # Full moon name
    # Find radius_km from MOON_SCALE_FACTORS (search all planets)
    radius_km = 1.0  # Default tiny
    for p_moons in MOON_SCALE_FACTORS.values():
        for m_data in p_moons:
            if m_data["name"].lower() == moon_name.lower():
                radius_km = m_data["radius_km"]
                break
        if radius_km > 1.0: break
    
    radius_scene = radius_km * SIZEKM_TO_SCENE if radius_scene == 0 else radius_scene
    if radius_km > 100:  # Large (e.g., Ganymede=2631 km)
        multiplier = 1.5
    elif radius_km > 10:  # Small (e.g., Phobos=11 km)
        multiplier = 3.0
    else:  # Tiny (e.g., S/2009 S 1 ~0.3 km effective)
        multiplier = max(10.0, 0.5 / radius_scene)  # Floor to avoid sub-pixel
    
    return multiplier * radius_scene


moon_colors = {
    "Moon": "#C0C0C0",  # Silver gray for Earth's Moon, dusty regolith
    "Phobos": "#A0522D",  # Sienna, reddish-brown, rocky with dust
    "Deimos": "#808080",  # Gray, rocky and dusty, slightly reddish
    "Metis": "#696969",  # Dim gray, small rocky inner moon of Jupiter
    "Adrastea": "#696969",  # Dim gray, small rocky inner moon
    "Amalthea": "#8B4513",  # Saddle brown, reddish rocky surface
    "Thebe": "#696969",  # Dim gray, small rocky moon
    "Io": "#FFD700",  # Gold, yellowish for volcanic sulfur-rich surface
    "Europa": "#F0F8FF",  # Alice blue, icy white with faint reddish streaks
    "Ganymede": "#A9A9A9",  # Dark gray, mix of icy and rocky terrain
    "Callisto": "#696969",  # Dim gray, heavily cratered dark surface
    "Themisto": "#CCCCCC",  # Light gray, small moon, limited data
    "Leda": "#CCCCCC",  # Light gray, small irregular moon
    "Ersa": "#CCCCCC",  # Light gray, small irregular moon
    "S/2018 J 2": "#CCCCCC",  # Light gray, small provisional moon
    "Himalia": "#A9A9A9",  # Dark gray, larger irregular rocky moon
    "Pandia": "#CCCCCC",  # Light gray, small irregular moon
    "Lysithea": "#CCCCCC",  # Light gray, small irregular moon
    "Elara": "#CCCCCC",  # Light gray, small irregular moon
    "S/2011 J 3": "#CCCCCC",  # Light gray, small provisional moon
    "Dia": "#CCCCCC",  # Light gray, small irregular moon
    "S/2018 J 4": "#CCCCCC",  # Light gray, small provisional moon
    "Carpo": "#CCCCCC",  # Light gray, small irregular moon
    "Valetudo": "#CCCCCC",  # Light gray, small irregular moon
    "Euporie": "#CCCCCC",  # Light gray, small irregular moon
    "S/2003 J 18": "#CCCCCC",  # Light gray, small provisional moon
    "Eupheme": "#CCCCCC",  # Light gray, small provisional moon
    "S/2021 J 3": "#CCCCCC",  # Light gray, small provisional moon
    "S/2010 J 2": "#CCCCCC",  # Light gray, small provisional moon
    "S/2016 J 1": "#CCCCCC",  # Light gray, small provisional moon
    "Mneme": "#CCCCCC",  # Light gray, small irregular moon
    "Euanthe": "#CCCCCC",  # Light gray, small irregular moon
    "S/2003 J 16": "#CCCCCC",  # Light gray, small provisional moon
    "S/2009 S 1": "#D3D3D3",  # Light gray, small ring moon of Saturn
    "Pan": "#A9A9A9",  # Dark gray, small rocky moon in Saturn’s rings
    "Daphnis": "#A9A9A9",  # Dark gray, small ring moon
    "Atlas": "#A9A9A9",  # Dark gray, small rocky moon
    "Prometheus": "#B0B0B0",  # Gray, small rocky moon with dusty surface
    "Pandora": "#B0B0B0",  # Gray, small rocky moon with craters
    "Epimetheus": "#A9A9A9",  # Dark gray, rocky co-orbital moon
    "Janus": "#A9A9A9",  # Dark gray, rocky co-orbital moon
    "Aegaeon": "#CCCCCC",  # Light gray, tiny ring moon, limited data
    "Mimas": "#D3D3D3",  # Light gray, icy with large craters
    "Methone": "#F0F8FF",  # Alice blue, smooth icy surface
    "Anthe": "#CCCCCC",  # Light gray, small icy moon
    "Pallene": "#CCCCCC",  # Light gray, small icy moon
    "Enceladus": "#FFFFFF",  # White, bright icy surface with geysers
    "Tethys": "#E0E0E0",  # Light gray, icy with large rift
    "Telesto": "#E0E0E0",  # Light gray, small icy Trojan moon
    "Calypso": "#E0E0E0",  # Light gray, small icy Trojan moon
    "Helene": "#D3D3D3",  # Light gray, small icy Trojan moon
    "Polydeuces": "#CCCCCC",  # Light gray, small icy Trojan moon
    "Dione": "#C0C0C0",  # Silver gray, icy with wispy streaks
    "Rhea": "#B0B0B0",  # Gray, icy with craters
    "Titan": "#FFA500",  # Orange, hazy methane-rich atmosphere
    "Hyperion": "#A9A9A9",  # Dark gray, porous rocky surface
    "Iapetus": "#808080",  # Gray, two-toned (dark and bright regions)
    "S/2023 S 1": "#CCCCCC",  # Light gray, provisional moon, limited data
    "S/2019 S 1": "#CCCCCC",  # Light gray, provisional moon, limited data
    "Cordelia": "#A9A9A9",  # Dark gray, small rocky inner moon of Uranus
    "Ophelia": "#A9A9A9",  # Dark gray, small rocky inner moon
    "S/2025 U 1": "#CCCCCC",  # Light gray, provisional moon, limited data
    "Bianca": "#A9A9A9",  # Dark gray, small rocky inner moon
    "Cressida": "#A9A9A9",  # Dark gray, small rocky inner moon
    "Desdemona": "#A9A9A9",  # Dark gray, small rocky inner moon
    "Juliet": "#A9A9A9",  # Dark gray, small rocky inner moon
    "Portia": "#A9A9A9",  # Dark gray, small rocky inner moon
    "Rosalind": "#A9A9A9",  # Dark gray, small rocky inner moon
    "Belinda": "#A9A9A9",  # Dark gray, small rocky inner moon
    "Puck": "#A9A9A9",  # Dark gray, larger rocky inner moon
    "Miranda": "#D3D3D3",  # Light gray, icy with dramatic cliffs
    "Ariel": "#E0E0E0",  # Light gray, bright icy surface
    "Umbriel": "#808080",  # Dark gray, dark icy surface
    "Titania": "#C0C0C0",  # Silver gray, icy with craters
    "Oberon": "#B0B0B0",  # Gray, icy with large craters
    "Francisco": "#CCCCCC",  # Light gray, small irregular moon
    "Caliban": "#696969",  # Dim gray, dark irregular moon
    "Stephano": "#CCCCCC",  # Light gray, small irregular moon
    "Trinculo": "#CCCCCC",  # Light gray, small irregular moon
    "Sycorax": "#696969",  # Dim gray, dark irregular moon
    "Margaret": "#CCCCCC",  # Light gray, small irregular moon
    "Prospero": "#CCCCCC",  # Light gray, small irregular moon
    "Setebos": "#CCCCCC",  # Light gray, small irregular moon
    "Ferdinand": "#CCCCCC",  # Light gray, small irregular moon
    "Naiad": "#A9A9A9",  # Dark gray, small rocky inner moon of Neptune
    "Thalassa": "#A9A9A9",  # Dark gray, small rocky inner moon
    "Despina": "#A9A9A9",  # Dark gray, small rocky inner moon
    "Galatea": "#A9A9A9",  # Dark gray, small rocky inner moon
    "Larissa": "#A9A9A9",  # Dark gray, small rocky inner moon
    "Hippocamp": "#CCCCCC",  # Light gray, small icy moon
    "Proteus": "#808080",  # Dark gray, large rocky moon
    "Triton": "#ADD8E6",  # Light blue, icy with nitrogen geysers
    "Nereid": "#C0C0C0",  # Silver gray, icy irregular moon
    "Halimede": "#CCCCCC",  # Light gray, small irregular moon
    "Sao": "#CCCCCC",  # Light gray, small irregular moon
    "Laomedeia": "#CCCCCC",  # Light gray, small irregular moon
    "Psamathe": "#CCCCCC",  # Light gray, small irregular moon
    "Neso": "#CCCCCC"  # Light gray, small irregular moon
}

if "MOON_SCALE_FACTORS" not in globals():
    raise ValueError("MOON_SCALE_FACTORS not found — please paste it above this line.")

def solve_kepler(M, e, tol=1e-10):
    E = M
    while True:
        dE = (E - e * math.sin(E) - M) / (1 - e * math.cos(E))
        E -= dE
        if abs(dE) < tol:
            break
    return E

def angular_speed_rad_per_frame_from_km(gm_km3_s2, r_km):
    # physical orbital angular speed (rad/s) for circular orbit: omega = sqrt(GM / r^3)
    omega_phys = np.sqrt(gm_km3_s2 / (r_km ** 3))
    # convert to simulated seconds per real second then to per frame
    return omega_phys * SECONDS_PER_FRAME

# ---------------------------- HELPER FUNCTIONS ----------------------------
def km_to_scene_default(km, planet_radius=SATURN_RADIUS):
    return (km / REAL_SATURN_MEAN_RADIUS_KM) * planet_radius

def calculate_keplerian_speed(moon_radius_km, gm):
    r_km = moon_radius_km  # Use real-world km directly
    omega = np.sqrt(gm / (r_km ** 3)) * SPEED_SCALE / FPS
    return omega
def create_moon_sphere(actor_radius, name="moon"):
    """
    Creates a realistic, physically solid moon with cratered surface texture
    using procedural 3D noise and lighting.
    """

    # Create high-resolution sphere (smooth geometry)
    sphere = pv.Sphere(radius=actor_radius, theta_resolution=96, phi_resolution=96)

    # Ensure points are float for filters (fixes your warning)
    sphere.points = sphere.points.astype(np.float32)

    # Convert Cartesian coordinates to spherical
    x, y, z = sphere.points[:, 0], sphere.points[:, 1], sphere.points[:, 2]
    r = np.sqrt(x**2 + y**2 + z**2)
    lon = np.arctan2(y, x)
    lat = np.arcsin(z / r)

    # --- Procedural texture generation (crater + surface roughness) ---
    noise = (0.3 * np.sin(3 * lon) * np.cos(5 * lat) +
             0.25 * np.sin(7 * lon + 3 * lat) +
             0.2 * np.sin(11 * lon - 5 * lat))
    noise += np.random.normal(0, 0.08, size=noise.shape)

    # Add crater-like depressions
    crater_mask = (np.sin(10 * lon)**2 + np.cos(15 * lat)**2)
    crater_depth = 0.15 * crater_mask * np.random.rand(len(lon))
    heightmap = 1 + 0.3 * noise - 0.1 * crater_depth

    # Perturb geometry to add real bumps
    sphere.points *= heightmap[:, None].astype(np.float32)

    # Compute color variation (grayscale + subtle tint)
    base_color = np.array([190, 190, 190]) / 255.0  # silver-gray regolith
    albedo = np.clip(0.8 + 0.4 * noise, 0, 1)
    colors = (base_color[None, :] * albedo[:, None])
    colors = np.clip(colors, 0, 1)
    colors = (colors * 255).astype(np.uint8)
    sphere.point_data["colors"] = colors

    # Add the moon actor with proper lighting
    actor = pl.add_mesh(
        sphere,
        scalars="colors",
        rgb=True,
        smooth_shading=True,
        ambient=0.25,
        diffuse=0.65,
        specular=0.1,
        name=name
    )

    return actor


def datetime_to_julian_day(dt):
    # dt must be timezone-aware UTC
    year = dt.year
    month = dt.month
    day = dt.day + (dt.hour + (dt.minute + dt.second / 60.0) / 60.0) / 24.0
    if month <= 2:
        year -= 1
        month += 12
    A = year // 100
    B = 2 - A + (A // 4)
    JD = int(365.25 * (year + 4716)) + int(30.6001 * (month + 1)) + day + B - 1524.5
    return JD

def gmst_from_jd(jd):
    """
    Compute Greenwich Mean Sidereal Time (GMST) in degrees from Julian Date.
    Formula from USNO/NASA: Accurate to ~0.1 arcsec for 1800–2200.
    """
    dj = jd - 2451545.0
    t = dj / 36525.0
    theta = (280.46061837 + 360.98564736629 * dj +
             0.000387933 * t**2 - t**3 / 38710000.0) % 360.0
    return theta

def get_current_jd():
    """
    Get current UTC JD (live, or fixed for testing).
    """
    now = datetime.now(timezone.utc)
    return datetime_to_julian_day(now)

def solve_kepler(M, e, tol=1e-10, max_it=50):
    # M in radians, returns E in radians
    if e < 0.8:
        E = M
    else:
        E = math.pi
    for _ in range(max_it):
        f = E - e*math.sin(E) - M
        fprime = 1 - e*math.cos(E)
        if abs(fprime) < 1e-12:
            break
        dE = -f / fprime
        E += dE
        if abs(dE) < tol:
            break
    return E

def kepler_to_state(a_km, e, i_deg, Omega_deg, omega_deg, M0_deg, epoch_jd, target_jd, mu=GM_SUN):
    """
    Compute position [km] and velocity vector [km/s] from osculating elements.
    Handles elliptic (e<1), parabolic (e~1), hyperbolic (e>1, a<0).
    Returns: x_km, y_km, z_km, vel_km_s (array [vx,vy,vz]), v_deg
    """
    delta_t_sec = (target_jd - epoch_jd) * 86400.0
    a_abs = abs(a_km)
    if a_abs == 0:
        return 0, 0, 0, np.array([0,0,0]), 0  # Degenerate fallback

    # Mean motion (rad/sec; positive for all)
    n = np.sqrt(mu / a_abs**3)
    M0_rad = np.radians(M0_deg)
    M = (M0_rad + n * delta_t_sec) % (2 * np.pi)  # Always add time; mod 2pi

    # Solve for true anomaly nu (radians)
    if abs(e - 1) < 1e-6:  # Parabolic approx
        # Barker's eq: D^3/3 + D - M = 0; solve for D = tan(nu/2)
        def barker_eq(D): return D**3 / 3 + D - M
        D = newton(barker_eq, M)  # Robust solver
        nu = 2 * np.arctan(D)
    elif e < 1:  # Elliptic
        def kepler_eq(E): return E - e * np.sin(E) - M
        E_guess = M if e < 0.8 else np.pi
        E = newton(kepler_eq, E_guess, tol=1e-12)
        nu = 2 * np.arctan2(np.sqrt(1 + e) * np.sin(E / 2), np.sqrt(1 - e) * np.cos(E / 2))
    else:  # Hyperbolic
        def hyp_kepler_eq(F): return e * np.sinh(F) - F - M
        # Better guess: asymptotic for large |M|
        F_guess = np.sign(M) * np.log(np.abs(2 * M / e) + 1.8) if np.abs(M) > 1e-6 else M
        F = newton(hyp_kepler_eq, F_guess, tol=1e-12)
        nu = 2 * np.arctan(np.sqrt((e + 1) / (e - 1)) * np.tanh(F / 2))

    v_deg = np.degrees(nu)
    cos_nu, sin_nu = np.cos(nu), np.sin(nu)

    # Radial distance r (km)
    if e < 1:
        r = a_abs * (1 - e**2) / (1 + e * cos_nu)
    else:  # Hyperbolic: a negative, but |a| * (e^2 - 1)
        r = a_abs * (e**2 - 1) / (1 + e * cos_nu)
    r = max(r, 1e-6)  # Min r guard

    # Perifocal position/velocity
    pos_peri = r * np.array([cos_nu, sin_nu, 0.0])
    v_mag = np.sqrt(mu * (2 / r - 1 / a_km))  # Vis-viva (valid for all conics)
    vel_peri = v_mag * np.array([-sin_nu, (e + cos_nu) / (1 + e * cos_nu) * (1 + e * cos_nu), 0.0])  # Corrected angular vel

    # ECI transformation (standard rotation)
    O, inc, w = map(np.radians, [Omega_deg, i_deg, omega_deg])
    Rz_w = np.array([[np.cos(w), -np.sin(w), 0], [np.sin(w), np.cos(w), 0], [0, 0, 1]])
    Rx_inc = np.array([[1, 0, 0], [0, np.cos(inc), -np.sin(inc)], [0, np.sin(inc), np.cos(inc)]])
    Rz_O = np.array([[np.cos(O), -np.sin(O), 0], [np.sin(O), np.cos(O), 0], [0, 0, 1]])
    R = np.dot(Rz_O, np.dot(Rx_inc, Rz_w))
    pos_eci = np.dot(R, pos_peri)
    vel_eci = np.dot(R, vel_peri)

    # Finite guard
    if not np.all(np.isfinite(pos_eci)) or not np.all(np.isfinite(vel_eci)):
        pos_eci = np.zeros(3)
        vel_eci = np.zeros(3)

    return pos_eci[0], pos_eci[1], pos_eci[2], vel_eci, v_deg


def solve_kepler_hyperbolic(M, e, tol=1e-12, max_it=100):
    # M in radians, returns H (hyperbolic anomaly) in radians
    sign_M = 1 if M > 0 else -1
    H = sign_M * math.log(2 * abs(M) / e + 1.8)  # Initial guess
    for _ in range(max_it):
        f = e * math.sinh(H) - H - M
        fprime = e * math.cosh(H) - 1
        if abs(fprime) < 1e-12:
            break
        dH = -f / fprime
        H -= dH
        if abs(dH) < tol:
            break
    return H

def set_body_position_from_elements(name, actor, elements, target_jd, mu=GM_SUN):
    x_km, y_km, z_km, _, _ = kepler_to_state(
        elements["a_km"],  # <-- Change to this (direct use of "a_km")
        elements["e"],
        elements["i_deg"],
        elements["Omega_deg"],
        elements["omega_deg"],
        elements["M0_deg"],
        elements["epoch_jd"],
        target_jd,
        mu=mu
    )
    pos_scene = np.array([x_km, y_km, z_km]) * KM_TO_SCENE
    actor.SetPosition(*pos_scene)
    
orbital_elements = {
    "mercury": {
        "a_km": 0.38709927 * AU_KM,
        "e": 0.20563593,
        "i_deg": 7.00497902,
        "Omega_deg": 48.33076593,
        "omega_deg": 29.12703035,
        "M0_deg": 174.79252722,
        "epoch_jd": J2000,
    },
    "venus": {
        "a_km": 0.72333566 * AU_KM,
        "e": 0.00677672,
        "i_deg": 3.39467605,
        "Omega_deg": 76.67984255,
        "omega_deg": 54.92262463,
        "M0_deg": 50.37663232,
        "epoch_jd": J2000,
    },
    "earth": {  # Earth-Moon Barycenter
        "a_km": 1.00000261 * AU_KM,
        "e": 0.01671123,
        "i_deg": -0.00001531,
        "Omega_deg": 0.0,
        "omega_deg": 102.93768193,
        "M0_deg": 357.52688973,
        "epoch_jd": J2000,
    },
    "mars": {
        "a_km": 1.52371034 * AU_KM,
        "e": 0.09339410,
        "i_deg": 1.84969142,
        "Omega_deg": 49.55953891,
        "omega_deg": 286.4968315,
        "M0_deg": 19.39019754,
        "epoch_jd": J2000,
    },
    "jupiter": {
        "a_km": 5.20288700 * AU_KM,
        "e": 0.04838624,
        "i_deg": 1.30439695,
        "Omega_deg": 100.47390909,
        "omega_deg": 274.25457074,
        "M0_deg": 19.66796068,
        "epoch_jd": J2000,
    },
    "saturn": {
        "a_km": 9.53667594 * AU_KM,
        "e": 0.05386179,
        "i_deg": 2.48599187,
        "Omega_deg": 113.66242448,
        "omega_deg": 338.93645383,
        "M0_deg": 317.35536592,
        "epoch_jd": J2000,
    },
    "uranus": {
        "a_km": 19.18916464 * AU_KM,
        "e": 0.04725744,
        "i_deg": 0.77263783,
        "Omega_deg": 74.01692503,
        "omega_deg": 96.93735127,
        "M0_deg": 142.28382821,
        "epoch_jd": J2000,
    },
    "neptune": {
        "a_km": 30.06992276 * AU_KM,
        "e": 0.00859048,
        "i_deg": 1.77004347,
        "Omega_deg": 131.78422574,
        "omega_deg": 273.18053653,
        "M0_deg": 259.91520804,
        "epoch_jd": J2000,
    },
    "moon": {  # Geocentric
        "a_km": 384400.0,
        "e": 0.0554,
        "i_deg": 5.16,
        "Omega_deg": 125.08,
        "omega_deg": 318.15,
        "M0_deg": 135.27,
        "epoch_jd": J2000,
    },
    "pluto": {
        "a_km": PLUTO_AU * AU_KM,
        "e": 0.2488,
        "i_deg": 17.16,
        "Omega_deg": 110.30,
        "omega_deg": 113.78,
        "M0_deg": 14.53,
        "epoch_jd": J2000,
    },
    "eris": {
        "a_km": ERIS_AU * AU_KM,
        "e": 0.4407,
        "i_deg": 44.04,
        "Omega_deg": 35.19,
        "omega_deg": 54.07,
        "M0_deg": 330.00,
        "epoch_jd": J2000,
    },
    "haumea": {
        "a_km": HAUMEA_AU * AU_KM,
        "e": 0.1951,
        "i_deg": 28.22,
        "Omega_deg": 338.07,
        "omega_deg": 240.85,
        "M0_deg": 359.00,
        "epoch_jd": J2000,
    },
    "makemake": {
        "a_km": MAKEMAKE_AU * AU_KM,
        "e": 0.159,
        "i_deg": 29.01,
        "Omega_deg": 150.35,
        "omega_deg": 148.72,
        "M0_deg": 289.00,
        "epoch_jd": J2000,
    },
    "ceres": {
        "a_km": CERES_AU * AU_KM,
        "e": 0.0758,
        "i_deg": 10.59,
        "Omega_deg": 80.33,
        "omega_deg": 73.60,
        "M0_deg": 95.00,
        "epoch_jd": J2000,
    },
}

def procedural_moon_actor(sphere, base_hex, name, ambient=0.2, diffuse=0.6, specular=0.1, emissive=False):
    # Use predefined color if available for the moon name
    if name in moon_colors:
        base_hex = moon_colors[name]
   
    if not (isinstance(base_hex, str) and base_hex.startswith("#") and len(base_hex) == 7):
        base_hex = "#cccccc"
    try:
        base_rgb = np.array([int(base_hex[i:i+2], 16) for i in (1, 3, 5)]) / 255.0
    except Exception:
        base_rgb = np.array([0.8, 0.8, 0.8])
    pts = sphere.points.copy()
    norms = np.linalg.norm(pts, axis=1, keepdims=True) + 1e-9
    x, y, z = pts[:, 0]/norms[:, 0], pts[:, 1]/norms[:, 0], pts[:, 2]/norms[:, 0]
    lon = np.arctan2(y, x)
    lat = np.arcsin(z)
    noise = (0.35*np.sin(3*lon)*np.cos(5*lat) +
             0.25*np.sin(7*lon + 3*lat) +
             0.15*np.sin(11*lon - 5*lat))
    rng = np.random.default_rng(abs(hash(name)) % (2**32))
    noise += rng.normal(0, 0.05, size=noise.shape)
    noise = (noise - noise.min()) / (noise.max() - noise.min() + 1e-9)
    albedo = 0.7 + 0.3 * noise
    lname = name.lower()
    # Initialize final_rgb with default value
    final_rgb = np.clip(base_rgb * albedo[:, None], 0, 1)
    if "io" in lname:
        base_rgb = np.array([0.93, 0.74, 0.24])
        albedo *= 0.9 + 0.2*np.sin(6*lon + 3*lat)
        final_rgb = np.clip(base_rgb * albedo[:, None], 0, 1)
    elif "europa" in lname:
        base_rgb = np.array([0.88, 0.93, 1.0])
        albedo *= 0.85 + 0.3*np.sin(8*lon)*np.cos(4*lat)
        final_rgb = np.clip(base_rgb * albedo[:, None], 0, 1)
    elif "ganymede" in lname:
        base_rgb = np.array([0.72, 0.72, 0.72])
        albedo *= 0.8 + 0.2*np.sin(5*lon - 2*lat)
        final_rgb = np.clip(base_rgb * albedo[:, None], 0, 1)
    elif "callisto" in lname:
        base_rgb = np.array([0.55, 0.47, 0.41])
        albedo *= 0.75 + 0.25*np.sin(4*lon + 3*lat)
        final_rgb = np.clip(base_rgb * albedo[:, None], 0, 1)
    elif "triton" in lname:
        base_rgb = np.array([0.95, 0.90, 0.93])
        tint_dark = np.array([0.55, 0.45, 0.50])
        noise_main = (0.35*np.sin(5*lon + 2*lat) +
                      0.25*np.sin(9*lat - 3*lon) +
                      0.15*np.sin(12*lon + 7*lat))
        geyser = np.sin(10*lon) * np.cos(4*lat)
        dark_mask = (geyser > 0.85).astype(float)
        noise_main = (noise_main - noise_main.min()) / (noise_main.max() - noise_main.min() + 1e-9)
        color_mix = base_rgb * (0.85 + 0.15 * noise_main[:, None])
        color_mix = np.where(dark_mask[:, None] > 0.5, tint_dark, color_mix)
        color_mix *= (0.85 + 0.15 * np.cos(lat)[:, None])
        ambient, diffuse, specular = 0.30, 0.30, 0.00
        final_rgb = np.clip(color_mix, 0, 1)
    elif "luna" in lname or "moon" in lname: # Earth's Moon
        base_rgb = np.array([0.75, 0.75, 0.75]) # Grayish for Luna
        albedo *= 0.8 + 0.2*np.sin(5*lon)*np.cos(3*lat) # Crater-like texture
        final_rgb = np.clip(base_rgb * albedo[:, None], 0, 1)
    # Add customs for Saturn moons if desired (e.g., Titan hazy orange)
    elif "titan" in lname:
        base_rgb = np.array([0.95, 0.75, 0.45]) # Orangish for hazy atmosphere
        albedo *= 0.7 + 0.3*np.sin(4*lon)*np.cos(6*lat) # Surface features
        final_rgb = np.clip(base_rgb * albedo[:, None], 0, 1)
    elif "enceladus" in lname:
        base_rgb = np.array([0.95, 0.95, 1.0]) # Icy white
        albedo *= 0.9 + 0.1*np.sin(10*lon + 5*lat) # Tiger stripes
        final_rgb = np.clip(base_rgb * albedo[:, None], 0, 1)
    elif "iapetus" in lname:
        # Two-toned: dark leading hemisphere, bright trailing
        dark_rgb = np.array([0.3, 0.3, 0.3]) # Dark side
        bright_rgb = np.array([0.9, 0.9, 0.9]) # Bright side
        # Mask for leading hemisphere (dark)
        dark_mask = np.cos(lon) > 0 # Approximate leading side
        albedo *= 0.7 + 0.3 * np.sin(4*lon + 2*lat) # Cratered texture
        final_rgb = np.where(dark_mask[:, None], dark_rgb * albedo[:, None], bright_rgb * albedo[:, None])
        final_rgb = np.clip(final_rgb, 0, 1)
        specular = 0.05 # Low specular for rough surface
    elif "rhea" in lname:
        base_rgb = np.array([0.7, 0.7, 0.75]) # Grayish icy
        albedo *= 0.85 + 0.15 * np.sin(6*lon) * np.cos(3*lat) # Craters and rays
        final_rgb = np.clip(base_rgb * albedo[:, None], 0, 1)
    elif "dione" in lname:
        base_rgb = np.array([0.75, 0.75, 0.8]) # Silvery icy
        albedo *= 0.8 + 0.2 * np.sin(7*lon - 4*lat) # Wispy terrain approximation
        final_rgb = np.clip(base_rgb * albedo[:, None], 0, 1)
    elif "mimas" in lname:
        base_rgb = np.array([0.85, 0.85, 0.9]) # Icy white
        # Approximate Herschel crater with a large low-frequency feature
        crater = 0.5 * np.sin(2*lon) * np.cos(lat) # Simple polar crater sim
        albedo *= 0.9 + 0.1 * crater + 0.2 * noise
        final_rgb = np.clip(base_rgb * albedo[:, None], 0, 1)
    elif "phobos" in lname:
        base_rgb = np.array([0.65, 0.2, 0.05]) # Reddish
        albedo *= 0.6 + 0.4 * np.sin(4*lon + lat) # Cratered, irregular
        final_rgb = np.clip(base_rgb * albedo[:, None], 0, 1)
    elif "deimos" in lname:
        base_rgb = np.array([0.5, 0.5, 0.5]) # Gray
        albedo *= 0.7 + 0.3 * np.sin(5*lon - 2*lat) # Dusty craters
        final_rgb = np.clip(base_rgb * albedo[:, None], 0, 1)
    elif "miranda" in lname:
        base_rgb = np.array([0.8, 0.8, 0.85]) # Light gray icy
        # Chaotic terrain: high-frequency noise
        chaos_noise = 0.4 * np.sin(10*lon) * np.cos(6*lat)
        albedo *= 0.75 + 0.25 * chaos_noise + 0.2 * noise
        final_rgb = np.clip(base_rgb * albedo[:, None], 0, 1)
    elif "titania" in lname or "oberon" in lname:
        base_rgb = np.array([0.75, 0.75, 0.8]) # Silvery gray
        albedo *= 0.8 + 0.2 * np.sin(5*lon) * np.cos(4*lat) # Cratered icy
        final_rgb = np.clip(base_rgb * albedo[:, None], 0, 1)
    elif "proteus" in lname:
        base_rgb = np.array([0.5, 0.5, 0.55]) # Dark gray rocky
        albedo *= 0.65 + 0.35 * np.sin(3*lon + 2*lat) # Irregular craters
        final_rgb = np.clip(base_rgb * albedo[:, None], 0, 1)
    # Pluto moons
    elif "charon" in lname:
        base_rgb = np.array([0.65, 0.65, 0.7])  # Grayish-red icy
        albedo *= 0.75 + 0.25 * np.sin(6*lon) * np.cos(3*lat)  # Craters, reddish tint
        final_rgb = np.clip(base_rgb * albedo[:, None], 0, 1)
    elif "styx" in lname:
        base_rgb = np.array([0.4, 0.4, 0.45])  # Dark gray rocky
        albedo *= 0.6 + 0.4 * np.sin(4*lon + lat)  # Cratered irregular
        final_rgb = np.clip(base_rgb * albedo[:, None], 0, 1)
    elif "nix" in lname:
        base_rgb = np.array([0.55, 0.55, 0.6])  # Medium gray icy
        albedo *= 0.7 + 0.3 * np.sin(5*lon - 2*lat)  # Dusty surface
        final_rgb = np.clip(base_rgb * albedo[:, None], 0, 1)
    elif "kerberos" in lname:
        base_rgb = np.array([0.45, 0.45, 0.5])  # Dark gray
        albedo *= 0.65 + 0.35 * np.sin(3*lon + 2*lat)  # Rough craters
        final_rgb = np.clip(base_rgb * albedo[:, None], 0, 1)
    elif "hydra" in lname:
        base_rgb = np.array([0.6, 0.6, 0.65])  # Light gray icy
        albedo *= 0.75 + 0.25 * np.sin(7*lon) * np.cos(4*lat)  # Elongated craters
        final_rgb = np.clip(base_rgb * albedo[:, None], 0, 1)
    # For small/provisional moons, use default procedural noise with low variation
    elif any(s in lname for s in ["/", "s/"]):
        albedo *= 0.8 + 0.2 * noise # Minimal texture
        final_rgb = np.clip(base_rgb * albedo[:, None], 0, 1)
        specular = 0.05 # Matte
    # Add more for other major moons like Rhea, Iapetus if needed
    sphere.point_data["colors"] = (final_rgb * 255).astype(np.uint8)
    moon_actor = pl.add_mesh(
        sphere, scalars="colors", rgb=True, smooth_shading=True,
        ambient=ambient, diffuse=diffuse, specular=specular, specular_power=10, name=name, lighting=True
    )
    if emissive:
        prop = moon_actor.GetProperty()
        prop.SetAmbient(min(1.0, ambient + 0.15))
        prop.SetDiffuse(diffuse * 0.9)
    return moon_actor

def generate_moons(planet_name):
    moons = []
    raw_moons = MOON_SCALE_FACTORS.get(planet_name, [])
    mu = GM_PLANET.get(planet_name, GM_EARTH)  # km³/s²

    for moon in raw_moons:
        a = moon.get("a_km", 0)
        if a <= 0:
            continue

        # Mean motion n = sqrt(μ / a³) in rad/s
        n = math.sqrt(mu / (a ** 3))
        M0_rad = math.radians(moon.get("M0_deg", 0))

        enhanced_moon = moon.copy()
        enhanced_moon.update({
            "omega": n,
            "M": M0_rad,
            "tilt_applied": False
        })
        moons.append(enhanced_moon)

    return moons


# ---------------------------- ASTEROID BELT GENERATION ----------------------------
def update_asteroid_opacity(caller, event):
    try:
        cam = pl.camera
        cam_pos = np.array(cam.position)
        cam_dir = np.array(cam.direction)
        cam_dir /= np.linalg.norm(cam_dir)

        # Belt roughly centered near 2.7 AU
        belt_center = np.array([AU_SCALE * 2.7, 0.0, 0.0])
        dist = np.linalg.norm(cam_pos - belt_center)

        # --- Distance-based fade parameters ---
        visibility_threshold = 0.8 * AU_SCALE
        fade_zone = visibility_threshold * 0.1

        if dist <= visibility_threshold - fade_zone:
            dist_opacity = 1.0
        elif dist >= visibility_threshold + fade_zone:
            dist_opacity = 0.0
        else:
            t = (dist - (visibility_threshold - fade_zone)) / (2 * fade_zone)
            dist_opacity = 1.0 - t

        # --- Angular visibility setup ---
        vectors_to_asteroids = positions - cam_pos
        dists = np.linalg.norm(vectors_to_asteroids, axis=1)
        unit_vectors = vectors_to_asteroids / dists[:, None]

        cos_angles = np.clip(np.sum(unit_vectors * cam_dir, axis=1), -1.0, 1.0)
        angles = np.arccos(cos_angles)

        visible_angle = np.pi / 3   # ±60° visibility cone
        fade_zone_angle = np.radians(10)  # 10° fade zone at edges

        angular_opacity = np.ones_like(angles)
        fade_mask = (angles > (visible_angle - fade_zone_angle)) & (angles < visible_angle)
        outside = angles >= visible_angle

        # Smooth angular fade near edge
        t = (angles[fade_mask] - (visible_angle - fade_zone_angle)) / fade_zone_angle
        angular_opacity[fade_mask] = 1.0 - t
        angular_opacity[outside] = 0.0

        # Combine both effects (distance + angle)
        final_opacity = dist_opacity * angular_opacity.mean()  # global belt fade strength

        # Apply to belt actor
        asteroid_actor.GetProperty().SetOpacity(float(final_opacity))

    except Exception as e:
        print(f"[WARN] Asteroid opacity update error: {e}")



ASTEROID_COUNT = 200
ASTEROID_INNER_AU = 2.10
ASTEROID_OUTER_AU = 3.30

ASTEROID_INNER_R = AU_SCALE * ASTEROID_INNER_AU
ASTEROID_OUTER_R = AU_SCALE * ASTEROID_OUTER_AU
rng = np.random.default_rng(42)
positions = []
colors = []
sizes = []
for _ in range(ASTEROID_COUNT):
    r = rng.uniform(ASTEROID_INNER_R, ASTEROID_OUTER_R)
    theta = rng.uniform(0, 2 * np.pi)
    inc = rng.uniform(-np.radians(5), np.radians(5))
    z = r * np.sin(inc) * 0.04
    x = r * np.cos(theta)
    y = r * np.sin(theta)
    positions.append([x, y, z])
    reflect_factor = 1.0 - (r - ASTEROID_INNER_R) / (ASTEROID_OUTER_R - ASTEROID_INNER_R)
    reflect_factor = np.clip(reflect_factor, 0.3, 1.0)  # clamp to avoid extremes
    g = rng.uniform(0.55, 0.85) * reflect_factor
    colors.append([g, g * 0.92, g * 0.84])
    sizes.append(rng.uniform(0.5, 1))
positions = np.array(positions)
colors = (np.array(colors) * 255).astype(np.uint8)
sizes = np.array(sizes )
asteroid_cloud = pv.PolyData(positions)
asteroid_cloud.point_data["rgb"] = colors
asteroid_cloud.point_data["size"] = sizes 
asteroid_actor = pl.add_points(
    asteroid_cloud, scalars="rgb", rgb=True, point_size=2, style='points',  # Set point_size to 5.0
    emissive=True, ambient=0, diffuse=0.4, name="AsteroidBelt", pickable=False
)
#asteroid_actor.GetMapper().SetVBOShiftScaleMethod(False)
asteroid_actor._belt_angle = 0.0
asteroid_actor._belt_speed = 0.03


# --- Dynamic fade for asteroid belt based on camera distance ---
if hasattr(pl.render_window, "GetInteractor"):
    interactor = pl.render_window.GetInteractor()
    if hasattr(interactor, "AddObserver"):
        interactor.AddObserver("RenderEvent", update_asteroid_opacity)
# ---------------------------- SETUP PLOTTER ----------------------------

pv.set_new_attribute(pl, 'pickpoint', np.array([0.0, 0.0, 0.0]))  # Fixed: Use set_new_attribute

# Background
space_bg = examples.download_cubemap_space_4k()
pl.set_environment_texture(space_bg, is_srgb=True, resample=True)
skybox = space_bg.to_skybox()
skybox.SetScale(1e6)
skybox.GetProperty().LightingOff()
skybox.GetProperty().SetOpacity(1.0)
skybox.PickableOff()
pl.add_actor(skybox)

def update_skybox_orientation(caller, event):
    cam = pl.camera
    try:
        orient = cam.GetOrientation()
        skybox.SetOrientation(orient)
    except Exception:
        pass

interactor = pl.render_window.GetInteractor()
if hasattr(interactor, "AddObserver"):
    interactor.AddObserver("RenderEvent", update_skybox_orientation)
elif hasattr(interactor, "add_observer"):
    interactor.add_observer("RenderEvent", update_skybox_orientation)

# Lighting
planet_light = pv.Light(
    light_type='scene light', position=(0.0, 0.0, 0.0), color=(1.0, 0.95, 0.8),
    intensity=0.06, positional=True, cone_angle=180, attenuation_values=(0.0, 1e-6, 1e-30)
)
pl.add_light(planet_light)
ambient_fill = pv.Light(light_type='headlight', color=(0.12, 0.12, 0.14), intensity=0.09)
pl.add_light(ambient_fill)

prograde_sign = -1

# Sun
sun_mesh = examples.planets.load_sun(radius=SUN_RADIUS)
sun_texture = examples.planets.download_sun_surface(texture=True)
sun_actor = pl.add_mesh(
    sun_mesh, texture=sun_texture, emissive=True, lighting=False,
    ambient=1.0, specular=0.0, name='sun'
)
sun_actor.RotateX(AXIAL_TILTS["sun"])

# Planets
jupiter_mesh = examples.planets.load_jupiter(radius=JUPITER_RADIUS)
jupiter_texture = examples.planets.download_jupiter_surface(texture=True)
jupiter_actor = pl.add_mesh(
    jupiter_mesh, texture=jupiter_texture, smooth_shading=True,
    ambient=0.1, diffuse=0.1, specular=0.01, name='jupiter'
)
jupiter_actor.RotateX(prograde_sign * AXIAL_TILTS["jupiter"])

saturn_mesh = examples.planets.load_saturn(radius=SATURN_RADIUS)
saturn_texture = examples.planets.download_saturn_surface(texture=True)
saturn_actor = pl.add_mesh(
    saturn_mesh, texture=saturn_texture, smooth_shading=True,
    ambient=0.12, diffuse=0.16, specular=0.01, name='saturn'
)
saturn_actor.RotateX(prograde_sign * AXIAL_TILTS["saturn"])

neptune_mesh = examples.planets.load_neptune(radius=NEPTUNE_RADIUS)
neptune_texture = examples.planets.download_neptune_surface(texture=True)
neptune_actor = pl.add_mesh(
    neptune_mesh, texture=neptune_texture, smooth_shading=True,
    ambient=0.1, diffuse=0.4, specular=0, name='neptune'
)
neptune_actor.RotateX(prograde_sign * AXIAL_TILTS["neptune"])

uranus_mesh = examples.planets.load_uranus(radius=URANUS_RADIUS)
uranus_texture = examples.planets.download_uranus_surface(texture=True)
uranus_actor = pl.add_mesh(
    uranus_mesh, texture=uranus_texture, smooth_shading=True,
    ambient=0.1, diffuse=0.3, specular=0, name='uranus'
)
uranus_actor.RotateX(prograde_sign * AXIAL_TILTS["uranus"])


mars_mesh = examples.planets.load_mars(radius=MARS_RADIUS)
mars_texture = examples.planets.download_mars_surface(texture=True)
mars_actor = pl.add_mesh(
    mars_mesh, texture=mars_texture, smooth_shading=True,
    ambient=0.15, diffuse=0.05, specular=0.0, name='mars'
)
mars_actor.RotateX(prograde_sign * AXIAL_TILTS["mars"])

earth_mesh = examples.planets.load_earth(radius=EARTH_RADIUS)
earth_texture = examples.load_globe_texture()
earth_actor = pl.add_mesh(
    earth_mesh, texture=earth_texture, smooth_shading=True,
    ambient=0.15, diffuse=0.07, specular=0.0, name='earth'
)

earth_actor.RotateX(prograde_sign * AXIAL_TILTS["earth"])

earth_actor.SetPosition(EARTH_ORBIT_RADIUS, 0.0, 0.0)

# === MERCURY ===
mercury_mesh = examples.planets.load_mercury(radius=MERCURY_RADIUS)
mercury_texture = examples.planets.download_mercury_surface(texture=True)
mercury_actor = pl.add_mesh(

    mercury_mesh,
    texture=mercury_texture,
    smooth_shading=True,
    ambient=0.1,
    diffuse=0.01,
    specular=0.0,
    name='mercury'
)
mercury_actor.RotateX(prograde_sign * AXIAL_TILTS["mercury"])
mercury_actor.SetPosition(MERCURY_ORBIT_RADIUS, 0.0, 0.0)


# Venus
venus_mesh = examples.planets.load_venus(radius=VENUS_RADIUS)
venus_texture = examples.planets.download_venus_surface(texture=True)
venus_actor = pl.add_mesh(
    venus_mesh, texture=venus_texture, smooth_shading=True,
    ambient=0.15, diffuse=0.014, specular=0.0, name='venus'
)
venus_actor.RotateX(prograde_sign * AXIAL_TILTS["venus"])
venus_actor.SetPosition(VENUS_ORBIT_RADIUS, 0.0, 0.0)

# Saturn Rings
ring_meshes = []
ring_definitions = [
    (RING_C_INNER, RING_C_OUTER, "ring_c", (0.25, 0.45)),
    (RING_B_INNER, RING_B_OUTER, "ring_b", (0.65, 0.9)),
    (RING_CASSINI_INNER, RING_CASSINI_OUTER, "ring_cassini", (0.05, 0.15)),
    (RING_A_INNER, RING_A_OUTER, "ring_a", (0.45, 0.7))
]
sun_dir = np.array([1.0, 0.2, 0.3]); sun_dir /= np.linalg.norm(sun_dir)
for inner_km, outer_km, name, opacity_range in ring_definitions:
    inner_r = km_to_scene_default(inner_km)
    outer_r = km_to_scene_default(outer_km)
    ring = pv.Disc(center=(0, 0, 0), inner=inner_r, outer=outer_r, r_res=2, c_res=1600)
    rpoints = ring.points.copy()
    rdist = np.linalg.norm(rpoints[:, :2], axis=1)
    rnorm = (rdist - rdist.min()) / max(1e-9, (rdist.max() - rdist.min()))
    rng = np.random.default_rng(abs(hash(name)) % (2**32))
    fine_grain = 0.15 * np.sin(100 * rnorm + rng.random()) + 0.1 * np.sin(400 * rnorm + rng.random()) + 0.05 * rng.normal(0, 1, len(rnorm))
    coarse_band = 0.3 * np.sin(10 * rnorm + rng.random() * 5)
    grain = 0.5 * fine_grain + 0.5 * coarse_band
    grain = np.clip(grain, -0.25, 0.25)
    base_color = np.stack([
        0.95 - 0.3 * rnorm + 0.4 * grain,
        0.92 - 0.35 * rnorm + 0.3 * grain,
        0.87 - 0.4 * rnorm + 0.2 * grain
    ], axis=1)
    base_color = np.clip(base_color, 0, 1)
    normals = np.tile(np.array([0, 0, 1]), (len(rpoints), 1))
    lambert = np.clip(0.5 + 0.5 * np.dot(normals, sun_dir), 0, 1).reshape(-1, 1)
    ring.point_data['colors'] = (base_color * lambert * 255).astype(np.uint8)
    ring_actor = pl.add_mesh(
        ring, scalars='colors', rgb=True, smooth_shading=True,
        ambient=0.8, diffuse=1, specular=1, specular_power=60,
        name=f'saturn_{name}', opacity=float(opacity_range[0]), lighting=True
    )
    ring_actor.SetOrientation(26.7, 0.0, 0.0)
    ring_meshes.append(ring_actor)

# Uranus Rings
uranus_ring_actors = []
uranus_real_rings = [
    ("ε", 51100, 52600, "#d1d4d8", 0.60),
    ("α", 44700, 45000, "#999999", 0.50),
    ("β", 45600, 46000, "#999999", 0.50),
]
URANUS_REF_KM = 25559.0
for name, inner_km, outer_km, color, opacity in uranus_real_rings:
    inner_r = (inner_km / URANUS_REF_KM) * URANUS_RADIUS
    outer_r = (outer_km / URANUS_REF_KM) * URANUS_RADIUS
    ring = pv.Disc(inner=inner_r, outer=outer_r, r_res=64, c_res=512)
    ring.rotate_x(98, inplace=True)
    pts = ring.points
    dist = np.linalg.norm(pts[:, :2], axis=1)
    norm_dist = (dist - dist.min()) / (dist.max() - dist.min() + 1e-9)
    rng_local = np.random.default_rng(abs(hash(name)) % (2**32))
    grain = 0.08 * np.sin(150 * norm_dist) + 0.03 * rng_local.normal(0, 1, len(norm_dist))
    base = np.array([int(color[1:3], 16), int(color[3:5], 16), int(color[5:7], 16)]) / 255.0
    rgb = (base * (0.8 + 0.2 * grain[:, None])).clip(0, 1)
    ring.point_data["colors"] = (rgb * 255).astype(np.uint8)
    actor = pl.add_mesh(
        ring, scalars="colors", rgb=True, opacity=opacity, smooth_shading=True,
        ambient=0.35, diffuse=0.8, specular=1, specular_power=80, name=f"uranus_ring_{name}"
    )
    uranus_ring_actors.append(actor)

# Dwarf Planets (with Pluto texture; others procedural)
dwarf_planets = {
    "pluto": None,  # Special handling below
    "eris": (ERIS_RADIUS, "#B0C4DE"),  # Light steel blue
    "haumea": (HAUMEA_RADIUS, "#F0E68C"),  # Khaki for rocky/icy
    "makemake": (MAKEMAKE_RADIUS, "#8B4513"),  # Saddle brown
    "ceres": (CERES_RADIUS, "#D2B48C"),  # Tan for carbonaceous
}

dwarf_actors = {}
# Special Pluto with texture
pluto_mesh = examples.planets.load_pluto()
pluto_texture = examples.planets.download_pluto_surface(texture=True)
pluto_actor = pl.add_mesh(
    pluto_mesh,
    texture=pluto_texture,
    smooth_shading=True,
    ambient=0.1,
    diffuse=0.2,
    specular=0.0,
    name="pluto"
)
pluto_actor.RotateX(AXIAL_TILTS["pluto"])
pluto_actor.SetPosition(PLUTO_ORBIT_RADIUS, 0.0, 0.0)  # Initial position
dwarf_actors["pluto"] = pluto_actor

# Others procedural
for name, (radius, color) in {k: v for k, v in dwarf_planets.items() if k != "pluto"}.items():
    mesh = pv.Sphere(radius=radius, theta_resolution=32, phi_resolution=32)
    actor = pl.add_mesh(
        mesh,
        color=color,
        smooth_shading=True,
        ambient=0.1,
        diffuse=0.2,
        specular=0.0,
        name=name
    )
    actor.RotateX(AXIAL_TILTS[name])
    actor.SetPosition(globals()[f"{name.upper()}_ORBIT_RADIUS"], 0.0, 0.0)  # Initial position
    dwarf_actors[name] = actor

# Ceres starts in asteroid belt; adjust if needed
ceres_actor = dwarf_actors["ceres"]
ceres_actor.SetPosition(CERES_ORBIT_RADIUS, 0.0, 0.0)




# Moon Actors
moon_key_map = {}
def create_point_actor_safe(plotter, planet_center, orbit_r, name, point_size=POINT_DISPLAY_SIZE):
    initial_pos = np.array([planet_center[0] + orbit_r, planet_center[1], planet_center[2]])
    pd = pv.PolyData([initial_pos])
    actor = plotter.add_points(pd, render_points_as_spheres=True, point_size=point_size, name=name, pickable=False, color="#ffffff")
    actor.name = name
    return actor

# -------------------- MOON ACTOR CREATION --------------------
planet_moons = {
    "earth": generate_moons("earth"),
    "mars": generate_moons("mars"),
    "jupiter": generate_moons("jupiter"),
    "saturn": generate_moons("saturn"),
    "uranus": generate_moons("uranus"),
    "neptune": generate_moons("neptune"),
    "pluto": generate_moons("pluto"),
    "eris": generate_moons("eris"),
    "haumea": generate_moons("haumea"),
    "makemake": generate_moons("makemake"),
}

# Build moon_actors dict 
moon_actors = {}
for planet in ["earth", "mars", "jupiter", "saturn", "uranus", "neptune", "pluto", "eris", "haumea", "makemake", "ceres"]:
    moon_data_list = generate_moons(planet)
    moon_actors[planet] = []
    planet_name = planet  # For consistent use in conditionals/naming

    for idx, moon_data in enumerate(moon_data_list):  # Use idx for key_map
        radius_km = moon_data["radius_km"]
        radius_scene = radius_km * SIZEKM_TO_SCENE * GLOBAL_SCALE_MULTIPLIER

        # LOD: Tiny moons as points (per your thresholds, e.g., POINT_LOD_THRESHOLD)
        if radius_scene < POINT_LOD_THRESHOLD:
            # Fast point proxy for tiny moons
            point_actor = pl.add_mesh(
                pv.Sphere(radius=radius_scene * 0.1, center=[0, 0, 0]),
                color='lightgray',
                point_size=POINT_DISPLAY_SIZE,
                name=f"{planet_name}_{moon_data['name']}_point",
                render_points_as_spheres=True,
                pickable=False
            )
            actor = point_actor
            smooth_shading = False
            theta_res = phi_res = 8  # Unused for points
        else:
            # Resolution / shading logic (your original)
            if radius_scene >= 0.1:  # Large/medium moons
                theta_res = 32 if radius_km > 100 else 24
                phi_res = 32 if radius_km > 100 else 24
                smooth_shading = True
            else:  # Small moons
                theta_res = 16 if radius_km > 10 else 8
                phi_res = 16 if radius_km > 10 else 8
                smooth_shading = False

            full_name = f"{planet_name}_{moon_data['name'].lower().replace(' ', '_')}"
            short_name = moon_data["name"]

            # SPECIAL CASE: EARTH'S MOON (real mesh + texture)
            if planet_name == "earth" and short_name == "Moon":
                try:
                    # Real Moon mesh (scaled to scene)
                    luna_mesh = examples.planets.load_moon(radius=radius_scene)
                    # Real texture (cached by PyVista)
                    luna_tex = examples.planets.download_moon_surface(texture=True)
                    # Exact same call as your planets (lighting, etc.)
                    luna_actor = pl.add_mesh(
                        luna_mesh,
                        texture=luna_tex,
                        smooth_shading=True,
                        ambient=0.15,
                        diffuse=0.03,
                        specular=0.0,
                        name=full_name,
                        lighting=True
                    )
                    # Optional invisible fallback (for focus/zoom code)
                    fallback = pv.Sphere(radius=radius_scene * 0.05, theta_resolution=8, phi_resolution=8)
                    pl.add_mesh(
                        fallback,
                        color="#ffffff",
                        opacity=0.0,
                        name=f"{full_name}_debug",
                        pickable=False
                    )
                    actor = luna_actor  # Use textured actor
                except Exception as e:
                    print(f"[WARN] Failed to load Moon: {e}. Using fallback sphere.")
                    # Fallback to procedural sphere
                    actor = pl.add_mesh(
                        pv.Sphere(
                            radius=radius_scene,
                            theta_resolution=theta_res,
                            phi_resolution=phi_res,
                            center=[0, 0, 0]
                        ),
                        color='lightblue',
                        smooth_shading=smooth_shading,
                        name=full_name,
                        lighting=True
                    )
            else:
                sphere = pv.Sphere(
                    radius=radius_scene,
                    theta_resolution=theta_res,
                    phi_resolution=phi_res,
                    center=[0, 0, 0]
                )
                base_hex = moon_colors.get(short_name, "#C0C0C0")  # Default gray
                # Call your function (replace with exact if different)
                try:
                    actor = procedural_moon_actor(sphere, base_hex, short_name)
                except NameError:
                    # Mock: Simple colored sphere if function missing
                    print(f"[WARN] procedural_moon_actor not found for {short_name}; using basic sphere.")
                    actor = pl.add_mesh(
                        sphere,
                        color=base_hex,
                        smooth_shading=smooth_shading,
                        ambient=0.3,
                        diffuse=0.7,
                        name=full_name,
                        lighting=True
                    )

        actor.SetPosition(0, 0, 0)
        key = f"{planet_name}_{idx}"
        moon_key_map[key] = actor
        # Reliable name for zoom-limiter/focus
        actor.name = full_name
        # Store (data, actor) tuple exactly like original
        moon_actors[planet_name].append((moon_data, actor))

        # DEBUG: Print once for Earth's Moon
        if planet_name == "earth" and short_name == "Moon":
            print(f"[INFO] Real Moon added → {full_name} (radius {radius_scene:.4f} scene units)")

print(f"Created {sum(len(moons) for moons in moon_actors.values())} moon actors across {len(moon_actors)} planets.")
def get_orbit_points(elements, mu=GM_SUN, num_points=200):
    a_km = elements["a_km"]
    e = elements["e"]
    i_deg = elements["i_deg"]
    Omega_deg = elements["Omega_deg"]
    omega_deg = elements["omega_deg"]
    points = []
    Ms = np.linspace(0, 2 * math.pi, num_points)
    for M in Ms:
        E = solve_kepler(M, e)
        sin_v = math.sqrt(1.0 - e**2) * math.sin(E) / (1.0 - e * math.cos(E))
        cos_v = (math.cos(E) - e) / (1.0 - e * math.cos(E))
        v = math.atan2(sin_v, cos_v)
        r = a_km * (1.0 - e * math.cos(E))
        x_orb = r * math.cos(v)
        y_orb = r * math.sin(v)
        z_orb = 0.0
        i = math.radians(i_deg)
        Omega = math.radians(Omega_deg)
        omega = math.radians(omega_deg)
        x1 = x_orb * math.cos(omega) - y_orb * math.sin(omega)
        y1 = x_orb * math.sin(omega) + y_orb * math.cos(omega)
        z1 = 0.0
        x2 = x1
        y2 = y1 * math.cos(i) - z1 * math.sin(i)
        z2 = y1 * math.sin(i) + z1 * math.cos(i)
        x3 = x2 * math.cos(Omega) - y2 * math.sin(Omega)
        y3 = x2 * math.sin(Omega) + y2 * math.cos(Omega)
        z3 = z2
        points.append([x3, y3, z3])
    return np.array(points) * KM_TO_SCENE




#  REALISTIC ATMOSPHERE HALO FOR EVERY MOON 
# --------------------------------------------------------------
# Assume moon_colors dict is in scope from procedural_moon_actor (e.g., {"io": "#FFFF00", ...})
# If not, define a basic one here or fallback to hex-to-rgb conversion
def hex_to_rgb(hex_str):
    """Convert hex color to RGB array (0-1)."""
    hex_str = hex_str.lstrip('#')
    return np.array([int(hex_str[i:i+2], 16)/255.0 for i in (0, 2, 4)])


def get_moon_base_color(name, moon_colors=None):
    """Get base RGB from moon_colors dict or default palette."""
    name_lower = name.lower()
    if moon_colors and name_lower in moon_colors:
        return hex_to_rgb(moon_colors[name_lower])
    
    # Extended palette (fallback/defaults)
    palette = {
        "titan": [1.0, 0.6, 0.2],      # orange haze
        "io": [1.0, 0.9, 0.4],          # sulfur yellow
        "europa": [0.8, 0.9, 1.0],      # icy blue-white
        "enceladus": [0.9, 0.95, 1.0],  # water-vapor white
        "triton": [0.2, 0.4, 0.8],      # blue nitrogen ice
        "ganymede": [0.6, 0.6, 0.6],    # gray cratered
        "callisto": [0.4, 0.4, 0.4],    # dark gray icy
        "rhea": [0.8, 0.8, 0.8],        # bright icy
        "dione": [0.7, 0.75, 0.8],      # pale icy
        "hyperion": [0.5, 0.5, 0.5],    # spongy gray
        "iapetus": [0.3, 0.3, 0.3],     # dark leading side
        "moon": [0.9, 0.9, 0.8],        # lunar gray
        "phobos": [0.4, 0.3, 0.2],      # rusty regolith
        "deimos": [0.5, 0.4, 0.3],      # dusty brown
        "charon": [0.6, 0.7, 0.8],      # reddish ice
        "proteus": [0.5, 0.5, 0.5],     # irregular gray
        "nereid": [0.7, 0.6, 0.5],      # reddish
    }
    return np.array(palette.get(name_lower, [0.9, 0.95, 1.0]))  # default faint white-blue for tiny/unknown

def add_atmosphere_halo(moon_actor, moon_data, parent_radius_km=None, moon_colors=None):
    """
    Adds a glowing atmospheric shell around a moon actor.
    - moon_actor : the PyVista actor of the moon
    - moon_data  : dict from MOON_SCALE_FACTORS (must contain radius_km)
    - parent_radius_km : optional – for very tiny moons we cap the halo size
    - moon_colors : optional dict for procedural colors (integrates with existing)
    """
    global pl  # Assume global plotter in scope

    # 1. Base radius (scene units)
    radius_km = moon_data["radius_km"]
    radius_scene = radius_km * SIZEKM_TO_SCENE * GLOBAL_SCALE_MULTIPLIER
    is_tiny = radius_km < 10.0  # Flag for special handling

    # 2. Halo thickness: 8% of radius, but boost for tiny moons; clamped 0.02–3.0 scene units
    base_thickness = radius_scene * 0.08
    if is_tiny:
        base_thickness = max(base_thickness * 2.0, radius_scene * 0.2)  # Larger relative for visibility
    thickness = max(min(base_thickness, 3.0), 0.02)

    # 3. Create outer sphere (low-res for perf; even lower for tiny)
    theta_res = 16 if is_tiny else 24
    phi_res = 16 if is_tiny else 24
    halo = pv.Sphere(
        radius=radius_scene + thickness,
        theta_resolution=theta_res,
        phi_resolution=phi_res,
        center=[0, 0, 0]
    )

    # 4. Limb-brightening gradient (brighter at edge/limb)
    normals = halo.point_normals
    # Fake view dir toward "sun" (static; could be dynamic in update if needed)
    view_dir = np.array([1.0, 0.3, 0.2])
    view_dir /= np.linalg.norm(view_dir)
    limb_factor = np.abs(np.dot(normals, view_dir))
    limb_factor = np.clip(limb_factor, 0.0, 1.0)  # 0 at back, 1 at front (but abs for rim glow)

    # 5. Base color from palette or procedural (integrates with moon_colors)
    name = moon_data["name"]
    base_rgb = get_moon_base_color(name, moon_colors)

    # 6. Apply gradient + micro-noise for realism
    rng = np.random.default_rng(abs(hash(name)) % (2**32))
    noise = rng.normal(0.0, 0.05 if is_tiny else 0.07, size=limb_factor.shape)
    intensity = 0.3 + 0.7 * limb_factor + noise  # Slightly lower base for subtlety
    if is_tiny:
        intensity *= 1.5  # Brighter for tiny moons to "feel" atmospheric
    intensity = np.clip(intensity, 0.0, 1.0)

    # 7. Alpha: Stronger at limb, fade to near-transparent inward; boost for tiny
    alpha = (limb_factor * 0.5 + 0.05)  # 0.05-0.55 base
    if is_tiny:
        alpha *= 1.2  # More opaque for visibility
    alpha = np.clip(alpha, 0.05, 0.8)

    # 8. Combine to RGBA (per-vertex colors with alpha for proper blending)
    rgba = np.column_stack([base_rgb * intensity[:, None], alpha])
    halo.point_data["colors"] = (rgba * 255).astype(np.uint8)

    # 9. Add to plotter: Use rgb=True for RGBA, opacity=1.0 to use alpha channel, no lighting for glow
    halo_actor = pl.add_mesh(
        halo,
        scalars="colors",
        rgb=True,  # Enables RGBA interpretation (alpha from scalars)
        opacity=1.0,  # Full opacity; alpha handled by scalars
        lighting=False,  # No lighting for uniform glow
        name=f"{moon_actor.name}_halo",
        pickable=False,
        ambient=0.8,
        diffuse=0.2,
        specular=0.0
    )

    # 10. Enable scalar visibility for colors/alpha (VTK mapper tweaks for blending)
    mapper = halo_actor.GetMapper()
    mapper.ScalarVisibilityOn()
    mapper.SetScalarModeToUsePointFieldData()
    mapper.SelectColorArray("colors")

    # 11. Property tweaks for translucent/glow effect (VTK-level)
    prop = halo_actor.GetProperty()
    prop.SetOpacity(1.0)  # Use alpha from scalars
    prop.SetInterpolationToPhong()  # Smooth shading
    # For additive-like glow in translucent mode (fallback if no shader)
    prop.ShadingOn()

    # 12. Parent-child sync (for animation loop)
    halo_actor.SetPosition(moon_actor.GetPosition())
    halo_actor.SetOrientation(moon_actor.GetOrientation())


    #print(f"[INFO] Added halo to {name} (r={radius_km:.1f}km, tiny={is_tiny}, thick={thickness:.3f})")
    return halo_actor


# ------------------------------------------------------------------
#  APPLY TO EVERY MOON 
# ------------------------------------------------------------------
halo_actors = {}
halo_parents = {}

for planet, moons in moon_actors.items():
    for moon_data, moon_actor in moons:
        name = moon_data["name"]
        full_key = f"{planet}_{name.lower().replace(' ', '_')}"
        # Skip if no radius (edge case)
        if "radius_km" not in moon_data or moon_data["radius_km"] <= 0:
            continue
        halo = add_atmosphere_halo(moon_actor, moon_data, moon_colors=moon_colors)  # Pass procedural colors
        if halo is not None:
            halo_actors[full_key] = halo
            halo_parents[full_key] = moon_actor  # NEW: Map full_key to parent actor





def _high_res_sphere(radius, theta_res=360, phi_res=360):
    return pv.Sphere(radius=radius,
                     theta_resolution=theta_res,
                     phi_resolution=phi_res)

# --------------------------------------------------------------
# Helper: layered noise (low-freq warp + mid-freq craters + high-freq dust)
# --------------------------------------------------------------
def _layered_noise(pts, dirs, rng, warp, crater, noise, radius):
    # 1. Low-frequency potato warp
    A, B, C = rng.uniform(-warp, warp, 3)
    distortion = A*dirs[:,0] + B*dirs[:,1] + C*dirs[:,2]
    pts = pts + dirs * (distortion * radius * 0.35)[:, None]

    # 2. Mid-frequency crater field (many small craters)
    n_craters = rng.integers(8, 15)
    for _ in range(n_craters):
        center = rng.normal(0, 1, 3)
        center /= np.linalg.norm(center)
        dot = np.dot(dirs, center)
        mask = dot > rng.uniform(0.85, 0.94)
        depth = (dot[mask] - 0.85) / 0.09 * crater * radius
        pts[mask] -= dirs[mask] * depth[:, None]

    # 3. High-frequency surface dust/roughness
    noise_vec = rng.normal(0, noise, pts.shape[0])
    pts += dirs * (noise_vec * radius * 0.04)[:, None]

    return pts

# --------------------------------------------------------------
# MAIN FUNCTION
# --------------------------------------------------------------
# Final updated create_comet_actors() function
def km_to_scene(v_km):
    """Convert km or km/s vector to scene units. Safe to call anywhere."""
    return np.asarray(v_km) * KM_TO_SCENE
def create_comet_actors(pl, elements):
    """
    Create all actors for comets: nucleus (high-res noisy sphere), coma (persistent shell-biased particle cloud),
    dust tail (initially empty), ion tail (initially empty).
    
    FIXED: Stores relative_coma_pts as 9th tuple item for stable updates.
    """
    actors = {}
    # ------------------------------------
    # SAFE OPACITY UPDATER FOR VOLUME COMA (KEPT FOR FUTURE USE)
    # ------------------------------------
    class OpacityProxy:
        def __init__(self, actor):
            self.actor = actor
        def set(self, new_op):
            try:
                if np.ndim(new_op) == 0:
                    arr = np.linspace(0, float(new_op), 256)
                else:
                    temp = np.asarray(new_op).ravel()
                    x_old = np.linspace(0, 1, temp.size)
                    x_new = np.linspace(0, 1, 256)
                    arr = np.interp(x_new, x_old, temp)
                lut = getattr(self.actor.mapper, "lookup_table", None)
                if lut is not None:
                    lut.apply_opacity(arr)
            except Exception:
                return
    # ------------------------------------
    # COMET VISUAL PROFILES (UNCHANGED)
    # ------------------------------------
    profiles = {
        "halley": {"albedo": "#4a4a4a", "rough": 0.92, "metal": 0.0,
                   "coma_op": 0.06, "tail_dust": "#d4a574", "tail_ion": "#87CEEB"},
        "halebopp": {"albedo": "#5f5a54", "rough": 0.88, "metal": 0.0,
                     "coma_op": 0.08, "tail_dust": "#e6c9a8", "tail_ion": "#87CEEB"},
        "enscke": {"albedo": "#3c3c3c", "rough": 0.94, "metal": 0.0,
                  "coma_op": 0.04, "tail_dust": "#d4a574", "tail_ion": "#87CEEB"},
        "lovejoy": {"albedo": "#6b6b6b", "rough": 0.90, "metal": 0.0,
                    "coma_op": 0.09, "tail_dust": "#e6c9a8", "tail_ion": "#87CEEB"},
    }
    # ------------------------------------
    # MAIN LOOP
    # ------------------------------------
    for name, el in elements.items():
        # Use lowercase name for profile lookup
        profile_name = name.lower()
        p = profiles.get(profile_name, profiles["halley"])
        radius = el["radius_km"] * SIZEKM_TO_SCENE
        rng = np.random.default_rng(hash(name) % (2**32))
      
        # 1) High res sphere (UNCHANGED)
        sphere = _high_res_sphere(radius)
        pts = sphere.points.copy()
        dirs = pts / np.linalg.norm(pts, axis=1)[:, None]
      
        # 2) Apply procedural noise (UNCHANGED)
        pts = _layered_noise(
            pts, dirs, rng,
            warp=0.42, crater=0.26, noise=0.22,
            radius=radius
        )
      
        # 3) Build nucleus mesh (UNCHANGED)
        mesh = pv.PolyData(pts, sphere.faces)
        mesh = mesh.triangulate()
        mesh.clean(tolerance=1e-6, inplace=True)
        mesh.compute_normals(
            cell_normals=False,
            point_normals=True,
            consistent_normals=True,
            inplace=True
        )
        mesh = mesh.smooth(n_iter=60, relaxation_factor=0.18)
      
        # 4) Roughness / subsurface data (UNCHANGED)
        rough = np.ones(mesh.n_points) * p["rough"]
        crater_mask = rng.normal(0, 0.2, mesh.n_points)
        rough[crater_mask > 0.6] = 0.65
        subsurface = np.full(mesh.n_points, 0.12)
        subsurface[crater_mask < -0.4] = 0.04
        mesh["roughness"] = rough
        mesh["subsurface"] = subsurface
      
        # 5) Random rotation (UNCHANGED)
        r = rng.uniform(0, 360, 3)
        mesh.rotate_x(r[0], inplace=True)
        mesh.rotate_y(r[1], inplace=True)
        mesh.rotate_z(r[2], inplace=True)
      
        # 6) Compute initial comet position (UNCHANGED)
        result = kepler_to_state(
            el["a_km"], el["e"], el["i_deg"], el["Omega_deg"],
            el["omega_deg"], el["M0_deg"], el["epoch_jd"], el["epoch_jd"],
            mu=GM_SUN
        )
      
        mx_km, my_km, mz_km, vel_km_s, v_deg = result
      
        pos_km = np.array([mx_km, my_km, mz_km])
        if not np.all(np.isfinite(pos_km)):
            print(f"[WARN] Invalid initial pos for comet {name}, fallback to safe position.")
            pos_km = np.array([0, 0, NEPTUNE_AU * AU_KM])
        pos_scene = km_to_scene(pos_km)
        if not np.all(np.isfinite(vel_km_s)):
            vel_km_s = np.zeros(3)
        vel_scene = vel_km_s * KM_TO_SCENE
      
        # 7) Add Nucleus with proper position (UNCHANGED, full opacity)
        nucleus_actor = pl.add_mesh(
            mesh,
            color=p["albedo"],
            pbr=True,
            metallic=p["metal"],
            roughness=p["rough"],
            opacity=1.0,
            smooth_shading=True,
            specular=0.03,
            name=f"nucleus_{name}"
        )
        nucleus_actor.SetPosition(*pos_scene)
      
    
        # FIXED: 9-item tuple with relative_pts last
        actors[name] = (
            nucleus_actor,

            el["Tp_jd"],
        )   
    # Cleanup (UNCHANGED)
    if hasattr(pl, "scalar_bars") and pl.scalar_bars:
        try:
            for key in list(pl.scalar_bars.keys()):
                pl.remove_scalar_bar(title=key)
        except Exception:
            pass
        pl.scalar_bars.clear()
    return actors





# --- Draw orbital paths for all planets around the Sun ---
##def add_orbit_paths(plotter):
##    """Draw thin white elliptical orbit lines for each planet."""
##    planet_names = ["mercury", "venus", "earth", "mars", "jupiter", "saturn", "uranus", "neptune"]
##    for planet_name in planet_names:
##        elements = orbital_elements[planet_name]
##        orbit_points = get_orbit_points(elements)
##        orbit_line = pv.Spline(orbit_points, 200)
##        plotter.add_mesh(
##            orbit_line,
##            color="white",
##            line_width=1.0,
##            opacity=0.4,
##            name=f"{planet_name}_orbit",
##            render_lines_as_tubes=False,
##            pickable=False
##        )
##
##add_orbit_paths(pl)

# ----------------------------------------------------------------------
#  FULL UPDATED SolarSystemAnimator CLASS
# ----------------------------------------------------------------------
class SolarSystemAnimator:
    def __init__(
        self,
        jactor,
        sactor,
        nactor,
        uactor,
        mactor,
        eactor,
        vactor,
        mercury_actor,  # Moved here for logical order
        moon_actors,
        plotter,
        rings,
        uranus_rings,
        asteroids,
        sun_actor,
        plutoactor=None,  # Optional kwargs-style for dwarfs
        erisactor=None,
        haumeaactor=None,
        makemakeactor=None,
        ceresactor=None,
        comet_actors=None,
    ):
        # ----- store references ------------------------------------------------
        self.jupiter_actor = jactor
        self.saturn_actor = sactor
        self.neptune_actor = nactor
        self.uranus_actor = uactor
        self.mars_actor = mactor
        self.earth_actor = eactor
        self.venus_actor = vactor
        self.mercury_actor = mercury_actor
        self.moon_actors = moon_actors
        self.plotter = plotter
        self.ring_actors = rings
        self.uranus_ring_actors = uranus_rings
        self.asteroid_actor = asteroids
        self.sun_actor = sun_actor
        self.comet_actors = comet_actors or {}
        self.prev_comet_pos = {}
        self.total_sim_time = 0.0


        
       
        # Dwarf actors (handle None for optional)
        self.pluto_actor = plutoactor if plutoactor else None
        self.eris_actor = erisactor if erisactor else None
        self.haumea_actor = haumeaactor if haumeaactor else None
        self.makemake_actor = makemakeactor if makemakeactor else None
        self.ceres_actor = ceresactor if ceresactor else None
       
        init_jd = get_current_jd()
        gmst_deg = gmst_from_jd(init_jd)
        self.earth_actor.RotateZ(gmst_deg + 180)
          
        # ----- animation state -------------------------------------------------
        self.last_time = time.time()
        self.current_focus_planet = None
        self.moon_number_buffer = ""
        self.current_focus_actor = None
        self.current_min_dist = MIN_ZOOM_DISTANCE
        self.plotter._in_update = False  # Flag for double-render prevention
        self.accumulator = 0.0
        self.comets_visible = True

       
        # Add a corner annotation for the simulated time display (upper left)
        self.time_text = self.plotter.add_text(
            "Date: Initializing...",
            position="upper_left",
            font_size=12,
            color="white",
            shadow=True,
            name="sim_time_display"
        )
       
        self.saturn_actor.RotateX(AXIAL_TILTS["saturn"])
        for ractor in self.ring_actors:
            ractor.RotateX(AXIAL_TILTS["saturn"])
        self.uranus_actor.RotateX(AXIAL_TILTS["uranus"])
        for ractor in self.uranus_ring_actors:
            ractor.RotateX(AXIAL_TILTS["uranus"])
       
        for name in ["sun", "mercury", "venus", "earth", "mars", "jupiter", "saturn", "uranus", "neptune"]:
            if hasattr(self, f"{name}_actor"):
                actor = getattr(self, f"{name}_actor")
                actor.RotateX(AXIAL_TILTS[name])
       
        # Dwarfs (if present)
        for name in ["pluto", "eris", "haumea", "makemake", "ceres"]:
            if hasattr(self, f"{name}_actor") and getattr(self, f"{name}_actor") is not None:
                actor = getattr(self, f"{name}_actor")
                actor.RotateX(AXIAL_TILTS[name])
               
        # ----- key bindings ----------------------------------------------------
        self.plotter.add_key_event('j', lambda: self.set_planet_focus('jupiter', self.jupiter_actor))
        self.plotter.add_key_event('s', lambda: self.set_planet_focus('saturn', self.saturn_actor))
        self.plotter.add_key_event('n', lambda: self.set_planet_focus('neptune',self.neptune_actor))
        self.plotter.add_key_event('u', lambda: self.set_planet_focus('uranus', self.uranus_actor))
        self.plotter.add_key_event('m', lambda: self.set_planet_focus('mars', self.mars_actor))
        self.plotter.add_key_event('g', lambda: self.set_planet_focus('earth', self.earth_actor))
        self.plotter.add_key_event('v', lambda: self.set_planet_focus('venus', self.venus_actor))
        self.plotter.add_key_event('y', lambda: self.set_planet_focus('mercury',self.mercury_actor))
        self.plotter.add_key_event('a', self.focus_asteroid_belt)
        self.plotter.add_key_event('t', self.focus_sun)


        self.plotter.add_key_event('z', lambda: self.focus_on_comet('halley'))
        self.plotter.add_key_event('x', lambda: self.focus_on_comet('halebopp'))
        self.plotter.add_key_event('l', lambda: self.focus_on_comet('enscke'))
        self.plotter.add_key_event('d', lambda: self.focus_on_comet('lovejoy'))

        
        # Dwarf keys (guarded to prevent crashes if actor is None)
        if self.pluto_actor is not None:
            self.plotter.add_key_event('o', lambda: self.set_planet_focus('pluto', self.pluto_actor))
        if self.eris_actor is not None:
            self.plotter.add_key_event('i', lambda: self.set_planet_focus('eris', self.eris_actor))  # Avoid conflict with Earth 'g'
        if self.haumea_actor is not None:
            self.plotter.add_key_event('h', lambda: self.set_planet_focus('haumea', self.haumea_actor))
        if self.makemake_actor is not None:
            self.plotter.add_key_event('k', lambda: self.set_planet_focus('makemake', self.makemake_actor))
        if self.ceres_actor is not None:
            self.plotter.add_key_event('c', lambda: self.set_planet_focus('ceres', self.ceres_actor))
        
        # Numeric keys for moons (0-9)
        for num in range(10):
            self.plotter.add_key_event(str(num), lambda n=num: self.handle_numeric_key(n))
        if hasattr(self, "mars_actor"):
            mars_gmst = gmst_deg * (ROTATION_PERIOD["earth"] / ROTATION_PERIOD["mars"]) % 360.0
            self.mars_actor.RotateX(AXIAL_TILTS["mars"])
            self.mars_actor.RotateZ(mars_gmst)
        if hasattr(self, "venus_actor"):
            venus_gmst = gmst_deg * (ROTATION_PERIOD["earth"] / abs(ROTATION_PERIOD["venus"])) % 360.0
            self.venus_actor.RotateX(AXIAL_TILTS["venus"])
            self.venus_actor.RotateZ(-venus_gmst)
        if hasattr(self, "mercury_actor"):
            merc_gmst = gmst_deg * (ROTATION_PERIOD["earth"] / abs(ROTATION_PERIOD["mercury"])) % 360.0
            self.mercury_actor.RotateX(AXIAL_TILTS["mercury"])
            self.mercury_actor.RotateZ(-merc_gmst)
           
        # ----- Consolidated Moon Initial Setup (single loop, no duplicates) -----
        for planet_name, moons in self.moon_actors.items():
            mu = GM_PLANET.get(planet_name, GM_EARTH)
            for i, (moon_dict, actor) in enumerate(moons):
                required = {"a_km", "e", "i_deg", "Omega_deg", "omega_deg", "M0_deg", "epoch_jd"}
                if not required.issubset(moon_dict.keys()):
                    continue
                # True anomaly at start-up
                _, _, _, _, v_deg = kepler_to_state(
                    moon_dict["a_km"], moon_dict["e"], moon_dict["i_deg"],
                    moon_dict["Omega_deg"], moon_dict["omega_deg"],
                    moon_dict["M0_deg"], moon_dict["epoch_jd"], init_jd, mu=mu
                )
                # Axial tilt (most moons ≈0°)
                moon_name_lower = moon_dict["name"].lower()
                tilt = AXIAL_TILTS.get(moon_name_lower, 0.0)
                prograde_moons = ["moon", "io", "europa", "ganymede", "callisto", "phobos", "deimos"]
                sign = -1 if moon_name_lower in prograde_moons else 1
                actor.RotateX(sign * tilt)
                # Lock near side to planet
                actor.RotateZ(v_deg)
                moon_dict["tilt_applied"] = True  # Prevent double in update()

        # ----- mouse-wheel zoom (clamped) --------------------------------------
        if hasattr(self.plotter, 'iren') and self.plotter.iren is not None:
            self.plotter.iren.add_observer("MouseWheelForwardEvent", self.limit_zoom)
            self.plotter.iren.add_observer("MouseWheelBackwardEvent", self.limit_zoom)

    # ----------------------------------------------------------------------
    # ZOOM LIMITER – uses auto-clipping, no double render
    # ----------------------------------------------------------------------
    def limit_zoom(self, caller=None, event=None):
        try:
            cam = self.plotter.camera
            pos = np.array(cam.GetPosition())
            focal = np.array(cam.GetFocalPoint())
            direction = focal - pos
            dist = np.linalg.norm(direction)
            if dist < 1e-9:
                return
            direction /= dist
            # forward = zoom-in, backward = zoom-out
            zoom_factor = 0.9 if event and 'Forward' in str(event) else 1.1
            new_dist = dist * zoom_factor
            new_dist = max(new_dist, self.current_min_dist)
            new_pos = focal - direction * new_dist
            cam.SetPosition(*new_pos)
            # AUTO CLIPPING – adjusts based on visible actors, fixes disappearing objects
            self.plotter.reset_camera_clipping_range()
            # Render only if NOT in update loop (prevents flicker)
            if not self.plotter._in_update:
                self.plotter.render()
        except Exception as e:
            print(f"[WARN] Zoom limiter error: {e}")
    # ----------------------------------------------------------------------
    # PLANET FOCUS – with auto clipping
    # ----------------------------------------------------------------------
    def set_planet_focus(self, planet_name, actor):
        self.current_focus_planet = planet_name
        self.moon_number_buffer = ""
        self.current_focus_actor = actor
        self.current_min_dist = get_object_min_dist(actor.name)
        self.update_camera_focus(actor)
        num_moons = len(self.moon_actors.get(planet_name, []))
        if num_moons > 0:
            print(f"Camera focus set to {planet_name.upper()}. Press 1–{num_moons} for moons. Min zoom: {self.current_min_dist:.2f} units")
        else:
            print(f"Camera focus set to {planet_name.upper()}. (No moons) Min zoom: {self.current_min_dist:.2f} units")
    # ----------------------------------------------------------------------
    # HANDLE MOON KEY – with auto clipping
    # ----------------------------------------------------------------------
    def handle_numeric_key(self, num):
        if self.current_focus_planet is None:
            return
        self.moon_number_buffer += str(num)
        try:
            moon_index = int(self.moon_number_buffer) - 1 # 0-based
        except ValueError:
            self.moon_number_buffer = ""
            return
        planet = self.current_focus_planet
        moons = self.moon_actors.get(planet, [])
        if moon_index < 0 or moon_index >= len(moons):
            print(f"Invalid moon index for {planet} (1-{len(moons)})")
            self.moon_number_buffer = ""
            return
        moon_data, moon_actor = moons[moon_index]
        moon_name = moon_data["name"]
        self.current_focus_actor = moon_actor
        self.current_min_dist = get_object_min_dist(moon_actor.name, moon_data.get("radius_km", 0) * SIZEKM_TO_SCENE)
        self.update_camera_focus(moon_actor)
        print(f"→ Focused on {planet.capitalize()}'s {moon_name}. Min zoom: {self.current_min_dist:.2f} units")
        self.moon_number_buffer = ""
    # ----------------------------------------------------------------------
    # UPDATE CAMERA FOCUS – with auto clipping
    # ----------------------------------------------------------------------
    def update_camera_focus(self, actor_or_pos):
        if isinstance(actor_or_pos, (list, tuple, np.ndarray)):
            actor_pos = np.array(actor_or_pos, dtype=float)
            actor_name = ""
        else:
            actor_pos = np.array(actor_or_pos.GetPosition(), dtype=float) if actor_or_pos else np.zeros(3)
            actor_name = getattr(actor_or_pos, "name", "").lower()
        cam = self.plotter.camera
        old_focal = np.array(cam.GetFocalPoint(), dtype=float)
        old_pos = np.array(cam.GetPosition(), dtype=float)
        view_vec = old_pos - old_focal
        vnorm = np.linalg.norm(view_vec)
        view_dir = view_vec / vnorm if vnorm > 1e-9 else np.array([0.0, 0.0, 1.0])
        base_radius = {
            "mercury": MERCURY_RADIUS,
            "venus": VENUS_RADIUS,
            "earth": EARTH_RADIUS,
            "mars": MARS_RADIUS,
            "jupiter": JUPITER_RADIUS,
            "saturn": SATURN_RADIUS,
            "uranus": URANUS_RADIUS,
            "neptune": NEPTUNE_RADIUS,
            "pluto": PLUTO_RADIUS,
            "eris": ERIS_RADIUS,
            "haumea": HAUMEA_RADIUS,
            "makemake": MAKEMAKE_RADIUS,
            "ceres": CERES_RADIUS,
            "sun": SUN_RADIUS,
        }.get(actor_name.split("_")[0] if "_" in actor_name else actor_name, JUPITER_RADIUS * 0.01)
        orbit_scale = {
            "mercury": MERCURY_ORBIT_RADIUS,
            "venus": VENUS_ORBIT_RADIUS,
            "earth": EARTH_ORBIT_RADIUS,
            "mars": MARS_ORBIT_RADIUS,
            "jupiter": JUPITER_ORBIT_RADIUS,
            "saturn": SATURN_ORBIT_RADIUS,
            "uranus": URANUS_ORBIT_RADIUS,
            "neptune": NEPTUNE_ORBIT_RADIUS,
            "pluto": PLUTO_ORBIT_RADIUS,
            "eris": ERIS_ORBIT_RADIUS,
            "haumea": HAUMEA_ORBIT_RADIUS,
            "makemake": MAKEMAKE_ORBIT_RADIUS,
            "ceres": CERES_ORBIT_RADIUS,
        }
        desired_distance = max(orbit_scale.get(actor_name, base_radius * 20) * 0.25, self.current_min_dist)
        new_focal = actor_pos
        new_position = new_focal + view_dir * desired_distance
        cam.SetFocalPoint(*new_focal)
        cam.SetPosition(*new_position)
        cam.SetViewUp(0, 0, 1)
        # AUTO CLIPPING – adjusts to visible bounds
        self.plotter.reset_camera_clipping_range()
        self.plotter.render()
    # ----------------------------------------------------------------------
    # MOVE CAMERA (DOLLY) – with auto clipping
    # ----------------------------------------------------------------------
    def move_camera_along_view(self, distance):
        cam = self.plotter.camera
        pos = np.array(cam.GetPosition())
        focal = np.array(cam.GetFocalPoint())
        direction = focal - pos
        dist = np.linalg.norm(direction)
        if dist < 1e-9:
            return
        direction /= dist
        new_pos = pos + direction * distance
        new_dist = np.linalg.norm(new_pos - focal)
        if new_dist < self.current_min_dist:
            new_pos = focal - (focal - new_pos) / new_dist * self.current_min_dist
            new_dist = self.current_min_dist
        cam.SetPosition(*new_pos)
        # AUTO CLIPPING
        self.plotter.reset_camera_clipping_range()
    # ----------------------------------------------------------------------
    # ASTEROID BELT FOCUS – with auto clipping
    # ----------------------------------------------------------------------
    def focus_asteroid_belt(self):
        cam = self.plotter.camera
        belt_center = np.array([AU_SCALE * 2.7, 0.0, 0.0])
        cam.SetFocalPoint(*belt_center)
        cam.SetPosition(belt_center[0] + 0.015 * AU_SCALE, belt_center[1] + 0.005 * AU_SCALE, belt_center[2] + 0.005 * AU_SCALE)
        cam.SetViewUp(0, 0, 1)
        # AUTO CLIPPING
        self.plotter.reset_camera_clipping_range()
        self.asteroid_actor.GetProperty().SetOpacity(1.0)
        self.current_focus_planet = None
        self.current_focus_actor = self.asteroid_actor
        self.current_min_dist = AU_SCALE * 0.5
        self.plotter.render()
    # ----------------------------------------------------------------------
    # SUN FOCUS – with auto clipping
    # ----------------------------------------------------------------------
    def focus_sun(self):
        cam = self.plotter.camera
        cam.SetFocalPoint(0, 0, 0)
        cam.SetPosition(AU_SCALE * 5, AU_SCALE * 3, AU_SCALE * 2)
        cam.SetViewUp(0, 0, 1)
        # AUTO CLIPPING
        self.plotter.reset_camera_clipping_range()
        self.plotter.render()
        print("Camera focus reset to SUN. Min zoom: {:.2f} units".format(get_object_min_dist("sun")))
        self.current_focus_planet = None
        self.current_focus_actor = None
        self.current_min_dist = get_object_min_dist("sun")

    def focus_on_comet(self, comet_name):
        if comet_name not in self.comet_actors:
            print(f"Unknown comet: {comet_name}")
            return
        nucleus_actor = self.comet_actors[comet_name][0]
        comet_pos = np.array(nucleus_actor.GetPosition())
        dist = np.linalg.norm(comet_pos) * 2 # 2x distance for good view
        if dist < 1: dist = 5 # Min zoom for inner
        # Set camera: from offset pos, look at comet, up vector
        offset_pos = comet_pos + np.array([dist, 0, 0]) # Simple X-offset; improve with direction if needed
        self.plotter.camera_position = [(offset_pos[0], offset_pos[1], offset_pos[2]), tuple(comet_pos), (0, 0, 1)]
        self.plotter.camera.view_angle = 20 # Tighter zoom
        self.plotter.reset_camera_clipping_range()
        print(f"Focused on {comet_name} at {comet_pos}")


        # ----------------------------------------------------------------------
    # RESOLVE ACTOR POSITION (helper)
    # ----------------------------------------------------------------------
    def resolve_actor_position(self, actor_or_key):
        if hasattr(actor_or_key, "GetPosition"):
            pos = np.array(actor_or_key.GetPosition(), dtype=float)
            if np.all(np.isfinite(pos)):
                return pos
        return np.zeros(3, dtype=float)





    
    #  MAIN UPDATE LOOP
    
    def update(self, obj=None, event=None):
        if self.plotter._in_update:
            return
        self.plotter._in_update = True
        current_time = time.time()
        dt_real = min(current_time - self.last_time, 0.1)  # Cap at 100ms to handle lag spikes
        self.last_time = current_time
        dt_sim = dt_real * SIM_SECONDS_PER_REAL_SECOND
        self.total_sim_time += dt_sim
        global SIM_TIME
        SIM_TIME += dt_sim
        sim_dt = datetime.fromtimestamp(SIM_TIME, tz=timezone.utc)
        formatted_time = sim_dt.strftime("%Y-%m-%d %H:%M:%S UTC")
        self.time_text.SetText(2, f"Simulated Date: {formatted_time}")
        target_jd = datetime_to_julian_day(sim_dt)

        # ----- Update planet positions -----
        set_body_position_from_elements("mercury", self.mercury_actor, orbital_elements["mercury"], target_jd)
        set_body_position_from_elements("venus", self.venus_actor, orbital_elements["venus"], target_jd)
        set_body_position_from_elements("earth", self.earth_actor, orbital_elements["earth"], target_jd)
        set_body_position_from_elements("mars", self.mars_actor, orbital_elements["mars"], target_jd)
        set_body_position_from_elements("jupiter", self.jupiter_actor, orbital_elements["jupiter"], target_jd)
        set_body_position_from_elements("saturn", self.saturn_actor, orbital_elements["saturn"], target_jd)
        set_body_position_from_elements("uranus", self.uranus_actor, orbital_elements["uranus"], target_jd)
        set_body_position_from_elements("neptune", self.neptune_actor, orbital_elements["neptune"], target_jd)
        set_body_position_from_elements("pluto", self.pluto_actor, orbital_elements["pluto"], target_jd)
        set_body_position_from_elements("eris", self.eris_actor, orbital_elements["eris"], target_jd)
        set_body_position_from_elements("haumea", self.haumea_actor, orbital_elements["haumea"], target_jd)
        set_body_position_from_elements("makemake", self.makemake_actor, orbital_elements["makemake"], target_jd)
        set_body_position_from_elements("ceres", self.ceres_actor, orbital_elements["ceres"], target_jd)

        # ----- Sun barycentric wobble -----
        jup_mass = GM_JUPITER / GM_SUN
        sat_mass = GM_SATURN / GM_SUN
        jup_pos = np.array(self.jupiter_actor.GetPosition())
        sat_pos = np.array(self.saturn_actor.GetPosition())
        sun_wobble = (jup_mass * jup_pos + sat_mass * sat_pos) / (1 + jup_mass + sat_mass)
        self.sun_actor.SetPosition(*sun_wobble)

        # ----- Planet + sun rotations -----
        for name, actor in [
            ("sun", self.sun_actor),
            ("mercury", self.mercury_actor),
            ("venus", self.venus_actor),
            ("earth", self.earth_actor),
            ("mars", self.mars_actor),
            ("jupiter", self.jupiter_actor),
            ("saturn", self.saturn_actor),
            ("uranus", self.uranus_actor),
            ("neptune", self.neptune_actor),
            ("pluto", self.pluto_actor),
            ("eris", self.eris_actor),
            ("haumea", self.haumea_actor),
            ("makemake", self.makemake_actor),
            ("ceres", self.ceres_actor),
        ]:
            deg = ROTATION_DEG_PER_SIMSEC[name] * dt_sim
            actor.RotateZ(deg % 360.0)

        # Ring follow
        for ractor in self.ring_actors:
            ractor.SetPosition(*self.saturn_actor.GetPosition())
        for ractor in self.uranus_ring_actors:
            ractor.SetPosition(*self.uranus_actor.GetPosition())
            
        sun_pos = np.array(self.sun_actor.GetPosition())
        sun_pos_km = sun_pos / KM_TO_SCENE  # Convert sun position to km for distance calculations

        # ----- Update halo positions/orientations (planets + comets) -----
        for full_key, halo in halo_actors.items():
            if full_key in halo_parents:  # Safe check
                parent = halo_parents[full_key]
                halo.SetPosition(*parent.GetPosition())
                halo.SetOrientation(*parent.GetOrientation())
                
                # Comet-specific tweaks (if key indicates comet)
                if 'comet' in full_key.lower() or full_key.startswith('comet_'):
                    # Get comet name from key (e.g., 'comet_halley')
                    comet_name = full_key.split('_')[-1] if '_' in full_key else full_key
                    r_km = np.linalg.norm(np.array(parent.GetPosition()) / KM_TO_SCENE - sun_pos_km)
                    
                    # Only show halo if within sublimation range
                    if r_km < SUBLIMATION_DISTANCE:
                        halo.SetVisibility(True)
                        # Fade opacity based on distance (brighter closer to sun)
                        halo_opacity = min(0.3 * (AU_KM / max(r_km, 0.1 * AU_KM)), 0.8)
                        halo.GetProperty().SetOpacity(halo_opacity)
                        halo.GetProperty().SetColor(0.7, 0.9, 1.0)  # Cyan for comets
                    else:
                        halo.SetVisibility(False)
                        # Optional: Set opacity to 0 even if visible (backup)
                        halo.GetProperty().SetOpacity(0.0)


        # ----- Update comet positions, halo & volumetric body -----
        if self.comets_visible:
            for name, (nucleus_actor, _) in self.comet_actors.items():
                try:
                    # ---- 1. Keplerian position ----
                    target_jd = datetime_to_julian_day(sim_dt)
                    mx_km, my_km, mz_km, _, _ = kepler_to_state(
                        COMET_ELEMENTS[name]["a_km"], COMET_ELEMENTS[name]["e"],
                        COMET_ELEMENTS[name]["i_deg"], COMET_ELEMENTS[name]["Omega_deg"],
                        COMET_ELEMENTS[name]["omega_deg"], COMET_ELEMENTS[name]["M0_deg"],
                        COMET_ELEMENTS[name]["epoch_jd"], target_jd, mu=GM_SUN
                    )
                    pos_km = np.array([mx_km, my_km, mz_km])
                    if not np.all(np.isfinite(pos_km)):
                        continue

                    # ---- 2. Scale & cap ----
                    pos_scene = km_to_scene(pos_km)

                    # ---- 3. Nucleus ----
                    current_pos = np.array(nucleus_actor.GetPosition())
                    new_pos = pos_scene
                    smoothed_pos = 0.7 * current_pos + 0.7 * new_pos  # Smooth motion
                    nucleus_actor.SetPosition(*smoothed_pos)
                    cam_pos = np.array(self.plotter.camera.GetPosition())
                    view_dist = np.linalg.norm(cam_pos - smoothed_pos)
                    point_size = max(1.0, 6.0 * (AU_SCALE / max(view_dist, AU_SCALE * 0.01)))
                    nucleus_actor.GetProperty().SetPointSize(point_size)
                    nucleus_actor.SetVisibility(True)

                    # ---- 4. Halo (coma glow) ----
                    halo_key = f"comet_{name.lower()}"
                    if halo_key in halo_actors:
                        halo = halo_actors[halo_key]
                        halo.SetPosition(*pos_scene)
                        halo.SetOrientation(*nucleus_actor.GetOrientation())
                        sun_pos_km = np.array(self.sun_actor.GetPosition()) / KM_TO_SCENE
                        r_km = np.linalg.norm(pos_km - sun_pos_km)
                        if r_km < SUBLIMATION_DISTANCE:
                            halo.SetVisibility(True)
                            halo_opacity = min(0.3 * (AU_KM / max(r_km, 0.1 * AU_KM)), 0.8)
                            halo.GetProperty().SetOpacity(halo_opacity)
                            halo.GetProperty().SetColor(0.7, 0.9, 1.0)  # Cyan glow
                        else:
                            halo.SetVisibility(False)
                            halo.GetProperty().SetOpacity(0.0)

                        

                    # ---- 5. Volumetric Body ----

                    if halo_key not in self.prev_comet_pos:
                        self.prev_comet_pos[halo_key] = pos_km.copy()

                    # Prevent division by zero if dt_sim is ~0
                    vel_km_s = np.zeros(3)
                    if dt_sim > 1e-9:
                        vel_km_s = (pos_km - self.prev_comet_pos[halo_key]) / dt_sim
                    
                    self.prev_comet_pos[halo_key] = pos_km.copy()
                    vel_scene_s = vel_km_s * KM_TO_SCENE

                    halo_radius_km = COMET_ELEMENTS[name]["radius_km"] * 120
                    head_scale = halo_radius_km * SIZEKM_TO_SCENE / COMET_GRID_HALF

                    # Generate density
                    new_density = make_volumetric_comet_body(
                        center=np.array([0.0, 0.0, 0.0]),
                        velocity=vel_scene_s,
                        t=self.total_sim_time,
                        head_scale=head_scale
                    )

                    # Update grid
                    vol_grid = vol_grids[halo_key]
                    vol_grid.cell_data["comet"][:] = new_density.flatten(order="F")

                    # REMOVE RECTANGULAR BLOCK: Zero out low density
                    vol_grid.cell_data["comet"][vol_grid.cell_data["comet"] < 0.01] = 0.0
                    vol_grid.Modified()

                    # Position & orient
                    vol_actor = vol_actors[halo_key] # Get the actor reference
                    vol_actor.SetPosition(*smoothed_pos)
                    if np.linalg.norm(vel_scene_s) > 1e-6:
                        tail_dir = -vel_scene_s / np.linalg.norm(vel_scene_s)
                        rot = R.align_vectors([tail_dir], [[1.0, 0.0, 0.0]])[0]
                        euler_deg = rot.as_euler('xyz', degrees=True)
                        vol_actor.SetOrientation(*euler_deg)

                    vol_actor.SetVisibility(int(r_km < SUBLIMATION_DISTANCE))

                    
                except Exception as e:
                    print(f"[COMET] {name} update error: {e}")
                    import traceback
                    traceback.print_exc()
                    continue





                # ----- MOONS -----
        planet_positions = {}
        for name in ["mercury", "venus", "earth", "mars", "jupiter", "saturn", "uranus", "neptune"]:
            actor = getattr(self, f"{name}_actor")
            if actor is not None:
                planet_positions[name] = np.array(actor.GetPosition())
        if self.pluto_actor is not None:
            planet_positions["pluto"] = np.array(self.pluto_actor.GetPosition())
        for name in ["eris", "haumea", "makemake", "ceres"]:
            actor = getattr(self, f"{name}_actor", None)
            if actor is not None:
                planet_positions[name] = np.array(actor.GetPosition())
        for planet_name, moons in self.moon_actors.items():
            if planet_name not in planet_positions:
                continue
            px, py, pz = planet_positions[planet_name]
            mu = GM_PLANET.get(planet_name, GM_EARTH)
            for moon, actor in moons:
                required = ["a_km", "e", "i_deg", "Omega_deg", "omega_deg", "M0_deg", "epoch_jd"]
                if not all(k in moon for k in required):
                    continue
                epoch_jd = moon["epoch_jd"]
                delta_jd = target_jd - epoch_jd
                TINY_MOON_THRESHOLD = 150000  # km
                TINY_MOON_SCALE = 0.1
                orbit_scale = TINY_MOON_SCALE if moon.get("a_km", 0) < TINY_MOON_THRESHOLD else 1.0
                effective_target_jd = epoch_jd + delta_jd * orbit_scale
                mx_km, my_km, mz_km, _, _ = kepler_to_state(
                    moon["a_km"], moon["e"], moon["i_deg"], moon["Omega_deg"],
                    moon["omega_deg"], moon["M0_deg"], epoch_jd, effective_target_jd, mu=mu
                )
                moon_offset = np.array([mx_km, my_km, mz_km]) * KM_TO_SCENE * (1 / 0.02)
                norm = np.linalg.norm(moon_offset)
                if norm < JUPITER_TINY_THRESHOLD and norm > 0:
                    moon_offset *= (JUPITER_TINY_THRESHOLD / norm)
                current_pos = np.array(actor.GetPosition())
                new_pos = np.array([px + moon_offset[0], py + moon_offset[1], pz + moon_offset[2]])
                smoothed_pos = 0.3 * current_pos + 0.7 * new_pos
                actor.SetPosition(*smoothed_pos)

                moon_name = moon["name"].lower()
                if "tilt_applied" not in moon:
                    tilt = AXIAL_TILTS.get(moon_name, 0.0)
                    prograde_moons = ["moon", "io", "europa", "ganymede", "callisto", "phobos", "deimos"]
                    sign = -1 if moon_name in prograde_moons else 1
                    actor.RotateX(sign * tilt)
                    moon["tilt_applied"] = True

                period = ROTATION_PERIOD.get(moon_name)
                if period is None:
                    a_km = moon["a_km"]
                    period = 2 * math.pi * math.sqrt(a_km**3 / mu) if a_km > 0 and mu > 0 else 0
                if period != 0:
                    spin_scale = orbit_scale if TINY_MOON_SPIN_SCALE_MATCH else 1.0
                    effective_period = period / spin_scale
                    deg = (360.0 / effective_period) * dt_sim
                    actor.RotateZ(deg % 360.0)

        # ----- Asteroid belt -----
        if self.asteroid_actor:
            try:
                cloud = self.asteroid_actor.GetMapper().GetInput()
                if cloud.GetNumberOfPoints() > 0:
                    points = cloud.points.copy()
                    x, y, _ = points[:, 0], points[:, 1], points[:, 2]
                    r_scene = np.sqrt(x**2 + y**2)
                    r_km = r_scene / KM_TO_SCENE
                    omega = np.sqrt(GM_SUN / r_km**3)
                    dtheta = omega * dt_sim
                    c, s = np.cos(dtheta), np.sin(dtheta)
                    new_x = x * c - y * s
                    new_y = x * s + y * c
                    cloud.points[:, 0], cloud.points[:, 1] = new_x, new_y
                    cloud.Modified()
                    self.asteroid_actor.GetMapper().Modified()
            except:
                pass


        # Render
        self.plotter.reset_camera_clipping_range()
        self.plotter.render()
        self.plotter._in_update = False

# ---------------------------- START SIMULATION ----------------------------
cam_distance = NEPTUNE_ORBIT_RADIUS * 0.6
pl.camera_position = [(cam_distance, cam_distance, cam_distance / 2), (0, 0, 0), (0, 0, 1)]
pl.camera.view_angle = 30
pl.reset_camera_clipping_range()  # Use auto clipping for initial view

comet_actors = create_comet_actors(pl, COMET_ELEMENTS)

# ------------------------------------------------------------------
# ADD VOLUMETRIC COMET BODY (DENSER & SMALLER)
# ------------------------------------------------------------------
COMET_TAIL_LENGTH_FACTOR = 5.0
COMET_TAIL_WIDTH_FACTOR  = 1.0

# === SHRINK THE GRID (makes comet smaller) ===
COMET_GRID_SIZE = 128
COMET_GRID_HALF = 1.0   # ← WAS 3.0 → now 60% smaller (1.8 / 3.0 = 0.6x size)
COMET_GRID = np.linspace(-COMET_GRID_HALF, COMET_GRID_HALF, COMET_GRID_SIZE)
Xg, Yg, Zg = np.meshgrid(COMET_GRID, COMET_GRID, COMET_GRID, indexing="ij")

# === DENSITY BOOST (preserves total "particles") ===
SIZE_REDUCTION_FACTOR = 3.0 / COMET_GRID_HALF  # ~1.67x
DENSITY_SCALE = SIZE_REDUCTION_FACTOR ** 3      # Volume scales with cube → ~4.6x denser

# === REAL-WORLD COMET COLORS (cmap for volumetric coma/tail) ===
COMET_CMAPS = {
    "halley": "Greens",      # Greenish coma from C2 gas
    "halebopp": "Blues",     # Blue ion tail, white dust
    "enscke": "YlOrBr",      # Yellowish dust tail
    "lovejoy": "Greens",     # Famous green coma from diatomic carbon
    "mcnaught": "Blues",     # Bright white dust, blue ion
    # Add more if your COMET_ELEMENTS has others
}

vol_grids = {}
vol_actors = {}
print(f"[INFO] Adding DENSER volumetric comet body (size ×{1/SIZE_REDUCTION_FACTOR:.2f}, density ×{DENSITY_SCALE:.1f})...")

def make_volumetric_comet_body(center, velocity, t, head_scale=1.0):
    safe_head_scale = max(head_scale, 1e-6)
    falloff_constant = math.log(1.5 / 0.1) / COMET_GRID_HALF / safe_head_scale

    vel_norm = np.linalg.norm(velocity)
    tail_dir = -velocity / vel_norm if vel_norm > 0 else np.array([1.0, 0.0, 0.0])

    # Slight forward shift to center density on nucleus
    center = center - tail_dir * 0.025

    # ---------- HEAD (coma) ----------
    rhead = np.sqrt((Xg-center[0])**2 + (Yg-center[1])**2 + (Zg-center[2])**2)
    head = np.exp(-rhead * falloff_constant)
    mask_radius = COMET_GRID_HALF * 0.8
    head_mask = rhead <= mask_radius
    head = head * head_mask.astype(float)

    # ---------- TAIL ----------
    tail_axis = (Xg-center[0])*tail_dir[0] + (Yg-center[1])*tail_dir[1] + (Zg-center[2])*tail_dir[2]
    tail_axis_stretched = tail_axis * COMET_TAIL_LENGTH_FACTOR
    tailmask = tail_axis_stretched > 0
    rcyl = np.sqrt(((Yg-center[1]) - tail_dir[1]*tail_axis_stretched)**2 +
                   ((Zg-center[2]) - tail_dir[2]*tail_axis_stretched)**2)
    width = (0.25 + 0.13*np.sin(t*1.5 + Yg*2.2 + Zg*2.2 + t*0.3)) * COMET_TAIL_WIDTH_FACTOR
    turbulence = 0.8 + 0.2*np.sin(2.5*t + Zg*1.7 + Xg*1.3 + t*0.5)
    dist = np.abs(center[0] - Xg) * turbulence
    vortex = 1.0 + 0.15 * np.sin(t*2.5 + Yg*3.0 + t*0.4) * \
            np.sin(3*np.arctan2(Zg-center[2], Yg-center[1]) + t*2 + t*0.6)
    tail = tailmask * np.exp(-dist*1.25) * np.exp(-rcyl**2/width**2) * vortex
    flicker = 0.8 + 0.2 * np.sin(t * 0.25)

    density = head * 1.5 + tail * flicker
    density = np.clip(density, 0, None)
    density = gaussian_filter(density, sigma=1.0)

    # === BOOST DENSITY TO PRESERVE TOTAL "MASS" ===
    density *= DENSITY_SCALE

    return density

for name, (nucleus_actor, _) in comet_actors.items():
    key = f"comet_{name.lower()}"
    halo_radius_km = COMET_ELEMENTS[name]["radius_km"] * 120
    head_scale = halo_radius_km * SIZEKM_TO_SCENE / COMET_GRID_HALF
    density = make_volumetric_comet_body(
        center=np.array([0.0, 0.0, 0.0]),
        velocity=np.array([0.0, 0.0, 0.0]),
        t=0.0,
        head_scale=head_scale
    )
    grid = pv.ImageData()
    grid.dimensions = np.array(density.shape) + 1
    step = COMET_GRID[1] - COMET_GRID[0]
    grid.origin = (COMET_GRID.min() - step / 2,) * 3
    grid.spacing = (step,) * 3
    density_flat = density.flatten(order="F")
    density_flat[density_flat < 0.1] = 0.0
    grid.cell_data["comet"] = density_flat
    # Use real-world specific cmap for this comet
    comet_cmap = COMET_CMAPS.get(name.lower(), "Blues")  # Default to Blues if not found
    vol_actor = pl.add_volume(
        grid,
        scalars="comet",
        cmap=comet_cmap,  # Per-comet color map
        opacity=[0, 0, 0.1, 0.4, 0.7, 0.95],
        blending="maximum",
        shade=True,
        name=f"volbody_{name}"
    )
    vol_grids[key] = grid
    vol_actors[key] = vol_actor
print("[INFO] Denser volumetric comet body added with real-world colors.")



# ------------------------------------------------------------------
#  APPLY HALOS TO COMETS (must be done after comet_actors creation)
# ------------------------------------------------------------------
print(f"[INFO] Creating halos for {len(comet_actors)} comets...")
for comet_name, comet_tuple in comet_actors.items():
    nucleus_actor = comet_tuple[0]
    comet_data = {
        "name": comet_name,
        "radius_km": COMET_ELEMENTS[comet_name]["radius_km"]
    }
      
    # OLD halo (keep for coma)
    try:
        halo = add_atmosphere_halo(nucleus_actor, comet_data, moon_colors=moon_colors)
        if halo is not None:
            full_key = f"comet_{comet_name.lower()}"
            halo_actors[full_key] = halo
            halo_parents[full_key] = nucleus_actor
            print(f"[INFO] Added halo for {comet_name}")
    except Exception as e:
        print(f"[ERROR] Failed to create halo for comet {comet_name}: {e}")

animator = SolarSystemAnimator(
    jupiter_actor, saturn_actor, neptune_actor, uranus_actor, mars_actor,
    earth_actor, venus_actor, mercury_actor,
    moon_actors,           # unified moon dictionary
    pl, ring_meshes, uranus_ring_actors, asteroid_actor, sun_actor,
    plutoactor=dwarf_actors["pluto"],
    erisactor=dwarf_actors["eris"],
    haumeaactor=dwarf_actors["haumea"],
    makemakeactor=dwarf_actors["makemake"],
    ceresactor=dwarf_actors["ceres"],
    comet_actors=comet_actors,
)

print("-" * 80)
print(f"Loaded system with {sum(len(v) for v in moon_actors.values())} total moons across {len(moon_actors)} planets.")
print("Controls: 'j' Jupiter, 's' Saturn, 'n' Neptune, 'u' Uranus, 'm' Mars, 'g' Earth, 'v' Venus, 'y' Mercury, 'a' Asteroid Belt, 't' Sun, 'c' Ceres, 'o' Pluto, 'i' Eris, 'h' Haumea, 'k' Makemake")
print(f"System scaled by {GLOBAL_SCALE_MULTIPLIER}x")
print("-" * 80)

pl.show(title=f"Solar System ({GLOBAL_SCALE_MULTIPLIER}x Scale)", auto_close=False, interactive=False)
interactor = pl.render_window.GetInteractor()
style = interactor.GetInteractorStyle()
style.SetMotionFactor(3.5)

try:
    if hasattr(interactor, "AddObserver"):
        interactor.AddObserver("TimerEvent", animator.update)
    elif hasattr(interactor, "add_observer"):
        interactor.add_observer("TimerEvent", animator.update)
except Exception as e:
    print(f"[WARN] Could not bind timer event: {e}")

timer_interval = int(1000 / FPS)
try:
    timer_id = interactor.CreateRepeatingTimer(timer_interval)
    interactor.Start()
    interactor.DestroyTimer(timer_id)
except Exception as e:
    print(f"[WARN] Timer handling error: {e}")

pl.close()
