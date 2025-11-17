A fully dynamic, physics-driven 3D Solar System simulation built using PyVista, VTK, NumPy, and real astronomical parameters.
This project features:
* Accurate orbital mechanics (Keplerian elements from JPL Horizons)
* True axial tilts & rotation periods for planets and moons
* 150+ moons with individual orbital elements
* Dynamic zoom limits per object
* Realistic irregular potato-shaped bodies
* Saturn’s ring system
* Comets with real orbital paths
* Scaled radii & distances for stable visualization
* High-performance rendering with level-of-detail logic
<img width="1358" height="663" alt="image" src="https://github.com/user-attachments/assets/5dbf978b-78ad-4503-ace0-a582004380af" />
Full Controls Overview
Navigate the simulation interactively with these keyboard shortcuts (case-sensitive; press once to activate). Focus on a planet to unlock moon viewing (numeric keys). Mouse: Left-drag to rotate view, right-drag to pan, wheel to zoom (clamped for clipping prevention). Arrow keys: Pan/orbit camera (use for gradual full-system zoom-out to reveal tiny outer moons without clamping issues—wheel zoom may obscure distant objects).
Key,Action,Details
t,Focus Sun,Centers camera on Sun (min zoom: ~5 units).
y,Focus Mercury,Zooms to innermost planet (no moons).
v,Focus Venus,Targets cloudy world (no moons).
g,Focus Earth,Homes in on home planet; press 1 for Luna.
m,Focus Mars,"Views Red Planet; press 1 for Phobos, 2 for Deimos."
j,Focus Jupiter,Gas giant view; press 1–30+ for 95 moons (inner: Metis=1 to Callisto=8; outer: tiny irregulars).
s,Focus Saturn,"Ringed wonder; press 1–83 for moons (e.g., 1=Mimas, 20=Titan)."
u,Focus Uranus,"Sideways spinner; press 1–27 for moons (e.g., 1=Cordelia, 15=Titania)."
n,Focus Neptune,"Windy blue; press 1–14 for moons (e.g., 1=Naiad, 8=Triton)."
o,Focus Pluto,Dwarf at edge; no moons here (Charon via separate logic).
i,Focus Eris,Distant dwarf (if enabled).
h,Focus Haumea,Elongated Kuiper object.
k,Focus Makemake,Icy world.
c,Focus Ceres,Asteroid belt dwarf.
a,Focus Asteroid Belt,Pans to main belt (~2.7 AU); reveals thousands of points.
z,Focus Halley Comet,Tracks famous periodic comet.
x,Focus Hale-Bopp,Great Comet of 1997.
l,Focus Encke,Short-period comet.
d,Focus Lovejoy,Recent sungrazer.
0–9 (Numerics),Select Moon,"After planet focus (e.g., j then 123 for Jupiter's moon #123); buffers multi-digit."
Arrow Keys,Pan/Orbit Camera,Smooth full-system navigation; ideal for unveiling tiny moons (avoids wheel clamp).
ESC,Pause/Exit,Halts animation or closes viewer.

Pro Tip: For crisp views of faint outer moons (e.g., Jupiter's 95+), start at Sun focus (t), then arrow-key zoom out to ~30 AU—wheel clamping keeps close-ups sharp but hides peripherals.

SUN
<img width="1364" height="665" alt="image" src="https://github.com/user-attachments/assets/1a4be5fb-58ce-4a51-8cf3-6805d13dae0d" />
The Sun serves as the central gravitational anchor and visual centerpiece of this solar system simulation, rendered as a glowing, rotating sphere using PyVista. Its implementation draws from real astronomical data (e.g., NASA/JPL sources) for authenticity, while incorporating dynamic effects like barycentric wobble due to planetary influences. Below is a breakdown of its key features, parameters, and behaviors.
* Radius: 696,340 km (scaled to ~0.02x relative to orbital distances for visualization; SUN_RADIUS = REAL_SUN_RADIUS_KM * SIZEKM_TO_SCENE).
* Gravitational Parameter (GM): 1.327 × 10¹¹ km³/s² (GM_SUN), used for orbital calculations of all bodies.
* Rotation Period: ~25 Earth days at the equator (differential rotation approximated; ROTATION_PERIOD["sun"] = 25.0 * 86400 seconds), converted to ~0.000166° per simulated second (ROTATION_DEG_PER_SIMSEC["sun"]).
* Axial Tilt: 7.25° relative to the ecliptic plane (AXIAL_TILTS["sun"]), applied via actor.RotateX(AXIAL_TILTS["sun"]) for realistic obliquity.
These values ensure Keplerian orbits and rotations align with J2000 epoch standards, with positions computed dynamically via kepler_to_state functions.

* Rendering: Created as a pv.Sphere with a bright yellow/orange color (color="yellow", opacity=0.8) and optional PBR shading for a photospheric glow (pbr=True, roughness=0.1). Textures (e.g., from examples.load_globe_texture()) can be layered for surface details like sunspots.
* Barycentric Wobble: The Sun isn't static—it "wobbles" slightly due to Jupiter and Saturn's pull, calculated in the update() loop.
  
* Rotation: Continuously rotates around its Z-axis in the animation loop (actor.RotateZ(deg % 360.0)), simulating solar day-night cycles (though simplified, as the Sun lacks a solid surface).
* Role in Simulation: Acts as the reference point for all distances (e.g., comet sublimation at 3 AU via SUBLIMATION_DISTANCE), with its position converted to km-scale for heliocentric calcs (sun_pos_km = sun_pos / KM_TO_SCENE).
* Focus on Sun: Press t to center the camera on the Sun (focus_sun() method), with automatic clipping range reset for smooth zooming (min distance: ~5 units).

This implementation balances performance (e.g., no heavy volumetric effects on the Sun) with educational value, making it easy to visualize solar dominance in the system. For enhancements, consider adding solar flares via particle systems or real-time data integration from NASA's Solar Dynamics Observatory (SDO).



MERCURY
<img width="1365" height="659" alt="image" src="https://github.com/user-attachments/assets/b926ef78-9f44-4bff-a8f6-a5d48f29ce81" />
Mercury, the smallest and innermost planet in our solar system, orbits closest to the Sun at an average distance of 0.387 AU. As a rocky terrestrial world, it boasts a heavily cratered, Moon-like surface scarred by ancient impacts, with vast plains of cooled lava and steep escarpments. Notably, Mercury lacks a true atmosphere; instead, it maintains a razor-thin exosphere—a diffuse envelope of scattered atoms primarily consisting of oxygen (O₂, ~42%), sodium (Na, ~29%), hydrogen (H₂, ~22%), helium (He, ~6%), and potassium (K)—continuously replenished by solar wind bombardment and micrometeorite strikes. This exosphere provides negligible insulation, exposing the planet to brutal temperature swings from a scorching 427°C (800°F) during the day to a frigid -173°C (-280°F) at night, and offers scant protection from solar radiation. NASA's Messenger mission (2008–2015) revealed these details, highlighting Mercury's iron-rich core (making up ~75% of its mass) and 3:2 spin-orbit resonance, where it rotates three times for every two orbits.
In this simulation, Mercury is modeled with high fidelity using Keplerian orbital elements and real-time dynamics, emphasizing its slow rotation and proximity to the Sun for educational visualization.
<img width="1365" height="680" alt="image" src="https://github.com/user-attachments/assets/c7dd0d4c-5178-4a62-91aa-606eb0f6fb2f" />

* Radius: 2,440 km (scaled to MERCURY_RADIUS = REAL_MERCURY_RADIUS_KM * SIZEKM_TO_SCENE for scene fitting).
* Gravitational Parameter (GM): 2.2032 × 10⁴ km³/s² (GM_MERCURY), influencing potential moon simulations (though Mercury has none).
* Orbital Distance: 0.387 AU (MERCURY_AU), yielding MERCURY_ORBIT_RADIUS = MERCURY_AU * AU_SCALE.
* Rotation Period: 58.646 Earth days (sidereal; ROTATION_PERIOD["mercury"] = 58.646 * 86400 seconds), resulting in ~0.000066° per simulated second (ROTATION_DEG_PER_SIMSEC["mercury"])—reflecting its 3:2 resonance (not fully simulated here for simplicity).
* Axial Tilt: 0.03° (AXIAL_TILTS["mercury"]), nearly equatorial, applied via actor.RotateX(AXIAL_TILTS["mercury"]).
* Rendering: Depicted as a grayish sphere (pv.Sphere with color="gray", opacity=1.0) textured to mimic its regolith-covered, cratered terrain (e.g., via procedural noise or Messenger imagery). No atmospheric effects are added due to its exosphere's sparsity, but subtle glow halos can be toggled for solar proximity visualization.
  <img width="1364" height="660" alt="image" src="https://github.com/user-attachments/assets/32576a2e-842c-42b2-bdce-ad1f4e486b1b" />

* Orbital Motion: Updated dynamically in the update() loop using osculating elements: set_body_position_from_elements("mercury", self.mercury_actor, orbital_elements["mercury"], target_jd). This computes heliocentric position in km, scaled to scene units (km_to_scene), with smoothing for fluid animation.
* Rotation: Slowly spins prograde around its Z-axis (actor.RotateZ(deg % 360.0)), with initial orientation adjusted for Greenwich Mean Sidereal Time (GMST) alignment:merc_gmst = gmst_deg * (ROTATION_PERIOD["earth"] / abs(ROTATION_PERIOD["mercury"])) % 360.0
self.mercury_actor.RotateZ(-merc_gmst).
* Focus on Mercury: Press y to center the camera (set_planet_focus('mercury', self.mercury_actor)), auto-adjusting clipping for close-up crater views (min distance: ~MERCURY_RADIUS * 20 units).

This module captures Mercury's harsh, airless isolation, ideal for demonstrating extreme orbital mechanics. Future updates could integrate real-time solar wind data from NASA's ACE satellite for dynamic exosphere effects.

VENUS
<img width="1364" height="671" alt="image" src="https://github.com/user-attachments/assets/cb6360c0-ec50-423f-ac80-74dd55660d15" />
Venus, often dubbed Earth's "evil twin" due to its similar size and rocky composition, is the second planet from the Sun at an average distance of 0.723 AU. Blanketed by a dense, toxic atmosphere composed mainly of carbon dioxide (96.5%) with traces of nitrogen (3.5%) and sulfuric acid clouds, it experiences a runaway greenhouse effect that traps heat, making it the hottest planet in the solar system with surface temperatures averaging 464°C (867°F)—hot enough to melt lead. Atmospheric pressure at the surface is a crushing 92 times that of Earth, rendering the landscape a hellish vista of volcanic plains, massive lava domes, and continent-sized highlands like Ishtar Terra. Its rotation is retrograde and excruciatingly slow (243 Earth days per spin, longer than its 225-day orbit), resulting in a "day" longer than its year, while radar mapping from NASA's Magellan mission (1989–1994) unveiled a surface dominated by volcanism and few impact craters. Venus's thick cloud cover reflects 70% of sunlight, giving it a brilliant white appearance from Earth, but hides a dynamic world probed by ongoing missions like Japan's Akatsuki (2015–present).
In this simulation, Venus is rendered with emphasis on its retrograde spin and atmospheric opacity, using Keplerian dynamics for precise orbital tracking and simplified cloud visuals for educational impact.
<img width="1365" height="670" alt="image" src="https://github.com/user-attachments/assets/f72758b6-b9e5-495d-a0e6-843c116c978f" />
* Radius: 6,052 km (scaled to VENUS_RADIUS = REAL_VENUS_RADIUS_KM * SIZEKM_TO_SCENE for relative sizing).
* Gravitational Parameter (GM): 3.249 × 10⁵ km³/s² (GM_VENUS), for potential future enhancements like surface simulations.
* Orbital Distance: 0.723 AU (VENUS_AU), yielding VENUS_ORBIT_RADIUS = VENUS_AU * AU_SCALE.
* Rotation Period: -243.025 Earth days (retrograde; ROTATION_PERIOD["venus"] = -243.025 * 86400 seconds), yielding ~ -0.000059° per simulated second (ROTATION_DEG_PER_SIMSEC["venus"])—the negative sign enforces backward spin.
* Axial Tilt: 177.4° (AXIAL_TILTS["venus"]), effectively a 2.6° tilt in the opposite direction due to retrograde motion, applied via actor.RotateX(AXIAL_TILTS["venus"]).
* Rendering: Modeled as a creamy-yellow sphere (pv.Sphere with color="yellow" or "orange" for cloud reflection, opacity=0.9) with optional procedural textures simulating sulfuric acid haze (e.g., layered noise for cloud bands). A faint atmospheric halo (add_atmosphere_halo) can represent the upper cloud deck, with opacity scaled to mimic 99% light scattering.
* Orbital Motion: Dynamically positioned in the update() loop via osculating elements: set_body_position_from_elements("venus", self.venus_actor, orbital_elements["venus"], target_jd). Converts heliocentric km coordinates to scene units (km_to_scene), with velocity smoothing for seamless animation.
<img width="1365" height="663" alt="image" src="https://github.com/user-attachments/assets/e0243889-53f5-4ce4-8626-079cfb084e44" />

* Rotation: Retrograde spin around the Z-axis (actor.RotateZ(-deg % 360.0) to account for negative period), initialized with GMST-adjusted orientation: venus_gmst = gmst_deg * (ROTATION_PERIOD["earth"] / abs(ROTATION_PERIOD["venus"])) % 360.0
self.venus_actor.RotateZ(-venus_gmst).
* Atmospheric Effects: While not volumetrically simulated (for performance), proximity to the Sun triggers subtle glow enhancements in halo logic (r_km = np.linalg.norm(pos_km - sun_pos_km)), evoking Venus's phase changes and superior conjunction visibility.
* Venus's implementation highlights its paradoxical "slow dance" with the Sun, keeping compute overhead low without moons.
* Focus on Venus: Press v to zoom in (set_planet_focus('venus', self.venus_actor)), with clipping auto-reset for detailed cloud layer views (min distance: ~VENUS_RADIUS * 20 units).

This captures Venus's veiled, volcanic fury, perfect for illustrating greenhouse extremes. Enhancements could include dynamic cloud rotation from Akatsuki data.

EARTH
<img width="1358" height="648" alt="image" src="https://github.com/user-attachments/assets/0c10b204-4eeb-4b8b-965c-41cad178cbc9" />
Earth, our blue marble and the only known cradle of life in the solar system, orbits at a comfortable 1 AU from the Sun, enabling a temperate climate that supports diverse ecosystems across seven continents, vast oceans covering 71% of its surface, and a protective atmosphere of nitrogen (78%) and oxygen (21%) that shields against cosmic rays while driving dynamic weather patterns like hurricanes and auroras. This gaseous envelope, laced with water vapor and trace gases, fosters the greenhouse effect that keeps global averages at a life-sustaining 15°C (59°F), though human-induced changes are accelerating warming. Earth's iron-nickel core generates a magnetic field deflecting solar wind, while its active geology—plate tectonics, volcanoes, and earthquakes—recycles its crust. Iconic missions like Apollo (1969–1972) and ongoing observations from satellites such as Terra and Aqua have mapped its ever-changing visage, revealing seasonal ice caps, swirling cloud vortices, and bioluminescent oceans under the night side.
In this simulation, Earth is the interactive hub, complete with axial precession, lunar orbit, and atmospheric glow, leveraging Keplerian elements for synchronized day-night cycles and tidal locking visuals.
<img width="1358" height="652" alt="image" src="https://github.com/user-attachments/assets/79858c6e-650d-4f77-b4e1-07a6432c35f7" />
A standout feature is Luna (Latin for Moon), Earth's sole natural satellite and the fifth-largest in the solar system, orbiting at an average 384,400 km with a diameter of 3,475 km—about a quarter of Earth's. Formed ~4.5 billion years ago from debris of a Mars-sized impactor (Theia), Luna stabilizes Earth's axial tilt for consistent seasons, drives tides influencing marine life and coastal erosion, and once hosted a molten magma ocean now solidified into basaltic maria (dark "seas") and rugged highlands pocked by craters like Tycho. Its thin exosphere (mostly argon and helium) and lack of atmosphere expose it to micrometeorites, creating perpetual "moonquakes" from tidal stresses. NASA's Artemis program (2020s) aims to return humans, building on Apollo's legacy of 382 kg of lunar samples.
<img width="1365" height="666" alt="image" src="https://github.com/user-attachments/assets/6db37606-291d-47f1-8c4b-81e51f6c1cac" />

* Radius: 6,371 km (scaled to EARTH_RADIUS = REAL_EARTH_RADIUS_KM * SIZEKM_TO_SCENE for balanced proportions).
* Gravitational Parameter (GM): 3.986 × 10⁵ km³/s² (GM_EARTH), central to moon orbital mechanics.
* Orbital Distance: 1.0 AU (EARTH_AU), defining EARTH_ORBIT_RADIUS = EARTH_AU * AU_SCALE as the baseline unit.
* Rotation Period: 0.99726968 Earth days (sidereal; ROTATION_PERIOD["earth"] = 0.99726968 * 86400 seconds), equating to ~0.004° per simulated second (ROTATION_DEG_PER_SIMSEC["earth"]) for diurnal motion.
* Axial Tilt: 23.44° (AXIAL_TILTS["earth"]), driving seasons via actor.RotateX(AXIAL_TILTS["earth"]).
* Rendering: A vibrant blue-green pv.Sphere (color="blue" with landmasses in green/brown via texture mapping, e.g., NASA's Blue Marble composite), overlaid with a translucent atmospheric halo (add_atmosphere_halo) simulating ozone layer scattering (opacity ~0.2 for Rayleigh blue skies). Clouds could be procedurally added as semi-transparent meshes for weather dynamism.
  <img width="1352" height="656" alt="image" src="https://github.com/user-attachments/assets/b175ed37-4b97-47b5-96fc-86d496c182d8" />

* Orbital Motion: Real-time heliocentric positioning in update(): set_body_position_from_elements("earth", self.earth_actor, orbital_elements["earth"], target_jd). Scaled from km to scene coordinates (km_to_scene), with inertial frame smoothing to mimic elliptical path.
* Rotation: Prograde spin synchronized to sidereal day (actor.RotateZ(deg % 360.0)), initialized with GMST for accurate longitude.
* Focus on Earth: Press g to home in (set_planet_focus('earth', self.earth_actor)), enabling moon navigation via numeric keys (1 for Luna; min distance: ~EARTH_RADIUS * 20 units).

MARS
<img width="1365" height="655" alt="image" src="https://github.com/user-attachments/assets/52da038b-17bc-46a7-9a22-c22756a929a1" />
Mars, the enigmatic Red Planet, orbits at 1.524 AU from the Sun, its rusty hue stemming from iron oxide (rust) dust blanketing vast deserts, towering canyons like Valles Marineris (four times Grand Canyon's length), and the solar system's largest volcano, Olympus Mons (22 km high). A thin atmosphere—mostly carbon dioxide (95.3%) with nitrogen (2.7%) and argon (1.6%)—creates wispy clouds, global dust storms that can engulf the planet for months, and seasonal polar ice caps of water and CO₂ frost, hinting at ancient rivers and potential subsurface oceans. Surface temperatures swing from -60°C ( -76°F) average to 20°C (68°F) highs, with low pressure (0.6% of Earth's) making liquid water unstable. Pioneering missions like Viking (1976), Pathfinder (1997), and NASA's Perseverance rover (2021–present) have uncovered organic molecules and methane plumes, fueling dreams of past microbial life and future human exploration.
In this simulation, Mars is vividly animated with its two potato-shaped moons, dust-red terrain, and orbital eccentricity, using Keplerian mechanics to showcase its elliptical path and synodic oppositions.
<img width="1351" height="659" alt="image" src="https://github.com/user-attachments/assets/783afb5c-0600-4da8-aae3-02345e188ed3" />
* Radius: 3,390 km (scaled to MARS_RADIUS = REAL_MARS_RADIUS_KM * SIZEKM_TO_SCENE for compact visualization).
* Gravitational Parameter (GM): 4.282 × 10⁴ km³/s² (GM_MARS), governing Phobos and Deimos orbits.
* Orbital Distance: 1.524 AU (MARS_AU), producing MARS_ORBIT_RADIUS = MARS_AU * AU_SCALE.
* Rotation Period: 1.025957 Earth days (sidereal; ROTATION_PERIOD["mars"] = 1.025957 * 86400 seconds), translating to ~0.0039° per simulated second (ROTATION_DEG_PER_SIMSEC["mars"]) for sol-like days.
* Axial Tilt: 25.19° (AXIAL_TILTS["mars"]), similar to Earth's for comparable seasons, via actor.RotateX(AXIAL_TILTS["mars"]).
* Rendering: A reddish pv.Sphere (color="red" or "maroon" with procedural craters and polar caps via noise textures, opacity=1.0), accented by a sparse atmospheric halo (add_atmosphere_halo) for dust storm hints—opacity modulated by solar distance (r_km = np.linalg.norm(pos_km - sun_pos_km)).
* Orbital Motion: Updated in the update() loop with osculating elements: set_body_position_from_elements("mars", self.mars_actor, orbital_elements["mars"], target_jd).  Scaled heliocentric positions (km_to_scene) with eccentricity-driven speed variations for realistic perihelion rushes.
* Rotation: Prograde axial spin (actor.RotateZ(deg % 360.0)), GMST-initialized for areographic accuracy.
  <img width="1353" height="657" alt="image" src="https://github.com/user-attachments/assets/b114dc97-0adc-4dea-819d-dc959788122d" />

* Moons and Atmosphere: Phobos (press 1) and Deimos (press 2) orbit in moon_actors["mars"], tidally locked with faint exosphere glows; optional dust animation could swirl via particle emitters during opposition.
