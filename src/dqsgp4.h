#pragma once

/**
 * @file dqsgp4.h
 * @brief Umbrella header: the library's public surface in one include.
 *
 * The library is built around the dual-quaternion numerical propagator: the
 * SE(3) state is integrated with pluggable force models, every value carries
 * a tracked error budget, and the numeric type T is a template parameter.
 * The analytical SGP4/SDP4 implementation is part of the library but is not
 * exported here; it serves as the verification reference and recovers epoch
 * states from TLEs (a TLE's elements are defined in terms of the SGP4 model,
 * so the TLE adapters below evaluate it internally at t = 0). Include
 * sgp4/sgp4_propagator.h explicitly to use the analytical model directly;
 * a test-suite check keeps the numerical propagator's dependence on it
 * confined to the two adapter headers.
 *
 * The map (each header names its derivation note in design/derivations/):
 *
 *   THE PROPAGATOR
 *     dynamics/dq_sgp4_propagator.h  the TLE-driven facade: propagate(minutes),
 *                                 authentic/boosted modes, DqForceOptions
 *                                 presets (lunisolar / drag / SRP), extra_forces
 *     dynamics/propagator.h       the force-list + integrator engine
 *                                 (propagate_to fixed-step, propagate_adaptive
 *                                 RKF7(8)); takes Cartesian states directly
 *     dynamics/propagatable.h     the propagate() concept shared by both
 *                                 propagator families
 *
 *   I/O & STATE
 *     tle/tle_parser.h, tle/omm_parser.h, tle/omm_xml_parser.h
 *     dynamics/state_from_tle.h   TLE -> epoch State (evaluates the SGP4 model
 *                                 at t = 0)
 *     dynamics/state_conversion.h State <-> StateVector adapters (m <-> km)
 *
 *   FORCES (each verified against an independent reference; see the headers)
 *     forces/geopotential.h       monopole+zonal+tesseral in one Cunningham pass
 *     forces/drag.h               the DensityModel interface + make_drag
 *     atmosphere/exponential_table.h  the Vallado 8-4 static atmosphere
 *     forces/third_body.h         Sun/Moon tidal force (verified against DE430)
 *     forces/srp.h                cannonball radiation pressure + shadow
 *
 *   EPHEMERIS & FRAMES
 *     ephemeris/sun_meeus.h, ephemeris/moon_meeus.h, ephemeris/body_position_gcrs.h
 *     astronomy/epoch.h           TimeScale + two-part Epoch
 *     astronomy/frames.h          rotation primitives + IAU2006 precession chain
 *     astronomy/earth_rotation.h  ERA/GMST06/GAST, polar motion, GCRS->ITRS
 *     astronomy/sidereal_time.h   the Aoki-82 GMST (the TEME convention)
 *
 *   FOUNDATIONS
 *     constants/constants_provider.h  wgs72/wgs84/grs80 with error budgets
 *     math/tracked_value.h        the three-error type everything computes in
 *
 * Not pulled in here: the analytical SGP4/SDP4 headers (sgp4/*.h, include
 * explicitly), the standalone integrator family (integrators/*.h — the facade
 * already selects RK4/RKF78), and the earlier low-precision ephemeris
 * instances (ephemeris/solar_ephemeris.h etc., superseded by the Meeus
 * instances).
 */

#include "astronomy/earth_rotation.h"
#include "astronomy/epoch.h"
#include "astronomy/frames.h"
#include "astronomy/sidereal_time.h"
#include "atmosphere/exponential_table.h"
#include "constants/constants_provider.h"
#include "dynamics/dq_sgp4_propagator.h"
#include "dynamics/propagatable.h"
#include "dynamics/propagator.h"
#include "dynamics/state_conversion.h"
#include "dynamics/state_from_tle.h"
#include "ephemeris/body_position_gcrs.h"
#include "ephemeris/moon_meeus.h"
#include "ephemeris/sun_meeus.h"
#include "forces/drag.h"
#include "forces/geopotential.h"
#include "forces/srp.h"
#include "forces/third_body.h"
#include "math/tracked_value.h"
#include "tle/omm_parser.h"
#include "tle/omm_xml_parser.h"
#include "tle/tle_parser.h"
