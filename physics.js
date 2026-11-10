/*
  NAVUMSHA physics engine (physics.js)
  ------------------------------------
  Used by lab.html. No libraries needed. Also runs in Node for testing:  node physics.js

  Coordinates (all in SI units: metres, seconds, kilograms, newtons):
    x = distance along the launch direction (+x is "forward")
    y = height above the y = 0 reference (up is positive)
    z = sideways drift (wind direction 90 degrees and sidespin push towards +z)

  What is solved
    Vacuum   : x = v0 cos(a) t,   y = h0 + v0 sin(a) t - g t^2 / 2
    With air : a = g + drag + Magnus, integrated with RK4 (4th-order Runge-Kutta)
      drag   : F = 1/2 rho Cd A |v_rel| v_rel (opposes the velocity relative to the air, v_rel = v - wind)
      Magnus : F = 1/2 rho A C_L |v_rel|^2 in the direction of (spin x v_rel),
               with spin parameter S = r |spin x v_rel| / |v_rel|^2 and C_L = S / (1 + 2 S)
               (a simple, widely used teaching model; real balls differ with surface and speed)
      Backspin (+) pushes up, topspin (-) pushes down, positive sidespin pushes towards +z.
      Spin is assumed constant during the flight (no spin decay).
    Landing  : the first time the object comes DOWN through the landing level. The exact moment is
               found by re-stepping the last RK4 step with bisection, so there is no step-size error.
    Time step: small and fixed for each flight (flight time / 4000, between 1 ms and 20 ms). In air,
               it is also shortened when drag is so strong that a long step could become unstable.
*/
(function (root) {
  'use strict';

  var PI = Math.PI;
  var LIM = {
    v0: [0, 200], angle: [0, 90], h0: [0, 500], ly: [-100, 200], mass: [0.001, 200], radius: [0.005, 0.5],
    cd: [0, 2], rho: [0, 100], windSpeed: [0, 50], windDir: [0, 360], backspin: [-6000, 6000], sidespin: [-6000, 6000], g: [0.1, 30]
  };
  var DEF = { v0: 30, angle: 45, h0: 0, ly: 0, mass: 0.43, radius: 0.11, cd: 0.25, rho: 1.225, windSpeed: 0, windDir: 0, backspin: 0, sidespin: 0, g: 9.81 };
  var SCALE_HEIGHT = 8500;     // metres: air density falls by a factor e every 8.5 km (halves about every 5.9 km)
  var MAX_TIME = 1800;         // seconds: flights longer than this are stopped (30 minutes)
  var MAX_STEPS = 800000;
  var MAX_SAMPLES = 12000;

  function clamp(v, a, b) { return Math.min(b, Math.max(a, v)); }
  function cosd(d) { var r = ((d % 360) + 360) % 360; if (r === 90 || r === 270) return 0; if (r === 0) return 1; if (r === 180) return -1; return Math.cos(d * PI / 180); }
  function sind(d) { var r = ((d % 360) + 360) % 360; if (r === 0 || r === 180) return 0; if (r === 90) return 1; if (r === 270) return -1; return Math.sin(d * PI / 180); }

  /* ---------- input checking: no negatives, no NaN, no divide-by-zero ---------- */
  function sanitize(p) {
    var q = {}, notes = [], k;
    p = p || {};
    for (k in LIM) {
      var v = Number(p[k]);
      if (p[k] === null || p[k] === undefined || p[k] === '' || !isFinite(v)) { v = DEF[k]; notes.push(k + ' was not a number, so ' + DEF[k] + ' was used'); }
      if (k === 'windDir') v = ((v % 360) + 360) % 360;
      var c = clamp(v, LIM[k][0], LIM[k][1]);
      if (c !== v) notes.push(k + ' was limited to ' + c);
      q[k] = c;
    }
    q.airOn = !!p.airOn; q.altDensity = !!p.altDensity;
    return { p: q, notes: notes };
  }

  /* ---------- air forces (shared by the integrator and the arrow display) ---------- */
  function makeAero(p) {
    var air = !!(p.airOn && p.rho > 0), A = PI * p.radius * p.radius, m = p.mass, r = p.radius;
    var wx = p.windSpeed * cosd(p.windDir), wz = p.windSpeed * sind(p.windDir);
    var K = PI / 30, omY = -p.sidespin * K, omZ = p.backspin * K;   // rpm -> rad/s; spin vector = (0, -sidespin, backspin)
    return {
      air: air,
      /* writes accelerations (m/s^2): dx,dy,dz = drag, mx,my,mz = Magnus, rate = a stiffness estimate */
      parts: function (vx, vy, vz, y, o) {
        o.dx = o.dy = o.dz = o.mx = o.my = o.mz = 0; o.rate = 0;
        if (!air) return;
        var rho = p.rho; if (p.altDensity) rho *= Math.exp(-y / SCALE_HEIGHT);
        var rx = vx - wx, ry = vy, rz = vz - wz, sp = Math.sqrt(rx * rx + ry * ry + rz * rz);
        if (sp < 1e-9) return;
        var kd = 0.5 * rho * p.cd * A / m * sp;
        o.dx = -kd * rx; o.dy = -kd * ry; o.dz = -kd * rz;
        var cx = omY * rz - omZ * ry, cy = omZ * rx, cz = -omY * rx, cm = Math.sqrt(cx * cx + cy * cy + cz * cz);
        if (cm > 1e-12) {
          var S = r * cm / (sp * sp), CL = S / (1 + 2 * S), f = 0.5 * rho * A * CL * sp * sp / m / cm;
          o.mx = f * cx; o.my = f * cy; o.mz = f * cz;
        }
        o.rate = 0.5 * rho * A / m * sp * (p.cd + 0.5);
      }
    };
  }
  /* forces in newtons for the arrows in the lab */
  function forcesAt(st, p0) {
    var p = sanitize(p0).p, a = makeAero(p), o = {};
    a.parts(st.vx, st.vy, st.vz || 0, st.y, o);
    var m = p.mass;
    return { drag: [o.dx * m, o.dy * m], magnus: [o.mx * m, o.my * m], dragZ: o.dz * m, magnusZ: o.mz * m, weight: m * p.g };
  }

  /* ---------- closed-form vacuum formulas (the textbook answers) ---------- */
  function vacuumFormula(p0) {
    var p = sanitize(p0).p, th = p.angle, vx0 = p.v0 * cosd(th), vy0 = p.v0 * sind(th), g = p.g, dh = p.h0 - p.ly;
    var disc = vy0 * vy0 + 2 * g * dh;
    if (disc < 0) return { reaches: false, range: 0, time: 0, maxHeight: p.h0 + (vy0 > 0 ? vy0 * vy0 / (2 * g) : 0), impactSpeed: 0, impactAngle: 0 };
    var t = (vy0 + Math.sqrt(disc)) / g; if (t < 0) t = 0;
    var vyEnd = vy0 - g * t, speed = Math.sqrt(vx0 * vx0 + vyEnd * vyEnd);
    return {
      reaches: true, range: vx0 * t, time: t, maxHeight: p.h0 + (vy0 > 0 ? vy0 * vy0 / (2 * g) : 0),
      impactSpeed: speed, impactAngle: Math.atan2(-vyEnd, Math.abs(vx0)) * 180 / PI
    };
  }
  /* best angle (degrees) with no air: arctan(v / sqrt(v^2 + 2 g dh)), dh = launch height - landing level */
  function optimalAngleFormula(v0, g, dh) {
    var rad = v0 * v0 + 2 * g * dh;
    if (!(v0 > 0) || rad < 0) return null;
    if (rad === 0) return 90;
    return Math.atan(v0 / Math.sqrt(rad)) * 180 / PI;
  }
  /* longest range with no air: (v / g) sqrt(v^2 + 2 g dh) */
  function maxRangeFormula(v0, g, dh) {
    var rad = v0 * v0 + 2 * g * dh;
    return rad < 0 ? 0 : v0 * Math.sqrt(rad) / g;
  }

  /* ---------- RK4 integrator ---------- */
  function simulate(p0, opt) {
    opt = opt || {};
    var p = sanitize(p0).p, g = p.g, L = p.ly, H0 = p.h0;
    var aero = makeAero(p), air = aero.air, o = {};
    var vx0 = p.v0 * cosd(p.angle), vy0 = p.v0 * sind(p.angle);

    function acc(vx, vy, vz, y, out) {
      aero.parts(vx, vy, vz, y, o);
      out[0] = o.dx + o.mx; out[1] = -g + o.dy + o.my; out[2] = o.dz + o.mz;
    }
    var a1 = [0, 0, 0], a2 = [0, 0, 0], a3 = [0, 0, 0], a4 = [0, 0, 0];
    function step(s, h, r) {
      var x = s[0], y = s[1], z = s[2], vx = s[3], vy = s[4], vz = s[5], h2 = h / 2;
      acc(vx, vy, vz, y, a1);
      var v2x = vx + h2 * a1[0], v2y = vy + h2 * a1[1], v2z = vz + h2 * a1[2];
      acc(v2x, v2y, v2z, y + h2 * vy, a2);
      var v3x = vx + h2 * a2[0], v3y = vy + h2 * a2[1], v3z = vz + h2 * a2[2];
      acc(v3x, v3y, v3z, y + h2 * v2y, a3);
      var v4x = vx + h * a3[0], v4y = vy + h * a3[1], v4z = vz + h * a3[2];
      acc(v4x, v4y, v4z, y + h * v3y, a4);
      var h6 = h / 6;
      r[0] = x + h6 * (vx + 2 * v2x + 2 * v3x + v4x);
      r[1] = y + h6 * (vy + 2 * v2y + 2 * v3y + v4y);
      r[2] = z + h6 * (vz + 2 * v2z + 2 * v3z + v4z);
      r[3] = vx + h6 * (a1[0] + 2 * a2[0] + 2 * a3[0] + a4[0]);
      r[4] = vy + h6 * (a1[1] + 2 * a2[1] + 2 * a3[1] + a4[1]);
      r[5] = vz + h6 * (a1[2] + 2 * a2[2] + 2 * a3[2] + a4[2]);
    }

    // choose the base step from a vacuum estimate of the flight time
    var disc = vy0 * vy0 + 2 * g * (H0 - L), Test = disc > 0 ? (vy0 + Math.sqrt(disc)) / g : (vy0 > 0 ? 2 * vy0 / g : 1);
    Test = Math.max(Test, 0.05);
    var dtBase = clamp(Test / 4000, 0.001, 0.02); if (opt.coarse) dtBase = Math.min(0.05, dtBase * 4);

    var S = [0, H0, 0, vx0, vy0, 0], S2 = [0, 0, 0, 0, 0, 0], T1 = [0, 0, 0, 0, 0, 0], T2 = [0, 0, 0, 0, 0, 0];
    var t = 0, steps = 0, rec = [0, 0, H0, 0, vx0, vy0, 0];
    var minX = 0, maxX = 0, minY = H0, maxY = H0, minZ = 0, maxZ = 0;
    var apex = { t: 0, x: 0, y: H0, z: 0 }, status = null, hitLevel = L, final = null, finalT = 0;

    while (true) {
      var h = dtBase;
      if (air) {
        aero.parts(S[3], S[4], S[5], S[1], o);
        if (o.rate > 0) h = Math.max(1e-5, Math.min(h, 0.2 / o.rate));
      }
      step(S, h, S2);
      steps++;
      var yP = S[1], yN = S2[1], hit = null, tauL = h, i, lo, hi, mid;

      if (yP >= L && yN < L) hit = L;
      else if (H0 < L && yP >= H0 && yN < H0) hit = H0;      // never reached the higher landing level

      if (hit !== null) {                                    // find the exact landing moment inside the step
        lo = 0; hi = h;
        for (i = 0; i < 60; i++) { mid = 0.5 * (lo + hi); step(S, mid, T1); if (T1[1] >= hit) lo = mid; else hi = mid; }
        tauL = 0.5 * (lo + hi);
      }
      if (S[4] > 0 && S2[4] <= 0) {                          // apex inside this step?
        lo = 0; hi = h;
        for (i = 0; i < 60; i++) { mid = 0.5 * (lo + hi); step(S, mid, T2); if (T2[4] > 0) lo = mid; else hi = mid; }
        var ta = 0.5 * (lo + hi);
        if (hit === null || ta <= tauL) { step(S, ta, T2); if (T2[1] > apex.y) apex = { t: t + ta, x: T2[0], y: T2[1], z: T2[2] }; }
      }
      if (hit !== null) {
        step(S, tauL, T1); T1[1] = hit;
        final = T1.slice(); finalT = t + tauL; status = (hit === L) ? 'landed' : 'short'; hitLevel = hit;
        rec.push(finalT, final[0], final[1], final[2], final[3], final[4], final[5]);
        break;
      }
      t += h;
      for (i = 0; i < 6; i++) S[i] = S2[i];
      rec.push(t, S[0], S[1], S[2], S[3], S[4], S[5]);
      if (S[0] < minX) minX = S[0]; if (S[0] > maxX) maxX = S[0];
      if (S[1] < minY) minY = S[1]; if (S[1] > maxY) maxY = S[1];
      if (S[2] < minZ) minZ = S[2]; if (S[2] > maxZ) maxZ = S[2];
      if (t >= MAX_TIME || steps >= MAX_STEPS) { status = 'timeout'; final = S.slice(); finalT = t; break; }
    }
    if (final[0] < minX) minX = final[0]; if (final[0] > maxX) maxX = final[0];
    if (final[2] < minZ) minZ = final[2]; if (final[2] > maxZ) maxZ = final[2];
    if (final[1] < minY) minY = final[1];
    if (apex.y > maxY) maxY = apex.y;

    // keep the sample list a sensible size (always keep the first and last point)
    var N = rec.length / 7, stride = Math.max(1, Math.ceil(N / MAX_SAMPLES)), keep = [], j;
    for (j = 0; j < N; j += stride) keep.push(j);
    if (keep[keep.length - 1] !== N - 1) keep.push(N - 1);
    var samples = new Float64Array(keep.length * 7), q;
    for (j = 0; j < keep.length; j++) for (q = 0; q < 7; q++) samples[j * 7 + q] = rec[keep[j] * 7 + q];

    var vh = Math.sqrt(final[3] * final[3] + final[5] * final[5]), speed = Math.sqrt(vh * vh + final[4] * final[4]);
    return {
      status: status, p: p, samples: samples, n: keep.length, T: finalT, dt: dtBase,
      landing: { x: final[0], y: status === 'timeout' ? final[1] : hitLevel, z: final[2], t: finalT, vx: final[3], vy: final[4], vz: final[5], speed: speed,
        angle: Math.atan2(-final[4], vh) * 180 / PI, dist: Math.sqrt(final[0] * final[0] + final[2] * final[2]) },
      apex: apex, minX: minX, maxX: maxX, minY: minY, maxY: maxY, minZ: minZ, maxZ: maxZ
    };
  }

  /* state at time t (cubic Hermite for position, linear for velocity) */
  function sampleAt(res, t) {
    var s = res.samples, n = res.n, T = s[(n - 1) * 7];
    t = clamp(isFinite(t) ? t : 0, s[0], T);
    var lo = 0, hi = n - 1, mid;
    while (hi - lo > 1) { mid = (lo + hi) >> 1; if (s[mid * 7] <= t) lo = mid; else hi = mid; }
    var b0 = lo * 7, b1 = Math.min(lo + 1, n - 1) * 7, t0 = s[b0], t1 = s[b1], h = t1 - t0, u = h > 0 ? (t - t0) / h : 0;
    var u2 = u * u, u3 = u2 * u, h00 = 2 * u3 - 3 * u2 + 1, h10 = u3 - 2 * u2 + u, h01 = -2 * u3 + 3 * u2, h11 = u3 - u2;
    function pos(k) { return h00 * s[b0 + k] + h10 * h * s[b0 + k + 3] + h01 * s[b1 + k] + h11 * h * s[b1 + k + 3]; }
    function vel(k) { return s[b0 + k] + (s[b1 + k] - s[b0 + k]) * u; }
    return { t: t, x: pos(1), y: pos(2), z: pos(3), vx: vel(4), vy: vel(5), vz: vel(6) };
  }

  /* scan launch angles for the longest range (x) under the current conditions */
  function findBestAngle(p0, opt) {
    opt = opt || {};
    var base = sanitize(p0).p, best = null;
    function rangeAt(a) {
      var r = simulate(Object.assign({}, base, { angle: a }), { coarse: !!opt.coarse });
      return r.status === 'landed' ? r.landing.x : -Infinity;
    }
    function scan(from, to, stepDeg) {
      var a, v;
      for (a = from; a <= to + 1e-9; a += stepDeg) {
        var aa = clamp(a, 0, 90); v = rangeAt(aa);
        if (best === null || v > best.dist) best = { angle: aa, dist: v };
      }
    }
    scan(0, 90, opt.coarse ? 3 : 1);
    if (best === null || best.dist === -Infinity) return null;
    if (!opt.coarse) {
      var c = best.angle; scan(Math.max(0, c - 1), Math.min(90, c + 1), 0.1);
      c = best.angle; scan(Math.max(0, c - 0.1), Math.min(90, c + 0.1), 0.01);
    }
    best.angle = Math.round(best.angle * 100) / 100;
    return best;
  }

  /* ---------- built-in tests (run on page load, and with: node physics.js) ---------- */
  function selfTest() {
    var results = [];
    function check(name, got, expected, tol) {
      var pass = isFinite(got) && Math.abs(got - expected) <= tol;
      results.push({ name: name, got: got, expected: expected, tol: tol, pass: pass });
    }
    function checkTrue(name, cond) { results.push({ name: name, got: cond ? 1 : 0, expected: 1, tol: 0, pass: !!cond }); }
    var vac = { airOn: false, g: 9.81, h0: 0, ly: 0 };
    function run(o) { return simulate(Object.assign({}, vac, o)); }

    // 1. The headline check: 45 degrees, 20 m/s, flat ground, Earth, no air
    var r = run({ v0: 20, angle: 45 });
    check('Range, 20 m/s at 45 degrees (v squared / g = 40.77 m)', r.landing.x, 400 / 9.81, 1e-6);
    check('Time of flight 2 v sin(a) / g', r.landing.t, 2 * 20 * Math.sin(PI / 4) / 9.81, 1e-6);
    check('Maximum height v squared sin squared(a) / 2g', r.apex.y, 400 * 0.5 / (2 * 9.81), 1e-6);
    check('Impact speed equals launch speed on level ground', r.landing.speed, 20, 1e-6);
    check('Impact angle equals launch angle on level ground (degrees)', r.landing.angle, 45, 1e-6);
    // 2. 30 and 60 degrees land together
    check('30 degrees and 60 degrees give the same range', run({ v0: 20, angle: 30 }).landing.x, run({ v0: 20, angle: 60 }).landing.x, 1e-6);
    // 3. Cliff: fall time depends on height only
    var c1 = run({ v0: 5, angle: 0, h0: 80 }), c2 = run({ v0: 50, angle: 0, h0: 80 });
    check('Cliff fall time sqrt(2h/g), slow throw', c1.landing.t, Math.sqrt(160 / 9.81), 1e-6);
    check('Cliff fall time sqrt(2h/g), fast throw', c2.landing.t, Math.sqrt(160 / 9.81), 1e-6);
    // 4. Best angle from a height, against the formula
    var hb = findBestAngle(Object.assign({}, vac, { v0: 20, h0: 10 }));
    check('Best angle from 10 m up (formula arctan(v / sqrt(v squared + 2gh)))', hb.angle, optimalAngleFormula(20, 9.81, 10), 0.05);
    checkTrue('Best angle from a height is below 45 degrees', hb.angle < 45);
    check('Longest range from a height (formula)', hb.dist, maxRangeFormula(20, 9.81, 10), 0.01);
    // 5. Landing higher needs more than 45 degrees
    var hu = findBestAngle(Object.assign({}, vac, { v0: 30, ly: 15 }));
    check('Best angle onto a 15 m higher level (formula)', hu.angle, optimalAngleFormula(30, 9.81, -15), 0.05);
    checkTrue('Best angle onto a higher level is above 45 degrees', hu.angle > 45);
    var rf = run({ v0: 30, angle: 60, ly: 15 }), ff = vacuumFormula(Object.assign({}, vac, { v0: 30, angle: 60, ly: 15 }));
    check('Range onto a higher level matches the formula', rf.landing.x, ff.range, 1e-6);
    // 6. Energy: landing speed from a height
    check('Landing speed from 30 m (energy: sqrt(v squared + 2 g h))', run({ v0: 25, angle: 35, h0: 30 }).landing.speed, Math.sqrt(625 + 2 * 9.81 * 30), 1e-6);
    // 7. Other planets: range scales as 1 / g
    check('Moon range / Earth range = 9.81 / 1.62', run({ v0: 20, angle: 45, g: 1.62 }).landing.x / r.landing.x, 9.81 / 1.62, 1e-6);
    // 8. Edge cases
    var e90 = run({ v0: 20, angle: 90 }); check('Straight up (90 degrees): range 0', e90.landing.x, 0, 1e-9); check('Straight up: height v squared / 2g', e90.apex.y, 400 / (2 * 9.81), 1e-6);
    check('Zero degrees from the ground: lands at once', run({ v0: 20, angle: 0 }).landing.t, 0, 1e-6);
    check('Zero speed from 10 m: drops straight down in sqrt(2h/g)', run({ v0: 0, h0: 10 }).landing.t, Math.sqrt(20 / 9.81), 1e-6);
    check('Huge speed (200 m/s at 45 degrees) range', run({ v0: 200, angle: 45 }).landing.x, 40000 / 9.81, 1e-4);
    var sh = run({ v0: 5, angle: 45, ly: 100 }); checkTrue('Too slow to reach a 100 m higher level: reported as short', sh.status === 'short');
    // 9. Air: terminal speed of a football dropped from 500 m: sqrt(2 m g / (rho Cd A))
    var air = { airOn: true, mass: 0.43, radius: 0.11, cd: 0.25, rho: 1.225, g: 9.81, h0: 500, ly: 0, v0: 0, angle: 0 };
    var vt = Math.sqrt(2 * 0.43 * 9.81 / (1.225 * 0.25 * PI * 0.11 * 0.11));
    check('Terminal speed with drag (m/s)', simulate(air).landing.speed, vt, vt * 1e-3);
    // 10. Air: range shorter, steeper fall, best angle below 45 degrees
    var a0 = { airOn: true, v0: 40, angle: 35, mass: 0.16, radius: 0.036, cd: 0.4, rho: 1.225, g: 9.81 };
    var ra = simulate(a0), rv = simulate(Object.assign({}, a0, { airOn: false }));
    checkTrue('Drag shortens the range', ra.landing.x < rv.landing.x);
    checkTrue('With drag the fall is steeper than the launch', ra.landing.angle > 35);
    checkTrue('With drag the impact speed is lower than the launch speed', ra.landing.speed < 40);
    checkTrue('With drag the best angle is below 45 degrees', findBestAngle(Object.assign({}, a0, { v0: 50 }), {}).angle < 45);
    // 11. Wind: tailwind helps, headwind hurts, crosswind pushes sideways
    var wbase = Object.assign({}, a0, { windSpeed: 10 });
    checkTrue('Tailwind increases range', simulate(Object.assign({}, wbase, { windDir: 0 })).landing.x > ra.landing.x);
    checkTrue('Headwind decreases range', simulate(Object.assign({}, wbase, { windDir: 180 })).landing.x < ra.landing.x);
    checkTrue('Crosswind (90 degrees) pushes towards +z', simulate(Object.assign({}, wbase, { windDir: 90 })).landing.z > 0.1);
    // 12. Spin: backspin lifts, topspin pushes down, sidespin curves to +z
    var sb = Object.assign({}, a0, { v0: 60, angle: 12, mass: 0.0459, radius: 0.02135, cd: 0.25 });
    checkTrue('Backspin gives more range than no spin (golf shot)', simulate(Object.assign({}, sb, { backspin: 3000 })).landing.x > simulate(sb).landing.x);
    checkTrue('Topspin gives less range than no spin', simulate(Object.assign({}, sb, { backspin: -3000 })).landing.x < simulate(sb).landing.x);
    var f1 = forcesAt({ x: 0, y: 5, z: 0, vx: 40, vy: 0, vz: 0 }, Object.assign({}, sb, { backspin: 3000 }));
    var f2 = forcesAt({ x: 0, y: 5, z: 0, vx: 40, vy: 0, vz: 0 }, Object.assign({}, sb, { backspin: -3000 }));
    checkTrue('Backspin Magnus force points up', f1.magnus[1] > 0); checkTrue('Topspin Magnus force points down', f2.magnus[1] < 0);
    checkTrue('Sidespin curves the ball towards +z', simulate(Object.assign({}, sb, { sidespin: 3000 })).landing.z > 0.1);
    // 13. Thin air at height helps range
    var hi = Object.assign({}, a0, { v0: 80, angle: 40, h0: 3000 });
    checkTrue('Thinner air with height gives more range', simulate(Object.assign({}, hi, { altDensity: true })).landing.x > simulate(hi).landing.x);
    // 14. Bad input is repaired
    var bad = sanitize({ v0: -5, angle: 'x', mass: 0, radius: -1, g: 0 }).p;
    checkTrue('Bad inputs are repaired (no negatives, no zero mass or gravity)', bad.v0 === 0 && bad.angle === 45 && bad.mass > 0 && bad.radius > 0 && bad.g > 0);
    var allOk = results.every(function (x) { return x.pass; });
    return { ok: allOk, results: results };
  }

  var api = { sanitize: sanitize, simulate: simulate, sampleAt: sampleAt, vacuumFormula: vacuumFormula, optimalAngleFormula: optimalAngleFormula,
    maxRangeFormula: maxRangeFormula, findBestAngle: findBestAngle, forcesAt: forcesAt, selfTest: selfTest, LIMITS: LIM };
  root.NavPhysics = api;
  if (typeof module !== 'undefined' && module.exports) {
    module.exports = api;
    if (typeof require !== 'undefined' && require.main === module) {
      var t0 = Date.now(), st = selfTest();
      st.results.forEach(function (x) { console.log((x.pass ? 'PASS ' : 'FAIL ') + x.name + '  got=' + x.got + ' expected=' + x.expected); });
      console.log(st.ok ? '\nALL ' + st.results.length + ' CHECKS PASSED (' + (Date.now() - t0) + ' ms)' : '\nSOME CHECKS FAILED');
      process.exit(st.ok ? 0 : 1);
    }
  }
})(typeof window !== 'undefined' ? window : globalThis);
