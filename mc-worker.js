// mc-worker.js
"use strict";

// ====== Constantes y utilidades compartidas ======
const ALIVE = 1, DEAD = 0;
const threshold = 1e-4, chance = 0.1;

// StepRand = 0:0.001:1
const StepRand = Array.from({ length: 1001 }, (_, k) => k / 1000);
function randFromStep() {
  const idx = Math.floor(Math.random() * StepRand.length);
  return StepRand[idx];
}
function randi10Over10() {
  const k = Math.floor(Math.random() * 10) + 1;
  return k / 10;
}

// --------- Geometría en cm ---------
function insideDomainInclusive(x, y, z, R, z_top) {
  return Math.abs(x) <= R && Math.abs(y) <= R && z >= 0 && z <= z_top;
}
function pickLayerByZInclusive(z, z_bound) {
  for (let i = 0; i < z_bound.length - 1; i++) {
    if (z_bound[i] <= z && z < z_bound[i + 1]) return i;
  }
  return z_bound.length - 2; // seguridad
}
function launchDirectionIsotropic() {
  const ct = 2 * Math.random() - 1;
  const st = Math.sqrt(1 - ct * ct);
  const psi = 2 * Math.PI * Math.random();
  return [st * Math.cos(psi), st * Math.sin(psi), ct];
}
function sampleHenyeyGreenstein(g) {
  const R = Math.random();
  if (g === 0) return 2 * R - 1;
  const t = (1 - g * g) / (1 - g + 2 * g * R);
  let ct = (1 + g * g - t * t) / (2 * g);
  if (ct > 1) ct = 1;
  if (ct < -1) ct = -1;
  return ct;
}

// --------- Esferas (overlay) ---------
function radialDistance(x, y, z, c) {
  const dx = x - c[0], dy = y - c[1], dz = z - c[2];
  return Math.hypot(dx, dy, dz);
}
// Devuelve el índice global de la primera esfera cuyo radio >= r (anidadas)
function pickSphereOverride(r, spheres) {
  for (let k = 0; k < spheres.length; k++) {
    if (r <= spheres[k].radius) return spheres[k].index;
  }
  return null;
}

// --------- Núcleo de simulación (VERSIÓN PURA, sin DOM) ---------
function simulate(N, mua, mus, g, albedo, z_bound_slab, Rxy, overlayCfg, laser) {
  const Nc = mua.length; // incluye todas las filas (slabs + esferas)
  const energy = new Float64Array(Nc);

  const slabIdx = overlayCfg?.slabIdx || [];
  const spheres = overlayCfg?.spheres || []; // [{index, radius, tipo}] ascendente
  const center  = overlayCfg?.center  || [0, 0, 0.15]; // (cm)
  const z_top   = z_bound_slab[z_bound_slab.length - 1];

  // LÁSER desde argumentos
  const laserX = Number(laser?.x) || 0;
  const laserY = Number(laser?.y) || 0;
  const laserZ = Number(laser?.z) || 0;

  for (let p = 0; p < N; p++) {
    let W = 1.0;
    let x = laserX;
    let y = laserY;
    let z = laserZ;
    let [ux, uy, uz] = launchDirectionIsotropic();

    // capa inicial por z (entre las slabs)
    let i = slabIdx[pickLayerByZInclusive(0, z_bound_slab)];

    while (true) {
      const rU = randFromStep();
      const s = (rU === 0) ? 1e9 : -Math.log(rU) / (mua[i] + mus[i]);

      x += s * ux;
      y += s * uy;
      z += s * uz;

      if (!insideDomainInclusive(x, y, z, Rxy, z_top)) { W = 0; break; }

      // 1) slab por z
      i = slabIdx[pickLayerByZInclusive(z, z_bound_slab)];

      // 2) overlay esférico (si cae dentro de alguna esfera, sobrescribe i)
      if (spheres.length) {
        const r = radialDistance(x, y, z, center);
        const iSph = pickSphereOverride(r, spheres);
        if (iSph !== null) i = iSph;
      }

      // absorción
      const absorb = W * (1 - albedo[i]);
      W -= absorb;
      energy[i] += absorb;

      // nueva dirección (HG)
      const ct = sampleHenyeyGreenstein(g[i]);
      const st = Math.sqrt(1 - ct * ct);
      const psi = 2 * Math.PI * Math.random();
      const cp = Math.cos(psi);
      const sp = (psi < Math.PI ? 1 : -1) * Math.sqrt(1 - cp * cp);

      let uxx, uyy, uzz;
      if (1 - Math.abs(uz) <= 1e-12) {
        uxx = st * cp;
        uyy = st * sp;
        uzz = ct * (uz >= 0 ? 1 : -1);
      } else {
        const t = Math.sqrt(1 - uz * uz);
        uxx = (st * ((ux * uz * cp) - (uy * sp)) / t) + (ux * ct);
        uyy = (st * ((uy * uz * cp) + (ux * sp)) / t) + (uy * ct);
        uzz = (-st * cp * t) + (uz * ct);
      }
      ux = uxx; uy = uyy; uz = uzz;

      // ruleta rusa
      if (W < threshold) {
        if (randi10Over10() <= chance) W /= chance;
        else break;
      }
    }
  }
  return energy;
}

// --------- Mensajería del worker ---------
self.onmessage = (ev) => {
  const {
    N, iters,
    mua, mus, g, albedo,
    z_bound_slab, Rxy,
    overlayCfg, laser
  } = ev.data;

  const Nc = mua.length;
  const energySum = new Float64Array(Nc);
  const losses = new Float64Array(iters); // (N - energyTit)/N por iter

  for (let j = 0; j < iters; j++) {
    const sim = simulate(N, mua, mus, g, albedo, z_bound_slab, Rxy, overlayCfg, laser);

    let energyTit = 0.0;
    for (let i = 0; i < Nc; i++) {
      energySum[i] += sim[i];
      energyTit += sim[i];
    }
    losses[j] = (N - energyTit) / N;

    // 🔔 Progreso: una simulación (iteración) terminada en este worker
    // (mandamos un mensajito pequeño sin transferibles)
    self.postMessage({ type: 'progress', inc: 1 });
  }

  // ✅ Resultado final del worker
  self.postMessage(
    { type: 'result', energySum, losses },
    [energySum.buffer, losses.buffer]
  );
};