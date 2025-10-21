/**************************************
 * Monte Carlo con slabs + esferas (cm)
 * NOTA: Todos los valores de longitud
 *       (Rxy, z_end, radios) están en
 *       CENTÍMETROS (cm).
 **************************************/

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
function isSphericalRow(t) {
  return /esfera|sphere|tumor/i.test(String(t || ""));
}
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

// --------- Núcleo de simulación ---------
function simulate(N, mua, mus, g, albedo, z_bound_slab, Rxy, overlayCfg) {
  const Nc = mua.length; // incluye todas las filas de la tabla (slabs + esferas)
  const energy = new Float64Array(Nc);

  const slabIdx = overlayCfg?.slabIdx || [];
  const spheres = overlayCfg?.spheres || []; // [{index, radius}] ascendente
  const center  = overlayCfg?.center  || [0, 0, 0.15]; // (cm) por defecto
  const z_top   = z_bound_slab[z_bound_slab.length - 1];

  for (let p = 0; p < N; p++) {
    let W = 1.0;
    let x = 0, y = 0, z = 0;
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

// ================== DIBUJO A ESCALA (cm) ==================
function parseGeometryFromTable(tipo, z_end, centerOverrideCm) {
  const isSphere = t => /esfera|sphere|tumor/i.test(String(t || ""));
  const slabIdx = [], spheres = [];

  for (let i = 0; i < tipo.length; i++) {
    if (isSphere(tipo[i])) spheres.push({ index: i, radius: Number(z_end[i]) || 0, tipo: String(tipo[i]) });
    else slabIdx.push(i);
  }

  // z-bounds de slabs (acumulados)
  const zBound = [0];
  for (const i of slabIdx) zBound.push(Number(z_end[i]) || 0);
  const zMax = zBound[zBound.length - 1];

  // radios ordenados (anidadas)
  spheres.sort((a, b) => a.radius - b.radius);

  // Centro del tumor: lo define el usuario (cm). Si no viene, 0.15 cm
  const tumorCenterZ = Number.isFinite(centerOverrideCm) ? centerOverrideCm : 0.15;

  return { slabIdx, zBound, zMax, spheres, tumorCenterZ };
}

function getLayerColor(tipo) {
  const t = String(tipo || "").toLowerCase();

  if (/tumor/.test(t))    return "#B71C1C";  // rojo intenso (solo para contorno de la esfera)
  if (/musculo|músculo/.test(t)) return "#E57373"; // músculo
  if (/grasa|adiposo/.test(t))   return "#FFF2B0"; // grasa
  if (/dermis|piel/.test(t))     return "#F8CFCF"; // piel/dermis
  if (/tejido/.test(t))          return "#F8E0C0"; // tejido genérico

  // por defecto (slab genérico)
  return "#F6E6D8";
}

// color de las esferas (contornos)
function getSphereStroke(tipo) {
  return /tumor/i.test(String(tipo)) ? "#B71C1C" : "#7FC5D6";
}



// Estado global (por canvas) para la animación
window.__corteAnim = window.__corteAnim || new Map();

// ======= Estado global por-canvas para la animación =======
window.__corteAnim = window.__corteAnim || new Map();

// ======= Estado global por-canvas para la animación =======
window.__corteAnim = window.__corteAnim || new Map();

function dibujarCorteEscala({ Rxy, zBound, zMax, spheres, tumorCenterZ, slabIdx, tipo, energy }) {
  const canvas = document.getElementById("canvas3D");
  if (!canvas) return;
  const ctx = canvas.getContext("2d");

  // ---------- helpers internos ----------
  const getLayerColor = (t) => {
    const s = String(t || "").toLowerCase();
    if (/musculo|músculo/.test(s)) return "#E57373"; // músculo
    if (/grasa|adiposo/.test(s))   return "#FFF2B0"; // grasa
    if (/dermis|piel/.test(s))     return "#F8CFCF"; // piel/dermis
    if (/tejido/.test(s))          return "#F8E0C0"; // genérico
    return "#F6E6D8";                                  // por defecto
  };
  const getSphereStroke = (t) => (/tumor/i.test(String(t)) ? "#B71C1C" : "#7FC5D6");

  // ---------- HiDPI / tamaño real de canvas ----------
  const DPR = Math.max(1, window.devicePixelRatio || 1);
  const cssW = canvas.clientWidth  || canvas.width;
  const cssH = canvas.clientHeight || canvas.height;
  if (canvas.width  !== Math.round(cssW * DPR)) canvas.width  = Math.round(cssW * DPR);
  if (canvas.height !== Math.round(cssH * DPR)) canvas.height = Math.round(cssH * DPR);
  ctx.setTransform(DPR, 0, 0, DPR, 0, 0);

  const W = cssW, H = cssH;
  ctx.clearRect(0, 0, W, H);

  // ---------- medir leyenda para margen derecho dinámico ----------
  ctx.save();
  ctx.font = "12px system-ui, Arial";
  const legendItems = spheres.slice().sort((a,b)=>a.radius-b.radius).map(s=>{
    const isTumor = /tumor/i.test(s.tipo);
    const E = Number(energy?.[s.index]) || 0;
    return `${isTumor?"Tumor":"Esfera"} • E=${E.toExponential(3)} • r=${(Number(s.radius)||0).toFixed(4)} cm`;
  });
  let maxLegendW = ctx.measureText("Capas esféricas:").width;
  for (const t of legendItems) maxLegendW = Math.max(maxLegendW, ctx.measureText(t).width);
  ctx.restore();

  // ---------- márgenes ----------
  const mL = 100, mT = 60, mB = 76;
  const mR = Math.max(240, 40 + maxLegendW); // aire + texto
  const plotW = W - mL - mR, plotH = H - mT - mB;

  // ---------- escalas ----------
  const widthCm  = 2 * Rxy;
  const heightCm = zMax || 1;
  const sY = plotH / heightCm;                   // escala física vertical (mantiene círculos)
  const widthPxPhysical = widthCm * sY;
  const widthPxForRect  = Math.max(widthPxPhysical, 360); // ensancha visual si Rxy es pequeño

  const x0 = mL + (plotW - widthPxForRect) / 2;
  const y0 = mT + (plotH - heightCm * sY) / 2;

  const scaleX = widthPxForRect / widthCm;       // escala horizontal (visual)
  const X = (cm)=> x0 + (cm + Rxy) * scaleX;     // x en [-Rxy,+Rxy]
  const Y = (cm)=> y0 + (cm) * sY;               // z en [0,zMax]

  // ---------- fondo / borde del dominio ----------
  ctx.strokeStyle = "#6b6b6b";
  ctx.lineWidth = 1;
  ctx.strokeRect(X(-Rxy), Y(0), widthPxForRect, heightCm * sY);

  // ---------- slabs coloreados ----------
  for (let k = 0; k < slabIdx.length; k++) {
    const idx = slabIdx[k];
    const y1 = Y(zBound[k]);
    const y2 = Y(zBound[k + 1]);
    const h  = y2 - y1;

    ctx.fillStyle = getLayerColor(tipo[idx]);
    ctx.globalAlpha = 0.85;
    ctx.fillRect(X(-Rxy), y1, widthPxForRect, h);
    ctx.globalAlpha = 1;

    // divisoria sutil
    if (k > 0) {
      ctx.strokeStyle = "#d3cfc9";
      ctx.beginPath(); ctx.moveTo(X(-Rxy), y1); ctx.lineTo(X(+Rxy), y1); ctx.stroke();
    }
  }
  // divisoria inferior
  ctx.strokeStyle = "#d3cfc9";
  ctx.beginPath(); ctx.moveTo(X(-Rxy), Y(zMax)); ctx.lineTo(X(+Rxy), Y(zMax)); ctx.stroke();

  // ---------- esferas contorno del modelo (no animadas) ----------
  const cx = X(0), cz = Y(tumorCenterZ);
  spheres.forEach((sph) => {
    const rpx = (Number(sph.radius) || 0) * sY;
    if (rpx <= 0) return;
    ctx.beginPath();
    ctx.strokeStyle = getSphereStroke(sph.tipo);
    ctx.lineWidth = /tumor/i.test(sph.tipo) ? 2 : 1.5;
    ctx.arc(cx, cz, rpx, 0, 2 * Math.PI);
    ctx.stroke();
  });

  // ---------- LASER ----------
  const lambda = Number(document.getElementById("lambda").value);
  let color_laser = 0x000;
  if (lambda <= 430) {
    color_laser = 0x8000FF;
  } else if (lambda > 430 && lambda <= 480) {
    color_laser = 0x0037FF;
  } else if (lambda > 480 && lambda <= 560) {
    color_laser = 0x00DB21;
  } else if (lambda > 560 && lambda <= 590) {
    color_laser = 0xDBD400;
  } else if (lambda > 590 && lambda <= 620) {
    color_laser = 0xDB9900;
  } else if (lambda > 620) {
    color_laser = 0xDB0000;
  }
  const color = `#${color_laser.toString(16).padStart(6, "0")}`;

  ctx.strokeStyle = color; 
  ctx.fillStyle = color;
  ctx.lineWidth = 6;

  ctx.beginPath();
  ctx.moveTo(X(0), Y(0) - 25);
  ctx.lineTo(X(0), Y(zMax) + 25);
  ctx.stroke();

  ctx.beginPath();
  ctx.lineWidth = 5;
  ctx.moveTo(X(0), Y(0) - 25);
  ctx.lineTo(X(0)-20, Y(0)-5);
  ctx.moveTo(X(0), Y(0) - 25);
  ctx.lineTo(X(0)+20, Y(0)-5);
  ctx.stroke();

  ctx.font = "17px Arial";
  ctx.fillStyle = color;
  ctx.fillText("Laser origin", X(0)-45, Y(0)-30);

  // ---------- etiquetas en slabs ----------
  ctx.font = "12px system-ui, Arial";
  for (let k = 0; k < slabIdx.length; k++) {
    const idx = slabIdx[k];
    const z1 = zBound[k], z2 = zBound[k+1];
    const zc = (z1 + z2) / 2;
    const h  = z2 - z1;
    const E  = Number(energy?.[idx]) || 0;

    const label = `${tipo[idx]} • E=${E.toExponential(3)} • h=${h.toFixed(3)} cm`;
    const padX = 6;
    const textW = ctx.measureText(label).width;
    const boxW = Math.min(textW + padX*2, widthPxForRect - 12);
    const boxH = 18;

    const bx = X(-Rxy) + (widthPxForRect - boxW)/2;
    const by = Y(zc) - boxH/2;

    ctx.fillStyle = "rgba(255,255,255,0.92)";
    ctx.strokeStyle = "rgba(0,0,0,0.12)";
    ctx.lineWidth = 1;
    ctx.beginPath(); ctx.rect(bx, by, boxW, boxH); ctx.fill(); ctx.stroke();

    ctx.fillStyle = "#222";
    ctx.fillText(label, bx + padX, by + boxH - 5);
  }

  // ---------- leyenda ----------
  const legendMaxW = mR - 24;
  const wrap = (text) => {
    const words = text.split(" ");
    let line = "", out = [];
    for (const w of words) {
      const test = line ? line + " " + w : w;
      if (ctx.measureText(test).width <= legendMaxW - 14) line = test; else { out.push(line); line = w; }
    }
    if (line) out.push(line);
    return out;
  };

  const rectRight = X(+Rxy);
  const legendX = Math.max(rectRight + 20, W - mR + 12);
  let legendY = y0 + 6;

  var linesPreview = ["Capas esféricas:"].concat(legendItems.flatMap(wrap));

  if (legendItems.length > 0) {
    const panelH = 24 + linesPreview.length * 18;
    const panelW = Math.min(mR - 16, 380);
    ctx.fillStyle = "rgba(255,255,255,0.96)";
    ctx.strokeStyle = "rgba(0,0,0,0.12)";
    ctx.lineWidth = 1;
    ctx.beginPath(); ctx.rect(legendX - 10, legendY - 18, panelW + 20, panelH); ctx.fill(); ctx.stroke();
    ctx.fillStyle = "#111";
    ctx.font = "16px system-ui, Arial";
    ctx.fillText("Capas esféricas:", legendX, legendY);
    legendY += 22;
  }

  ctx.font = "12px system-ui, Arial";
  spheres.slice().sort((a,b)=>a.radius-b.radius).forEach(sph => {
    const idx = sph.index;
    const isTumor = /tumor/i.test(sph.tipo);
    const E = Number(energy?.[idx]) || 0;
    const text = `${isTumor?"Tumor":"Esfera"} • E=${E.toExponential(3)} • r=${(Number(sph.radius)||0).toFixed(4)} cm`;

    ctx.fillStyle = getSphereStroke(sph.tipo);
    ctx.fillRect(legendX, legendY - 9, 10, 10);

    ctx.fillStyle = "#222";
    const lines = wrap(text);
    ctx.fillText(lines[0], legendX + 14, legendY);
    for (let i = 1; i < lines.length; i++) {
      legendY += 16;
      ctx.fillText(lines[i], legendX + 14, legendY);
    }
    legendY += 18;
  });

  // ---------- cotas globales ----------
  ctx.strokeStyle = "#6b6b6b"; ctx.fillStyle = "#444"; ctx.lineWidth = 1;

  // Alto (izquierda)
  ctx.beginPath(); ctx.moveTo(X(-Rxy) - 40, Y(0)); ctx.lineTo(X(-Rxy) - 40, Y(zMax)); ctx.stroke();
  ctx.beginPath();
  ctx.moveTo(X(-Rxy)-40, Y(0));   ctx.lineTo(X(-Rxy)-35, Y(0)+7);   ctx.lineTo(X(-Rxy)-45, Y(0)+7);   ctx.closePath(); ctx.fill();
  ctx.beginPath();
  ctx.moveTo(X(-Rxy)-40, Y(zMax)); ctx.lineTo(X(-Rxy)-35, Y(zMax)-7); ctx.lineTo(X(-Rxy)-45, Y(zMax)-7); ctx.closePath(); ctx.fill();
  ctx.fillText(`${heightCm.toFixed(2)} cm`, X(-Rxy)-70, Y(0)+(heightCm*sY)/2 + 4);

  // Ancho (abajo)
  ctx.beginPath(); ctx.moveTo(X(-Rxy), Y(zMax) + 28); ctx.lineTo(X(+Rxy), Y(zMax) + 28); ctx.stroke();
  ctx.beginPath();
  ctx.moveTo(X(-Rxy), Y(zMax)+28); ctx.lineTo(X(-Rxy)+7, Y(zMax)+23); ctx.lineTo(X(-Rxy)+7, Y(zMax)+33); ctx.closePath(); ctx.fill();
  ctx.beginPath();
  ctx.moveTo(X(+Rxy), Y(zMax)+28); ctx.lineTo(X(+Rxy)-7, Y(zMax)+23); ctx.lineTo(X(+Rxy)-7, Y(zMax)+33); ctx.closePath(); ctx.fill();
  ctx.fillText(`${widthCm.toFixed(2)} cm`, X(-Rxy)+(widthPxForRect)/2 - 18, Y(zMax) + 46);


  // ===================== ANIMACIÓN DE ESFERAS (overlay) =====================

  // ----- utilidades -----
  const rand = (min, max) => Math.random() * (max - min) + min;
  const hexToRgba = (hex, a=1) => {
    let h = hex.replace("#", "");
    if (h.length === 3) h = h.split("").map(c=>c+c).join("");
    const r = parseInt(h.slice(0,2), 16);
    const g = parseInt(h.slice(2,4), 16);
    const b = parseInt(h.slice(4,6), 16);
    return `rgba(${r},${g},${b},${a})`;
  };

  // CDF por energía
  function buildEnergyCDF(energy, slabIdx) {
    const weights = slabIdx.map((idx) => Math.max(0, Number(energy?.[idx]) || 0));
    const total = weights.reduce((a, b) => a + b, 0);
    const probs = (total > 0)
      ? weights.map(w => w / total)
      : weights.map(() => 1 / weights.length);
    let acc = 0;
    const cdf = probs.map(p => (acc += p));
    cdf[cdf.length - 1] = 1; // cierre numérico
    return cdf;
  }
  function pickLayerByCDF(cdf) {
    const u = Math.random();
    let lo = 0, hi = cdf.length - 1;
    while (lo < hi) {
      const mid = (lo + hi) >> 1;
      if (u <= cdf[mid]) hi = mid; else lo = mid + 1;
    }
    return lo;
  }

  // ---- estado por-canvas ----
  let state = window.__corteAnim.get(canvas);
  if (!state) {
    state = {
      particles: [],
      last: performance.now(),
      running: true,
      spawnAcc: 0,
      MAX: 120,     // máximo de esferas simultáneas
      RATE: 40,     // esferas por segundo aprox.
      advanceOnNextDraw: false,
      _started: false
    };
    window.__corteAnim.set(canvas, state);

    // Pausa/Reanuda con 'p'
    window.addEventListener("keydown", (e) => {
      if (e.key.toLowerCase() === "p") {
        state.running = !state.running;
        if (state.running) {
          state.last = performance.now();
          requestAnimationFrame(tick);
        }
      }
    }, { passive: true });
  }

  // Guardar referencias para el loop y actualizar CDF
  state.lastArgs   = { Rxy, zBound, zMax, spheres, tumorCenterZ, slabIdx, tipo, energy };
  state.energyCDF  = buildEnergyCDF(energy, slabIdx);
  state.zBoundRef  = zBound;
  state.RxyRef     = Rxy;
  state.zMaxRef    = zMax;
  state.colorLaser = color;  // para usar mismo color
  state.sY         = sY;     // para diámetro en cm si se desea

  // (opcional) escalar tasa con energía máxima
  {
    const maxE = Math.max(...slabIdx.map(i => Number(energy?.[i]) || 0), 0);
    const kRate = 1 + Math.log10(1 + maxE); // factor suave
    state.RATE = 40 * kRate;
  }

  // ---- SPWAN con sesgo por energía ----
  function spawnParticle() {
  // ===== DIÁMETRO FIJO EN CM =====
  const DIAMETER_CM = 0.01;                    // <-- pon aquí tu valor fijo
  let r = (DIAMETER_CM / 2) * state.sY;        // convierte a píxeles con la escala actual

  // (opcional) Asegurar visibilidad mínima y evitar burbujas gigantes:
  r = Math.max(0.8, Math.min(1.2, r));            // clamp en píxeles (min 3px, máx 28px)

  // ... el resto de tu spawn (elección de capa por CDF, x/y, velocidades, vida, etc.)
  const k = state.energyCDF?.length ? pickLayerByCDF(state.energyCDF) : 0;
  const z1 = state.zBoundRef[k];
  const z2 = state.zBoundRef[k + 1];
  const left = X(-state.RxyRef), right = X(+state.RxyRef);
  const xMin = left + r,  xMax = right - r;
  const yMin = Y(z1) + r, yMax = Y(z2) - r;

  state.particles.push({
      x: Math.random() * Math.max(0, (xMax - xMin)) + xMin,
      y: Math.random() * Math.max(0, (yMax - yMin)) + yMin,
      r,
      t: 0,
      life: rand(0.8, 2.2),
      fadeIn: rand(0.2, 0.5),
      fadeOut: rand(0.2, 0.5),
      vx: rand(-20, 20),
      vy: rand(-20, 20)
    });
  }


  // Avanzar física si corresponde
  if (state.advanceOnNextDraw && state.running) {
    const now = performance.now();
    const dt = Math.min(0.05, (now - state.last) / 1000);
    state.last = now;

    // Spawning controlado
    state.spawnAcc += state.RATE * dt;
    while (state.spawnAcc >= 1 && state.particles.length < state.MAX) {
      spawnParticle();
      state.spawnAcc -= 1;
    }

    // Actualizar partículas
    for (let i = state.particles.length - 1; i >= 0; i--) {
      const p = state.particles[i];
      p.t += dt;
      p.x += p.vx * dt;
      p.y += p.vy * dt;

      // Rebote en límites del dominio
      const left = X(-state.RxyRef), right = X(+state.RxyRef);
      const top  = Y(0), bottom = Y(state.zMaxRef);
      if (p.x < left + p.r)   { p.x = left + p.r;   p.vx *= -1; }
      if (p.x > right - p.r)  { p.x = right - p.r;  p.vx *= -1; }
      if (p.y < top + p.r)    { p.y = top + p.r;    p.vy *= -1; }
      if (p.y > bottom - p.r) { p.y = bottom - p.r; p.vy *= -1; }

      if (p.t >= p.life) state.particles.splice(i, 1);
    }

    state.advanceOnNextDraw = false;
  }

  // Dibujo de partículas (mismo color que el láser, con opacidad animada)
  for (const p of state.particles) {
    let a = 1;
    if (p.t < p.fadeIn) {
      const u = p.t / p.fadeIn; a = 1 - Math.pow(1 - u, 3);   // ease-out
    } else if (p.life - p.t < p.fadeOut) {
      const u = 1 - (p.life - p.t) / p.fadeOut; a = 1 - Math.pow(u, 3); // ease-in
    }
    ctx.beginPath();
    ctx.arc(p.x, p.y, p.r, 0, Math.PI * 2);
    // Si quisieras color por partícula: hexToRgba(p.color || state.colorLaser, a)
    ctx.fillStyle = hexToRgba(state.colorLaser, a);
    ctx.fill();
  }

  // Loop de animación
  function tick() {
    if (!state.running) return;
    state.advanceOnNextDraw = true;
    dibujarCorteEscala(state.lastArgs); // re-dibuja todo y hace el step de física
    requestAnimationFrame(tick);
  }
  if (!state._started) {
    state._started = true;
    state.last = performance.now();
    requestAnimationFrame(tick);
  }
}

// ----- Render auxiliar (sin cambios relevantes) -----
function renderEsquemaDesdeTabla(Rxy, tipo, z_end, centerCm, energyArr) {
  const geom = parseGeometryFromTable(tipo, z_end, centerCm);
  if (!geom.zBound.length) return;
  dibujarCorteEscala({
    Rxy,
    zBound: geom.zBound,
    zMax: geom.zMax,
    spheres: geom.spheres,
    tumorCenterZ: geom.tumorCenterZ,
    slabIdx: geom.slabIdx,
    tipo,
    energy: energyArr || new Float64Array(tipo.length) // por si viene vacío
  });
}


// (tu renderEsquemaDesdeTabla permanece igual)
function renderEsquemaDesdeTabla(Rxy, tipo, z_end, centerCm, energyArr) {
  const geom = parseGeometryFromTable(tipo, z_end, centerCm);
  if (!geom.zBound.length) return;
  dibujarCorteEscala({
    Rxy,
    zBound: geom.zBound,
    zMax: geom.zMax,
    spheres: geom.spheres,
    tumorCenterZ: geom.tumorCenterZ,
    slabIdx: geom.slabIdx,
    tipo,
    energy: energyArr || new Float64Array(tipo.length)
  });
}


// --------- UI ---------
document.addEventListener("DOMContentLoaded", function () {
  const NphotonsEl = document.getElementById("Nphotons");
  const RxyEl = document.getElementById("Rxy");
  const NCEl = document.getElementById("NC");
  const initBtn = document.getElementById("init");
  const runBtn = document.getElementById("run");
  const tabla = document.getElementById("tabla");
  const tbody = tabla.querySelector("tbody");
  const salida = document.getElementById("salida");
  const theightEl = document.getElementById("theight");

  initBtn.addEventListener("click", function () {
    const Nc = Math.max(1, parseInt(NCEl.value, 10) || 1);
    tbody.innerHTML = "";
    for (let i = 0; i < Nc; i++) {
      const tr = document.createElement("tr");
      tr.innerHTML = `
        <td>${i + 1}</td>
        <td><input type="text" step="any" placeholder="tipo (slab/sphere/tumor)" required></td>
        <td><input type="number" step="any" placeholder="mua" required></td>
        <td><input type="number" step="any" placeholder="mus" required></td>
        <td><input type="number" step="any" placeholder="g"   required></td>
        <td><input type="number" step="any" placeholder="nt"  required></td>
        <td><input type="number" step="any" placeholder="z_end (cm) o radio (cm)" required></td>
      `;
      tbody.appendChild(tr);
    }
  });

  runBtn.addEventListener("click", function () {
    const N = Math.max(1, parseInt(NphotonsEl.value, 10) || 1000000);

    // Rxy en cm (ignora el "(m)" del label)
    const Rxy = parseFloat(RxyEl.value);
    if (!isFinite(Rxy) || Rxy <= 0) {
      salida.textContent = "Error: Rxy debe ser un número en centímetros.";
      return;
    }

    const rows = Array.from(tbody.querySelectorAll("tr"));
    if (rows.length === 0) {
      salida.textContent = "Primero crea la tabla de capas.";
      return;
    }

    const Nc = rows.length;
    const tipo = new Array(Nc);
    const mua = new Float64Array(Nc);
    const mus = new Float64Array(Nc);
    const g   = new Float64Array(Nc);
    const nt  = new Float64Array(Nc);
    const z_end = new Float64Array(Nc); // slab: límite z acumulado (cm); esfera/tumor: radio (cm)

    for (let i = 0; i < Nc; i++) {
      const inputs = rows[i].querySelectorAll("input");
      tipo[i]  = inputs[0].value;
      mua[i]   = parseFloat(inputs[1].value);
      mus[i]   = parseFloat(inputs[2].value);
      g[i]     = parseFloat(inputs[3].value);
      nt[i]    = parseFloat(inputs[4].value);
      z_end[i] = parseFloat(inputs[5].value);

      if (![mua[i], mus[i], g[i], nt[i], z_end[i]].every(Number.isFinite)) {
        salida.textContent = `Error: revisa los valores numéricos en la fila ${i + 1}.`;
        return;
      }
    }

    // Particionamos en slabs (por z) y esferas (radio)
    const slabIdx = [];
    const spheres = []; // {index, radius, tipo}
    for (let i = 0; i < Nc; i++) {
      if (isSphericalRow(tipo[i])) {
        spheres.push({ index: i, radius: z_end[i], tipo: tipo[i] });
      } else {
        slabIdx.push(i);
      }
    }
    if (slabIdx.length === 0) {
      salida.textContent = "Error: define al menos una capa 'slab' (tipo sin 'esfera'/'tumor').";
      return;
    }

    // z_end estrictamente creciente solo entre slabs
    for (let k = 1; k < slabIdx.length; k++) {
      const a = z_end[slabIdx[k - 1]], b = z_end[slabIdx[k]];
      if (!(b > a)) {
        salida.textContent = "Error: z_end (cm) debe ser estrictamente creciente para las capas slab.";
        return;
      }
    }

    // Límite por z (cm) usando solo slabs
    const z_bound_slab = new Float64Array(slabIdx.length + 1);
    z_bound_slab[0] = 0.0;
    for (let k = 0; k < slabIdx.length; k++) {
      z_bound_slab[k + 1] = z_end[slabIdx[k]];
    }

    // Centro del tumor (cm) definido por el usuario
    let tumorCenterZ = Number(theightEl?.value);
    if (!Number.isFinite(tumorCenterZ)) tumorCenterZ = 0.15; // default MATLAB
    const zMax = z_bound_slab[z_bound_slab.length - 1];
    tumorCenterZ = Math.max(0, Math.min(zMax, tumorCenterZ)); // clamp al dominio

    // Ordenamos esferas por radio ascendente (anidadas)
    spheres.sort((a, b) => a.radius - b.radius);

    // Albedo para todas las filas (slabs + esferas)
    const albedo = new Float64Array(Nc);
    for (let i = 0; i < Nc; i++) albedo[i] = mus[i] / (mus[i] + mua[i]);

    // --- Simulación ---
    const t0 = performance.now();
    const energy = simulate(
      N, mua, mus, g, albedo, z_bound_slab, Rxy,
      {
        center: [0, 0, tumorCenterZ], // centro del tumor según usuario
        spheres,
        slabIdx
      }
    );
    const t1 = performance.now();

    // Resultados
    let energyT = 0.0;
    for (let i = 0; i < energy.length; i++) energyT += energy[i];
    const lost = (N - energyT) / N;
    const elapsed = (t1 - t0) / 1000;

    let energyStr = "<ul>";
    for (let i = 0; i < energy.length; i++) {
      energyStr += "<li>";
      energyStr += `<b>${tipo[i]}:</b> ${energy[i]}\n`;
      energyStr += "</li>";
    }
    energyStr += "</ul>";

    console.log(tumorCenterZ)
    const threePayload = {
      Rxy,
      zBound: Array.from(z_bound_slab),
      slabIdx: Array.from(slabIdx),
      tipo: Array.from(tipo),
      spheres: spheres.map(s => ({ index: s.index, radius: s.radius, tipo: s.tipo })),
      tumorCenterZ,
      energy: Array.from(energy),
      tumorColor:0x005AD9,
      tumorCenterZ:tumorCenterZ,
      energyloss:lost.toFixed(6),
      lambda:parseInt(document.getElementById("lambda").value)
    };
    
    try{
      open3DViewer(threePayload);
    }catch{
      alert("Tiempo de espera excesivo, no se pudo abrir la simulación 3D");
    }


    salida.innerHTML =
      "<b>Energy per layer: </b>\n" + energyStr +
      "<b>Total absorbed energy: </b>" + energyT.toFixed(6) + "\n\n" +
      "<b>Total energy loss: </b><i style='color:red'>" + lost.toFixed(6) + "</i>\n\n" +
      "<b>Time: </b>" + elapsed.toFixed(3) + " s";

      // --- Dibujo (corte a escala) ---
      renderEsquemaDesdeTabla(Rxy, tipo, z_end, tumorCenterZ, energy);
    });
});
// --- Cargar CSV y poblar tabla ---
(function attachCsvLoader(){
  const fileInput = document.getElementById("fileCSV");
  if (!fileInput) return;

  fileInput.addEventListener("change", async function () {
    const file = this.files && this.files[0];
    if (!file) return;

    try {
      const text = await file.text();
      // Parseo CSV simple (asume coma como separador)
      const lines = text.split(/\r?\n/).filter(l => l.trim().length);
      if (lines.length === 0) throw new Error("Archivo vacío.");

      // Detectar encabezado
      const header = lines[0].split(",").map(h => h.trim().toLowerCase());
      const map = {
        tipo: header.indexOf("tipo"),
        mua: header.indexOf("mua"),
        mus: header.indexOf("mus"),
        g: header.indexOf("g"),
        nt: header.indexOf("nt"),
        zend: header.indexOf("z_end"),
      };
      if (Object.values(map).some(idx => idx === -1)) {
        throw new Error('Encabezado esperado: "Tipo,mua,mus,g,nt,z_end"');
      }

      // Construir arreglo de filas
      const rows = [];
      for (let i = 1; i < lines.length; i++) {
        const cols = lines[i].split(",").map(c => c.trim());
        if (cols.length < header.length) continue; // línea incompleta
        rows.push({
          tipo: cols[map.tipo],
          mua: parseFloat(cols[map.mua]),
          mus: parseFloat(cols[map.mus]),
          g:   parseFloat(cols[map.g]),
          nt:  parseFloat(cols[map.nt]),
          z_end: parseFloat(cols[map.zend]),
        });
      }
      if (!rows.length) throw new Error("No se encontraron filas válidas.");

      fillTableFromData(rows);
    } catch (err) {
      alert("Error al cargar CSV: " + (err?.message || err));
    } finally {
      this.value = ""; // reset para poder cargar el mismo archivo de nuevo si quieres
    }
  });

  function fillTableFromData(rows){
    const NCEl = document.getElementById("NC");
    const initBtn = document.getElementById("init");
    const tbody = document.querySelector("#tabla tbody");
    if (!NCEl || !initBtn || !tbody) {
      alert("No se encontró la tabla en el DOM.");
      return;
    }
    // Crear tabla con el número de filas del archivo
    NCEl.value = String(rows.length);
    initBtn.click();

    // Rellenar
    const trs = Array.from(tbody.querySelectorAll("tr"));
    trs.forEach((tr, i) => {
      const r = rows[i];
      const inputs = tr.querySelectorAll("input");
      inputs[0].value = r.tipo ?? "";
      inputs[1].value = isFinite(r.mua) ? r.mua : "";
      inputs[2].value = isFinite(r.mus) ? r.mus : "";
      inputs[3].value = isFinite(r.g)   ? r.g   : "";
      inputs[4].value = isFinite(r.nt)  ? r.nt  : "";
      inputs[5].value = isFinite(r.z_end) ? r.z_end : "";
    });
  }
})();
function open3DViewer(payload) {
  const viewer = window.open("viewer.html", "_blank");

  // Envía los datos al visor cuando esté listo
  const sendData = () => {
    if (viewer) {
      viewer.postMessage({ type: "mc3d", payload }, "*");
    }
  };

  // Espera a que cargue el visor antes de enviar
  viewer.addEventListener("load", sendData);
  setTimeout(sendData, 800); // respaldo por si load no dispara
}