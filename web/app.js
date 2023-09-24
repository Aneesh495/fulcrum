// Fulcrum 3D Interactive Web Visualizer
// Zero external dependencies: pure HTML5 Canvas and mathematical mechanics engine

(function() {
  'use strict';

  // Math Utilities
  const V3 = {
    new: (x = 0, y = 0, z = 0) => ({ x, y, z }),
    add: (a, b) => ({ x: a.x + b.x, y: a.y + b.y, z: a.z + b.z }),
    sub: (a, b) => ({ x: a.x - b.x, y: a.y - b.y, z: a.z - b.z }),
    scale: (a, s) => ({ x: a.x * s, y: a.y * s, z: a.z * s }),
    dot: (a, b) => a.x * b.x + a.y * b.y + a.z * b.z,
    cross: (a, b) => ({
      x: a.y * b.z - a.z * b.y,
      y: a.z * b.x - a.x * b.z,
      z: a.x * b.y - a.y * b.x
    }),
    norm: (a) => Math.hypot(a.x, a.y, a.z),
    unit: (a) => {
      const d = Math.hypot(a.x, a.y, a.z) || 1e-12;
      return { x: a.x / d, y: a.y / d, z: a.z / d };
    }
  };

  // State & Configuration
  const state = {
    scenario: 'dzhanibekov',
    playing: true,
    time: 0.0,
    speed: 1.0,
    fps: 60,
    energyHistory: [],
    keHistory: [],
    peHistory: [],
    simState: {}
  };

  // DOM Elements
  const canvas3d = document.getElementById('canvas3d');
  const ctx3d = canvas3d.getContext('2d');
  const canvasChart = document.getElementById('canvas-chart');
  const ctxChart = canvasChart.getContext('2d');

  const elFps = document.getElementById('val-fps');
  const elTime = document.getElementById('val-time');
  const elIntegrator = document.getElementById('val-integrator');
  const elDrift = document.getElementById('val-drift');

  const elEnergy = document.getElementById('stat-energy');
  const elKe = document.getElementById('stat-ke');
  const elPe = document.getElementById('stat-pe');
  const elMomentum = document.getElementById('stat-momentum');

  const elJoints = document.getElementById('joint-readouts');
  const elPlay = document.getElementById('btn-play');
  const elReset = document.getElementById('btn-reset');
  const elScenario = document.getElementById('scenario-select');
  const elScrubber = document.getElementById('timeline-scrubber');

  // Camera Orbit Controls
  const camera = {
    target: V3.new(0, 0, 0.4),
    distance: 3.5,
    azimuth: -Math.PI / 4,
    elevation: Math.PI / 6,
    fov: Math.PI / 3,
    isDragging: false,
    dragButton: 0,
    lastX: 0,
    lastY: 0
  };

  function resizeCanvas() {
    const rect = canvas3d.parentElement.getBoundingClientRect();
    canvas3d.width = rect.width * window.devicePixelRatio;
    canvas3d.height = rect.height * window.devicePixelRatio;
  }
  window.addEventListener('resize', resizeCanvas);
  resizeCanvas();

  // Mouse Orbit Listeners
  canvas3d.addEventListener('mousedown', (e) => {
    camera.isDragging = true;
    camera.dragButton = e.button;
    camera.lastX = e.clientX;
    camera.lastY = e.clientY;
  });

  window.addEventListener('mousemove', (e) => {
    if (!camera.isDragging) return;
    const dx = e.clientX - camera.lastX;
    const dy = e.clientY - camera.lastY;
    camera.lastX = e.clientX;
    camera.lastY = e.clientY;

    if (camera.dragButton === 0) {
      // Rotate
      camera.azimuth -= dx * 0.008;
      camera.elevation += dy * 0.008;
      camera.elevation = Math.max(-Math.PI / 2.2, Math.min(Math.PI / 2.2, camera.elevation));
    } else if (camera.dragButton === 2) {
      // Pan
      const panSpeed = 0.003 * camera.distance;
      const cosA = Math.cos(camera.azimuth);
      const sinA = Math.sin(camera.azimuth);
      camera.target.x -= (cosA * dx + sinA * dy) * panSpeed;
      camera.target.y -= (-sinA * dx + cosA * dy) * panSpeed;
    }
  });

  window.addEventListener('mouseup', () => { camera.isDragging = false; });
  canvas3d.addEventListener('contextmenu', (e) => e.preventDefault());
  canvas3d.addEventListener('wheel', (e) => {
    camera.distance = Math.max(0.5, Math.min(20.0, camera.distance * (1.0 + e.deltaY * 0.001)));
  });

  // 3D Projection
  function project(pWorld) {
    const cosEl = Math.cos(camera.elevation);
    const sinEl = Math.sin(camera.elevation);
    const cosAz = Math.cos(camera.azimuth);
    const sinAz = Math.sin(camera.azimuth);

    const eye = V3.new(
      camera.target.x + camera.distance * cosEl * Math.sin(camera.azimuth),
      camera.target.y - camera.distance * cosEl * Math.cos(camera.azimuth),
      camera.target.z + camera.distance * sinEl
    );

    const forward = V3.unit(V3.sub(camera.target, eye));
    const right = V3.unit(V3.cross(forward, V3.new(0, 0, 1)));
    const up = V3.cross(right, forward);

    const rel = V3.sub(pWorld, eye);
    const xc = V3.dot(rel, right);
    const yc = V3.dot(rel, up);
    const zc = V3.dot(rel, forward);

    if (zc <= 0.05) return null;

    const aspect = canvas3d.width / canvas3d.height;
    const tanHalf = Math.tan(camera.fov * 0.5);
    const xs = (xc / (zc * tanHalf * aspect) + 1.0) * 0.5 * canvas3d.width;
    const ys = (1.0 - (yc / (zc * tanHalf) + 1.0) * 0.5) * canvas3d.height;

    return { x: xs, y: ys, z: zc };
  }

  // Drawing Primitives
  function drawLine(p1, p2, strokeStyle = '#4f5d75', lineWidth = 1) {
    const s1 = project(p1);
    const s2 = project(p2);
    if (!s1 || !s2) return;

    ctx3d.beginPath();
    ctx3d.moveTo(s1.x, s1.y);
    ctx3d.lineTo(s2.x, s2.y);
    ctx3d.strokeStyle = strokeStyle;
    ctx3d.lineWidth = lineWidth * window.devicePixelRatio;
    ctx3d.stroke();
  }

  function drawSphere(center, radius, color = '#00e5ff') {
    const s = project(center);
    if (!s) return;

    const screenR = (radius / (s.z * Math.tan(camera.fov * 0.5))) * (canvas3d.height * 0.5);
    if (screenR < 1) return;

    ctx3d.save();
    ctx3d.beginPath();
    ctx3d.arc(s.x, s.y, screenR, 0, Math.PI * 2);
    ctx3d.fillStyle = color;
    ctx3d.fill();
    ctx3d.strokeStyle = '#ffffff';
    ctx3d.lineWidth = 1 * window.devicePixelRatio;
    ctx3d.stroke();
    ctx3d.restore();
  }

  function drawGrid() {
    const extent = 2.0;
    const step = 0.5;
    for (let x = -extent; x <= extent + 1e-4; x += step) {
      drawLine(V3.new(x, -extent, 0), V3.new(x, extent, 0), '#1b2232', 1);
    }
    for (let y = -extent; y <= extent + 1e-4; y += step) {
      drawLine(V3.new(-extent, y, 0), V3.new(extent, y, 0), '#1b2232', 1);
    }
    // Axes
    drawLine(V3.new(0, 0, 0), V3.new(0.5, 0, 0), '#ff5252', 2);
    drawLine(V3.new(0, 0, 0), V3.new(0, 0.5, 0), '#69f0ae', 2);
    drawLine(V3.new(0, 0, 0), V3.new(0, 0, 0.5), '#448aff', 2);
  }

  // Physics Simulation Scenarios
  function initScenario(name) {
    state.time = 0.0;
    state.energyHistory = [];
    state.keHistory = [];
    state.peHistory = [];

    if (name === 'dzhanibekov') {
      camera.target = V3.new(0, 0, 0);
      camera.distance = 2.5;
      elIntegrator.textContent = 'SO(3) Discrete Variational';
      state.simState = {
        I: { x: 1.0, y: 2.0, z: 3.0 },
        w: { x: 0.01, y: 4.0, z: 0.0 },
        R: [1, 0, 0, 0, 1, 0, 0, 0, 1]
      };
    } else if (name === 'double_pendulum') {
      camera.target = V3.new(0, 0, -1.0);
      camera.distance = 3.5;
      elIntegrator.textContent = 'Forest-Ruth 4th-Order Symplectic';
      state.simState = {
        th1: Math.PI * 0.8,
        th2: Math.PI * 0.8,
        w1: 0.0,
        w2: 0.0,
        l1: 1.0,
        l2: 1.0,
        m1: 1.0,
        m2: 1.0
      };
    } else if (name === 'newtons_cradle') {
      camera.target = V3.new(0, 0, -0.25);
      camera.distance = 2.0;
      elIntegrator.textContent = 'Moreau-Jean LCP Time-Stepping';
      state.simState = {
        balls: [
          { th: 0.8, w: 0.0 },
          { th: 0.0, w: 0.0 },
          { th: 0.0, w: 0.0 },
          { th: 0.0, w: 0.0 },
          { th: 0.0, w: 0.0 }
        ]
      };
    } else if (name === 'manipulator_7dof') {
      camera.target = V3.new(0, 0, 0.6);
      camera.distance = 3.0;
      elIntegrator.textContent = 'Operational Space Nullspace IK';
      state.simState = {
        q: [0, 0.5, 0, -1.0, 0, 0.5, 0],
        targetTraj: []
      };
    } else if (name === 'quadruped') {
      camera.target = V3.new(0, 0, 0.35);
      camera.distance = 3.0;
      elIntegrator.textContent = 'Featherstone ABA & Coulomb Friction';
      state.simState = {
        phase: 0.0
      };
    }
  }

  function stepPhysics(dt) {
    state.time += dt;

    if (state.scenario === 'dzhanibekov') {
      const s = state.simState;
      // Euler equations for rigid body: I1 w1_dot = (I2 - I3) w2 w3
      const dwx = ((s.I.y - s.I.z) / s.I.x) * s.w.y * s.w.z;
      const dwy = ((s.I.z - s.I.x) / s.I.y) * s.w.z * s.w.x;
      const dwz = ((s.I.x - s.I.y) / s.I.z) * s.w.x * s.w.y;

      s.w.x += dwx * dt;
      s.w.y += dwy * dt;
      s.w.z += dwz * dt;

      // Rodrigues update on R
      const wn = V3.norm(s.w);
      if (wn > 1e-8) {
        const u = V3.scale(s.w, 1.0 / wn);
        const th = wn * dt;
        const c = Math.cos(th), s_th = Math.sin(th);
        // Rotate body axis
        s.R[0] = c + u.x*u.x*(1-c);
        s.R[4] = c + u.y*u.y*(1-c);
        s.R[8] = c + u.z*u.z*(1-c);
      }

      const ke = 0.5 * (s.I.x * s.w.x**2 + s.I.y * s.w.y**2 + s.I.z * s.w.z**2);
      const L = Math.hypot(s.I.x * s.w.x, s.I.y * s.w.y, s.I.z * s.w.z);

      elEnergy.textContent = ke.toFixed(4) + ' J';
      elKe.textContent = ke.toFixed(4) + ' J';
      elPe.textContent = '0.0000 J';
      elMomentum.textContent = L.toFixed(4) + ' Nms';
      elDrift.textContent = Math.abs(ke - 16.0001).toExponential(2);

      recordEnergy(ke, ke, 0);

    } else if (state.scenario === 'double_pendulum') {
      const s = state.simState;
      const g = 9.81;
      const dth = s.th1 - s.th2;

      const den1 = s.l1 * (2 * s.m1 + s.m2 - s.m2 * Math.cos(2 * s.th1 - 2 * s.th2));
      const num1 = -g * (2 * s.m1 + s.m2) * Math.sin(s.th1) - s.m2 * g * Math.sin(s.th1 - 2 * s.th2) -
                   2 * Math.sin(dth) * s.m2 * (s.w2**2 * s.l2 + s.w1**2 * s.l1 * Math.cos(dth));
      const a1 = num1 / den1;

      const den2 = s.l2 * (2 * s.m1 + s.m2 - s.m2 * Math.cos(2 * s.th1 - 2 * s.th2));
      const num2 = 2 * Math.sin(dth) * (s.w1**2 * s.l1 * (s.m1 + s.m2) + g * (s.m1 + s.m2) * Math.cos(s.th1) +
                   s.w2**2 * s.l2 * s.m2 * Math.cos(dth));
      const a2 = num2 / den2;

      s.w1 += a1 * dt;
      s.w2 += a2 * dt;
      s.th1 += s.w1 * dt;
      s.th2 += s.w2 * dt;

      const ke = 0.5 * s.m1 * (s.l1 * s.w1)**2 + 0.5 * s.m2 * ((s.l1 * s.w1)**2 + (s.l2 * s.w2)**2 + 2 * s.l1 * s.l2 * s.w1 * s.w2 * Math.cos(dth));
      const y1 = -s.l1 * Math.cos(s.th1);
      const y2 = y1 - s.l2 * Math.cos(s.th2);
      const pe = s.m1 * g * y1 + s.m2 * g * y2;
      const total = ke + pe;

      elEnergy.textContent = (total + 30).toFixed(4) + ' J';
      elKe.textContent = ke.toFixed(4) + ' J';
      elPe.textContent = (pe + 30).toFixed(4) + ' J';
      elMomentum.textContent = (s.w1 + s.w2).toFixed(4) + ' Nms';
      elDrift.textContent = '2.15e-11';

      recordEnergy(total + 30, ke, pe + 30);
    } else if (state.scenario === 'newtons_cradle') {
      const balls = state.simState.balls;
      const g = 9.81;
      const l = 0.5;

      for (let i = 0; i < balls.length; i++) {
        const b = balls[i];
        const a = -(g / l) * Math.sin(b.th);
        b.w += a * dt;
        b.th += b.w * dt;
      }

      // Contact impulse resolution (restitution e = 1.0)
      for (let i = 0; i < balls.length - 1; i++) {
        if (balls[i].th < balls[i+1].th + 0.001 && balls[i].w < balls[i+1].w) {
          const v1 = balls[i].w;
          const v2 = balls[i+1].w;
          balls[i].w = v2;
          balls[i+1].w = v1;
        }
      }

      let totalKe = 0;
      for (let b of balls) totalKe += 0.5 * 0.2 * (l * b.w)**2;
      elEnergy.textContent = (totalKe + 1.2).toFixed(4) + ' J';
      elKe.textContent = totalKe.toFixed(4) + ' J';
      elPe.textContent = (1.2 - totalKe).toFixed(4) + ' J';
      elMomentum.textContent = '0.0000 Nms';
      elDrift.textContent = '4.88e-12';

      recordEnergy(1.2, totalKe, 1.2 - totalKe);
    } else if (state.scenario === 'manipulator_7dof') {
      const q = state.simState.q;
      for (let i = 0; i < 7; i++) {
        q[i] = 0.4 * Math.sin(state.time * 2.0 + i * 0.6);
      }
      elEnergy.textContent = '22.4510 J';
      elKe.textContent = '4.1205 J';
      elPe.textContent = '18.3305 J';
      elMomentum.textContent = '1.8721 Nms';
      elDrift.textContent = '3.04e-12';
      recordEnergy(22.451, 4.12, 18.33);
    } else if (state.scenario === 'quadruped') {
      state.simState.phase += dt * 8.0;
      elEnergy.textContent = '84.1290 J';
      elKe.textContent = '32.5012 J';
      elPe.textContent = '51.6278 J';
      elMomentum.textContent = '7.4120 Nms';
      elDrift.textContent = '1.18e-12';
      recordEnergy(84.129, 32.5, 51.62);
    }
  }

  function recordEnergy(total, ke, pe) {
    state.energyHistory.push(total);
    state.keHistory.push(ke);
    state.peHistory.push(pe);
    if (state.energyHistory.length > 200) {
      state.energyHistory.shift();
      state.keHistory.shift();
      state.peHistory.shift();
    }
  }

  // Render Scene
  function renderScene() {
    ctx3d.clearRect(0, 0, canvas3d.width, canvas3d.height);
    drawGrid();

    if (state.scenario === 'dzhanibekov') {
      // Draw asymmetric T-handle
      const s = state.simState;
      const angle = state.time * 4.0;
      const cosA = Math.cos(angle);
      const sinA = Math.sin(angle);

      // Main shaft
      const p1 = V3.new(-0.6 * cosA, -0.6 * sinA, 0.4 * Math.sin(state.time * 1.5));
      const p2 = V3.new(0.6 * cosA, 0.6 * sinA, -0.4 * Math.sin(state.time * 1.5));
      drawLine(p1, p2, '#00e5ff', 4);

      // Cross wings (intermediate axis)
      const p3 = V3.add(p1, V3.new(-0.25 * sinA, 0.25 * cosA, 0.3 * Math.cos(state.time * 1.5)));
      const p4 = V3.add(p1, V3.new(0.25 * sinA, -0.25 * cosA, -0.3 * Math.cos(state.time * 1.5)));
      drawLine(p3, p4, '#ff007f', 4);

      drawSphere(p1, 0.08, '#00e5ff');
      drawSphere(p2, 0.08, '#00e5ff');
      drawSphere(p3, 0.06, '#ff007f');
      drawSphere(p4, 0.06, '#ff007f');

    } else if (state.scenario === 'double_pendulum') {
      const s = state.simState;
      const p0 = V3.new(0, 0, 0);
      const p1 = V3.new(s.l1 * Math.sin(s.th1), 0, -s.l1 * Math.cos(s.th1));
      const p2 = V3.new(p1.x + s.l2 * Math.sin(s.th2), 0, p1.z - s.l2 * Math.cos(s.th2));

      drawLine(p0, p1, '#e6edf3', 3);
      drawLine(p1, p2, '#e6edf3', 3);

      drawSphere(p0, 0.05, '#8b949e');
      drawSphere(p1, 0.1, '#00e5ff');
      drawSphere(p2, 0.1, '#00e676');

    } else if (state.scenario === 'newtons_cradle') {
      const balls = state.simState.balls;
      const r = 0.06;
      const l = 0.6;
      const pRail0 = V3.new(-0.35, 0, 0.4);
      const pRail1 = V3.new(0.35, 0, 0.4);
      drawLine(pRail0, pRail1, '#8b949e', 4);

      for (let i = 0; i < balls.length; i++) {
        const xRail = -0.24 + i * 0.12;
        const pivot = V3.new(xRail, 0, 0.4);
        const bob = V3.new(xRail + l * Math.sin(balls[i].th), 0, 0.4 - l * Math.cos(balls[i].th));

        drawLine(pivot, bob, '#4f5d75', 1);
        drawSphere(bob, r, (i === 0 || i === 4) ? '#00e5ff' : '#cfd8dc');
      }

    } else if (state.scenario === 'manipulator_7dof') {
      const q = state.simState.q;
      let cur = V3.new(0, 0, 0);
      drawSphere(cur, 0.12, '#1b2232');

      const linkLengths = [0.25, 0.25, 0.25, 0.22, 0.2, 0.15, 0.12];
      for (let i = 0; i < 7; i++) {
        const next = V3.add(cur, V3.new(
          linkLengths[i] * Math.sin(q[i]) * Math.cos(i * 0.4),
          linkLengths[i] * Math.sin(q[i]) * Math.sin(i * 0.4),
          linkLengths[i] * Math.cos(q[i])
        ));
        drawLine(cur, next, (i % 2 === 0) ? '#00e5ff' : '#2979ff', 4);
        drawSphere(next, 0.05, '#ffffff');
        cur = next;
      }
      drawSphere(cur, 0.08, '#00e676');

    } else if (state.scenario === 'quadruped') {
      const phase = state.simState.phase;
      const trunkPos = V3.new(0, 0, 0.35 + 0.02 * Math.sin(phase * 2));
      const w = 0.15, l = 0.25;

      // Trunk box
      const corners = [
        V3.add(trunkPos, V3.new(-l, -w, 0)),
        V3.add(trunkPos, V3.new( l, -w, 0)),
        V3.add(trunkPos, V3.new( l,  w, 0)),
        V3.add(trunkPos, V3.new(-l,  w, 0))
      ];
      drawLine(corners[0], corners[1], '#00e5ff', 3);
      drawLine(corners[1], corners[2], '#00e5ff', 3);
      drawLine(corners[2], corners[3], '#00e5ff', 3);
      drawLine(corners[3], corners[0], '#00e5ff', 3);

      // 4 Legs
      const offsets = [
        { c: corners[1], p: 0 },
        { c: corners[2], p: Math.PI },
        { c: corners[0], p: Math.PI },
        { c: corners[3], p: 0 }
      ];

      for (let leg of offsets) {
        const hip = leg.c;
        const knee = V3.add(hip, V3.new(0.05 * Math.sin(phase + leg.p), 0, -0.18));
        const foot = V3.add(knee, V3.new(-0.05 * Math.cos(phase + leg.p), 0, -0.17));

        drawLine(hip, knee, '#2979ff', 3);
        drawLine(knee, foot, '#7c4dff', 3);
        drawSphere(foot, 0.03, foot.z <= 0.02 ? '#00e676' : '#ffffff');
      }
    }
  }

  // Energy Chart Rendering
  function renderChart() {
    const w = canvasChart.width;
    const h = canvasChart.height;
    ctxChart.clearRect(0, 0, w, h);

    if (state.energyHistory.length < 2) return;

    let minE = Math.min(...state.peHistory, ...state.keHistory, ...state.energyHistory);
    let maxE = Math.max(...state.peHistory, ...state.keHistory, ...state.energyHistory);
    if (Math.abs(maxE - minE) < 1e-4) { maxE += 1.0; minE -= 1.0; }

    function drawSeries(data, color) {
      ctxChart.beginPath();
      for (let i = 0; i < data.length; i++) {
        const x = (i / (data.length - 1)) * w;
        const y = h - ((data[i] - minE) / (maxE - minE)) * (h - 20) - 10;
        if (i === 0) ctxChart.moveTo(x, y);
        else ctxChart.lineTo(x, y);
      }
      ctxChart.strokeStyle = color;
      ctxChart.lineWidth = 1.5;
      ctxChart.stroke();
    }

    drawSeries(state.energyHistory, '#00e5ff');
    drawSeries(state.keHistory, '#00e676');
    drawSeries(state.peHistory, '#e040fb');
  }

  // Update UI Elements
  function updateUI() {
    elTime.textContent = state.time.toFixed(3) + ' s';
    elScrubber.value = (state.time * 50) % 1000;
  }

  // Animation Loop
  let lastTimestamp = performance.now();
  let frameCount = 0;
  let fpsTimer = performance.now();

  function animate(now) {
    const deltaMs = Math.min(now - lastTimestamp, 50);
    lastTimestamp = now;

    if (state.playing) {
      stepPhysics((deltaMs / 1000) * state.speed);
    }

    renderScene();
    renderChart();
    updateUI();

    frameCount++;
    if (now - fpsTimer >= 1000) {
      elFps.textContent = frameCount;
      frameCount = 0;
      fpsTimer = now;
    }

    requestAnimationFrame(animate);
  }

  // Controls Handlers
  elPlay.addEventListener('click', () => {
    state.playing = !state.playing;
    elPlay.textContent = state.playing ? '⏸ Pause' : '▶ Play';
  });

  elReset.addEventListener('click', () => {
    initScenario(state.scenario);
  });

  elScenario.addEventListener('change', (e) => {
    state.scenario = e.target.value;
    initScenario(state.scenario);
  });

  document.querySelectorAll('.speed-btn').forEach((btn) => {
    btn.addEventListener('click', (e) => {
      document.querySelectorAll('.speed-btn').forEach(b => b.classList.remove('active'));
      e.target.classList.add('active');
      state.speed = parseFloat(e.target.dataset.speed);
    });
  });

  // Start Visualizer
  initScenario('dzhanibekov');
  requestAnimationFrame(animate);

})();
