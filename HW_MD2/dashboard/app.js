/* Offline playback of saved data. This file never integrates equations of motion. */
'use strict';
const DATA = JSON.parse(document.getElementById('trajectory-data').textContent);
const $ = id => document.getElementById(id);
const ids = ['play', 'reset', 'frame', 'metric', 'particle', 'coordinates', 'window'];
const names = {displacement: '去质心位移 u', velocity: '去质心速度 v′', kinetic: '单粒子内部动能 K′'};
const units = {displacement: '位移', velocity: '速度', kinetic: '能量'};
const diverging = [[0, '#2166ac'], [0.5, '#f7f7f7'], [1, '#b2182b']];
const sequential = ['#000004', '#180f3d', '#440f76', '#721f81', '#9e2f7f',
  '#cd4071', '#f1605d', '#fd9668', '#feca8d', '#fcfdbf']
  .map((color, index, colors) => [index / (colors.length - 1), color]);
const state = {frame: 0, particle: 0, metric: 'velocity', coordinates: 'com10',
  range: [DATA.time[0], Math.min(DATA.time.at(-1), DATA.time[0] + 2)], playing: false};
const panels = [];
let rangeLock = 0;
let paintTask = null;
let dirty = false;
let playToken = 0;
let sizeObserver;
const plotConfig = {responsive: true, displaylogo: false, scrollZoom: false,
  modeBarButtonsToRemove: ['lasso2d', 'select2d', 'autoScale2d'],
  toImageButtonOptions: {format: 'png', scale: 2}};
const baseLayout = () => ({paper_bgcolor: 'white', plot_bgcolor: 'white',
  font: {family: 'Noto Sans CJK SC, Microsoft YaHei, system-ui, sans-serif', size: 11, color: '#405169'},
  margin: {l: 62, r: 65, t: 20, b: 42},
  xaxis: {title: {text: '时间（约化单位）'}, range: [...state.range], gridcolor: '#eef1f5', zeroline: false},
  yaxis: {gridcolor: '#eef1f5', zerolinecolor: '#dce2ea'},
  legend: {orientation: 'h', x: 0, y: 1.22, font: {size: 11}},
  uirevision: 'keep-view'});

function nearest(values, target) {
  let lo = 0, hi = values.length - 1;
  while (lo < hi) { const mid = (lo + hi) >> 1; if (values[mid] < target) lo = mid + 1; else hi = mid; }
  return lo > 0 && Math.abs(values[lo - 1] - target) <= Math.abs(values[lo] - target) ? lo - 1 : lo;
}
const transpose = matrix => matrix[0].map((_, j) => matrix.map(row => row[j]));
const fmt = value => Math.abs(value) < 1e-4 && value !== 0 ? value.toExponential(2) : value.toFixed(5);
function colorOptions(metric) {
  return {cmin: metric === 'kinetic' ? 0 : -DATA.limits[metric], cmax: DATA.limits[metric],
    colorscale: metric === 'kinetic' ? sequential : diverging};
}
function cursor(time) {
  return {type: 'line', xref: 'x', yref: 'paper', x0: time, x1: time, y0: 0, y1: 1,
    line: {color: '#233c58', width: 1.5}};
}
function fieldScale(metric) {
  const colors = colorOptions(metric);
  return {zmin: colors.cmin, zmax: colors.cmax, colorscale: colors.colorscale};
}
function displayPositions(item, frame) {
  if (state.coordinates === 'lab') return item.positions_unwrapped[frame];
  if (state.coordinates === 'wrapped') return item.positions[frame];
  const gain = state.coordinates === 'com10' ? 10 : 1;
  return item.positions_unwrapped[0].map((x, j) => x + gain * item.fields.displacement[frame][j]);
}
function chainRange() {
  let lo = Infinity, hi = -Infinity;
  for (const panel of panels) {
    for (let k = 0; k < panel.item.time.length; k++) {
      displayPositions(panel.item, k).forEach((x, j) => {
        const tip = x + panel.speedScale * panel.item.fields.velocity[k][j];
        lo = Math.min(lo, x, tip); hi = Math.max(hi, x, tip);
      });
    }
  }
  const padding = Math.max((hi - lo) * 0.07, 0.1);
  return [lo - padding, hi + padding];
}

async function makePanel(item, index) {
  const section = document.createElement('section');
  section.className = 'case';
  const label = item.boundary === 'free' ? '自由边界' : '周期边界';
  const n = item.positions[0].length;
  const config = item.metadata.config;
  section.innerHTML = `<div class="case-header"><h2>${label}</h2>
    <div class="case-meta">N = ${n} · Δt = ${config.dt} · rc = ${item.metadata.rcut.toFixed(2)} · NVE</div>
    <div class="sample-label" id="sample-${index}"></div></div>
    <p class="chart-title">粒子运动与速度方向</p><div class="chart chain" id="chain-${index}" role="img" aria-label="${label}粒子位置和速度箭头"></div>
    <p class="plot-note">灰圈为初始位置；箭头表示 v′ 的方向和相对大小。</p>
    <div class="readings" id="readings-${index}"></div>
    <p class="chart-title" id="heat-title-${index}">时空图</p><div class="chart heatmap" id="heat-${index}" role="img" aria-label="${label}粒子随时间变化的热图"></div>
    <p class="chart-title" id="history-title-${index}">所选粒子的运动</p><div class="chart history" id="history-${index}" role="img" aria-label="${label}所选粒子的位移速度历史"></div>
    <p class="chart-title">动势能交换与温度</p><div class="chart thermo" id="thermo-${index}" role="img" aria-label="${label}能量变化和温度历史"></div>`;
  $('cases').append(section);
  const panel = {item, index, n, section, heat: $(`heat-${index}`), chain: $(`chain-${index}`),
    history: $(`history-${index}`), thermo: $(`thermo-${index}`),
    matrices: Object.fromEntries(Object.entries(item.fields).map(([key, values]) => [key, transpose(values)])),
    speedScale: 0.6 * item.metadata.spacing / DATA.limits.velocity,
    lastParticle: null, lastMetric: null, lastCoordinates: null};
  panels.push(panel);
  const zeros = Array(n).fill(0);
  const particleIds = Array.from({length: n}, (_, j) => j + 1);
  const chainLayout = baseLayout();
  chainLayout.margin = {l: 40, r: 35, t: 26, b: 50};
  chainLayout.xaxis = {title: {text: '显示位置 x(0) + 10u'}, zeroline: false, showgrid: false};
  chainLayout.yaxis = {visible: false, range: [-0.35, 0.55], fixedrange: true};
  chainLayout.showlegend = false;
  await Plotly.newPlot(panel.chain, [
    {x: item.positions_unwrapped[0], y: zeros, mode: 'markers', marker: {symbol: 'circle-open', color: '#b5becb', size: 17}, hoverinfo: 'skip'},
    {x: [], y: [], mode: 'lines', line: {color: '#8997a8', width: 1.5}, hoverinfo: 'skip'},
    {x: [], y: [], mode: 'markers', marker: {color: '#8997a8', size: 7}, hoverinfo: 'skip'},
    {x: item.positions[0], y: zeros, mode: 'markers+text', text: particleIds, textposition: 'top center',
      marker: {size: 16, color: item.fields.velocity[0], ...colorOptions('velocity')},
      hovertemplate: '粒子 %{customdata[0]}<br>u=%{customdata[1]:.5f}<br>v′=%{customdata[2]:.5f}<br>K′=%{customdata[3]:.6f}<br>x 连续=%{customdata[4]:.5f}<extra></extra>'},
    {x: [], y: [0], mode: 'markers', marker: {symbol: 'circle-open', size: 24, color: '#172c45', line: {width: 2}}, hoverinfo: 'skip'}
  ], chainLayout, plotConfig);
  const heatLayout = baseLayout();
  heatLayout.margin = {l: 52, r: 92, t: 12, b: 42};
  heatLayout.yaxis = {title: {text: '粒子编号'}, range: [0.5, n + 0.5], dtick: n <= 16 ? 1 : undefined, fixedrange: true};
  await Plotly.newPlot(panel.heat, [{type: 'heatmap', x: item.time_edges, y: particleIds,
    z: panel.matrices.velocity, zsmooth: false, ...fieldScale('velocity'),
    customdata: particleIds.map(() => item.time),
    colorbar: {title: {text: '速度'}, thickness: 12, len: 0.9},
    hovertemplate: 't=%{customdata:.4f}<br>粒子 %{y}<br>值=%{z:.6f}<extra></extra>'}], heatLayout, plotConfig);
  const historyLayout = baseLayout();
  historyLayout.yaxis.title = {text: '位移 u'};
  historyLayout.yaxis.range = [-1.08 * DATA.limits.displacement, 1.08 * DATA.limits.displacement];
  historyLayout.yaxis2 = {title: {text: '速度 v′'}, overlaying: 'y', side: 'right', showgrid: false, zeroline: false,
    range: [-1.08 * DATA.limits.velocity, 1.08 * DATA.limits.velocity]};
  await Plotly.newPlot(panel.history, [
    {x: item.time, y: [], mode: 'lines', name: 'u', line: {color: '#2868a5', width: 1.3}},
    {x: item.time, y: [], mode: 'lines', name: 'v′', yaxis: 'y2', line: {color: '#b96026', width: 1.1}},
    {x: [], y: [], mode: 'markers', showlegend: false, marker: {color: '#2868a5', size: 7}},
    {x: [], y: [], mode: 'markers', yaxis: 'y2', showlegend: false, marker: {color: '#b96026', symbol: 'diamond', size: 7}}
  ], historyLayout, plotConfig);
  const thermoLayout = baseLayout();
  thermoLayout.yaxis.title = {text: '能量变化'};
  let energyMax = 1e-15, temperatureMax = 1e-15;
  for (const data of DATA.cases) for (const row of data.thermal) {
    energyMax = Math.max(energyMax, Math.abs(row[2] - data.thermal[0][2]), Math.abs(row[3] - data.thermal[0][3]));
    temperatureMax = Math.max(temperatureMax, row[5]);
  }
  thermoLayout.yaxis.range = [-1.08 * energyMax, 1.08 * energyMax];
  thermoLayout.yaxis2 = {title: {text: 'T*'}, overlaying: 'y', side: 'right', showgrid: false, zeroline: false,
    range: [0, 1.08 * temperatureMax]};
  await Plotly.newPlot(panel.thermo, [
    {x: item.time, y: item.thermal.map(row => row[2] - item.thermal[0][2]), name: 'Δ势能', mode: 'lines', line: {color: '#2868a5', width: 1.2}},
    {x: item.time, y: item.thermal.map(row => row[3] - item.thermal[0][3]), name: 'Δ动能', mode: 'lines', line: {color: '#b96026', width: 1.2}},
    {x: item.time, y: item.thermal.map(row => row[5]), name: 'T*（右轴）', yaxis: 'y2', mode: 'lines', line: {color: '#6d5796', width: 1, dash: 'dot'}}
  ], thermoLayout, plotConfig);
  panel.chain.on('plotly_click', event => {
    const point = event.points[0];
    if (point.curveNumber === 3) selectParticle(point.pointNumber);
  });
  panel.heat.on('plotly_click', event => {
    pause();
    const point = event.points[0];
    state.particle = Math.max(0, Math.min(Number($('particle').options.length) - 1, Math.round(point.y) - 1));
    // A heatmap x is the cell center; use its actual saved time from customdata.
    state.frame = boundedFrame(nearest(DATA.time, Number(point.customdata ?? point.x)));
    render();
  });
  for (const plot of [panel.heat, panel.history, panel.thermo]) {
    plot.on('plotly_relayout', event => {
      if (rangeLock) return;
      if (event['xaxis.autorange']) { applyRange(DATA.time[0], DATA.time.at(-1)); return; }
      const lo = event['xaxis.range[0]'] ?? event['xaxis.range']?.[0];
      const hi = event['xaxis.range[1]'] ?? event['xaxis.range']?.[1];
      if (lo !== undefined && hi !== undefined) applyRange(Number(lo), Number(hi));
    });
  }
  for (const plot of [panel.history, panel.thermo]) plot.on('plotly_click', event => {
    pause(); state.frame = boundedFrame(nearest(DATA.time, event.points[0].x)); render();
  });
}

function boundedFrame(frame) { return Math.max(Number($('frame').min), Math.min(Number($('frame').max), Math.round(frame))); }
function pause() { state.playing = false; playToken++; $('play').textContent = '播放'; }
async function selectParticle(particle) {
  state.particle = Math.max(0, Math.min($('particle').options.length - 1, Number(particle)));
  return render();
}
async function applyRange(lo, hi) {
  lo = Math.max(DATA.time[0], lo); hi = Math.min(DATA.time.at(-1), hi);
  if (!Number.isFinite(lo) || !Number.isFinite(hi) || hi <= lo) return;
  pause(); state.range = [lo, hi];
  const earlyEnd = Math.min(DATA.time[0] + 2, DATA.time.at(-1));
  $('window').value = Math.abs(lo - DATA.time[0]) < 1e-9 && Math.abs(hi - earlyEnd) < 1e-9 ? 'early' :
    Math.abs(lo - DATA.time[0]) < 1e-9 && Math.abs(hi - DATA.time.at(-1)) < 1e-9 ? 'full' : 'custom';
  let first = nearest(DATA.time, lo), last = nearest(DATA.time, hi);
  if (DATA.time[first] < lo && first < last) first++;
  if (DATA.time[last] > hi && last > first) last--;
  $('frame').min = first; $('frame').max = last;
  state.frame = boundedFrame(state.frame);
  rangeLock++;
  try {
    await Promise.all(panels.flatMap(panel => [panel.heat, panel.history, panel.thermo])
      .map(plot => Plotly.relayout(plot, {'xaxis.range': [...state.range], 'xaxis.autorange': false})));
  } finally { rangeLock--; }
  return render();
}

async function draw() {
  const target = DATA.time[state.frame];
  const metric = state.metric, particle = state.particle, coordinate = state.coordinates;
  $('frame').value = state.frame; $('particle').value = particle; $('metric').value = metric; $('coordinates').value = coordinate;
  $('time-label').textContent = `t = ${target.toFixed(4)} / ${state.range[1].toFixed(2)}`;
  const coordinateChanged = panels.some(panel => panel.lastCoordinates !== coordinate);
  const xrange = coordinateChanged ? chainRange() : null;
  await Promise.all(panels.map(async panel => {
    const item = panel.item, k = nearest(item.time, target), actual = item.time[k];
    panel.section.dataset.frame = k;
    const x = displayPositions(item, k), velocity = item.fields.velocity[k];
    const arrowX = [], arrowY = [], tips = [];
    for (let j = 0; j < panel.n; j++) {
      const tip = x[j] + panel.speedScale * velocity[j];
      arrowX.push(x[j], tip, null); arrowY.push(0, 0, null); tips.push(tip);
    }
    const updates = [
      Plotly.restyle(panel.chain, {x: [arrowX], y: [arrowY]}, [1]),
      Plotly.restyle(panel.chain, {x: [tips], y: [Array(panel.n).fill(0)],
        'marker.symbol': [velocity.map(v => v >= 0 ? 'triangle-right' : 'triangle-left')],
        'marker.size': [velocity.map(v => Math.abs(v) < 1e-10 ? 0 : 7)]}, [2]),
      Plotly.restyle(panel.chain, {x: [x], 'marker.color': [item.fields[metric][k]],
        'marker.cmin': colorOptions(metric).cmin, 'marker.cmax': colorOptions(metric).cmax,
        'marker.colorscale': [colorOptions(metric).colorscale],
        customdata: [x.map((_, j) => [j + 1, item.fields.displacement[k][j], velocity[j], item.fields.kinetic[k][j], item.positions_unwrapped[k][j]])]}, [3]),
      Plotly.restyle(panel.chain, {x: [[x[particle]]]}, [4]),
      Plotly.restyle(panel.history, {x: [[actual], [actual]],
        y: [[item.fields.displacement[k][particle]], [velocity[particle]]]}, [2, 3]),
      Plotly.relayout(panel.heat, {shapes: [cursor(actual), {type: 'rect', xref: 'paper', yref: 'y',
        x0: 0, x1: 1, y0: particle + .55, y1: particle + 1.45,
        fillcolor: 'rgba(0,0,0,0)', line: {color: '#233c58', width: 1}}]}),
      Plotly.relayout(panel.history, {shapes: [cursor(actual)]}),
      Plotly.relayout(panel.thermo, {shapes: [cursor(actual)]})
    ];
    if (panel.lastParticle !== particle) {
      updates.push(Plotly.restyle(panel.history, {y: [panel.matrices.displacement[particle], panel.matrices.velocity[particle]]}, [0, 1]));
      $(`history-title-${panel.index}`).textContent = `粒子 ${particle + 1} · 位移与速度`;
    }
    if (panel.lastMetric !== metric) {
      updates.push(Plotly.restyle(panel.heat, {z: [panel.matrices[metric]], ...fieldScale(metric),
        colorscale: [fieldScale(metric).colorscale], 'colorbar.title.text': units[metric]}, [0]));
      $(`heat-title-${panel.index}`).textContent = `时空图 · ${names[metric]}`;
    }
    if (coordinateChanged) {
      const title = coordinate === 'com10' ? '显示位置 x(0) + 10u（仅显示放大）' : coordinate === 'com' ? '去质心位置 x(0) + u' : coordinate === 'lab' ? '连续实验室坐标 x' : '位置 x（周期模式折回盒内）';
      const boundaries = coordinate === 'wrapped' && item.boundary === 'periodic' ? [0, item.metadata.box_length].map(value =>
        ({type: 'line', xref: 'x', yref: 'paper', x0: value, x1: value, y0: 0, y1: 1, line: {color: '#8090a4', dash: 'dot', width: 1}})) : [];
      updates.push(Plotly.relayout(panel.chain, {'xaxis.range': xrange, 'xaxis.title.text': title, shapes: boundaries}));
    }
    const thermal = item.thermal[k];
    $(`sample-${panel.index}`).textContent = `采样 t = ${actual.toFixed(4)} · step ${thermal[0]} · 箭头显示倍率 ${panel.speedScale.toFixed(2)}`;
    $(`readings-${panel.index}`).textContent = `粒子 ${particle + 1}：u = ${fmt(item.fields.displacement[k][particle])} · v′ = ${fmt(velocity[particle])} · K′ = ${fmt(item.fields.kinetic[k][particle])} | T* = ${fmt(thermal[5])} · ΔE = ${(thermal[4] - item.thermal[0][4]).toExponential(2)}`;
    await Promise.all(updates);
    panel.lastParticle = particle; panel.lastMetric = metric; panel.lastCoordinates = coordinate;
  }));
}

function render() {
  dirty = true;
  if (!paintTask) paintTask = (async () => {
    while (dirty) { dirty = false; await draw(); }
  })().finally(() => { paintTask = null; });
  return paintTask;
}

function play() {
  if (state.playing) { pause(); return; }
  if (state.frame >= Number($('frame').max)) state.frame = Number($('frame').min);
  state.playing = true; $('play').textContent = '暂停';
  const token = ++playToken, origin = DATA.time[state.frame], started = performance.now();
  const speed = Number($('speed').value), end = DATA.time[Number($('frame').max)];
  async function tick() {
    if (!state.playing || token !== playToken) return;
    const target = Math.min(end, origin + (performance.now() - started) * speed / 1000);
    state.frame = boundedFrame(nearest(DATA.time, target));
    await render();
    if (!state.playing || token !== playToken) return;
    if (target >= end) { pause(); return; }
    setTimeout(tick, 65);
  }
  setTimeout(tick, 0);
}

async function resizeCharts() {
  await Promise.all(panels.flatMap(panel => [panel.chain, panel.heat, panel.history, panel.thermo])
    .map(plot => Plotly.Plots.resize(plot)));
}

async function initialize() {
  const particleCount = Math.min(...DATA.cases.map(item => item.positions[0].length));
  for (let j = 0; j < particleCount; j++) $('particle').add(new Option(String(j + 1), String(j)));
  state.particle = Math.min(particleCount - 1, DATA.cases[0].metadata.initialization.kick_particle_index_zero_based);
  for (let index = 0; index < DATA.cases.length; index++) await makePanel(DATA.cases[index], index);
  // Adding the second card changes the first card from full width to half width.
  // A window resize event does not fire for this grid reflow.
  await resizeCharts();
  const widths = new WeakMap();
  sizeObserver = new ResizeObserver(entries => {
    let changed = false;
    for (const entry of entries) {
      const previous = widths.get(entry.target);
      widths.set(entry.target, entry.contentRect.width);
      changed ||= previous === undefined || Math.abs(previous - entry.contentRect.width) > 0.5;
    }
    if (changed) resizeCharts().catch(console.error);
  });
  panels.forEach(panel => sizeObserver.observe(panel.section));
  await applyRange(...state.range);
  ids.forEach(id => $(id).disabled = false);
  $('status').textContent = `已加载 ${DATA.cases.length} 组轨迹 · 所有图表共享时间与粒子选择 · 播放 1× = 每秒推进 1 个约化时间单位`;
  $('frame').addEventListener('input', () => { pause(); state.frame = boundedFrame(Number($('frame').value)); render(); });
  $('particle').addEventListener('change', () => selectParticle($('particle').value));
  $('metric').addEventListener('change', () => { state.metric = $('metric').value; render(); });
  $('coordinates').addEventListener('change', () => { state.coordinates = $('coordinates').value; render(); });
  $('window').addEventListener('change', () => applyRange(DATA.time[0], $('window').value === 'early' ? Math.min(DATA.time[0] + 2, DATA.time.at(-1)) : DATA.time.at(-1)));
  $('play').addEventListener('click', play);
  $('speed').addEventListener('change', () => { if (state.playing) { pause(); play(); } });
  $('reset').addEventListener('click', () => { pause(); state.frame = Number($('frame').min); render(); });
  document.addEventListener('visibilitychange', () => { if (document.hidden) pause(); });
  return true;
}

// Small inspection interface for automated end-to-end checks of the actual page.
window.mdDashboard = {
  getState: () => ({...state, range: [...state.range], caseFrames: panels.map(panel => Number(panel.section.dataset.frame))}),
  setFrame: async frame => { pause(); state.frame = boundedFrame(frame); await render(); },
  setParticle: selectParticle,
  setMetric: async metric => { if (!(metric in names)) throw Error('Unknown metric'); state.metric = metric; await render(); },
  setCoordinates: async value => { if (!['com10', 'com', 'lab', 'wrapped'].includes(value)) throw Error('Unknown coordinate mode'); state.coordinates = value; await render(); },
  setWindow: applyRange, play, pause,
  ready: initialize()
};
window.mdDashboard.ready.catch(error => {
  $('status').setAttribute('role', 'alert');
  $('status').textContent = `图表加载失败：${error.message}`;
  console.error(error);
});
