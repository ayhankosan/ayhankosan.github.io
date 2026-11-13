/***********************************************************************
 * Palisades Fire Burn Severity Analysis - Los Angeles, California (2025)
 *
 * Purpose:
 *   Maps burn severity with the differenced Normalized Burn Ratio (dNBR),
 *   offset-corrected against unburned surroundings, classified with the
 *   Key & Benson (2006) / USGS thresholds. Adds the Relativized Burn Ratio
 *   (RBR) and Sentinel-1 radar change layers for a multi-sensor view.
 *
 * Data:
 *   - Sentinel-2 SR Harmonized          COPERNICUS/S2_SR_HARMONIZED
 *   - Cloud Score+ for Sentinel-2        GOOGLE/CLOUD_SCORE_PLUS/V1/S2_HARMONIZED
 *   - Sentinel-1 GRD, IW, VV+VH (dB)     COPERNICUS/S1_GRD
 *   - JRC Global Surface Water v1.3      JRC/GSW1_3/GlobalSurfaceWater
 *   - Optional: fire perimeter imported as a FeatureCollection named
 *     `firePerimeter` (e.g. from the NIFC/WFIGS interagency perimeters)
 *
 * References:
 *   - Key, C.H. & Benson, N.C. (2006). Landscape Assessment (LA). FIREMON,
 *     USDA Forest Service RMRS-GTR-164-CD.
 *   - Parks, S.A. et al. (2014). A new metric for quantifying burn severity:
 *     the Relativized Burn Ratio. Remote Sensing 6(3), 1827-1844.
 *   - UN-SPIDER Recommended Practice: Burn Severity (NBR)
 *
 * Author: Ayhan Kosan
 ***********************************************************************/

// -----------------------------------------------------------------------
// 1. CONFIGURATION
// -----------------------------------------------------------------------
var CONFIG = {
  // Fire ignition 2025-01-07; containment 2025-01-31
  preFireStart: '2024-12-01',
  preFireEnd: '2025-01-07',      // filterDate end is exclusive: last day used is 2025-01-06
  postFireStart: '2025-02-01',
  postFireEnd: '2025-02-13',     // ends before the 2025-02-13/14 storm, which wets soils and alters radar backscatter

  csThreshold: 0.60,             // Cloud Score+ cs_cdf: pixels below this value are treated as cloud or shadow
  s1Pass: 'DESCENDING',          // pre and post radar scenes share orbit direction and relative orbit
  speckleRadius: 30,             // metres, focal mean applied in linear power

  opticalScale: 20,
  sarScale: 10,
  maxPixels: 1e10,
  tileScale: 4,
  exportFolder: 'GEE_Palisades_Outputs',

  waterOccurrencePct: 70,        // JRC occurrence at or above this value is masked as permanent water
  offsetRing: [500, 2500],       // metres outside the perimeter used as the unburned reference ring

  // dNBR class breaks (scaled by 1000), Key & Benson (2006)
  burnThresholds: {
    regrowthHigh: -250,
    regrowthLow: -100,
    unburnedUpper: 100,
    lowUpper: 270,
    moderateLowUpper: 440,
    moderateHighUpper: 660
  }
};

// Study area: drawn geometry when present, otherwise a box around the burn
var AOI = (typeof geometry !== 'undefined')
  ? geometry
  : ee.Geometry.Rectangle([-118.70, 34.00, -118.48, 34.13]);

// Statistics region: fire perimeter when imported, otherwise the AOI
var HAS_PERIMETER = (typeof firePerimeter !== 'undefined');
var STATS_REGION = HAS_PERIMETER ? firePerimeter.geometry() : AOI;


// -----------------------------------------------------------------------
// 2. MASKS AND HELPERS
// -----------------------------------------------------------------------
// Permanent inland/coastal water from JRC; open ocean is removed per scene via SCL class 6
var permanentWater = ee.Image('JRC/GSW1_3/GlobalSurfaceWater').select('occurrence')
  .unmask(0).gte(CONFIG.waterOccurrencePct).rename('permanent_water');
var landMask = permanentWater.not();

// Sentinel-2: Cloud Score+ plus Scene Classification Layer (no data, saturated,
// shadow, water, cloud, cirrus, snow); reflectance scaled to 0-1
function maskS2(image) {
  var scl = image.select('SCL');
  var sclOk = scl.neq(0).and(scl.neq(1)).and(scl.neq(3)).and(scl.neq(6))
    .and(scl.neq(8)).and(scl.neq(9)).and(scl.neq(10)).and(scl.neq(11));
  var clear = image.select('cs_cdf').gte(CONFIG.csThreshold);
  return image.select(['B2', 'B3', 'B4', 'B8', 'B12'])
    .divide(10000)
    .updateMask(sclOk.and(clear))
    .copyProperties(image, ['system:time_start']);
}

function addIndices(image) {
  var nbr = image.normalizedDifference(['B8', 'B12']).rename('NBR');
  var ndvi = image.normalizedDifference(['B8', 'B4']).rename('NDVI');
  var evi = image.expression('2.5 * ((NIR - RED) / (NIR + 6 * RED - 7.5 * BLUE + 1))', {
    NIR: image.select('B8'), RED: image.select('B4'), BLUE: image.select('B2')
  }).rename('EVI');
  return image.addBands([nbr, ndvi, evi]);
}

function s2Composite(start, end) {
  var cs = ee.ImageCollection('GOOGLE/CLOUD_SCORE_PLUS/V1/S2_HARMONIZED');
  var col = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
    .filterBounds(AOI)
    .filterDate(start, end)
    .linkCollection(cs, ['cs_cdf'])
    .map(maskS2)
    .map(addIndices);
  print('Sentinel-2 scenes ' + start + ' to ' + end + ':', col.size());
  return col.median().clip(AOI).updateMask(landMask);
}

// Sentinel-1 GRD in Earth Engine is already in dB. Speckle filtering and
// averaging run in linear power, then convert back to dB.
function toLinear(db) { return ee.Image(10).pow(db.divide(10)); }
function toDb(lin) { return lin.log10().multiply(10); }

function s1Collection(start, end) {
  return ee.ImageCollection('COPERNICUS/S1_GRD')
    .filterBounds(AOI)
    .filterDate(start, end)
    .filter(ee.Filter.eq('instrumentMode', 'IW'))
    .filter(ee.Filter.eq('orbitProperties_pass', CONFIG.s1Pass))
    .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
    .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))
    .select(['VV', 'VH']);
}

function s1Composite(col) {
  var lin = col.map(function (img) {
    return toLinear(img).focalMean({radius: CONFIG.speckleRadius, units: 'meters'})
      .copyProperties(img, ['system:time_start']);
  }).mean();
  var vv = toDb(lin.select('VV')).rename('VV');
  var vh = toDb(lin.select('VH')).rename('VH');
  var ratio = vh.subtract(vv).rename('VH_VV_dB');                          // ratio in dB is a difference
  var rvi = lin.expression('4 * VH / (VV + VH)', {VV: lin.select('VV'), VH: lin.select('VH')})
    .rename('RVI');                                                         // radar vegetation index, linear power
  return ee.Image.cat([vv, vh, ratio, rvi]).clip(AOI).updateMask(landMask);
}

// Area in hectares per class inside a region
function classAreaHa(classImage, region) {
  return ee.Image.pixelArea().divide(10000).addBands(classImage)
    .reduceRegion({
      reducer: ee.Reducer.sum().group({groupField: 1, groupName: 'class'}),
      geometry: region, scale: CONFIG.opticalScale,
      maxPixels: CONFIG.maxPixels, tileScale: CONFIG.tileScale
    });
}

// -----------------------------------------------------------------------
// 3. OPTICAL: NBR, dNBR (offset-corrected), RBR, NDVI, EVI
// -----------------------------------------------------------------------
var pre = s2Composite(CONFIG.preFireStart, CONFIG.preFireEnd);
var post = s2Composite(CONFIG.postFireStart, CONFIG.postFireEnd);

var dNBRraw = pre.select('NBR').subtract(post.select('NBR')).multiply(1000).rename('dNBR_raw');

// Offset: median dNBR of an unburned ring outside the perimeter (Key & Benson 2006).
// Removes phenology and illumination differences between the two composites.
var offset = ee.Number(0);
if (HAS_PERIMETER) {
  var ring = firePerimeter.geometry().buffer(CONFIG.offsetRing[1])
    .difference(firePerimeter.geometry().buffer(CONFIG.offsetRing[0]));
  offset = ee.Number(dNBRraw.reduceRegion({
    reducer: ee.Reducer.median(), geometry: ring, scale: CONFIG.opticalScale,
    maxPixels: CONFIG.maxPixels, tileScale: CONFIG.tileScale
  }).get('dNBR_raw'));
  print('dNBR offset from unburned ring (scaled):', offset);
}
var dNBR = dNBRraw.subtract(offset).rename('dNBR');

// Relativized Burn Ratio: dNBR relative to pre-fire NBR (Parks et al. 2014)
var RBR = dNBR.divide(pre.select('NBR').add(1.001)).rename('RBR');

var ndviDiff = pre.select('NDVI').subtract(post.select('NDVI')).rename('NDVI_loss');
var eviDiff = pre.select('EVI').subtract(post.select('EVI')).rename('EVI_loss');

// -----------------------------------------------------------------------
// 4. RADAR: backscatter and RVI change, same orbit geometry pre and post
// -----------------------------------------------------------------------
var postS1Raw = s1Collection(CONFIG.postFireStart, CONFIG.postFireEnd);
var relOrbit = ee.Number(postS1Raw.first().get('relativeOrbitNumber_start'));
var preS1Raw = s1Collection(CONFIG.preFireStart, CONFIG.preFireEnd)
  .filter(ee.Filter.eq('relativeOrbitNumber_start', relOrbit));
postS1Raw = postS1Raw.filter(ee.Filter.eq('relativeOrbitNumber_start', relOrbit));
print('Sentinel-1 relative orbit used:', relOrbit);
print('Sentinel-1 scenes pre / post:', preS1Raw.size(), postS1Raw.size());

var preSar = s1Composite(preS1Raw);
var postSar = s1Composite(postS1Raw);
// Sign convention: post minus pre, so negative values mean backscatter or RVI loss
var sarChange = postSar.subtract(preSar)
  .rename(['VV_change_dB', 'VH_change_dB', 'VH_VV_change_dB', 'RVI_change']);

// -----------------------------------------------------------------------
// 5. BURN SEVERITY CLASSES (Key & Benson 2006)
// -----------------------------------------------------------------------
// 0 Enhanced regrowth high  (< -250)
// 1 Enhanced regrowth low   (-250 to -101)
// 2 Unburned                (-100 to 99)
// 3 Low                     (100 to 269)
// 4 Moderate-low            (270 to 439)
// 5 Moderate-high           (440 to 659)
// 6 High                    (>= 660)
var t = CONFIG.burnThresholds;
var burnClass = ee.Image(2)
  .where(dNBR.lt(t.regrowthHigh), 0)
  .where(dNBR.gte(t.regrowthHigh).and(dNBR.lt(t.regrowthLow)), 1)
  .where(dNBR.gte(t.unburnedUpper).and(dNBR.lt(t.lowUpper)), 3)
  .where(dNBR.gte(t.lowUpper).and(dNBR.lt(t.moderateLowUpper)), 4)
  .where(dNBR.gte(t.moderateLowUpper).and(dNBR.lt(t.moderateHighUpper)), 5)
  .where(dNBR.gte(t.moderateHighUpper), 6)
  .updateMask(dNBR.mask())
  .rename('burn_class')
  .toByte();
var CLASS_NAMES = ['Regrowth high', 'Regrowth low', 'Unburned', 'Low', 'Moderate-low', 'Moderate-high', 'High'];

var moderateToHigh = burnClass.gte(4).selfMask().rename('moderate_to_high');

// -----------------------------------------------------------------------
// 6. STATISTICS
// -----------------------------------------------------------------------
var areaStats = classAreaHa(burnClass, STATS_REGION);
print('Area by burn class (ha), ' + (HAS_PERIMETER ? 'inside fire perimeter' : 'inside AOI') + ':', areaStats);

var meanChange = ee.Image.cat([dNBR, RBR, ndviDiff, eviDiff]).reduceRegion({
  reducer: ee.Reducer.mean(), geometry: STATS_REGION, scale: CONFIG.opticalScale,
  maxPixels: CONFIG.maxPixels, tileScale: CONFIG.tileScale
});
var meanSar = sarChange.reduceRegion({
  reducer: ee.Reducer.mean(), geometry: STATS_REGION, scale: CONFIG.sarScale,
  maxPixels: CONFIG.maxPixels, tileScale: CONFIG.tileScale
});
print('Mean optical change:', meanChange);
print('Mean radar change (post minus pre):', meanSar);


// -----------------------------------------------------------------------
// 7. VECTORS AND EXPORTS (run from the Tasks tab)
// -----------------------------------------------------------------------
var burnPolygons = moderateToHigh.toByte().reduceToVectors({
  geometry: STATS_REGION, scale: CONFIG.opticalScale, geometryType: 'polygon',
  labelProperty: 'moderate_to_high', maxPixels: CONFIG.maxPixels, tileScale: CONFIG.tileScale
});

function exportImage(image, name, scale) {
  Export.image.toDrive({
    image: image, description: name, folder: CONFIG.exportFolder, fileNamePrefix: name.toLowerCase(),
    region: AOI, scale: scale, crs: 'EPSG:32611', maxPixels: CONFIG.maxPixels
  });
}
exportImage(ee.Image.cat([dNBRraw, dNBR, RBR]).toFloat(), 'Palisades_2025_dNBR_RBR', CONFIG.opticalScale);
exportImage(burnClass, 'Palisades_2025_burn_class', CONFIG.opticalScale);
exportImage(ee.Image.cat([ndviDiff, eviDiff]).toFloat(), 'Palisades_2025_vegetation_loss', CONFIG.opticalScale);
exportImage(sarChange.toFloat(), 'Palisades_2025_radar_change', CONFIG.sarScale);
Export.table.toDrive({
  collection: burnPolygons, description: 'Palisades_2025_moderate_to_high_polygons',
  folder: CONFIG.exportFolder, fileFormat: 'SHP'
});
Export.table.toDrive({
  collection: ee.FeatureCollection(ee.List(areaStats.get('groups')).map(function (g) {
    g = ee.Dictionary(g);
    return ee.Feature(null, {class_id: g.get('class'), class_name: ee.List(CLASS_NAMES).get(g.get('class')), area_ha: g.get('sum')});
  })),
  description: 'Palisades_2025_class_area', folder: CONFIG.exportFolder, fileFormat: 'CSV'
});

// -----------------------------------------------------------------------
// 8. APP: control panel, analysis map, before/after swipe, inspector
// -----------------------------------------------------------------------
// Colour-blind-safe palettes (ColorBrewer). Diverging scales are centred on zero.
var CLASS_COLORS = ['#2166ac', '#92c5de', '#ececec', '#fed976', '#fd8d3c', '#e31a1c', '#800026'];
var DIV = ['#2166ac', '#4393c3', '#92c5de', '#d1e5f0', '#f7f7f7', '#fddbc7', '#f4a582', '#d6604d', '#b2182b'];
var DIV_REV = DIV.slice().reverse();   // radar: post minus pre, losses in red

// Light, label-only basemap so the analysis layers carry the colour
var LIGHT_STYLE = [
  {elementType: 'geometry', stylers: [{color: '#f4f3ef'}]},
  {elementType: 'labels.icon', stylers: [{visibility: 'off'}]},
  {elementType: 'labels.text.fill', stylers: [{color: '#5f5f5f'}]},
  {elementType: 'labels.text.stroke', stylers: [{color: '#f4f3ef'}]},
  {featureType: 'poi', stylers: [{visibility: 'off'}]},
  {featureType: 'transit', stylers: [{visibility: 'off'}]},
  {featureType: 'road', elementType: 'geometry', stylers: [{color: '#ffffff'}]},
  {featureType: 'road.highway', elementType: 'geometry', stylers: [{color: '#e3e1da'}]},
  {featureType: 'water', elementType: 'geometry', stylers: [{color: '#cfdbe4'}]}
];

// Terrain context: 1/3 arc-second (~10 m) USGS 3DEP hillshade
var hillshade = ee.Terrain.hillshade(ee.Image('USGS/3DEP/10m').select('elevation')).clip(AOI);

var LAYERS = {
  'Burn severity (classes)': {image: burnClass, vis: {min: 0, max: 6, palette: CLASS_COLORS}, legend: 'classes'},
  'dNBR, offset-corrected': {image: dNBR, vis: {min: -800, max: 800, palette: DIV}, legend: 'bar', units: 'dNBR × 1000'},
  'Relativized Burn Ratio (RBR)': {image: RBR, vis: {min: -600, max: 600, palette: DIV}, legend: 'bar', units: 'RBR × 1000'},
  'NDVI loss': {image: ndviDiff, vis: {min: -0.5, max: 0.5, palette: DIV}, legend: 'bar', units: 'pre minus post'},
  'EVI loss': {image: eviDiff, vis: {min: -0.4, max: 0.4, palette: DIV}, legend: 'bar', units: 'pre minus post'},
  'Radar VH change': {image: sarChange.select('VH_change_dB'), vis: {min: -3, max: 3, palette: DIV_REV}, legend: 'bar', units: 'dB, post minus pre'},
  'Radar vegetation index change': {image: sarChange.select('RVI_change'), vis: {min: -0.3, max: 0.3, palette: DIV_REV}, legend: 'bar', units: 'post minus pre'},
  'Post-fire SWIR false color': {image: post, vis: {bands: ['B12', 'B8', 'B4'], min: 0, max: 0.4, gamma: 1.1}, legend: 'none'}
};
var COMPOSITES = {
  'True color': {bands: ['B4', 'B3', 'B2'], min: 0, max: 0.3, gamma: 1.3},
  'SWIR false color (burn scars dark red)': {bands: ['B12', 'B8', 'B4'], min: 0, max: 0.4, gamma: 1.1}
};

var state = {layer: 'Burn severity (classes)', opacity: 0.85, hideUnburned: true, terrain: true, composite: 'True color'};

// ---- Maps
function styledMap() {
  var m = ui.Map();
  m.setOptions('Light', {Light: LIGHT_STYLE}, ['Light', 'SATELLITE', 'HYBRID']);
  m.setControlVisibility({layerList: false, drawingToolsControl: false, fullscreenControl: false, mapTypeControl: true});
  return m;
}
var analysisMap = styledMap();
var preMap = styledMap();
var postMap = styledMap();
ui.Map.Linker([analysisMap, preMap, postMap]);
analysisMap.style().set('cursor', 'crosshair');

var perimeterLayer = HAS_PERIMETER
  ? ee.FeatureCollection(firePerimeter).style({color: '111111', fillColor: '00000000', width: 2})
  : ee.FeatureCollection([ee.Feature(AOI)]).style({color: '111111', fillColor: '00000000', width: 1, lineType: 'dashed'});

function mapLabel(text) {
  return ui.Label(text, {position: 'top-left', fontWeight: 'bold', fontSize: '13px', padding: '4px 8px', backgroundColor: 'rgba(255,255,255,0.9)'});
}
preMap.add(mapLabel('Before · Dec 1, 2024 – Jan 6, 2025'));
postMap.add(mapLabel('After · Feb 1 – Feb 12, 2025'));
var swipe = ui.SplitPanel({firstPanel: preMap, secondPanel: postMap, wipe: true, style: {stretch: 'both'}});

function drawSwipe() {
  var vis = COMPOSITES[state.composite];
  preMap.layers().reset([ui.Map.Layer(pre, vis, 'Before'), ui.Map.Layer(perimeterLayer, {}, 'Perimeter')]);
  postMap.layers().reset([ui.Map.Layer(post, vis, 'After'), ui.Map.Layer(perimeterLayer, {}, 'Perimeter')]);
}

function drawAnalysis() {
  var def = LAYERS[state.layer];
  var img = def.image;
  if (state.layer === 'Burn severity (classes)' && state.hideUnburned) {
    img = burnClass.updateMask(burnClass.neq(2));
  }
  analysisMap.layers().reset([
    ui.Map.Layer(hillshade, {min: 120, max: 255}, 'Hillshade', state.terrain, 0.55),
    ui.Map.Layer(img, def.vis, state.layer, true, state.opacity),
    ui.Map.Layer(perimeterLayer, {}, 'Perimeter')
  ]);
  drawLegend();
}

// ---- Legend
var legendPanel = ui.Panel();
function colorBar(palette, min, max, units) {
  var bar = ui.Thumbnail({
    image: ee.Image.pixelLonLat().select(0),
    params: {bbox: [0, 0, 1, 0.1], dimensions: '300x12', format: 'png', min: 0, max: 1, palette: palette},
    style: {stretch: 'horizontal', margin: '0', maxHeight: '14px'}
  });
  var ticks = ui.Panel([
    ui.Label(String(min), {fontSize: '11px', margin: '2px 0'}),
    ui.Label(units || '', {fontSize: '11px', color: '#777', margin: '2px 0', textAlign: 'center', stretch: 'horizontal'}),
    ui.Label(String(max), {fontSize: '11px', margin: '2px 0'})
  ], ui.Panel.Layout.flow('horizontal'));
  return ui.Panel([bar, ticks]);
}
function drawLegend() {
  var def = LAYERS[state.layer];
  legendPanel.clear();
  if (def.legend === 'classes') {
    CLASS_NAMES.forEach(function (name, i) {
      if (state.hideUnburned && i === 2) { return; }
      legendPanel.add(ui.Panel([
        ui.Label('', {backgroundColor: CLASS_COLORS[i], padding: '7px', margin: '2px 8px 2px 0', border: '1px solid #bbb'}),
        ui.Label(name + '  ' + BREAKS[i], {fontSize: '12px', margin: '2px 0'})
      ], ui.Panel.Layout.flow('horizontal')));
    });
  } else if (def.legend === 'bar') {
    legendPanel.add(colorBar(def.vis.palette, def.vis.min, def.vis.max, def.units));
  } else {
    legendPanel.add(ui.Label('Bands 12-8-4: burned ground appears dark red-brown, vegetation green.', {fontSize: '12px', color: '#555'}));
  }
}
var BREAKS = ['< −250', '−250 to −101', '−100 to 99', '100 to 269', '270 to 439', '440 to 659', '≥ 660'];

// ---- Control panel
var H = {fontSize: '11px', fontWeight: 'bold', color: '#555', margin: '16px 0 6px 0'};
var panel = ui.Panel({style: {width: '360px', padding: '14px 16px'}});
panel.add(ui.Label('Palisades Fire 2025', {fontSize: '22px', fontWeight: 'bold', margin: '0 0 2px 0'}));
panel.add(ui.Label('Burn severity from Sentinel-2 and Sentinel-1', {fontSize: '14px', color: '#444', margin: '0 0 10px 0'}));
panel.add(ui.Label('Offset-corrected dNBR classified with Key & Benson (2006) thresholds, the Relativized Burn Ratio, vegetation index loss and radar change. ' +
  (HAS_PERIMETER ? 'Statistics inside the fire perimeter.' : 'Statistics inside the study area.'), {fontSize: '12px', color: '#555'}));

// KPI tiles
function kpi(label) {
  var value = ui.Label('…', {fontSize: '20px', fontWeight: 'bold', margin: '0'});
  var tile = ui.Panel([value, ui.Label(label, {fontSize: '11px', color: '#666', margin: '2px 0 0 0'})], null,
    {width: '100px', padding: '8px', margin: '0 6px 6px 0', border: '1px solid #ddd'});
  return {tile: tile, value: value};
}
var kBurned = kpi('ha burned, low to high');
var kMod = kpi('ha moderate to high');
var kHigh = kpi('ha high severity');
panel.add(ui.Panel([kBurned.tile, kMod.tile, kHigh.tile], ui.Panel.Layout.flow('horizontal')));

panel.add(ui.Label('VIEW', H));
var modeSelect = ui.Select({
  items: ['Analysis map', 'Before / after swipe'], value: 'Analysis map',
  onChange: function (v) { showMode(v); }, style: {stretch: 'horizontal'}
});
panel.add(modeSelect);

var analysisControls = ui.Panel();
analysisControls.add(ui.Label('LAYER', H));
analysisControls.add(ui.Select({
  items: Object.keys(LAYERS), value: state.layer, style: {stretch: 'horizontal'},
  onChange: function (v) { state.layer = v; drawAnalysis(); }
}));
analysisControls.add(ui.Panel([
  ui.Label('Opacity', {fontSize: '12px', margin: '8px 8px 0 0'}),
  ui.Slider({min: 0, max: 1, value: state.opacity, step: 0.05, style: {stretch: 'horizontal'},
    onChange: function (v) { state.opacity = v; analysisMap.layers().get(1).setOpacity(v); }})
], ui.Panel.Layout.flow('horizontal'), {stretch: 'horizontal'}));
analysisControls.add(ui.Checkbox({label: 'Hide unburned pixels', value: state.hideUnburned, style: {fontSize: '12px'},
  onChange: function (v) { state.hideUnburned = v; drawAnalysis(); }}));
analysisControls.add(ui.Checkbox({label: 'Terrain hillshade', value: state.terrain, style: {fontSize: '12px'},
  onChange: function (v) { state.terrain = v; analysisMap.layers().get(0).setShown(v); }}));
analysisControls.add(ui.Label('LEGEND', H));
analysisControls.add(legendPanel);
panel.add(analysisControls);

var swipeControls = ui.Panel([
  ui.Label('COMPOSITE', H),
  ui.Select({items: Object.keys(COMPOSITES), value: state.composite, style: {stretch: 'horizontal'},
    onChange: function (v) { state.composite = v; drawSwipe(); }}),
  ui.Label('Drag the divider to compare before and after.', {fontSize: '12px', color: '#555'})
]);
swipeControls.style().set('shown', false);
panel.add(swipeControls);

// Area chart (one bar per class, class colours)
panel.add(ui.Label('AREA BY CLASS', H));
var chartHolder = ui.Panel([ui.Label('Calculating…', {fontSize: '12px', color: '#777'})]);
panel.add(chartHolder);

// Point inspector
panel.add(ui.Label('INSPECTOR', H));
var inspector = ui.Label('Click the map to read values at a point.', {fontSize: '12px', color: '#555', whiteSpace: 'pre'});
panel.add(inspector);

panel.add(ui.Label('SOURCES', H));
panel.add(ui.Label('Sentinel-2 SR Harmonized with Cloud Score+; Sentinel-1 GRD; JRC Global Surface Water; USGS 3DEP 10 m.', {fontSize: '11px', color: '#666'}));
panel.add(ui.Label('UN-SPIDER Recommended Practice: Burn Severity', {fontSize: '11px'},
  'https://un-spider.org/advisory-support/recommended-practices/recommended-practice-burn-severity'));
panel.add(ui.Label('Analysis: Ayhan Kosan', {fontSize: '11px', color: '#666', margin: '10px 0 0 0'}));

// ---- Mode switching
var mapHolder = ui.Panel({widgets: [analysisMap], style: {stretch: 'both'}});
function showMode(mode) {
  var swipeMode = (mode === 'Before / after swipe');
  analysisControls.style().set('shown', !swipeMode);
  swipeControls.style().set('shown', swipeMode);
  if (swipeMode) { drawSwipe(); }
  mapHolder.widgets().reset([swipeMode ? swipe : analysisMap]);
}

// ---- Inspector
analysisMap.onClick(function (c) {
  var pt = ee.Geometry.Point([c.lon, c.lat]);
  var marker = ui.Map.Layer(pt, {color: '111111'}, 'Selected point');
  if (analysisMap.layers().length() > 3) { analysisMap.layers().set(3, marker); } else { analysisMap.layers().add(marker); }
  inspector.setValue('Reading values…');
  ee.Image.cat([burnClass, dNBR, RBR, ndviDiff, sarChange.select('VH_change_dB')])
    .reduceRegion({reducer: ee.Reducer.first(), geometry: pt, scale: 10})
    .evaluate(function (v) {
      if (!v || v.burn_class === null || v.burn_class === undefined) { inspector.setValue('No data at this point (water, cloud or outside the study area).'); return; }
      var f = function (x, d) { return (x === null || x === undefined) ? 'n/a' : Number(x).toFixed(d); };
      inspector.setValue(
        'Lat ' + c.lat.toFixed(5) + ', Lon ' + c.lon.toFixed(5) + '\n' +
        'Class: ' + CLASS_NAMES[v.burn_class] + '\n' +
        'dNBR: ' + f(v.dNBR, 0) + '   RBR: ' + f(v.RBR, 0) + '\n' +
        'NDVI loss: ' + f(v.NDVI_loss, 2) + '\n' +
        'Radar VH change: ' + f(v.VH_change_dB, 2) + ' dB');
    });
});

// ---- Statistics into KPIs and chart
areaStats.evaluate(function (res) {
  var groups = (res && res.groups) ? res.groups : [];
  var byClass = {};
  groups.forEach(function (g) { byClass[g['class']] = g.sum; });
  var sum = function (ids) { return ids.reduce(function (a, i) { return a + (byClass[i] || 0); }, 0); };
  var fmt = function (x) { return Math.round(x).toLocaleString('en-US'); };
  kBurned.value.setValue(fmt(sum([3, 4, 5, 6])));
  kMod.value.setValue(fmt(sum([4, 5, 6])));
  kHigh.value.setValue(fmt(sum([6])));
  var table = [['Class', 'Hectares', {role: 'style'}]];
  CLASS_NAMES.forEach(function (name, i) { table.push([name, byClass[i] || 0, CLASS_COLORS[i]]); });
  chartHolder.clear();
  chartHolder.add(ui.Chart(table, 'BarChart', {
    legend: {position: 'none'}, hAxis: {title: 'hectares', textStyle: {fontSize: 10}}, vAxis: {textStyle: {fontSize: 11}},
    chartArea: {left: 95, top: 6, width: '62%', height: '85%'}, height: 210, bar: {groupWidth: '75%'}
  }));
});

// ---- Layout
ui.root.clear();
ui.root.setLayout(ui.Panel.Layout.flow('horizontal'));
ui.root.add(panel);
ui.root.add(mapHolder);
drawAnalysis();
analysisMap.centerObject(STATS_REGION, 13);
