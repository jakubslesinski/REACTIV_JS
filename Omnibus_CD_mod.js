//Based on: https://github.com/mortcanty

//Define styling and determine the color of the shapefile  
var shp_styl = {
  color: 'yellow',
  fillColor: '00000000',
};

// var Czkalowsk_poligon = ee.FeatureCollection('users/slesinskijakub/Czkalowsk_poli');
// Map.addLayer(Czkalowsk_poligon.style(shp_styl));
var Belbek_poligon = ee.FeatureCollection('users/slesinskijakub/Belbek');
Map.addLayer(Belbek_poligon.style(shp_styl));

var geometry = MPS;

Map.setOptions('satellite');
Map.style().set('cursor', 'crosshair');

// ******************************************
// Front End for Sequential Omnibus Algorithm
// ******************************************
var omb = require('users/mortcanty/changedetection:omnibus');
var util = require('users/mortcanty/changedetection:utilities');
var jet = ['black','blue','cyan', 'yellow','red'];
// Input data 
var significance = 0.0001;   
var relorbitnumber = 'any'; 
var orbitpass = 'ASCENDING';
var median = true;
var startDate = '2019-07-15';
var endDate = '2019-11-01';
var assetExportId = 'users/mortcanty/omnibus/testxx';
var bnds = ee.List(geometry.bounds().coordinates().get(0));
var p0 = ee.Geometry.Point(bnds.get(0))
var p1 = ee.Geometry.Point(bnds.get(1))
var p2 = ee.Geometry.Point(bnds.get(2))
var p3 = ee.Geometry.Point(bnds.get(3))
var collection = ee.ImageCollection('COPERNICUS/S1_GRD') 
                   .filterBounds(p0) 
                   .filterBounds(p1)
                   .filterBounds(p2) 
                   .filterBounds(p3) 
                   .filterDate(ee.Date(startDate), ee.Date(endDate)) 
                   .filter(ee.Filter.eq('transmitterReceiverPolarisation', ['VV','VH'])) 
                   .filter(ee.Filter.eq('resolution_meters', 10)) 
                   .filter(ee.Filter.eq('instrumentMode', 'IW'))
                   .filter(ee.Filter.eq('orbitProperties_pass', orbitpass)); 
if (relorbitnumber != 'any'){
      var collection = collection.filter(ee.Filter.eq('relativeOrbitNumber_start', relorbitnumber));                        
} 
var collection = collection.sort('system:time_start');  
var acquisition_times = ee.List(collection.aggregate_array('system:time_start'));
var count = acquisition_times.length().getInfo();
if (count===0){ print('No images found')}
else{
  print('Timestamps');
  print(ee.List(acquisition_times.getInfo().map(function t(d){ return new Date(d)})));
  print('Relative Orbit Numbers');
  print(ee.List(collection.aggregate_array('relativeOrbitNumber_start'))); 
  var vis = {min:0, max:count, palette:jet};
  Map.add(util.makeLegend(vis));
// Create a list of clipped images
  var pList = collection.map(omb.get_vvvh).toList(count);
  var first = ee.Dictionary({imlist:ee.List([]),geom:geometry});
  var imList = ee.List(ee.Dictionary(pList.iterate(omb.clipList,first)).get('imlist')); 
// run the algorithm
  var result = ee.Dictionary(omb.omnibus(imList,significance,median));
// get change maps 
  var cmap = ee.Image(result.get('cmap')).byte();
  var smap = ee.Image(result.get('smap')).byte();
  var fmap = ee.Image(result.get('fmap')).byte();
  var bmap = ee.Image(result.get('bmap')).byte();
  var cnames = ['cmap','smap','fmap'];
  for (var i = 1; i < count; i++){
    if (i < 10) {var label = 'bmap0'} else {var label= 'bmap'}
    cnames = cnames.concat([label.concat(i.toString())]);
  }
  cnames = cnames.concat(['background']);
// background image for video export
  var background =collection.mean()
                            .select(0)
                            .multiply(ee.Image.constant(Math.log(10.0)/10.0)).exp();
  background = background.where(background.gte(1),1).clip(geometry);
// concatenate change maps and export
  var cmaps = ee.Image.cat(cmap,smap,fmap,bmap,background).rename(cnames);
  var exportTask = Export.image.toAsset({image:cmaps,assetId:assetExportId,scale:10,maxPixels:1e9});
// display change maps. NOTE: inaccurate because pyramid levels have different ENLs
  Map.centerObject(cmap,15);
  Map.addLayer(cmaps.select('cmap'),vis,'cmap', false, 0.6);
  Map.addLayer(cmaps.select('smap'),vis,'smap', true, 0.6);
  Map.addLayer(cmaps.select('fmap').multiply(2),vis,'fmap*2', false, 0.6); 
}