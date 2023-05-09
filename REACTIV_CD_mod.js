// Based on:
// Elise Colin Koeniguer et al, 
// Colored visualization of multitemporal SAR data for change detection: issues and methods
// EUSAR 2018
// 
// &
//
// Visualisation des changements sur séries temporelles radar : méthode REACTIV 
// évaluée à l’échelle mondiale sous Google Earth Engine
// Elise Colin Koeniguer et al. 
// CFPT 2018
// https://rfiap2018.ign.fr/sites/default/files/ARTICLES/CFPT2018/Oraux/CFPT2018_paper_koeniguer.pdf
// -------------------------------------------------------------


//Define styling and determine the color of the shapefile  
var shp_styl = {
  color: 'yellow',
  fillColor: '00000000',
};

var Czkalowsk_poligon = ee.FeatureCollection('users/slesinskijakub/Czkalowsk_poli');
Map.addLayer(Czkalowsk_poligon.style(shp_styl));

var Belbek_poligon = ee.FeatureCollection('users/slesinskijakub/Belbek');
Map.addLayer(Belbek_poligon.style(shp_styl));


// -------------------------------------------------------------
// In this version, the method is applied on one unique orbit.
// The chosen orbit is the most frequent one over the Center of the Map 

// ------------------------------------------------------------
// Parameters: DATES, ASCENDING OR DESCENDING, POLARISATION VH OR VV
var str2='2019-11-01';
var str1='2019-07-01';
var str='ASCENDING';
var polar='VV';

// ------------------------------------------------------------
// date selection
var date2 = ee.Date(str2);
var date1 = ee.Date(str1);
var ds = date2.difference(date1, 'day')

// Centering
var pos = Belbek;
var pos = Map.getCenter();
print('Coordinate of the Center of the Map',pos);
Map.setCenter (33.5800, 44.6933, 15);

// Load the Sentinel-1 ImageCollection centered on the location "pos"
// Necessity to have the stack centered on this location in order to find next the orbit numbers
var sentinel1_liste = ee.ImageCollection('COPERNICUS/S1_GRD')
  .filterDate(date1, date2)
  .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
  .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))
  .filter(ee.Filter.eq('instrumentMode', 'IW'))
  .filterBounds(pos)
  .filter(ee.Filter.eq('orbitProperties_pass', str));

// sentinel collection of the world without the restriction of the position
var sentinel1_liste2 = ee.ImageCollection('COPERNICUS/S1_GRD')
  .filterDate(date1, date2)
  .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
  .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))
  .filter(ee.Filter.eq('instrumentMode', 'IW'))
  .filter(ee.Filter.eq('orbitProperties_pass', str));



// a solution to get metadata value of images of a collection
var NbOrbit = sentinel1_liste.aggregate_count_distinct('relativeOrbitNumber_start');
print('Number of Orbits',NbOrbit);
var ListOrbits = sentinel1_liste.aggregate_array('relativeOrbitNumber_start');
print('ListOrbits:', ListOrbits);


// find orbit numbers and their frequency
var freq = ee.Dictionary(ListOrbits.reduce(ee.Reducer.frequencyHistogram()));  // sortuje orbity wg. czestosci
 print('freq',freq.values());
//print('test', freq.keys())
// var array=ee.Array([freq.keys()['0']*1,freq.keys()[1]*1, freq.values()])
var array = ee.Array([freq.keys().map(function(x) {return ee.Number.parse(x)}), freq.values()]);
 print('array',array);

// orbit choice : first, the one with the max frequency
var frequences = array.slice(0,-1);
var arraysort = array.sort(frequences);
var index = NbOrbit.add(-1);
var orbite = arraysort.get([0,index]);
print('Selected orbit=',orbite);


// find images with the choice orbit
var sentinel1 = sentinel1_liste2.filterMetadata('relativeOrbitNumber_start', 'equals', orbite);

var sentinel1_lista = sentinel1.toList(sentinel1.size());
print('sentinel1_lista: ', sentinel1_lista)

/////////////////////////////////////////
// Optional: Mask for oceans
//var elev = ee.Image('USGS/GMTED2010');
//var ocean = elev.lte(0);


// This function applies to each image the linear scale
var linear = function(image) {
  var imlin = image.expression(
    '10**(amplitude/20)', {
      'amplitude': image.select(polar)
  });
  return imlin; // conversion in linear, then compute mean: classical mean
};

var stdLinear = sentinel1.select(polar).map(linear).reduce(ee.Reducer.stdDev()); // ee.Reducer.stdDev() - returns a Reducer that computes the standard deviation of its inputs
var meanLinear = sentinel1.select(polar).map(linear).reduce(ee.Reducer.mean()); //ee.Reducer.mean() - returns a Reducer that computes the (weighted) arithmetic mean of its inputs. 
var CV = stdLinear.divide(meanLinear);

print(meanLinear.bandNames());

var imagemax = sentinel1.select(polar).max();
var imax = imagemax.expression(
    '10 ** (amplitude/20)', {  //x **= y  ==  x = x * y
      'amplitude': imagemax.select(polar)
});



// This function affects value of days for pixels where maximum is reached
var time = function(image) {
  var days = image.date().difference(date1, 'day').multiply(0.85).divide(ds); //divide by the period of time observed
  return image.where(image.lt(imagemax),0).where(image.gte(imagemax),days);
};
var days=sentinel1.select(polar).map(time).sum();  // sumowanie w funkcji do pojedynczego kanału dla całego szeregu czasowego
//print('Days= ',days);

// Images of Number of images: sizepile
var unit = function(image) {
  var imunit = image.multiply(0).add(1);
  return imunit; // conversion in linear, then compute mean: classical mean
};
var sizepile=sentinel1.select(polar).map(unit).sum(); 
print('sizepile= ',sizepile);

// Parameter for dynamics
var mu=0.2286; // Theoretical mean for Rayleigh Nakagam L=4.9
var stdmu=ee.Image(0.1616);
var stdmu=stdmu.divide(sizepile.sqrt()); // Theoretical std for Rayleigh Nakagami L=4.9
var CV_norm=CV.subtract(mu).divide(stdmu.multiply(10)).add(0.25).clamp(0,1); // ee.Image.clamp - Clamps the values in all bands of an image to all lie within the specified range.

Map.addLayer({
  eeObject: CV_norm,
  //visParams: {min: [0.75], max: [1], palette: palette}, 
  name: 'CV_norm detection Map', 
  shown: 0,
  opacity: 1
});

var testH=ee.String(polar)
var facteurA=(testH.compareTo('VV').eq(0).multiply(2.6).add(testH.compareTo('VH').eq(0).multiply(0.8))).add((testH.compareTo('HH').eq(0).multiply(1.4))); // to recast VV differently from VH

var rgb=ee.Image.cat(days,CV_norm,imax.subtract(0.0).pow(facteurA)).hsvToRgb();  //ee.Image.pow(image2) - Raises the first value to the power of the second for each matched pair
// dynamics are different for VV and VH polarization

var rgb_cut = rgb.clip(Belbek_poligon);

Map.addLayer({
  eeObject: rgb_cut,
  visParams: visparams,
  name: 'REACTIV Visualization', 
  opacity: 0.5,
  shown: 1
  });
  
//Map.addLayer(CV_norm, {min:0, max:1});


// TEMPORAL LEGEND
var vis = {min:0, max:0.85, palette:['FF0000','FF9900','CCFF00','33FF00','00FF66','00FFFF','0066FF','3300FF','CC00FF','FF0099']};//,'FF0000']};
// var palettes = require('users/gena/packages:palettes');
// var palette1 = palettes.colorbrewer.RdYlGn[9].reverse();
// var vis = {min:0, max:0.85, palette: palette1};


function makeLegend(vis) {
  var lon = ee.Image.pixelLonLat().select('longitude');
  var gradient = lon.multiply((vis.max-vis.min)/100.0).add(vis.min);
  var legendImage = gradient.visualize(vis);
  var thumb = ui.Thumbnail({
    image: legendImage, 
    params: {bbox:'0,0,100,4', dimensions:'285x20'},  
    style: {padding: '1px', position: 'bottom-center',backgroundColor:'black'}
  });
  

  var panel = ui.Panel({
    widgets: [
      ui.Label(str1),
      ui.Label(str),
      ui.Label(polar),
      ui.Label(str2),
    ],
    layout: ui.Panel.Layout.flow('horizontal'),
    style: {stretch: 'horizontal',backgroundColor:'white',color:'black'}
  });
  return ui.Panel({style: {backgroundColor: 'white'}}).add(panel).add(thumb);

}

Map.add(makeLegend(vis));

// Wczytanie podkładu google satellite
var palette = [
  '000000','298A08' ,'80FF00','FFFF00','FF8000','FF0000']; //% od najniższego do najwyzszego
var change = CV_norm.gt(0.8).and(imagemax.gt(-10));
Map.setOptions('satellite');

Map.addLayer({
  eeObject: CV_norm.updateMask(change),
  visParams: {min: [0.75], max: [1], palette: palette}, 
  name: 'Change detection Map', 
  shown: 0,
  opacity: 0.5
})


//Map.addLayer(Geilenkirchen_2017, Params);

var maskedW21 = function(image) {
  var RED = image.select(['red']);
  var GREEN = image.select(['green']);
  var BLUE = image.select(['blue']);
  
  // Define masks according to the desired thresholds
  var maskRED = RED.gt(0.5);
  var maskGREEN = GREEN.gt(0.5);
  var maskBLUE = BLUE.gt(0.5);
  
  // Make a mask that fulfills the three masks and rename the band
  var threshold_mask = maskRED.or(maskGREEN).or(maskBLUE)
  .rename('maskedW21');

  // Add the mask as a band to the image and update the mask for the image
  return image.addBands(threshold_mask).updateMask(threshold_mask);
};

var threshold_image = maskedW21(rgb).clip(Belbek_poligon);
Map.addLayer({
  eeObject: threshold_image,
  //visParams: {min: [0.75], max: [1], palette: palette}, 
  name: 'threshold image', 
  shown: 0,
  opacity: 0.5
})
print(threshold_image.bandNames());

// Create a panel to hold the chart.
Map.style().set('cursor', 'crosshair');

var panel = ui.Panel();
panel.style().set({
  width: '500px',
  position: 'bottom-left'
});

Map.add(panel);
// Register a function to draw a chart when a user clicks on the map.
Map.onClick(function(coords) {
  panel.clear();
  var point = ee.Geometry.Point(coords.lon, coords.lat);
  //panel.add(date);
  var chart2=ui.Chart.image.series(sentinel1.select('VV','VH'), point, null, 30)
  //.setChartType('ScatterChart')
      .setOptions({
      title: 'Zmiana intensywności w czasie',
      hAxis: {title: 'Data akwizycji'},
      vAxis: {title: 'Intensywność [dB]'},
      //pointSize: 5,
      //pointShape: 'circle'
      });
  panel.add(chart2);
});


// Create a threshold image.
//var threshold_rgb = rgb.gt(0.5);

// // Compute connected pixel counts; stop searching for connected pixels
// // once the size of the connected neighborhood reaches 30 pixels, and
// // use 8-connected rules.
// var conn = threshold_rgb.connectedPixelCount({
//   maxSize: 4,
//   eightConnected: true
// });

// // Make a binary image of small clusters.
// var smallClusters = conn.lt(30);

// //Map.addLayer(rgb, {min: 0, max: 1}, 'original');
// Map.addLayer(smallClusters.updateMask(smallClusters),
//         {min: 0, max: 1}, 'cc');



// Map.addLayer(threshold_rgb.updateMask(threshold_rgb),
//         {min: 0, max: 1}, 'bright');
         
// // Uniquely label the hotspot image objects.
// var objectId = threshold_rgb.connectedComponents({
//   //connectedness: ee.Kernel.plus(1),
//   connectedness: ee.Kernel.square(1),
//   maxSize: 128
// });

// // Display the uniquely ID'ed objects to the Map.
// Map.addLayer(objectId.randomVisualizer(), null, 'Objects');

// // Compute the number of pixels in each object defined by the "labels" band.
// var objectSize = objectId.select('labels')
//   .connectedPixelCount({
//     maxSize: 128, eightConnected: false
//   });

// // Display object pixel count to the Map.
// Map.addLayer(objectSize, null, 'Object n pixels');


// Display a series by region chart

//var points = ee.FeatureCollection(Geilenkirchen);
// var points_Geilenkirchen = ee.FeatureCollection('users/slesinskijakub/punkty');
// var points_Engels = ee.FeatureCollection('users/slesinskijakub/Engels_pkt');

// Map.addLayer(points_Engels.style({color: 'red', pointShape: 'star5'}));

// var chart0 = ui.Chart.image.seriesByRegion({
//   imageCollection: sentinel1,
//   regions: points_Geilenkirchen, 
//   band: polar,
//   reducer:ee.Reducer.mean(),  
//   scale:10,
//   seriesProperty: 'label'
// })
//   .setChartType('AreaChart')
//       .setOptions({
//       title: 'Zmiana intensywności w czasie',
//       hAxis: {title: 'Data akwizycji'},
//       vAxis: {title: 'Intensywność '+ polar + ' [dB]'},
//       });
// //print(chart0);

// var chart1 = ui.Chart.image.seriesByRegion({
//   imageCollection: sentinel1_liste,
//   regions: points_Geilenkirchen, 
//   band: polar,
//   reducer:ee.Reducer.mean(),  
//   scale:10,
//   seriesProperty: 'label'
// })
//   .setOptions({
//       title: 'Zmiana intensywności w czasie',
//       hAxis: {title: 'Data akwizycji'},
//       vAxis: {title: 'Intensywność '+ polar + ' [dB]'},
//       });
      
// //.setChartType('AreaChart');
// //print(chart1);

// // var chart2 = ui.Chart.image.seriesByRegion({
// //   imageCollection: sentinel1,
// //   regions: points, 
// //   band: polar,
// //   reducer:ee.Reducer.mean(),  
// //   scale:10,
// //   seriesProperty: 'label'
// // })
// // .setChartType('Column Chart');
// // print(chart2);


// // // Display a time-series chart
// // var chart3 = ui.Chart.image.series({
// //   imageCollection: sentinel1,
// //   region: geometry,
// //   reducer: ee.Reducer.mean(),
// //   scale: 20
// // }).setOptions({
// //       lineWidth: 1,
// //       title: 'Sentinel-1 Time Series',
// //       interpolateNulls: true,
// //       vAxis: {title: polar},
// //       hAxis: {title: '', format: 'YYYY-MM'}
// //     })
// //     .setChartType('AreaChart');
// // print(chart3);

// // var chart4 = ui.Chart.image.series({
// //   imageCollection: sentinel1.select('VH', 'VV'),
// //   region: geometry,
// //   reducer: ee.Reducer.mean(),
// //   scale: 50,
// //   xProperty: 'system:time_start'
// // }).setOptions({
// //       lineWidth: 1,
// //       title: 'Sentinel-1 Time Series',
// //       interpolateNulls: true,
// //       vAxis: {title: polar},
// //       hAxis: {title: '', format: 'MM-YYYY'}
// //     })
// //     .setChartType('AreaChart');
// // print(chart4);

// // Create a panel to hold charts.
// var chartpanel = ui.Panel();
// chartpanel.style().set('width', '500px');

// // panels to hold lon/lat values
// var lon = ui.Label();
// var lat = ui.Label();
// chartpanel.add(ui.Panel([lon, lat], ui.Panel.Layout.flow('horizontal')));

// chartpanel.widgets().set(2, chart0);
// chartpanel.widgets().set(3, chart1);

// //ui.root.insert(0, chartpanel);

// var Geilenkirchen_poligon = ee.FeatureCollection('users/slesinskijakub/Geilenkirchen_poli');


// // This function clip each image in collection
// function clp(img) {
//   return img.clip(Czkalowsk_poligon)
// }

// var clipped = sentinel1_liste.map(clp)
// print(clipped)

// // Import Tools
// var batch = require('users/fitoprincipe/geetools:batch');

// // Change Folder Name
// var folder = 'Czkalowsk_EXPORT';

// // Download S1 Collection
// batch.Download.ImageCollection.toDrive(clipped, folder, 
//                 {name: polar+'_'+str+'_'+'{id}',
//                   scale: 10,
//                   maxPixels: 1e13,
//                 region: Czkalowsk_poligon, 
//                 type: 'float'});

var rgb_clip = rgb.clip(Czkalowsk_poligon);

Export.image.toDrive({
  image: rgb_clip,
  description: 'Czkalowsk_rgb_'+str+'_'+polar+'_2021',
  scale: 10,
  region: Czkalowsk_poligon,
  fileFormat: 'GeoTIFF'
})