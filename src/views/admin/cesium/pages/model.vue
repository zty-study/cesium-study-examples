<template>
  <div>
    <MapInfo class="map-base-info z-2" />
  </div>
</template>

<script setup lang="ts">
import { MapInfo } from '@/packages/vue3-cesium-use'

import * as Cesium from 'cesium'
import { getViewer } from '@/packages/vue3-cesium-use'

const viewer = getViewer()
const token = '1b7e6d0d69b7c4ac48e13111aa6dbf81'
// 服务域名
const tdtUrl = 'https://t{s}.tianditu.gov.cn/'
// 服务负载子域
const subdomains = ['0', '1', '2', '3', '4', '5', '6', '7']

viewer.entities.removeAll()

const position = Cesium.Cartesian3.fromDegrees(119.04707, 23.85515, 20200)
// const heading = Cesium.Math.toRadians(135)
const heading = Cesium.Math.toRadians(80)
const pitch = Cesium.Math.toRadians(10)
const roll = Cesium.Math.toRadians(10)
const hpr = new Cesium.HeadingPitchRoll(heading, pitch, roll)
const orientation = Cesium.Transforms.headingPitchRollQuaternion(position, hpr)

const entity = viewer.entities.add({
  position: position,
  orientation: orientation,
  model: {
    uri: '/source/model/Cesium_Air.glb',
    // uri: '/public/source/model/douglas_xb19/scene.gltf',
    minimumPixelSize: 128,
    maximumScale: 20000
  }
})
// viewer.trackedEntity = entity

viewer.entities.add({
  position: Cesium.Cartesian3.fromDegrees(118.38719, 24.449, 0),
  //   orientation: orientation,
  model: {
    uri: '/public/source/model/satellite_ground_station/scene.gltf',
    // uri: '/public/source/model/douglas_xb19/scene.gltf',
    minimumPixelSize: 56,
    maximumScale: 200,
    color: Cesium.Color.fromCssColorString('#ff0000'),
    colorBlendMode: Cesium.ColorBlendMode.HIGHLIGHT
  }
})

// 叠加国界服务天地图
const iboMap = new Cesium.UrlTemplateImageryProvider({
  url: tdtUrl + 'DataServer?T=ibo_w&x={x}&y={y}&l={z}&tk=' + token,
  subdomains: subdomains,
  tilingScheme: new Cesium.WebMercatorTilingScheme(),
  maximumLevel: 10
})
viewer.imageryLayers.addImageryProvider(iboMap)

// 矢量注记
// viewer.imageryLayers.addImageryProvider(
//   new Cesium.WebMapTileServiceImageryProvider({
//     url:
//       'http://t0.tianditu.com/cva_w/wmts?service=wmts&request=GetTile&version=1.0.0&LAYER=cva&tileMatrixSet=w&TileMatrix={TileMatrix}&TileRow={TileRow}&TileCol={TileCol}&style=default.jpg&tk=' +
//       token,
//     layer: 'tdtAnnoLayer',
//     style: 'default',
//     format: 'image/jpeg',
//     tileMatrixSetID: 'GoogleMapsCompatible'
//   })
// )

// 影像注记
// viewer.imageryLayers.addImageryProvider(
//   new Cesium.WebMapTileServiceImageryProvider({
//     url:
//       'http://t0.tianditu.com/cia_w/wmts?service=wmts&request=GetTile&version=1.0.0&LAYER=cia&tileMatrixSet=w&TileMatrix={TileMatrix}&TileRow={TileRow}&TileCol={TileCol}&style=default.jpg&tk=' +
//       token,
//     layer: 'tdtAnnoLayer',
//     style: 'default',
//     format: 'image/jpeg',
//     tileMatrixSetID: 'GoogleMapsCompatible'
//   })
// )

onMounted(() => {
  viewer.camera.setView({
    // destination: Cesium.Cartesian3.fromDegrees(105.30018, 26.996108, 1500)
    destination: new Cesium.Cartesian3(-2596707.29925, 5292676.85157, 2575124.67888),
    orientation: {
      heading: Cesium.Math.toRadians(81.32),
      pitch: Cesium.Math.toRadians(-12.1),
      roll: Cesium.Math.toRadians(0.41)
    }
  })
})
</script>

<style lang="scss" scoped>
.map-base-info {
  left: 0;
  top: 0;
  position: absolute;
}
</style>
