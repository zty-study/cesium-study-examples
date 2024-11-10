<template>
  <div id="base-cesium-viewer" class="cesium-layout wh-100 relative">
    <router-view v-if="isMounted" />
  </div>
</template>

<script setup lang="ts">
import * as Cesium from 'cesium'
import { provideViewer } from '@/packages/vue3-cesium-use'

const token = '1b7e6d0d69b7c4ac48e13111aa6dbf81'
// 服务域名
const tdtUrl = 'https://t{s}.tianditu.gov.cn/'
// 服务负载子域
const subdomains = ['0', '1', '2', '3', '4', '5', '6', '7']

const { isMounted } = provideViewer(() => {
  const viewer = new Cesium.Viewer('base-cesium-viewer', {
    infoBox: false,
    baseLayerPicker: false,
    geocoder: false,
    homeButton: false,
    sceneModePicker: false,
    navigationHelpButton: false,
    // animation: false,
    shouldAnimate: true,
    timeline: false,
    fullscreenButton: false,
    selectionIndicator: false, // 取消绿色选框
    vrButton: false,
    sceneMode: Cesium.SceneMode.SCENE3D
  })

  Cesium.createWorldTerrainAsync({
    requestVertexNormals: true,
    requestWaterMask: true
  }).then((terrainProvider) => {
    viewer.terrainProvider = terrainProvider
  })

  // 抗锯齿
  viewer.scene.postProcessStages.fxaa.enabled = false
  // 水雾特效
  viewer.scene.globe.showGroundAtmosphere = true

  viewer.resolutionScale = window.devicePixelRatio

  // 叠加国界服务天地图
  // const iboMap = new Cesium.UrlTemplateImageryProvider({
  //   url: tdtUrl + 'DataServer?T=ibo_w&x={x}&y={y}&l={z}&tk=' + token,
  //   subdomains: subdomains,
  //   tilingScheme: new Cesium.WebMercatorTilingScheme(),
  //   maximumLevel: 10
  // })
  // viewer.imageryLayers.addImageryProvider(iboMap)

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

  return viewer
})
</script>

<style lang="scss">
.cesium-widget-credits {
  display: none !important;
  visibility: hide !important;
}
.cesium-viewer-animationContainer {
  display: none !important;
}
</style>
