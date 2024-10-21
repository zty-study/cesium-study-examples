<template>
  <div id="base-cesium-viewer" class="cesium-layout wh-100 relative">
    <router-view v-if="isMounted" />
  </div>
</template>

<script setup lang="ts">
import * as Cesium from 'cesium'
import { provideViewer } from '@/packages/vue3-cesium-use'

const { isMounted } = provideViewer(() => {
  const viewer = new Cesium.Viewer('base-cesium-viewer', {
    infoBox: false,
    baseLayerPicker: false,
    geocoder: false,
    homeButton: false,
    sceneModePicker: false,
    navigationHelpButton: false,
    // animation: false,
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
