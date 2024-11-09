<template>
  <div>
    <MapInfo class="map-base-info z-2" />
  </div>
</template>

<script setup lang="ts">
import { MapInfo } from '@/packages/vue3-cesium-use'
import * as turf from '@turf/turf'
import * as Cesium from 'cesium'
import { getViewer } from '@/packages/vue3-cesium-use'

const viewer = getViewer()

onMounted(() => {
  viewer.camera.setView({
    // destination: Cesium.Cartesian3.fromDegrees(105.30018, 26.996108, 1500)
    destination: new Cesium.Cartesian3(-2836751.16204, 6397045.27142, 1724458.85728),
    orientation: {
      heading: Cesium.Math.toRadians(19.23),
      pitch: Cesium.Math.toRadians(-41.54),
      roll: Cesium.Math.toRadians(0.54)
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
