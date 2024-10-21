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

viewer.entities.removeAll()

const position = Cesium.Cartesian3.fromDegrees(105.30018, 26.996108, 1500)
// const heading = Cesium.Math.toRadians(135)
const heading = Cesium.Math.toRadians(268.89)
const pitch = Cesium.Math.toRadians(-11.62)
const roll = 0
const hpr = new Cesium.HeadingPitchRoll(heading, pitch, roll)
const orientation = Cesium.Transforms.headingPitchRollQuaternion(position, hpr)

const entity = viewer.entities.add({
  position: position,
  orientation: orientation,
  model: {
    uri: '/source/model/Cesium_Air.glb',
    minimumPixelSize: 128,
    maximumScale: 20000
  }
})
viewer.trackedEntity = entity

onMounted(() => {
  viewer.camera.setView({
    destination: Cesium.Cartesian3.fromDegrees(105.30018, 26.996108, 1500)
    // orientation: {
    //   heading: 4.417425186487677,
    //   pitch: -0.3107427037808055,
    //   roll: 6.283161438815717
    // }
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
