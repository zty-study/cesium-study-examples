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

viewer.entities.removeAll()

var center = turf.point([118.38719, 24.449])
var radius = 100
var bearing1 = 95
var bearing2 = 185

var sector = turf.sector(center, radius, bearing1, bearing2)
const positions = sector.geometry.coordinates[0].map((coor) => {
  return Cesium.Cartesian3.fromDegrees(coor[0], coor[1])
})

console.log('positions', positions)

viewer.entities.add({
  polygon: {
    hierarchy: {
      positions: positions as Cesium.Cartesian3[]
    },
    material: Cesium.Color.fromCssColorString('#dcd062').withAlpha(0.7),
    outline: true,
    outlineColor: Cesium.Color.BLACK
  }
})

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
    colorBlendMode: Cesium.ColorBlendMode.HIGHLIGHT
  }
})

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
