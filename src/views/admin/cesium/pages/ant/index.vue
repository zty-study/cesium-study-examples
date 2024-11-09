<template>
  <div>
    <MapInfo class="map-base-info z-2" />
  </div>
</template>

<script setup lang="ts">
import { MapInfo } from '@/packages/vue3-cesium-use'
import * as turf from '@turf/turf'
import * as Cesium from 'cesium'
import { getViewer, useModel, useSkyBox } from '@/packages/vue3-cesium-use'
import { siteList } from './core/sites'

const viewer = getViewer()
const { modelList } = useModel()
const { setSkyBox } = useSkyBox()

setSkyBox({
  positiveX: '/source/skybox/3/tycho2t3_80_px.jpg',
  negativeX: '/source/skybox/3/tycho2t3_80_mx.jpg',
  positiveY: '/source/skybox/3/tycho2t3_80_py.jpg',
  negativeY: '/source/skybox/3/tycho2t3_80_my.jpg',
  positiveZ: '/source/skybox/3/tycho2t3_80_pz.jpg',
  negativeZ: '/source/skybox/3/tycho2t3_80_mz.jpg'
})

viewer.camera.setView({
  destination: new Cesium.Cartesian3(-2836751.16204, 6397045.27142, 1724458.85728),
  orientation: {
    heading: Cesium.Math.toRadians(19.23),
    pitch: Cesium.Math.toRadians(-41.54),
    roll: Cesium.Math.toRadians(0.54)
  }
})

modelList.value = siteList
</script>

<style lang="scss" scoped>
.map-base-info {
  left: 0;
  top: 0;
  position: absolute;
}
</style>
