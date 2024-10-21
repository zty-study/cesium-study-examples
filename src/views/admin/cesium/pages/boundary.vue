<template>
  <div></div>
</template>

<script setup lang="ts">
import * as Cesium from 'cesium'
import { getViewer } from '@/packages/vue3-cesium-use'

const viewer = getViewer()

Cesium.GeoJsonDataSource.load('https://geo.datav.aliyun.com/areas_v3/bound/100000.json', {
  stroke: Cesium.Color.YELLOW,
  strokeWidth: 100,
  fill: Cesium.Color.PALETURQUOISE.withAlpha(0) //填充区域设置为透明
}).then((dataSource) => {
  viewer.dataSources.add(dataSource)

  // viewer.zoomTo(dataSource)
})

// Cesium.Cesium3DTileset.fromUrl('/source/3dtiles/china/tileset.json').then((tileset) => {
//   viewer.scene.primitives.add(tileset)
// })

onMounted(() => {
  viewer.camera.setView({
    destination: Cesium.Cartesian3.fromDegrees(101.4728, 38.6736, 10000000)
    // orientation: {
    //   heading: 4.417425186487677,
    //   pitch: -0.3107427037808055,
    //   roll: 6.283161438815717
    // }
  })

  //   viewer.scene.primitives.add(await Cesium.createOsmBuildingsAsync())
})
</script>

<style lang="scss" scoped></style>
