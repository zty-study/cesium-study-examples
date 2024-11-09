<template>
  <div>
    <MapInfo class="map-base-info z-2" />
  </div>
</template>

<script setup lang="ts">
import { MapInfo } from '@/packages/vue3-cesium-use'
import * as Cesium from 'cesium'
import { getViewer, useModel, useSkyBox, toCartesian3 } from '@/packages/vue3-cesium-use'

const viewer = getViewer()
const { modelList } = useModel()
const { setSkyBox } = useSkyBox()

const models = ref([
  {
    id: '1',
    position: [118.38719, 24.449, 0],
    hpr: [-20, 0, 0]
  },
  {
    id: '2',
    position: [121.29509, 28.44508, 0],
    hpr: [-30, 0, 0]
  },
  {
    id: '3',
    position: [115.93484, 23.00697, 0],
    hpr: [-40, 0, 0]
  }
])

watchEffect(() => {
  modelList.value = models.value.map((it) => {
    const position = toCartesian3(it.position)
    const heading = Cesium.Math.toRadians(it.hpr[0])
    const pitch = Cesium.Math.toRadians(it.hpr[1])
    const roll = Cesium.Math.toRadians(it.hpr[2])
    const hpr = new Cesium.HeadingPitchRoll(heading, pitch, roll)
    const orientation = Cesium.Transforms.headingPitchRollQuaternion(position, hpr)
    return {
      id: it.id,
      position,
      orientation,
      model: {
        uri: '/source/model/satellite_ground_station/scene.gltf',
        minimumPixelSize: 80,
        maximumScale: 500,
        heightReference: Cesium.HeightReference.RELATIVE_TO_GROUND,
        colorBlendMode: Cesium.ColorBlendMode.HIGHLIGHT
      }
    }
  })
})

const methods = {
  // 初始化场景
  initSecen() {
    methods.setSceneSkyBox()
    methods.setCameraView()
  },

  // 设置场景天空盒
  setSceneSkyBox() {
    setSkyBox({
      positiveX: '/source/skybox/3/tycho2t3_80_px.jpg',
      negativeX: '/source/skybox/3/tycho2t3_80_mx.jpg',
      positiveY: '/source/skybox/3/tycho2t3_80_py.jpg',
      negativeY: '/source/skybox/3/tycho2t3_80_my.jpg',
      positiveZ: '/source/skybox/3/tycho2t3_80_pz.jpg',
      negativeZ: '/source/skybox/3/tycho2t3_80_mz.jpg'
    })
  },

  // 设置相机视角
  setCameraView() {
    viewer.camera.setView({
      destination: new Cesium.Cartesian3(-2836751.16204, 6397045.27142, 1724458.85728),
      orientation: {
        heading: Cesium.Math.toRadians(19.23),
        pitch: Cesium.Math.toRadians(-41.54),
        roll: Cesium.Math.toRadians(0.54)
      }
    })
  }
}

methods.initSecen()
</script>

<style lang="scss" scoped>
.map-base-info {
  left: 0;
  top: 0;
  position: absolute;
}
</style>
