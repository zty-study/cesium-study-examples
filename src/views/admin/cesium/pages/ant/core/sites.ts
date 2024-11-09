import * as Cesium from 'cesium'
import type { ModelOptions } from '@/packages/vue3-cesium-use'
import { toCartesian3 } from '@/packages/vue3-cesium-use'

const position = Cesium.Cartesian3.fromDegrees(119.04707, 23.85515, 20200)
// const heading = Cesium.Math.toRadians(135)

const models = [
  {
    id: '1',
    position: [118.38719, 24.449, 0],
    hpr: [0, 0, 0]
  },
  {
    id: '2',
    position: [121.29509, 28.44508, 0],
    hpr: [80, 0, 0]
  },
  {
    id: '3',
    position: [115.93484, 23.00697, 0],
    hpr: [80, 0, 0]
  }
]

export const siteList: ModelOptions[] = models.map((it) => {
  const heading = Cesium.Math.toRadians(80)
  const pitch = Cesium.Math.toRadians(10)
  const roll = Cesium.Math.toRadians(10)
  const hpr = new Cesium.HeadingPitchRoll(heading, pitch, roll)
  const orientation = Cesium.Transforms.headingPitchRollQuaternion(position, hpr)
  return {
    id: it.id,
    position: toCartesian3(it.position),
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
