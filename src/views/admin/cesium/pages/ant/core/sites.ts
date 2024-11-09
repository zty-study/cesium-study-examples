import * as Cesium from 'cesium'
import type { ModelOptions } from '@/packages/vue3-cesium-use'

export const siteList: ModelOptions[] = [
  {
    id: '1',
    position: Cesium.Cartesian3.fromDegrees(118.38719, 24.449, 0),
    //   orientation: orientation,
    model: {
      uri: '/source/model/satellite_ground_station/scene.gltf',
      // uri: '/public/source/model/douglas_xb19/scene.gltf',
      minimumPixelSize: 80,
      maximumScale: 500,
      heightReference: Cesium.HeightReference.RELATIVE_TO_GROUND,
      colorBlendMode: Cesium.ColorBlendMode.HIGHLIGHT
    }
  }
]
