import * as Cesium from 'cesium'
import { getViewer } from '../../useViewer'
import { onBeforeUnmount } from 'vue'

export const useDrawSensor = () => {
  const viewer = getViewer()

  // const sensorCollection = new Cesium.CustomDataSource('sensor')
  // viewer.dataSources.add(sensorCollection)
  const graphicLaayer = new xt3d.layer.GraphicLayer({
    viewer: viewer
  })

  const s1 = new xt3d.graphic.satellite.RectSensor({
    position: [115.93484, 23.00697, 700],
    style: {
      radius: 400000,
      pitch: 0,
      heading: -45,
      roll: 80,
      vAngle: 20,
      hAngle: 90,
      scanPlaneRate: 7,
      color: 'rgba(255, 1, 1, 0.2)'
    }
  })
  const s2 = new xt3d.graphic.satellite.RectSensor({
    position: [121.29509, 28.44508, 200],
    style: {
      radius: 500000,
      pitch: 0,
      heading: -45,
      roll: 80,
      vAngle: 20,
      hAngle: 90
    }
  })
  graphicLaayer.add(s1)
  graphicLaayer.add(s2)
  graphicLaayer.flyTo()

  onBeforeUnmount(() => {
    // graphicLaayer.remove(s)
    graphicLaayer.removeAll()
  })
}
