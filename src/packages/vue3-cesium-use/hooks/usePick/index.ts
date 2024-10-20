import * as Cesium from 'cesium'
import { ref } from 'vue'
import { getViewer } from '../useViewer'

export const usePick = () => {
  const viewer = getViewer()

  const pickItem = ref<any>()
  const pickList = ref<any[]>([])
  const mousePoint = ref<Cesium.Cartesian2>(new Cesium.Cartesian2(0, 0))
  const mouseCartesian = ref<Cesium.Cartesian3>(new Cesium.Cartesian3())
  const mouseCoordniates = ref<Cesium.Cartographic>(new Cesium.Cartographic())

  const handler = new Cesium.ScreenSpaceEventHandler(viewer.canvas)

  handler.setInputAction(({ endPosition }: Cesium.ScreenSpaceEventHandler.MotionEvent) => {
    mouseCartesian.value = viewer.scene.pickPosition(endPosition) || new Cesium.Cartesian3()
    mouseCoordniates.value =
      Cesium.Cartographic.fromCartesian(mouseCartesian.value)?.clone() || new Cesium.Cartographic()
    mousePoint.value = Cesium.Cartesian2.clone(endPosition)
    pickItem.value = viewer.scene.pick(endPosition)
    pickList.value = viewer.scene.drillPick(endPosition)
  }, Cesium.ScreenSpaceEventType.MOUSE_MOVE)

  return { pickItem, pickList, mousePoint, mouseCartesian, mouseCoordniates }
}
