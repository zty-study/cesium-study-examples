<template>
  <div class="cesium-map-base-info">
    <div class="cesium-map-base-info-line-wrap">
      <span class="cesium-map-base-info-title">【屏幕】</span>
      <span>X: {{ mousePoint.x.toFixed(2) }}</span>
      <span>Y: {{ mousePoint.y.toFixed(2) }}</span>
      <span></span>
    </div>
    <div class="cesium-map-base-info-line-wrap">
      <span class="cesium-map-base-info-title">【坐标】</span>

      <span>经度: {{ coordinates.longitude }}</span>
      <span>纬度: {{ coordinates.latitude }}</span>
      <span>高度: {{ coordinates.altitude }}</span>
    </div>
    <div class="cesium-map-base-info-line-wrap">
      <span class="cesium-map-base-info-title">【相机】</span>

      <span>方位: {{ (cameraInfo.heading || 0).toFixed(2) }}°</span>
      <span>俯仰: {{ (cameraInfo.pitch || 0).toFixed(2) }}°</span>
      <span>侧翻: {{ (cameraInfo.roll || 0).toFixed(2) }}°</span>
    </div>

    <div class="cesium-map-base-info-line-wrap">
      <span class="cesium-map-base-info-title">【相机】</span>
      <span>X: {{ (cameraInfo.position?.x || 0).toFixed(5) }}</span>
      <span>Y: {{ (cameraInfo.position?.y || 0).toFixed(5) }}</span>
      <span>Z: {{ (cameraInfo.position?.z || 0).toFixed(5) }}</span>
    </div>

    <div class="cesium-map-base-info-line-wrap">
      <span class="cesium-map-base-info-title">【笛卡尔】</span>

      <span>X: {{ mouseCartesian.x.toFixed(5) }}</span>
      <span>Y: {{ mouseCartesian.x.toFixed(5) }}</span>
      <span>Z: {{ mouseCartesian.x.toFixed(5) }}</span>
    </div>

    <div class="cesium-map-base-info-line-wrap">
      <span class="cesium-map-base-info-title">【性能】</span>
      <span>{{ performanceInfo.fps }}</span>
      <span>{{ performanceInfo.ms }}</span>
      <span></span>
    </div>

    <div class="cesium-map-base-info-line-wrap">
      <span class="cesium-map-base-info-title">【其他】</span>
      <span>视高: {{ (cameraInfo.height || 0).toFixed(2) }}m</span>
      <span>层级: {{ cameraInfo.level || 0 }}</span>
      <span></span>
    </div>
  </div>
</template>

<script setup lang="ts">
import { reactive } from 'vue'
import { usePick } from '../../hooks/usePick'
import { getViewer } from '../../hooks/useViewer'
import { useCamera } from '../../hooks/useCamera'
import { useThrottleFn } from '@vueuse/core'
import { prettifyCoordinates } from '../../utils/coordinates'

const viewer = getViewer()
viewer.scene.debugShowFramesPerSecond = true

const performanceInfo = reactive({
  fps: '0 FPS',
  ms: '0 MS'
})

const onScenePostRender = useThrottleFn((scene) => {
  performanceInfo.fps = scene._performanceDisplay?._fpsText.nodeValue
  performanceInfo.ms = scene._performanceDisplay?._msText.nodeValue
  scene._performanceDisplay._container.style.display = 'none'
}, 250)

viewer.scene.postRender.addEventListener(onScenePostRender)

const { cameraInfo } = useCamera()
const { mousePoint, mouseCartesian, mouseCoordniates } = usePick()

const coordinates = computed(() =>
  prettifyCoordinates(mouseCoordniates.value.longitude, mouseCoordniates.value.latitude, {
    height: mouseCoordniates.value.height
  })
)
</script>

<style lang="scss" scoped>
.cesium-map-base-info {
  font-size: 12px;
  color: #fff;
  padding: 5px;
  background-color: rgba(0, 0, 0, 0.6);
  min-width: 500px;
  width: fit-content;

  &-line-wrap {
    display: flex;
    align-items: center;
    justify-content: space-between;
    margin-bottom: 5px;
    &:last-child {
      margin-bottom: 0;
    }
    span {
      flex: 1;
      margin-right: 20px;
      white-space: nowrap;
      &:last-child {
        margin-right: 0;
      }

      &.cesium-map-base-info-title {
        flex: 0 0 80px;
        margin-right: 0px;
      }
    }
  }
}
</style>
