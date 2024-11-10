<template>
  <div class="cesium-use-compass-contain">
    <div
      class="cesium-use-compass-wrap"
      :class="positionState.classes.value"
      :style="positionState.style.value"
      @dblclick="compassState.handleDoubleClick"
      @mousedown="compassState.handleMouseDown"
      @mouseup="compassState.resetRotater"
      @touchend="compassState.resetRotater"
      @touchstart="compassState.handleMouseDown"
    >
      <img :src="outerImg" class="wh-100" :style="outerImgCss" />
      <img :src="innerImg" class="wh-100" />
      <div class="compass-rotation-marker absolute-center" :style="rotationMarkerStyle">
        <svg viewBox="0 0 1024 1024">
          <path
            d="M0 506.590189C0 226.82566 226.82566 0 506.590189 0v173.886792C322.849811 173.886792 173.886792 322.849811 173.886792 506.590189"
            fill="#00aefe"
            opacity=".8"
          ></path>
        </svg>
      </div>
    </div>
  </div>
</template>

<script setup lang="ts">
import outSvg from './sources/outerSvg.svg'
import innerSvg from './sources/innerSvg.svg'
import useCompass from './useCompass'
import { type CSSProperties, reactive, computed, onBeforeUnmount } from 'vue'
import { usePosition, type Position } from '../../shared/usePosition'

// const outSvg = require('./sources/outerSvg.svg')
// const innerSvg = require('./sources/innerSvg.svg')

const props = withDefaults(
  defineProps<{
    position?: Position
    offset?: [number, number]
    outerImg?: string
    innerImg?: string
  }>(),
  {
    position: 'top-right',
    offset: () => [0, 0],
    outerImg: outSvg,
    innerImg: innerSvg
  }
)
const compassState = useCompass()

// 导航盘位置
const positionState = computed(() => usePosition(props.position, props.offset))

const outerImgCss = computed(() => {
  return {
    transform: 'rotate(-' + compassState.heading.value + 'rad)',
    WebkitTransform: 'rotate(-' + compassState.heading.value + 'rad)'
  }
})

const rotationMarkerStyle = computed(() => {
  return {
    transform: 'rotate(-' + compassState.orbitCursorAngle.value + 'rad)',
    WebkitTransform: 'rotate(-' + compassState.orbitCursorAngle.value + 'rad)',
    opacity: compassState.orbitCursorOpacity.value
  }
})

compassState.load()

onBeforeUnmount(() => {
  compassState.unload()
})
</script>

<style lang="scss" scoped>
.cesium-use-compass-contain {
  position: absolute;
  top: 0;
  left: 0;
  right: 0;
  bottom: 0;
  display: flex;
  justify-content: center;
  align-items: center;
  z-index: 9999;
  pointer-events: none;
}
.cesium-use-compass-wrap {
  height: 120px;
  width: 120px;
  user-select: none;
  position: absolute;
  overflow: hidden;
  border-radius: 50%;
  background-color: rgba(0, 0, 0, 0.5);
  pointer-events: all;

  img {
    position: absolute;
    user-select: none;
    -webkit-user-drag: none;
  }
}
.wh-100 {
  width: 100%;
  height: 100%;
}
.absolute {
  position: absolute;
}
.absolute-top {
  top: 0;
}
.absolute-right {
  right: 0;
}
.absolute-bottom {
  bottom: 0;
}
.absolute-left {
  left: 0;
}
.absolute-top-left {
  top: 0;
  left: 0;
}
.absolute-top-right {
  top: 0;
  right: 0;
}
.absolute-bottom-left {
  bottom: 0;
  left: 0;
}
.absolute-bottom-right {
  bottom: 0;
  right: 0;
}
</style>
