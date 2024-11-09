<template>
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
    <img :src="outerImg" class="wh-100" :style="outerCircleStyle" />
    <img :src="innerImg" class="wh-100" />
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
const rootStyle = reactive<CSSProperties>({})
const compassState = useCompass()

// 导航盘位置
const positionState = computed(() => usePosition(props.position, props.offset))

const outerCircleStyle = computed(() => {
  return {
    transform: 'rotate(-' + compassState.heading.value + 'rad)',
    WebkitTransform: 'rotate(-' + compassState.heading.value + 'rad)'
  }
})

compassState.load()

onBeforeUnmount(() => {
  compassState.unload()
})
</script>

<style lang="scss" scoped>
.cesium-use-compass-wrap {
  z-index: 9999;
  height: 120px;
  width: 120px;
  user-select: none;
  position: absolute;
  overflow: hidden;
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
.fixed-top,
.absolute-top {
  top: 0;
  left: 0;
  right: 0;
}
.fixed-right,
.absolute-right {
  top: 0;
  right: 0;
  bottom: 0;
}
.fixed-bottom,
.absolute-bottom {
  right: 0;
  bottom: 0;
  left: 0;
}
.fixed-left,
.absolute-left {
  top: 0;
  bottom: 0;
  left: 0;
}
.fixed-top-left,
.absolute-top-left {
  top: 0;
  left: 0;
}
.fixed-top-right,
.absolute-top-right {
  top: 0;
  right: 0;
}
.fixed-bottom-left,
.absolute-bottom-left {
  bottom: 0;
  left: 0;
}
.fixed-bottom-right,
.absolute-bottom-right {
  bottom: 0;
  right: 0;
}
.absolute-center {
  top: 50%;
  left: 50%;
  transform: translate(-50%, -50%);
}
</style>
