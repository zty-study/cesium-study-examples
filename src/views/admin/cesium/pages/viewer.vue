<template>
  <div class="cesium-viewer-page">
    <MapInfo class="map-base-info z-2" />
  </div>
</template>

<script setup lang="ts">
import { MapInfo } from '@/packages/vue3-cesium-use'
</script>

<style lang="scss" scoped>
.map-base-info {
  left: 0;
  top: 0;
  position: absolute;
}
</style>
