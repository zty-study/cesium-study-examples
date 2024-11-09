<script setup lang="ts">
import type { Cartesian3 } from 'cesium'
import { type Component, ref } from 'vue'
import { type MaybeCoordinates } from '../../utils/coordinates'

import { type UseLocatedOptions, useLocated } from '../../hooks/useLocated'

type Nullable<T> = T | null | undefined

defineOptions({
  name: 'Located'
})

const props = withDefaults(
  defineProps<{
    coordinate: Nullable<Cartesian3 | MaybeCoordinates>

    placement?: UseLocatedOptions['placement']

    offset?: UseLocatedOptions['offset']

    as?: string | Component
  }>(),
  {
    as: 'div'
  }
)
const state = defineModel<boolean>({
  default: true
})

const el = ref()
const { style, show } = useLocated(el, {
  state,
  coordinate: () => props.coordinate,
  offset: props.offset,
  placement: props.placement
})
</script>

<template>
  <component
    :is="as"
    ref="el"
    :style="[
      style,
      {
        position: 'absolute',
        zIndex: 999,
        visibility: show ? 'visible' : 'hidden'
      }
    ]"
  >
    <slot />
  </component>
</template>
