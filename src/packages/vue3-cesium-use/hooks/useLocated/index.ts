import type { Cartesian3 } from 'cesium'
import { toRefs, useElementBounding } from '@vueuse/core'
import { type MaybeRefOrGetter, type Ref, computed, ref, toValue } from 'vue'
import { type MaybeCoordinates } from '../../utils/coordinates'
import { toCartesian3, useEventHandler, getViewer } from '../../index'

type Nullable<T> = T | null | undefined

type Placement =
  | 'topLeft'
  | 'top'
  | 'topRight'
  | 'right'
  | 'bottomRight'
  | 'bottom'
  | 'bottomLeft'
  | 'left'

export interface UseLocatedOptions {
  state: Ref<boolean>
  coordinate: MaybeRefOrGetter<Nullable<Cartesian3 | MaybeCoordinates>>

  /**
   * @default true
   */
  initialValue: boolean

  /**
   * @default 'bottomRight'
   */
  placement: Placement
  offset: {
    left?: number
    top?: number
  }
}

export function useLocated(
  el: MaybeRefOrGetter<Nullable<HTMLElement>>,
  options: Partial<UseLocatedOptions> = {}
) {
  const {
    state: _state,
    coordinate: _coordinate,
    initialValue = true,
    placement = 'bottomRight',
    offset = {}
  } = options

  const state = _state ?? ref(initialValue)
  const coordinate = _coordinate
    ? computed(() => toValue(_coordinate))
    : ref<Nullable<Cartesian3 | MaybeCoordinates>>()
  const __show = ref(true)
  const position = ref({ x: 0, y: 0 })

  const { width, height } = useElementBounding(el)

  const placementMap: Record<
    Placement,
    (rect: { width: number; height: number }) => { width: number; height: number }
  > = {
    topLeft: (rect) => ({
      width: -rect.width,
      height: -rect.height
    }),
    top: (rect) => ({ width: -rect.width / 2, height: -rect.height }),
    topRight: (rect) => ({ width: 0, height: -rect.height }),
    right: (rect) => ({ width: 0, height: -rect.height / 2 }),
    bottomRight: (_rect) => ({ width: 0, height: 0 }),
    bottom: (rect) => ({ width: -rect.width / 2, height: 0 }),
    bottomLeft: (rect) => ({ width: -rect.width, height: 0 }),
    left: (rect) => ({ width: -rect.width, height: -rect.height / 2 })
  }

  const viewer = getViewer()

  const onPostRender = useEventHandler(viewer.scene.postRender)
  onPostRender(() => {
    if (!state.value) return

    if (!coordinate.value) {
      __show.value = false
      return
    }

    const pos = viewer.scene.cartesianToCanvasCoordinates(toCartesian3(coordinate.value))

    if (!pos) {
      __show.value = false
      return
    }
    __show.value = true

    const offsetPlacement = placementMap[placement]({
      width: width.value,
      height: height.value
    })

    position.value = {
      x: pos.x + (offsetPlacement.width ?? 0) + (offset.left ?? 0),
      y: pos.y + (offsetPlacement.height ?? 0) + (offset.top ?? 0)
    }
  })

  return {
    ...toRefs(position),
    position,
    state,
    show: computed(() => state.value && __show.value),
    coordinate,
    style: computed(() => `left:${position.value.x}px;top:${position.value.y}px`)
  }
}
