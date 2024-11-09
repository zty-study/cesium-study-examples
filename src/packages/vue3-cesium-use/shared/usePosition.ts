import { computed, type CSSProperties, ref } from 'vue'
export type Position =
  | 'top-right'
  | 'top-left'
  | 'bottom-right'
  | 'bottom-left'
  | 'top'
  | 'right'
  | 'bottom'
  | 'left'

export const usePosition = (position: Position, offset?: [number, number]) => {
  const attach = computed(() => {
    return {
      top: position.indexOf('top') > -1,
      right: position.indexOf('right') > -1,
      bottom: position.indexOf('bottom') > -1,
      left: position.indexOf('left') > -1,
      vertical: position === 'top' || position === 'bottom',
      horizontal: position === 'left' || position === 'right'
    }
  })

  const top = ref(0)
  const right = ref(0)
  const left = ref(0)
  const bottom = ref(0)

  const style = computed(() => {
    let posX: number | string = 0
    let posY: number | string = 0

    const side = attach.value
    const dir = 1

    if (side.top === true && top.value !== 0) {
      posY = `${top.value}px`
    } else if (side.bottom === true && bottom.value !== 0) {
      posY = `${-bottom.value}px`
    }

    if (side.left === true && left.value !== 0) {
      posX = `${dir * left.value}px`
    } else if (side.right === true && right.value !== 0) {
      posX = `${-dir * right.value}px`
    }

    const css: CSSProperties = {
      transform: `translate(${posX}, ${posY})`
    }

    if (offset) {
      css.margin = `${offset[1]}px ${offset[0]}px`
    }

    if (side.vertical === true) {
      if (left.value !== 0) {
        css['right'] = `${left.value}px`
      }
      if (right.value !== 0) {
        css['left'] = `${right.value}px`
      }
    } else if (side.horizontal === true) {
      if (top.value !== 0) {
        css.top = `${top.value}px`
      }
      if (bottom.value !== 0) {
        css.bottom = `${bottom.value}px`
      }
    }

    return css
  })

  const classes = computed(() => `absolute absolute-${position}`)

  return {
    attach,
    style,
    classes
  }
}
