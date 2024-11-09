import { getViewer } from '../../index'
import * as Cesium from 'cesium'

/**
 * ## 天空盒格式
 */
export interface SkyBoxOptions {
  positiveX: string
  negativeX: string
  positiveY: string
  negativeY: string
  positiveZ: string
  negativeZ: string
}

export const useSkyBox = () => {
  const viewer = getViewer()

  /**
   * ## 设置天空盒
   * ---
   * @param sources 天空盒图片
   * @example
   * sources: {
   *     positiveX: 'right.jpg',    // 右侧图片
   *     negativeX: 'left.jpg',     // 左侧图片
   *     positiveY: 'front.jpg',    // 前方图片
   *     negativeY: 'back.jpg',     // 后方图片
   *     positiveZ: 'top.jpg',      // 顶部图片
   *     negativeZ: 'bottom.jpg'    // 底部图片
   * }
   */
  const setSkyBox = (sources: SkyBoxOptions) => {
    viewer.scene.skyBox = new Cesium.SkyBox({
      sources: {
        positiveX: sources.positiveX,
        negativeX: sources.negativeX,
        positiveY: sources.positiveY,
        negativeY: sources.negativeY,
        positiveZ: sources.positiveZ,
        negativeZ: sources.negativeZ
      }
    })
  }

  return {
    setSkyBox
  }
}
