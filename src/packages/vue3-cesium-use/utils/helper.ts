import * as Cesium from 'cesium'

// 计算高度对应的level
export function heightToLevel(altitude: number) {
  // 粗略计算
  const A = 40487.57
  const B = 0.00007096758
  const C = 91610.74
  const D = -40467.74

  return Math.round(D + (A - D) / (1 + Math.pow(altitude / C, B)))
}
