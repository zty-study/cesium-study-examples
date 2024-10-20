import * as Cesium from 'cesium'

export type MaybeCoordinates =
  | (number | string)[]
  | {
      longitude: number | string
      latitude: number | string
      height?: number | string
    }

export function normalizeCoordinates(source: MaybeCoordinates) {
  const value = Array.isArray(source)
    ? {
        longitude: Number(source[0]),
        latitude: Number(source[1]),
        height: source[2] ? Number(source[2]) : undefined
      }
    : {
        longitude: Number(source.longitude),
        latitude: Number(source.latitude),
        height: source.height ? Number(source.height) : undefined
      }

  if (
    Object.entries(value).some(([k, v]) => {
      if (k === 'height') return !v ? false : Number.isNaN(v)

      return Number.isNaN(v)
    })
  ) {
    throw new Error('Invalid value which cannot be parsed to number.')
  }

  return value
}

/**
 * 地理坐标格式化.
 *
 * @param {Number} longitude 经度
 * @param {Number} latitude 纬度
 * @param {Number} height 高度
 */
export function prettifyCoordinates(
  longitude: number,
  latitude: number,
  options?: {
    decimal?: number
    rangeType?: 0 | 1 | 2
    height?: number
    errorBar?: number
  }
) {
  {
    const result = {
      latitude: '-',
      longitude: '-',
      altitude: '-'
    }

    const { Math: CesiumMath } = Cesium
    const lat = CesiumMath.toDegrees(latitude)
    const lon = CesiumMath.toDegrees(longitude)

    const optionsDefaulted = Cesium.defaultValue(options, {})
    const decimal = Cesium.defaultValue(optionsDefaulted.decimal, 5)
    const rangeType = Cesium.defaultValue(optionsDefaulted.rangeType, 0)

    if (rangeType === 0) {
      result.latitude = Math.abs(lat).toFixed(decimal) + '°' + (lat < 0.0 ? 'S' : 'N')
      result.longitude = Math.abs(lon).toFixed(decimal) + '°' + (lon < 0.0 ? 'W' : 'E')
    } else if (rangeType === 1) {
      result.latitude = lat.toFixed(decimal) + '°'
      result.longitude = lon.toFixed(decimal) + '°'
    } else if (rangeType === 2) {
      result.latitude = lat.toFixed(decimal) + '°'
      result.longitude = (lon < 0 ? 360 + lon : lon).toFixed(decimal) + '°'
    }

    if (Cesium.defined(optionsDefaulted.height)) {
      result.altitude =
        (optionsDefaulted.height > 0 ? optionsDefaulted.height : 0).toFixed(2) +
        (Cesium.defined(optionsDefaulted.errorBar)
          ? '±' + Math.round(optionsDefaulted.errorBar)
          : '') +
        'm'
    } else {
      result.altitude = ''
    }

    return result
  }
}
