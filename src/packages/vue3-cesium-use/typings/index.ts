import * as Cesium from 'cesium'

export type MaybeCoordinates =
  | (number | string)[]
  | {
      longitude: number | string
      latitude: number | string
      height?: number | string
    }

export type AnyFunction<T = void> = (...args: any[]) => T

export interface AnyObject {
  [propName: string]: any
}

export type Nullable<T> = T | null

export type CesiumPosition =
  | Cesium.Cartesian3
  | Cesium.CompositePositionProperty
  | Cesium.ConstantPositionProperty
  | Cesium.SampledPositionProperty
  | Cesium.TimeIntervalCollectionPositionProperty
  | Cesium.CallbackProperty
