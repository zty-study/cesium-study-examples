export type MaybeCoordinates =
  | (number | string)[]
  | {
      longitude: number | string
      latitude: number | string
      height?: number | string
    }
