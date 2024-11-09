import { ref } from 'vue'
import * as Cesium from 'cesium'
import { getViewer } from '../useViewer'

export interface ModelOptions {
  id: string
  position: Cesium.Cartesian3
  model: Cesium.ModelGraphics | Cesium.ModelGraphics.ConstructorOptions
}

export const useModel = () => {
  const modelList = ref<ModelOptions[]>([])

  const viewer = getViewer()
  const modelCollection = new Cesium.CustomDataSource('models')
  viewer.dataSources.add(modelCollection)

  // 监听modelList的变化, 移除modelList中被删除的model，添加modelList中新增的model
  watch(modelList, (newVal, oldVal) => {
    const removedModels = oldVal.filter((model) => !newVal.includes(model))
    removedModels.forEach(({ id }) => {
      modelCollection.entities.removeById(id)
    })
    newVal.forEach((modelOptions) => {
      if (!modelList.value.find((m) => m.id === modelOptions.id)) return
      const model = new Cesium.Entity({
        id: modelOptions.id,
        position: modelOptions.position,
        model: modelOptions.model
      })
      modelCollection.entities.add(model)
    })
  })

  return { modelList }
}
