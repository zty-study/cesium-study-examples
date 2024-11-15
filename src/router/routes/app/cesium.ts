import type { RouteRecordRaw } from 'vue-router'
import CesiumLayout from '@/layouts/cesium/index.vue'

export const CESIUM_ROUTES: RouteRecordRaw[] = [
  {
    path: 'cesium',
    name: 'admin-cesium',
    meta: { title: 'cesium', icon: 'earth-fill', noKeepAlive: true },
    redirect: { name: 'admin-cesium-viewer' },
    component: CesiumLayout,
    children: [
      {
        path: 'viewer',
        name: 'admin-cesium-viewer',
        component: () => import('@/views/admin/cesium/pages/viewer.vue'),
        meta: { title: '场景视图' }
      },
      {
        path: 'boundary',
        name: 'admin-cesium-boundary',
        component: () => import('@/views/admin/cesium/pages/boundary.vue'),
        meta: { title: '国境边界' }
      },
      {
        path: 'model',
        name: 'admin-cesium-model',
        component: () => import('@/views/admin/cesium/pages/model.vue'),
        meta: { title: '模型加载' }
      },
      {
        path: 'ant',
        name: 'admin-cesium-ant',
        component: () => import('@/views/admin/cesium/pages/ant/index.vue'),
        meta: { title: '天线侦查' }
      }
    ]
  }
]
