import type { RouteRecordRaw } from 'vue-router'
import CesiumLayout from '@/layouts/cesium/index.vue'

export const CESIUM_ROUTES: RouteRecordRaw[] = [
  {
    path: 'cesium',
    name: 'admin-cesium',
    meta: { title: 'cesium', icon: 'earth-fill' },
    redirect: { name: 'admin-cesium-viewer' },
    component: CesiumLayout,
    children: [
      {
        path: 'viewer',
        name: 'admin-cesium-viewer',
        component: () => import('@/views/admin/cesium/pages/viewer.vue'),
        meta: { title: 'viewer' }
      }
    ]
  }
]
