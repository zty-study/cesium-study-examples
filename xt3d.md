## 波束雷达

`new xt3d.graphic.satellite.BeamRadar(options)`

```json
options: {
    position: [111,28,0], // 经纬度[x,y,z]
    style: {
        color: 'rgba(255,0,0,1)', // 颜色
        selectedColor: 'rgba(255,255,0,1)', // 选中颜色
        length: 2000, // 长度(米)
        radius: 100, // 半径(米)
        heading: 0, // 偏航角(角度)
        pitch: 0, // 俯仰角(角度)
        roll: 0, // 滚转角(角度)
    }
}
```

### 方法
1. `setTargetPosition(position, calLength)`: 设置目标位置,自动指向目标点
    - `position: [x,y,z]` // 目标点经纬度
    - `calLength: Boolean` // 是否自动计算长度


## 椎体雷达

`new xt3d.graphic.satellite.CylinderRadar(options)`

```json
options: {
    position: [111,28,0], // 经纬度[x,y,z]
    style: {
        color: 'rgba(255,0,0,1)', // 颜色
        selectedColor: 'rgba(255,255,0,1)', // 选中颜色
        angle: 22, // 张角(角度)
        radius: 100, // 半径(米)
        heading: 0, // 偏航角(角度)
        pitch: 0, // 俯仰角(角度)
        roll: 0, // 滚转角(角度)
        lineShow: true, // 是否显示线
        lineColor: 'rgba(228, 228, 228, 0.5)', // 线颜色
        topShow: true, // 是否显示顶盖
        topline: true, // 是否显示顶盖线
    }
}
```


## 探测雷达

`new xt3d.graphic.satellite.ProbeRadar(options)`

```json
options: {
    position: [111,28,0], // 经纬度[x,y,z]
    style: {
        color: 'rgba(255,0,0,1)', // 颜色
        selectedColor: 'rgba(255,255,0,1)', // 选中颜色
        length: 2, // 长度(米)
        radius: 100, // 半径(米)
        heading: 0, // 偏航角(角度)
        pitch: 0, // 俯仰角(角度)
        roll: 0, // 滚转角(角度)
        repeat: 20, // 纹理重复次数
        thickness: 0.3 // 纹理宽度
    }
}
```


## 范围雷达


`new xt3d.graphic.satellite.RangeRadar(options)`

```json
options: {
    position: [111,28,0], // 经纬度[x,y,z]
    style: {
        color: 'rgba(255,0,0,1)', // 颜色
        selectedColor: 'rgba(255,255,0,1)', // 选中颜色
        lineShow: true, // 是否显示线
        lineColor: 'rgba(228, 228, 228, 0.5)', // 线颜色
        topShow: true, // 是否显示顶盖
        topline: true, // 是否显示顶盖线
        heading: 0, // 偏航角(角度)
        pitch: 0, // 俯仰角(角度)
        roll: 0, // 滚转角(角度)
        outerRadius: 100, // 外半径(米)
        innerRadius: 50, // 内半径(米)
        startFovH: 0, // 水平起始角度(度)
        endFovH: 360, // 水平结束角度(度)
        startFovV: 0, // 垂直起始角度(度)
        endFovV: 90, // 垂直结束角度(度)
    }
}
```


## 雷达四棱锥体

`new xt3d.graphic.satellite.RectPyramid(options)`

```json
options: {
    position: [111,28,0], // 经纬度[x,y,z]
    style: {
        color: 'rgba(255,0,0,1)', // 颜色
        selectedColor: 'rgba(255,255,0,1)', // 选中颜色
        lineShow: true, // 是否显示线
        lineColor: 'rgba(228, 228, 228, 0.5)', // 线颜色
        topShow: true, // 是否显示顶盖
        topline: true, // 是否显示顶盖线
        heading: 0, // 偏航角(角度)
        pitch: 0, // 俯仰角(角度)
        roll: 0, // 滚转角(角度)
        radius: 100, // 半径(米)
        hAngle: 22, // 水平张角
        vAngle: 22, // 垂直张角
    }
}
```


## 相控雷达
`new xt3d.graphic.satellite.RectSensor(options)`

```json
options: {
    position: [111,28,0], // 经纬度[x,y,z]
    style: {
        color: 'rgba(255,0,0,1)', // 颜色
        selectedColor: 'rgba(255,255,0,1)', // 选中颜色
        lineShow: true, // 是否显示线
        lineColor: 'rgba(228, 228, 228, 0.5)', // 线颜色
        topShow: true, // 是否显示顶盖
        topline: true, // 是否显示顶盖线
        heading: 0, // 偏航角(角度)
        pitch: 0, // 俯仰角(角度)
        roll: 0, // 滚转角(角度)
        radius: 100, // 半径(米)
        hAngle: 22, // 水平张角
        vAngle: 22, // 垂直张角
        showScanPlane: true, // 是否显示扫描面
        scanPlaneColor: 'rgba(255, 0, 0, 1)', // 扫描面颜色
        scanPlaneRate: 10, // 扫描频率
    }
}
```