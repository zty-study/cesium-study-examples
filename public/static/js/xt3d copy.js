function a0_0x3b79(t, e) {
  const i = a0_0x270b()
  return (a0_0x3b79 = function (t, e) {
    return i[(t -= 113)]
  })(t, e)
}
function a0_0x270b() {
  const t = [
    'dynamicHeight',
    '_cartesian3',
    'updateMaterial',
    'label',
    'json',
    'CylinderGlowFlowWallType',
    '_measureEntities',
    'height: 10px; width: 10px; background: rgb(251, 11, 11); border-radius: 50%;border: 1px solid black;position:absolute;pointer-events: none;z-index:0',
    'ViewDomeMaterial',
    'theme',
    '_removeTooltip',
    'lineShow',
    '_pickCommand',
    '视椎体视屏图元',
    'position',
    '光晕线',
    '_addMoveAnchor',
    'spacing',
    'info-input1',
    'showDrivenLine',
    '_materialCache',
    '点击鼠标右键取消绘制',
    '_oldQuaternion',
    'log2',
    'upAxis',
    'fillColor',
    'mousemove',
    'clearRect',
    '_createOuterCurveCommand',
    "\n #ifdef GL_OES_standard_derivatives\n    #extension GL_OES_standard_derivatives : enable\n#endif\n\nuniform bool u_showIntersection;\nuniform bool u_showThroughEllipsoid;\n\nuniform float u_radius;\nuniform float u_xHalfAngle;\nuniform float u_yHalfAngle;\nuniform float u_normalDirection;\nuniform vec4 u_color;\n\nin vec3 v_position;\nin vec3 v_positionWC;\nin vec3 v_positionEC;\nin vec3 v_normalEC;\n\nvec4 getColor(float sensorRadius, vec3 pointEC)\n{\n    czm_materialInput materialInput;\n\n    vec3 pointMC = (czm_inverseModelView * vec4(pointEC, 1.0)).xyz;\n    materialInput.st = sensor2dTextureCoordinates(sensorRadius, pointMC);\n    materialInput.str = pointMC / sensorRadius;\n\n    vec3 positionToEyeEC = -v_positionEC;\n    materialInput.positionToEyeEC = positionToEyeEC;\n\n    vec3 normalEC = normalize(v_normalEC);\n    materialInput.normalEC = u_normalDirection * normalEC;\n\n    czm_material material = czm_getMaterial(materialInput);\n\n    material.diffuse = u_color.rgb;\n    material.alpha = u_color.a; \n    return mix(czm_phong(normalize(positionToEyeEC), material, czm_lightDirectionEC), vec4(material.diffuse, material.alpha), 0.4);\n\n}\n\nbool isOnBoundary(float value, float epsilon)\n{\n    float width = getIntersectionWidth();\n    float tolerance = width * epsilon;\n\n#ifdef GL_OES_standard_derivatives\n    float delta = max(abs(dFdx(value)), abs(dFdy(value)));\n    float pixels = width * delta;\n    float temp = abs(value); \n    return temp < tolerance && temp < pixels || (delta < 10.0 * tolerance && temp - delta < tolerance && temp < pixels);\n#else\n    return abs(value) < tolerance;\n#endif\n}\n\nvec4 shade(bool isOnBoundary)\n{\n    if (u_showIntersection && isOnBoundary)\n    {\n        return getIntersectionColor();\n    }\n    return getColor(u_radius, v_positionEC);\n}\n\nfloat ellipsoidSurfaceFunction(vec3 point)\n{\n    vec3 scaled = czm_ellipsoidInverseRadii * point;\n    return dot(scaled, scaled) - 1.0;\n}\n\nvoid main()\n{\n    vec3 sensorVertexWC = czm_model[3].xyz;      // (0.0, 0.0, 0.0) in model coordinates\n    vec3 sensorVertexEC = czm_modelView[3].xyz;  // (0.0, 0.0, 0.0) in model coordinates\n \n    float positionX = v_position.x;\n    float positionY = v_position.y;\n    float positionZ = v_position.z;\n\n    vec3 zDir = vec3(0.0, 0.0, 1.0);\n    vec3 lineX = vec3(positionX, 0 ,positionZ);\n    vec3 lineY = vec3(0, positionY, positionZ);\n    float resX = dot(normalize(lineX), zDir);\n    if(resX < cos(u_xHalfAngle) - 0.0001){\n        discard;\n    }\n    float resY = dot(normalize(lineY), zDir);\n    if(resY < cos(u_yHalfAngle)- 0.0001){\n        discard;\n    }\n\n\n    float ellipsoidValue = ellipsoidSurfaceFunction(v_positionWC);\n\n    // Occluded by the ellipsoid?\n\tif (!u_showThroughEllipsoid)\n\t{\n\t    // Discard if in the ellipsoid\n\t    // PERFORMANCE_IDEA: A coarse check for ellipsoid intersection could be done on the CPU first.\n\t    if (ellipsoidValue < 0.0)\n\t    {\n            discard;\n\t    }\n\n\t    // Discard if in the sensor's shadow\n\t    if (inSensorShadow(sensorVertexWC, v_positionWC))\n\t    {\n\t        discard;\n\t    }\n    }\n\n    // Notes: Each surface functions should have an associated tolerance based on the floating point error.\n    bool isOnEllipsoid = isOnBoundary(ellipsoidValue, czm_epsilon3);\n    out_FragColor = u_color; \n\n}",
    'createFsShader',
    '/A61BA73E80614F0ABD07D8682C673B61.png',
    'removeChild',
    '_updateGeometry',
    '_createDom',
    'MaterialSupport',
    '_definitionChanged',
    'offset',
    '_clearEntity',
    '\nin vec3 position;\nin vec3 normal;\nin vec2 st;\nin vec4 color;\nin float batchId;\nout vec3 v_positionEC;\nout vec3 v_normalEC;\nout vec2 v_st;\nout vec4 v_pickColor;\nvoid main()\n{\n    v_positionEC = (czm_modelView * vec4(position, 1.0)).xyz;       // position in eye coordinates\n    v_normalEC = czm_normal * normal;                               // normal in eye coordinates\n    v_st = st;\n    v_pickColor = color;\n    gl_Position = czm_modelViewProjection * vec4(position, 1.0);\n}\n',
    'PixelFormat',
    '_geoJsonDataSource',
    'drawStart',
    '_waterPlane',
    'minScaleHeight',
    'u_intersectionWidth',
    'primitiveType',
    '燕尾直箭头',
    'Cartographic',
    'closedCurve',
    '_applayTexture',
    'scanPlaneRate',
    '墙体流动材质',
    '_positionsTemp',
    '_maskTexture',
    'fourOindices',
    '_inverseMatrix4',
    '_setSelectedStyle',
    '_selectedGraphic',
    '_removeSelectedDom',
    '_computeSquareShape',
    '_pickTooltip',
    '_playing',
    'selectedEnable',
    '_distance',
    'EventType',
    'UrlTemplate',
    'func',
    'EllipsoidSurfaceAppearance',
    'height',
    'rotation',
    'u_hiddenColor',
    '_setPositionsHeight',
    'LineString',
    'NeonPointFS',
    '_circleRadius',
    'scaleHeight',
    'Plane',
    'alpha',
    'board',
    'divideByScalar',
    'topPpositions',
    'ComponentDatatype',
    'EllipseGeometry',
    '_getArrowHeadPoints',
    'loadFromUrl',
    '_positionProperty',
    'HorizontalOrigin',
    'context',
    'WallCoolType',
    'fill',
    'lockView',
    '7299174FrLZUL',
    'defined',
    'index',
    'minimumHeights',
    'positions',
    'outerRadius',
    ' \n            in vec3 position;\n            in vec3 normal;\n            // attribute vec2 st;\n            // attribute float batchId;\n            \n            out vec3 v_positionEC;\n            out vec3 v_normalEC;\n            // varying vec2 v_st;\n           \n            void main()\n            { \n                \n                v_positionEC = (czm_modelView * vec4(position, 1.0)).xyz;       // position in eye coordinates\n                v_normalEC = czm_normal * normal;                               // normal in eye coordinates\n                // v_st = st;\n                gl_Position = czm_modelViewProjection * vec4(position, 1.0);\n            }',
    '_graphicType',
    'vid',
    'load',
    'u_lineColor',
    'StencilOperation',
    'scale',
    '_initEvent',
    '_neckWidthFactor',
    'wallLight',
    'globalAlpha',
    'area',
    'normalMatrix',
    '着色器特效霓虹点',
    'segmentH',
    '_runing',
    '\nprecision highp float;\nprecision highp int;\nuniform float u_speed;\nuniform vec3 u_color;\nuniform float u_time;\nuniform float u_glow; \nin vec3 v_positionEC;\nin vec3 v_normalEC;\nin vec2 v_st;\nmat4 mat  = mat4 ( vec4 ( 1.0 , 0.0 , 0.0 , 0.0 ),\n      vec4 ( 0.0 , 1.0 , 0.0 , 0.0 ),\n      vec4 ( 0.0 , 0.0 , 1.0 , 0.0 ),\n      vec4 ( 0.0 , 0.0 , 0.0 , 1.0 ) );\n\nvec2 pos;\n\nvec4 col = vec4 ( 0., 0., 0., 1000.0 );\n  \nvoid Line2 ( vec2 a, vec2 b );\nvoid Line2 ( vec2 a, vec2 b ) {\n  float d = distance ( pos , a ) + distance ( pos , b ) - distance ( a , b ) + 1e-5;\n  col += max ( 1. - pow ( d * 14. , 0.1 ) , -0.01 );\n}\n\nvoid Line4 ( vec4 a, vec4 b );\nvoid Line4 ( vec4 a, vec4 b ) {\n  a = mat * a;\n  a.xyz /= 1.5 + a.w * 2.;\n  b = mat * b;\n  b.xyz /= 1.5 + b.w * 2.;\n  Line2 ( a.xy , b.xy );\n}\n\nvoid Point ( vec4 p );\nvoid Point ( vec4 p ) {\n  p = mat * p;\n  p.xyz /= 1.5 + p.w * 2.;\n  \n  float d = distance ( pos , p.xy );\n  \n  if ( d < .3 )\n  if ( p.z < col.a ) {\n    col.b += max ( 1.0 - pow ( d * 5.0 , .1 ) , 0.0 );\n  }\n}\n\nvoid Rotate ( float angle, float d1, float d2, float d3, float d4);\nvoid Rotate ( float angle, float d1, float d2, float d3, float d4) {\n  float c = cos (angle), s = sin (angle);\n  mat *= mat4 ( vec4 (  c*d1+(1.-d1),  s * d2 * d1 , -s * d3 * d1 ,  s * d4 * d1 ),\n          vec4 ( -s * d1 * d2 ,  c*d2+(1.-d2),  s * d3 * d2 , -s * d4 * d2 ),\n          vec4 (  s * d1 * d3 , -s * d2 * d3 ,  c*d3+(1.-d3),  s * d4 * d3 ),\n          vec4 ( -s * d1 * d4 ,  s * d2 * d4 , -s * d3 * d4 ,  c*d4+(1.-d4)) );\n}\n\nvoid main() {\n\n  vec3 positionToEyeEC = -v_positionEC;\n  vec3 normalEC = normalize(v_normalEC);\n  normalEC = faceforward(normalEC, vec3(0.0, 0.0, 1.0), -normalEC);\n\n  float time = czm_frameNumber / 60.0;\n  time = u_speed * time;\n  pos = v_st - 0.5;\n  float pi = 3.141592;\n  \n  \n\n  Rotate ( pi / 3.,      0.0, 1.0, 1.0, 0.0 );\n  Rotate ( time,      1.0, 1.0, 0.0, 0.0 );\n\n  vec4 point1 = vec4 ( 0., 0., .2, -.2 );\n  vec4 point2 = vec4 ( .2, .2,-.2,-.2 );\n  vec4 point3 = vec4 ( -.2, .2,-.2,-.2 );\n  vec4 point4 = vec4 ( .2,-.2,-.2,-.2 );\n  vec4 point5 = vec4 ( -.2,-.2,-.2,-.2 );\n\n  Line4 ( point1, point2 );\n  Line4 ( point1, point3 );\n  Line4 ( point1, point4 );\n  Line4 ( point1, point5 );\n\n  Line4 ( point2, point3 );\n  Line4 ( point2, point4 );\n  Line4 ( point3, point5 );\n  Line4 ( point4, point5 );\n  Point ( point1 );\n  Point ( point2 );\n  Point ( point3 );\n  Point ( point4 );\n  Point ( point5 );\n  \n  out_FragColor = vec4( col.xyz * u_glow, 1.0 );\n out_FragColor =vec4(out_FragColor.rgb + u_color, out_FragColor.b* out_FragColor.b);\n}',
    'imageLabel',
    'drivenTime',
    'WallLightMaterialProperty',
    'rgba(76, 16, 81, 1)',
    '着色器特效魔法球',
    '_videoPlane',
    '_headWidthFactor',
    'param',
    'DeveloperError',
    'background-image:linear-gradient(135deg, transparent 25px, ',
    'MeasureAreaVertex',
    'auto',
    '_postProcessStage',
    'coordinates',
    'graphicId',
    'TRIANGLES',
    'fixPointCount',
    '闪烁点',
    "\n// change this to get different explosions :)\n#define EXPLOSION_SEED 0.\n\n// uncomment this to get a cross section view\n//#define CROSS_SECTION\n\n// the bounding sphere of the explosion. this is less general but means that\n// ray cast is only performed for nearby pixels, and raycast can begin from the sphere\n// (instead of walking out from the camera)\nfloat expRadius;\nvec3 expCenter;\n\n//iq's LUT 3D noise\nfloat noise( in vec3 x )\n{\n    vec3 f = fract(x);\n    vec3 p = x - f; // this avoids the floor() but doesnt affect performance for me.\n    f = f*f*(3.0-2.0*f);\n    \n    vec2 uv = (p.xy+vec2(37.0,17.0)*p.z) + f.xy;\n    vec2 rg = texture( iChannel0, (uv+ 0.5)/256.0, 0.0 ).yx;\n    return mix( rg.x, rg.y, f.z );\n}\n\n// assign colour to the media\nvec3 computeColour( float density, float radius )\n{\n// these are almost identical to the values used by iq\n\n// colour based on density alone. gives impression of occlusion within\n// the media\nvec3 result = mix( 1.1*vec3(1.0,0.9,0.8), vec3(0.4,0.15,0.1), density );\n\n// colour added for explosion\nvec3 colBottom = 3.1*vec3(1.0,0.5,0.05);\nvec3 colTop = 2.*vec3(0.48,0.53,0.5);\nresult *= mix( colBottom, colTop, min( (radius+.5)/1.7, 1.0 ) );\n\nreturn result;\n}\n\n// maps 3d position to colour and density\nfloat densityFn( in vec3 p, in float r, out float rawDens, in float rayAlpha )\n{\n// density has dependency on mouse y coordinate (linear radial ramp)\nfloat mouseIn = 0.85;\nfloat mouseY = 1.0 - mouseIn;\n    float den = -0.1 - 1.5*r*(4.*mouseY+.5);\n    \n// offset noise based on seed\n    float t = EXPLOSION_SEED;\n    vec3 dir = vec3(0.,1.,0.);\n    \n    // participating media    \n    float f;\n    vec3 q = p - dir* t; f  = 0.50000*noise( q );\nq = q*2.02 - dir* t; f += 0.25000*noise( q );\nq = q*2.03 - dir* t; f += 0.12500*noise( q );\nq = q*2.01 - dir* t; f += 0.06250*noise( q );\nq = q*2.02 - dir* t; f += 0.03125*noise( q );\n\n// add in noise with scale factor\nrawDens = den + 4.0*f;\n\n    den = clamp( rawDens, 0.0, 1.0 );\n    \n// thin out the volume at the far extends of the bounding sphere to avoid\n// clipping with the bounding sphere\nden *= 1.-smoothstep(0.8,1.,r/expRadius);\n\n#ifdef CROSS_SECTION\nden *= smoothstep(.0,.1,-p.x);\n#endif\n\nreturn den;\n}\n\nvec4 raymarch( in vec3 rayo, in vec3 rayd, in float expInter)\n{\n    vec4 sum = vec4( 0.0 );\n    \n    float step = 0.075;\n    \n    // dither start pos to break up aliasing\nvec3 pos = rayo + rayd * (expInter + step*texture( iChannel0, fragCoord.xy/iResolution.x ).x);\n\n    for( int i=0; i<25; i++ )\n    {\n        if( sum.a > 0.99 ) continue;\n    \n    float radiusFromExpCenter = length(pos - expCenter);\n    \n    if( radiusFromExpCenter > expRadius+0.01 ) continue;\n    \n    float dens, rawDens;\n    \n        dens = densityFn( pos, radiusFromExpCenter, rawDens, sum.a );\n    \n    vec4 col = vec4( computeColour(dens,radiusFromExpCenter), dens );\n    \n    // uniform scale density\n    col.a *= 0.6;\n    \n    // colour by alpha\n    col.rgb *= col.a;\n    \n    // alpha blend in contribution\n    sum = sum + col*(1.0 - sum.a);  \n    \n    // take larger steps through negative densities.\n    // something like using the density function as a SDF.\n    float stepMult = 1. + 2.5*(1.-clamp(rawDens+1.,0.,1.));\n    \n    // step along ray\n    pos += rayd * step * stepMult;\n    }\n\n    return clamp( sum, 0.0, 1.0 );\n}\n\n// iq's sphere intersection\nfloat iSphere(in vec3 ro, in vec3 rd, in vec4 sph)\n{\n//sphere at origin has equation |xyz| = r\n//sp |xyz|^2 = r^2.\n//Since |xyz| = ro + t*rd (where t is the parameter to move along the ray),\n//we have ro^2 + 2*ro*rd*t + t^2 - r2. This is a quadratic equation, so:\nvec3 oc = ro - sph.xyz; //distance ray origin - sphere center\n\nfloat b = dot(oc, rd);\nfloat c = dot(oc, oc) - sph.w * sph.w; //sph.w is radius\nfloat h = b*b - c; // delta\nif(h < 0.0) \n    return -1.0;\nfloat t = (-b - sqrt(h)); //Again a = 1.\n\nreturn t;\n}\n\nvec3 computePixelRay( in vec2 p, out vec3 cameraPos )\n{\n    // camera orbits around explosion\n\n    float camRadius = 3.8;\n// use mouse x coord\nfloat a = iTime*20.;\nfloat theta = -(a-iResolution.x)/80.;\n    float xoff = camRadius * cos(theta);\n    float zoff = camRadius * sin(theta);\n    cameraPos = vec3(xoff,expCenter.y,zoff);\n    \n    // camera target\n    vec3 target = vec3(0.,expCenter.y,0.);\n    \n    // camera frame\n    vec3 fo = normalize(target-cameraPos);\n    vec3 ri = normalize(vec3(fo.z, 0., -fo.x ));\n    vec3 up = normalize(cross(fo,ri));\n    \n    // multiplier to emulate a fov control\n    float fov = .5;\n\n    // ray direction\n    vec3 rayDir;\n    if(u_param.y == 1.0){\n        rayDir = normalize(fo + fov*p.x*ri + fov*p.y*up);\n    }else{\n        rayDir = normalize(fo + fov*p.y*ri + fov*p.x*up);\n    }\n\n    return rayDir;\n}\n\nvoid main()\n{\n// get aspect corrected normalized pixel coordinate\nvec2 p = 2.*v_st.xy - vec2(1., 1.);\n    \nexpRadius = 1.75;\nexpCenter = vec3(0.,expRadius,0.);\n\nvec3 rayDir, cameraPos;\nrayDir = computePixelRay( p, cameraPos );\n\nvec4 col = vec4(0.);\n\n// does pixel ray intersect with exp bounding sphere?\nfloat boundingSphereInter = iSphere( cameraPos, rayDir, vec4(expCenter,expRadius) );\nif( boundingSphereInter > 0. )\n{\n    // yes, cast ray\n    col = raymarch( cameraPos, rayDir, boundingSphereInter ); \n}\n\n    // smoothstep final color to add contrast\n    fragColor.xyz = col.xyz*col.xyz*(3.0-2.0*col.xyz);\n}\n\n",
    '_drawEnd',
    '_scanUniforms',
    '_viewShadowMap',
    'MagicRingFS',
    'preload',
    '失败，数据为空！',
    'fromRotationX',
    '_update',
    '动画点',
    '#001e0f',
    'translucent',
    'strokeRect',
    '_updateVirtualCamera',
    'shaders',
    'ShieldFS',
    '_locationPoint',
    'vertexShaderSource',
    'PolylineFlow',
    '_removeHook',
    '_color',
    'clock',
    'editGraphic',
    'margin:5px 1px;',
    'pow',
    'bottomRadius',
    '#ffff00',
    '_colorFramebuffer',
    'shadowMapMatrix',
    'addHWAnchor',
    'zoy',
    '_far',
    '_registerEvents',
    'TRANSLUCENT',
    'semiMajorAxis',
    'px;',
    '着色器特效科技球',
    'fromCache',
    'crossorigin',
    '导入失败！类型错误，请确保是通过toGeoJson（）方法导出的数据',
    '拖拽平移',
    'MODIFY_MATERIAL',
    '_hpr',
    '_handleMoveEvent',
    'routeEnd',
    '扫描材质_1',
    'primitive',
    'reflectWater',
    'lightCamera',
    '_removeLabelDom',
    'webkitHidden',
    'code',
    '_computeRepeat',
    '_createVideoEle',
    'KEEP',
    'CircleGraphic',
    '_beforDiv',
    '_minimumHeights',
    'NeonLightFS',
    'dynamic-border-label-container',
    'true',
    'rgba(13, 234, 168, 1)',
    '_createPolyline',
    'tailedFineArrow',
    'fromElements',
    'CLAMPED',
    '不支持的文件类型！',
    '小数位：',
    '_minHeight',
    '模型图层',
    'svg图元',
    'FLOAT',
    '_initChart',
    '  \n        #define time iTime*1.25\n        #define p0 0.5, 0.5, 0.5,  0.5, 0.5, 0.5,  1.0, 1.0, 1.0,  0.0, 0.33, 0.67\t\n\n        const float numParticles = 25.;\n        const float numRings = 5.;\n        const float offsetMult = 30.;\n        const float tau = 6.23813;\n\n        vec3 palette( in float t, in float a0, in float a1, in float a2, in float b0, in float b1, in float b2,\n                    in float c0, in float c1, in float c2,in float d0, in float d1, in float d2)\n        {\n            return vec3(a0,a1,a2) + vec3(b0,b1,b2)*cos( tau*(vec3(c0,c1,c2)*t+vec3(d0,d1,d2)) );\n        }\n\n        vec3 particleColor(vec2 uv, float radius, float offset, float periodOffset)\n        {\n            vec3 color = palette(.4 + offset / 4., p0);\n            uv /= pow(periodOffset, .75) * sin(periodOffset * iTime) + sin(periodOffset + iTime);\n            vec2 pos = vec2(cos(offset * offsetMult + time + periodOffset),\n                    sin(offset * offsetMult + time * 5. + periodOffset * tau));\n            \n            float dist = radius / distance(uv, pos);\n            return color * pow(dist, 2.) * 1.75;\n        } \n\n        void main()\n        {\n            vec2 uv = 2.*v_st.xy - vec2(1., 1.);\n            uv *= 3.45;\n\n            fragColor.rgb = vec3(0.);\n            \n            for (float n = 0.; n <= numRings; n++)\n            {\n                for (float i = 0.; i <= numParticles; i++) {\n                fragColor.rgb += particleColor(uv, .03, i / numParticles, n / 2.);\n            }\n            }\n\n            fragColor.rgb = fragColor.rgb * u_glow;\n    }',
    '_shadowMap',
    'autoplay',
    '_normal',
    '_covertMaximumHeights',
    'moveHeight',
    'FeatureCollection',
    '_highlightAnchor',
    'LEFT',
    'minimumClusterSize',
    'classList',
    'ne-resize',
    '_setPlane',
    'handleMouseMove',
    '_geometryInstance',
    'fragmentShader',
    '_handleDeActivate',
    'height: 10px; width: 10px; background: yellow; border-radius: 50%;border: 1px solid black;position:absolute;pointer-events: none;z-index:0',
    'Cesium3DTilePassState',
    '_layer',
    'asin',
    'mousedown',
    'files',
    'labelField',
    'wall',
    '\n              void fragmentMain(FragmentInput fsInput, inout czm_modelMaterial material)\n              {\n                vec4 position = czm_inverseModelView * vec4(fsInput.attributes.positionEC,1); // 位置 \n                vec4  color = vec4(',
    '_positions',
    'u_alpha',
    'semiMinorAxis',
    'getBoundingClientRect',
    '_contentDoms',
    'YELLOW',
    'pre-topCard-list-item-line',
    'ns-resize',
    '_tooltipContainer',
    '_scanPlaneBackCommand',
    'secondsDifference',
    'intensitys',
    'circleWave',
    'polylineVolume',
    'multiplyByPoint',
    'out_FragColor = vec4(out_FragColor.rgb + u_color, out_FragColor.r + out_FragColor.r * 0.1);',
    '_addData',
    'endDraw',
    'EllipsoidGeometry',
    'fromUrl',
    'cullFaceType',
    'handler',
    '_computeCircleShape',
    'trunc',
    'vAngle',
    '动效边框',
    '/assets/',
    '_particle',
    '_addTooltip',
    'shaderEffetGlowPoint',
    '_backFaceRS',
    'matrixAxis',
    'sampleHeight',
    'ColorfulPointFS',
    'currentTime',
    '视频图元',
    '_computeLineInfo',
    'PostProcessStage',
    '_createText',
    'tilesLoaded',
    '着色器特效浮游点',
    'imageryLayer',
    'vite-plugin-css-injected-by-js',
    'map',
    'shaderEffetBase',
    'se-resize',
    'setPositions',
    '_createImage_1',
    'getDate',
    'color.rgb',
    'hAngle',
    '_editMode',
    '_sectorSegmentLineCommand',
    'rgba(255, 255, 255, 0.2)',
    'LINEAR',
    'xt3d_token',
    'IntersectionTests',
    '_removeAllChild',
    'fontFamily',
    'pause',
    'position:absolute;left:0px;top:0px;z-index:-9999;overflow:hidden;',
    '_vertexShader',
    'maxScaleWidth',
    '无效的graphicType！ ：',
    'xt3d-divgraphic',
    'uniform sampler2D u_map;\nuniform bool u_invert;\nuniform float u_alpha;\n\nin vec2 v_textureCoordinates;\nvoid main(){\n    vec2 uv = v_textureCoordinates;\n    vec4 color = texture( u_map , uv );\n    //计算灰度\n    float grayVal =( color.r + color.g + color.b ) / 3.;\n    //灰度反转\n    if(u_invert){\n        grayVal = 1. - grayVal;\n    }\n    //应用前景透明度，以便和背景混合\n    float alpha = color.a * u_alpha;\n    out_FragColor=vec4( vec3( grayVal ) , alpha );\n} ',
    'speed',
    '_sp',
    '_pointBias',
    '_createResultLabel',
    '_handleLeftClick',
    'precision',
    'DebugCameraPrimitive',
    '.xt3d-component-animate-marker_boder',
    'destination',
    '_createRightCrossSectionCommand',
    'computeNormal',
    '\n                    uniform float u_speed;\n                    uniform vec3 u_color;\n                    uniform float u_time;\n                    uniform float u_glow;\n                    in vec2 v_st;\n                    precision highp float;\n                    precision highp int;\n                    #define pi 3.1415926535\n                    #define PI2RAD 0.01745329252\n                    #define TWO_PI (2. * PI)\n\n\n                    float time;\n                    float rands(float p){\n                        return fract(sin(p) * 10000.0);\n                    }\n                    float noise(vec2 p){\n                        float t = time / 20000.0;\n                        if(t > 1.0) t -= floor(t);\n                        return rands(p.x * 14. + p.y * sin(t) * 0.5);\n                    }\n                    vec2 sw(vec2 p){\n                        return vec2(floor(p.x), floor(p.y));\n                    }\n                    vec2 se(vec2 p){\n                        return vec2(ceil(p.x), floor(p.y));\n                    }\n                    vec2 nw(vec2 p){\n                        return vec2(floor(p.x), ceil(p.y));\n                    }\n                    vec2 ne(vec2 p){\n                        return vec2(ceil(p.x), ceil(p.y));\n                    }\n                    float smoothNoise(vec2 p){\n                        vec2 inter = smoothstep(0.0, 1.0, fract(p));\n                        float s = mix(noise(sw(p)), noise(se(p)), inter.x);\n                        float n = mix(noise(nw(p)), noise(ne(p)), inter.x);\n                        return mix(s, n, inter.y);\n                    }\n                    float fbm(vec2 p){\n                        float z = 2.0;\n                        float rz = 0.0;\n                        vec2 bp = p;\n                        for(float i = 1.0; i < 6.0; i++){\n                            rz += abs((smoothNoise(p) - 0.5)* 2.0) / z;\n                            z *= 2.0;\n                            p *= 2.0;\n                        }\n                        return rz;\n                    }\n                    void main()\n                    {\n                        vec2 vUv = v_st;\n                        time = czm_frameNumber / 60.0;\n                        time = u_speed * time;\n                        vec2 uv = vUv;\n                        vec2 uv2 = vUv;\n                        if (uv.y < 0.5) {\n                            discard;\n                        }\n                        uv *= 4.;\n                        float rz = fbm(uv);\n                        uv /= exp(mod(time * 2.0, pi));\n                        rz *= pow(15., 1.0);\n                        vec4 color = mix(vec4(u_color, 1.0) / rz, vec4(u_color, 0.1), 0.5);\n                        if (uv2.x < 0.05) {\n                            color = mix(vec4(u_color, 0.1), color, uv2.x / 0.05);\n                        }\n                        if (uv2.x > 0.95){\n                            color = mix(color, vec4(u_color, 0.1), (uv2.x - 0.95) / 0.05);\n                        }\n\n                        float alpha = color.a * 1.5;\n                        vec3 diffuse = max(color.rgb + color.rgb * alpha, color.rgb);\n\n                        out_FragColor = vec4(diffuse * u_glow,alpha);\n                    }\n    ',
    '_createEntiy',
    'http://www.xt3d.online',
    '_isAdd',
    'amplitude',
    'u_normalDirection',
    'moveEnd',
    'push',
    'file',
    'LightingModel',
    'heightRatio',
    'globe',
    'WallGradientType',
    'yHalfAngle',
    '_infiniteProjection',
    '_setOuterCylinderAppearence',
    'Water',
    'face',
    'createUniformMap',
    ' 50.1%, transparent 50%)',
    '_sectorSegmentLineVA',
    'PerInstanceColorAppearance',
    'assaultDirection',
    'CallbackProperty',
    '_addBillboard',
    '_precision',
    '_headingPitchRange',
    '_scanPlaneFrontCommand',
    '_leftDown',
    'img',
    'keys',
    'computeModelMatrix',
    'getType',
    'fromCartesian',
    'shaderEffetColorfulPoint',
    'logarithmicDepthBuffer',
    'display',
    'blendingEnabled',
    'circle',
    '3123270lShmdb',
    '火焰粒子',
    'createPlayer',
    'height: 58px; width: 140px;',
    'CircleScanMaterialProperty_3',
    'invert',
    '个点，按下鼠标左键确定第',
    'RealFlameFS',
    '_bottomRing',
    'clusterEvent',
    'entities',
    'endScale',
    'halfFloatingPointTexture',
    '拖拽改变水平位置',
    '_getCirclePosition',
    'slice',
    'aspectRatio',
    'control',
    'timeSum',
    'getById',
    '_neckHeightFactor',
    'GradientRingFS',
    'tilesetPassState',
    'head',
    '图元不存在！',
    '_headHeightFactor',
    'pass',
    '视频投射',
    '_tilesetOpts',
    'endFrame',
    'geometry',
    'createPostProcessStage',
    '_outFragColorBody',
    '当前测量类型：空间距离',
    '.title',
    'bounceMarker',
    '_init',
    'startFovH',
    '_handleRightClick',
    '_clearEditObject',
    '_headAngle',
    'rgba(228, 228, 228, 0.5)',
    'circleScan_4',
    'bottom',
    'detachEvent',
    'geometryType',
    'execute',
    'tooltip',
    'MouseEvents',
    '_createCylinderInstance',
    'home',
    'gradient',
    'close',
    '数据格式有问题！',
    'lightPositionEC',
    'brightenDiv',
    'outWidth',
    'SvgGraphic',
    'CircleScanType_1',
    'hasChildNodes',
    'offsetWidth',
    '点击将内容复制到粘贴板',
    '\nczm_material czm_getMaterial(czm_materialInput materialInput)\n{\n    czm_material material = czm_getDefaultMaterial(materialInput);\n    vec2 st = materialInput.st;\n    vec4 colorImage = texture(image, vec2(fract(st.t), st.t));\n    material.alpha = colorImage.a;\n    material.diffuse = colorImage.rgb * 2.5 ;\n    return material;\n}',
    'thickness',
    'flowDir',
    'CylinderGlowCircle',
    'toString',
    '_createEditAnchors',
    'innerHTML',
    '_getSceneCenterPosition',
    'ShaderProgram',
    '#FFFF00',
    ' 25px)',
    'gatheringPlace',
    'drawEnd',
    'HeadingPitchRoll',
    'normalOffsetScale',
    'imageUrl',
    'hidden',
    'shaderProgram',
    '_start',
    'stroke',
    'graphicLayer',
    '着色器特效发光点',
    '_tempPositions',
    '\nczm_material czm_getMaterial(czm_materialInput materialInput) {\n    czm_material material = czm_getDefaultMaterial(materialInput);\n    material.diffuse = 1.5 * color.rgb;\n    vec2 st = materialInput.st;\n    vec3 str = materialInput.str;\n    float dis = distance(st, vec2(0.5, 0.5));\n    float per = fract(czm_frameNumber * speed / 1000.0);\n    if (abs(str.z) > 0.001) {\n      discard;\n    }\n    if (dis > 0.5) {\n      discard;\n    } else {\n      float perDis = 0.5 / count;\n      float disNum;\n      float bl = .0;\n      for (int i = 0; i <= 9; i++) {\n        if (float(i) <= count) {\n          disNum = perDis *float(i) - dis + per / count;\n          if (disNum > 0.0) {\n            if (disNum < perDis) {\n              bl = 1.0 - disNum / perDis;\n            } else if(disNum - perDis < perDis) {\n              bl = 1.0 - abs(1.0 - disNum / perDis);\n            }\n            material.alpha = pow(bl, gradient);\n          }\n        }\n      }\n    }\n    return material;\n  }',
    'Cartesian2',
    'particleSize',
    'acos',
    'rgba(227,108,9, 0.5)',
    '\n               in vec3 position3DHigh;\n               in vec3 position3DLow;\n               in vec3 normal;\n               in vec2 st;\n               in float batchId;\n               \n               out vec3 v_positionEC;\n               out vec3 v_normalEC;\n               out vec2 v_st;\n               \n               uniform mat4 reflectorProjectionMatrix;\n               uniform mat4 reflectorViewMatrix;\n               uniform mat4 reflectMatrix;\n               out vec4 v_worldPosition;  // 世界坐标\n               out vec4 v_uv;             // 纹理坐标\n               \n               \n               \n               void main()\n               {\n                   vec4 p = czm_computePosition();\n               \n                   v_positionEC = (czm_modelViewRelativeToEye * p).xyz;      // position in eye coordinates\n                   v_normalEC = czm_normal * normal;                         // normal in eye coordinates\n                   v_st = st;\n               \n               \n                   mat4 modelView = reflectorViewMatrix * reflectMatrix * czm_model;\n                   modelView[3][0] = 0.0;\n                   modelView[3][1] = 0.0;\n                   modelView[3][2] = 0.0;\n                   v_uv = reflectorProjectionMatrix * modelView * p;\n                   vec4 positionMC = vec4( position3DHigh + position3DLow, 1.0 );\n                   v_worldPosition = czm_model * positionMC;\n               \n                   gl_Position = czm_modelViewProjectionRelativeToEye * p;\n          }',
    'latitude',
    '_setEditMode',
    'parse',
    'u_color',
    'glowCylinder',
    'startDraw',
    'startTime',
    'maximumScale',
    '_createTopOutlineGeometry',
    'shadows',
    'LineGraphic',
    'PolyElevationContourMaterialProperty',
    'ConstellationChainFS',
    'slices',
    'faceColor',
    'defaultVS',
    'UNLIT',
    '_userData',
    'colorBlendMode',
    '点击鼠标左键确定位置，点击鼠标右键取消绘制',
    '_rotation',
    'shadowMapTexelSizeDepthBiasAndNormalShadingSmooth',
    'east',
    '_topGeometry',
    'modelMatrix',
    'materialType',
    'pixelSize',
    'normal',
    '_getTempPoint4',
    'u_speed',
    'multiplyByTranslation',
    '_setNodeStyle',
    'ShaderSource',
    '_image',
    'readAsText',
    '_createShadowMap',
    'fromTranslationQuaternionRotationScale',
    'TranslationRotationScale',
    'drawingBufferHeight',
    'styleType',
    'shaderSource',
    '_stop',
    '_getCenterPosition',
    'beginFrame',
    'shaderEffetPetal',
    'measureStart',
    '_viewer',
    'handleMouseMoveALT',
    '_createFlvVideoEle',
    '_windowPosition',
    'root',
    'assign',
    '\n            in vec3 position3DHigh;\n            in vec3 position3DLow;\n            in vec2 st;\n            in float batchId;\n            uniform sampler2D image_0; \n            uniform float maxHeight_1;\n            uniform float heightRatio_2;\n            out vec3 v_positionMC;\n            out vec3 v_positionEC;\n            out vec2 v_st; \n            in vec3 normal;\n            out vec3 v_normalEC;\n            void main(){\n                vec4 p = czm_computePosition();\n                v_positionMC = position3DHigh + position3DLow;\n                v_positionEC = (czm_modelViewRelativeToEye * p).xyz;\n                v_st = st;\n                v_normalEC = czm_normal * normal;  \n                vec4 color = texture(image_0, v_st); \n                vec3 upDir = normalize(v_positionMC.xyz); \n                vec3 disPos = upDir * color.r *  maxHeight_1 * heightRatio_2;\n                p += vec4(disPos, 0.0);\n                gl_Position = czm_modelViewProjectionRelativeToEye * p;\n            }',
    '_area',
    '\nvec2 rotate2D(vec2 _st, float _angle){\n    //_st-= 0.5;\n    _st =  mat2(cos(_angle),-sin(_angle),\n                sin(_angle),cos(_angle)) * _st;\n    _st += 0.5;\n    return _st;\n} \nczm_material czm_getMaterial(czm_materialInput materialInput){\n  czm_material material = czm_getDefaultMaterial(materialInput);\n  float time = czm_frameNumber * speed / 1000.0; \n  vec2 st = materialInput.st - 0.5;  \n  st/=fract(time);\n  st=rotate2D(st,0.);\n  vec4 colorImage = texture(image, st);\n  material.alpha =colorImage.a ;\n  material.diffuse = color.rgb; \n  return material;\n} ',
    'sunShiny',
    'showSkirts',
    'Module',
    '_depthStencilTexture',
    'minScaleWidth',
    'image',
    '浮动点',
    '_points',
    'defaultColor',
    '_drivenTime',
    '_cartesian3ArrayTemp',
    'route',
    'getUniform',
    '_sectorVA',
    'LabelStyle',
    'image3',
    '_primitiveBias',
    'TWO_PI',
    '_mouseTip',
    'WavePetalsFS',
    'BACK',
    'rgba(0,255,255,0.3)',
    'WGS84',
    'length',
    ';  \n                     }',
    '_graphicLayer',
    '_drivenPositionsTemp',
    'clientHeight',
    '_lineOpts',
    'clampAnimations',
    '_rightClickEvent',
    'gif',
    '着色器特效霓虹灯',
    '#FFFFFF',
    '点击鼠标左键新增点，点击鼠标右键结束',
    'drawCancel',
    'viewMatrix',
    '_orientation',
    'px;width:',
    'minPointCount',
    ',\n                z: ',
    'showGlowRing',
    '_getPoint',
    'href',
    'bottomWidth',
    'enabled',
    '_modelGeometryNeedsUpdate',
    'lightColor',
    'addContent',
    '_drivenPositions',
    'createCamera',
    'WallCoolMaterialProperty',
    'createOutlineGeometry',
    'SingleTile',
    'editStart',
    'getTranslation',
    'Cartesian3',
    'anchorImage',
    'UniformState',
    'uniform vec4 u_color;\nczm_material czm_getMaterial(czm_materialInput materialInput){\n    czm_material material = czm_getDefaultMaterial(materialInput);\n    vec2 st = materialInput.st;\n    float powerRatio = 1./(fract(czm_frameNumber / 30.0) +  1.) ;\n    float alpha = pow(1. - st.t,powerRatio);\n    vec4 color = vec4(u_color.rgb, alpha*u_color.a);\n    material.diffuse = color.rgb;\n    material.alpha = color.a;\n    return material;\n}',
    '_geometryType',
    '&xt3d_token=',
    'video',
    '_leftDownPosition',
    'out_FragColor',
    'divIndicator-drag',
    'rgb(1 255 244)',
    'screenSpaceEventHandler',
    'material',
    'PolylineTrail',
    '_minPointCount',
    'MillitaryGraphic',
    'circleColorful',
    'CustomShaderMode',
    'minimumPixelSize',
    'withAlpha',
    'activeVideoListener',
    'commandList',
    'multiplyTransformation',
    'horizontal',
    'updateFrustum',
    'dispose',
    '_createLayer',
    'highDynamicRange',
    '着色器特效雷达扫描',
    'SCENE3D',
    'DOUBLE',
    'uniform sampler2D u_map;\nuniform vec4 u_bgColor;\n    \nin vec2 v_textureCoordinates;\nvoid main(){\n    vec2 uv = v_textureCoordinates;\n    vec4 bgColor = u_bgColor;\n    vec4 color = texture( u_map , uv );\n    if( color.a > 0. ){\n        out_FragColor = bgColor;\n    }\n}',
    'getFrustumPositions',
    '_style',
    '_unRegisterEvents',
    '着色器特效发光十字架',
    'fixedFrameToEastNorthUpTransform',
    'substr',
    'far',
    'lineHeight',
    '_fixPointCount',
    'silhouetteSize',
    '总距离：',
    '_cameraEvent',
    'angleBetween',
    'styleChange',
    'PopupType',
    '_matrixState',
    'picth',
    'lineLength',
    '_container',
    'reflectMatrix',
    'tool',
    'toDegrees',
    'cascadesEnabled',
    'center',
    '_angle',
    'video-js',
    '\nuniform vec4 color;\nuniform float speed;\nuniform float percent;\nuniform float gradient;\n\nczm_material czm_getMaterial(czm_materialInput materialInput){\n    czm_material material = czm_getDefaultMaterial(materialInput);\n    vec2 st = materialInput.st;\n    float t =fract(czm_frameNumber * speed / 1000.0);\n    t *= (1.0 + percent);\n    float alpha = smoothstep(t- percent, t, st.s) * step(-t, -st.s);\n    alpha += gradient;\n    material.diffuse = color.rgb*2.5;\n    material.alpha = alpha;\n    return material;\n}',
    'BYTES_PER_ELEMENT',
    'loop',
    'e-resize',
    '_graphicEdit',
    'setData',
    '_element',
    '_frustumGeometry',
    'outlineWidth',
    'topOindices',
    'fadeFactor',
    'normalize',
    'headingPitchRange',
    'angularity',
    'PolyGradientType',
    'multiplyByMatrix3',
    'innerRadius',
    'getPrototypeOf',
    'reflectivity',
    '_pointOpts',
    'fragmentShaderText',
    'hasSelected',
    'floatHeight',
    '_convertData',
    '_createBottomRing',
    '_getStyle',
    '可拖拽 div',
    'viewer',
    '等高线',
    'minRadius',
    'extrudedHeight',
    'light',
    'rgba(247, 226, 6, 0.46)',
    'Framebuffer',
    'cursor',
    '_log2FarDepthFromNearPlusOne',
    '\n#extension GL_OES_standard_derivatives : enable \nprecision highp float; \nuniform float time;\t\nuniform vec2  click;\nuniform vec2  resolution;\nuniform vec2 position;\nuniform vec2 mouse;\nuniform float t;\nuniform float a;\nin vec2 v_st;\n#define PI 3.14\nmat2 rotate3d(float angle){\n\treturn mat2(cos(angle), -sin(angle), sin(angle), cos(angle));\n}\nvoid main( void ) { \n       vec2 resolution = czm_viewport.zw;\n       float time = czm_frameNumber / 1000.0;\n\t   vec2 p = (gl_FragCoord.xy * 2.0 - resolution) / min(resolution.x, resolution.y);\n       vec2  st= v_st * 2.0 - 1.0;\n\t   p = rotate3d((time * 10.0) * PI) * st * 2.0 ;\n       float t = 0.025 / abs(abs(sin(time)) - length(p));\n\n       vec4 col= vec4(vec3(t) * vec3(p.x,p.y,80.0), 1.0);\n       float alpha=0.5;\n       alpha= col.z / 1000.0; \n       out_FragColor =  vec4(col.xyz,alpha);\n    }\n',
    '拖拽缩放高度',
    'shaderType',
    '_getOpts',
    'exports',
    'downloadurl',
    '_removeAnimateEntity',
    '_sectorLineCommand',
    'fillRect',
    '_tempPoint4',
    'background-image: linear-gradient( 135deg, transparent 30px, ',
    '_setVisible',
    'time',
    '水柱粒子',
    'closePath',
    '_createPostcess',
    '_createRawCommand',
    'LEFT_DOUBLE_CLICK',
    '/FF4AD0069FD64074BE8B10BCE8B4C27C.png',
    '\nuniform vec4 color;\n        uniform float speed;\n        #define PI 3.14159265359\n\n        vec2 rotate2D (vec2 _st, float _angle) {\n        _st =  mat2(cos(_angle),-sin(_angle),  sin(_angle),cos(_angle)) * _st;\n        return _st;\n        }\n        czm_material czm_getMaterial(czm_materialInput materialInput){\n        czm_material material = czm_getDefaultMaterial(materialInput);\n        vec2 st = materialInput.st * 2.0 - 1.0;\n        st *= 1.6;\n        float time = czm_frameNumber * speed / 1000.0;\n        float r = length(st);\n        float w = .3;\n        st = rotate2D(st,(r*PI*6.-time*2.));\n        float a = smoothstep(-w,.2,st.x) * smoothstep(w,.2,st.x);\n        float b = abs(1./(sin(pow(r,2.)*2.-time*1.3)*6.))*.4;\n        material.alpha = a * b ;\n        material.diffuse = color.rgb * a * b  * 3.0;\n        return material;\n }',
    '.title-line',
    '_cullFaceType',
    '_pointLightRadius',
    'beamRadar',
    '视频地址是必须的！',
    'rgba(255,255,255,0.5)',
    '\n#ifdef GL_ES\nprecision mediump float;\n#endif\n\n#extension GL_OES_standard_derivatives : enable\n\n  float time;\n  vec2 mouse;\n  vec2 resolution;\n  in vec2 v_st;\nconst int num_balls = 280; // change this, but be careful \nconst float coordinate_scale = 1000.;\nconst float structure_size = 1.20; // Size from 0 to 1, max size = 1\nconst float glow_decay = 1.40; \nconst float trail_len = 151.0;\n#define pi 3.14159265358\nconst float speed = 0.189;\nconst float rot_speed = speed*10.;\nconst float starting_pt = -506.0; // This is a good number\n\nvec4 draw_ball(float i, float j, float size) {\n\tfloat balls = float(num_balls);\n\tfloat dt = starting_pt + time * speed;\n\t// Map coordinates to 0-1\n\tvec2 coord =  v_st;//gl_FragCoord.xy/resolution.xy;\n\t//map coordinates to -coord_scale -> +coord_scale\n\tcoord = coord *coordinate_scale-coordinate_scale/2.;\n\tcoord -= vec2(coord.x/2.,coord.y/2.);\n\t\n\t//Controls motion of balls\n\tfloat spacing = (2.*pi)/balls;\n\t\n\tfloat x =  (sin(dt*i*spacing)*100. - cos(dt*j*spacing)) * structure_size;\n\tfloat y =  (cos(dt*j*spacing)*100. - sin(dt*i*spacing)) * structure_size;\n\ty *= ((j - dt)/-dt) + sin(i*spacing + dt*sin(dt/100.));\n\tx *= ((i - dt)/dt)  + sin(i*spacing - dt*cos(dt/100.));\n\t//Correct aspect ratio\n\tcoord.x *= resolution.x/resolution.y; \n\tvec2 pos = vec2(x,y);\n\tmat2 rot = mat2(cos(dt*rot_speed), -sin(dt*rot_speed), sin(dt*rot_speed), cos(dt*rot_speed));\n\tpos *= rot;\n\tfloat dist = length(coord - pos);\n\t\n\t//Controls how quickly brightness falls off\n\tfloat intensity = pow(size/dist, glow_decay);\n\t\n\tvec4 color = vec4(vec3(1.0) * abs(sin(vec3(time*0.25,time/2.,time/3.))), 1.0);\n\treturn color * intensity;\n}\n\n\nvoid main( void ) {\n\t\n    mouse=vec2(1.0,1.0);\n    time = czm_frameNumber /100.0;\n    resolution = czm_viewport.zw; \n\n\t vec4 col = vec4(0.0);\n\t\n\tfor (int i = 0; i < num_balls; i++) {\n\t\tvec2 pt = vec2(float(i),float(i));\n\t\tcol += draw_ball(float(i),float(i), 1.2-distance(pt,vec2(0.))/coordinate_scale);\n\t}\n  float alpha=(col.x+ col.y + col.z) / 5.0; \n  if(alpha < 0.2){ \n      discard; \n  } \n  col*=2.0;\n  col.a=alpha;\n\tout_FragColor =col;// vec4( col,alpha );\n\t//out_FragColor = col;\n}',
    '_createContentDom',
    'textures',
    'ShadowMap',
    '__removeHook',
    'scaleByDistance',
    'cesiumWidget',
    '_handler',
    '_createPrimitive',
    'http://www.xt3d.online:9999/resource/?assetId=',
    '_endFovH',
    'model',
    'autoPoistion',
    '_addPrimitive',
    'graphic',
    '_renderRequested',
    '_getPolylinePositions',
    'xHalfAngle',
    'reflectorViewMatrix',
    'setOption',
    'fire',
    'mp4',
    '拖拽缩放宽度',
    'svg获取失败！',
    'pixelRange',
    '_innerCylinder',
    'rgba(255, 0, 0, 0.5)',
    '_drawGraphicOpts',
    'preMultiplyAlpha',
    'textarea',
    'addSample',
    '_createBottomCircle',
    '_container_fixed',
    '竖立文本标签',
    '_updateWallGeometry',
    'WallGraphic',
    'TRANSPARENT',
    '_createVertex',
    'isSelected',
    '_createImage_2',
    'addSeconds',
    'enableCursorStyle',
    'catch',
    'CylinderRadar',
    'removeById',
    'Math',
    '_rect',
    'CircleWaveType',
    '.b-t',
    'circleScan_3',
    'WallScrollType',
    'classificationType',
    'shaderEffetConstellationChain',
    '坐标值：',
    '_removeHeatContainer',
    'TextureMinificationFilter',
    'documentElement',
    '_covertMinimumHeights',
    '_inputC',
    '_outImage',
    'clustering',
    'fetchImage',
    'RectangleGeometry',
    'attachEvent',
    'scaleAll',
    'indices',
    '_computePlane',
    'wallCool',
    'start',
    '_createCircleImage',
    'TextureMagnificationFilter',
    '着色器特效病毒',
    '拖拽修改位置',
    '_clearSelectedObject',
    'textBaseline',
    'PetalFS',
    'selectedColor',
    'startScale',
    'bottomHeight',
    '_imageCache',
    'Microsoft YaHei',
    '20UaEEGw',
    '_dragDom_mousedown',
    '波纹扩散',
    '_createPolygonEntity',
    'RenderState',
    'transformWindowToDrawingBuffer',
    'out_FragColor =  vec4(out_FragColor.rgb + u_color * 2., out_FragColor.r * .3)',
    '_cameraQ',
    'stop',
    'hanging',
    '_domeBackCommand',
    'tilesetOpts',
    'topRadius',
    'GlowStarFS',
    'Event',
    '_isMeasuring',
    '，最少需要',
    '_typeName',
    'outlineStyle',
    'data:image/octet-stream;base64,iVBORw0KGgoAAAANSUhEUgAAAGQAAAABCAYAAAAo2wu9AAAAOklEQVQoU2P8//+/McMoIDsEfv/+zcLKysrOwMDAAcVg9t+/fzmYmZkxxP/9+8fBxMSEIc7AwAAWAwD/kwzHVTmPqQAAAABJRU5ErkJggg==',
    '分隔符：',
    '水滴球',
    'wallImageFlow',
    'once',
    'firstChild',
    'axisGraphic',
    'fillStyle',
    'assetsPath',
    'BLACK',
    '_valueType',
    'EditAnchor',
    'ALT',
    'feature',
    'setScale',
    '_frustumOutlineGeometry',
    '_handleActivate',
    'PixelDatatype',
    '_merge',
    'diffuseWall',
    'ReflectWater',
    'result',
    'rgba(0,183,239, 0.5)',
    'frequency',
    'flv',
    'restore',
    '_billboard',
    'save',
    'xt3d-divgraphic-hide',
    '_oneOverLog2FarDepthFromNearPlusOne',
    'apply',
    'pickMousePosition',
    'updateAllPositionByMatrix',
    'graphicOpts',
    '_fov',
    '_polygonOpts',
    'clusterStyleType',
    'initAnchorPositions',
    '_options',
    'CircleScanMaterialProperty_1',
    'translation',
    ';font-size:',
    'modelGraphic',
    'setColumn',
    'image/svg+xml',
    'width',
    'pickTooltip',
    '_attributeLocations',
    'WallImageFlowMaterialProperty',
    '_trs',
    '扫描材质_2',
    'rgba(255, 255, 0, 0.4)',
    'rgba(0, 255, 255, 0.44)',
    '点聚合图层',
    'circleLightRing',
    '_darkness',
    'rightWC',
    'ColorfulStellarChainFS',
    'WHITE',
    'setOpts',
    '_createImage_4',
    'LightGraphic',
    '_outFragColor',
    '   \nprecision mediump float;\nin vec2 v_st;\n\nvoid main(void){\nfloat time = czm_frameNumber /60.0;\nvec2 resolution=czm_viewport.zw;\n\nvec2  st=   v_st /2.0;// * 2.0 -1.0 ;\n\nvec3 destColor = vec3(0.52, 0.2, 0.1);\nvec2 p =  v_st * 2.0 - 1.0 ;//(st.xy * 2.0 - resolution) / min(resolution.x, resolution.y); \t\nfloat a = atan(p.y / p.x) * 2.0; // Instead of * 2.0, try * 26 or * 128 and LOWER\nfloat l = 0.05 / abs(length(p) - 0.8 + sin(a + time * 4.5) * 0.1);\ndestColor *= 1.9+ sin(a + time * 00.13) * 0.03;\n\nvec3 destColor2 = vec3(0.0, 0.2, 0.9);\nvec2 p2 = (st.xy * 3.0 - resolution) / min(resolution.x, resolution.y); \nfloat a2 = atan(p.y / p.x) * 3.0;\nfloat l2 = 0.05 / abs(length(p) + 0.1 - (tan(time/2.)+0.5) + sin(a + time * 13.5) * (0.1 * l));\ndestColor2 *= ( 0.5 + sin(a + time * 00.03) * 0.03 ) * 4.0;\n\nvec3 destColor3 = vec3(0.2, 0.9, 0.35);\nvec2 p3 = (st.xy * 2.0 - resolution) / min(resolution.x, resolution.y); \nfloat a3 = atan(p.y / p.x) * 10.0;\nfloat l3 = 0.05 / abs(length(p) - 0.4 + sin(a + time * 23.5) * (0.1 * l2));\ndestColor3 *= 0.5 + sin(a + time * 10.23) * 0.03;\n\nvec3 color=l*destColor   + l3*destColor3;\n// vec3 color=l*destColor + l2*destColor2 + l3*destColor3;\nfloat alpha=(color.x+color.y+color.z)/10.0;\nif(alpha<0.1){\nalpha=0.0;\n} \n// out_FragColor = vec4(l*destColor + l2*destColor2 + l3*destColor3, alpha);\nout_FragColor = vec4(l*destColor +   l3*destColor3, alpha);\n}',
    'xt3d外观',
    '_appearanceProxy',
    'HeatGraphic',
    'command',
    '\n #ifdef GL_OES_standard_derivatives\n    #extension GL_OES_standard_derivatives : enable\n#endif\n\nuniform bool u_showIntersection;\nuniform bool u_showThroughEllipsoid;\n\nuniform float u_radius;\nuniform float u_xHalfAngle;\nuniform float u_yHalfAngle;\nuniform float u_normalDirection;\nuniform float u_type;\nuniform vec4 u_color;\nuniform vec4 czm_pickColor;\n\n\nin vec3 v_position;\nin vec3 v_positionWC;\nin vec3 v_positionEC;\nin vec3 v_normalEC;\nvec4 getColor(float sensorRadius, vec3 pointEC)\n{\n    czm_materialInput materialInput;\n\n    vec3 pointMC = (czm_inverseModelView * vec4(pointEC, 1.0)).xyz;\n    materialInput.st = sensor2dTextureCoordinates(sensorRadius, pointMC);\n    materialInput.str = pointMC / sensorRadius;\n\n    vec3 positionToEyeEC = -v_positionEC;\n    materialInput.positionToEyeEC = positionToEyeEC;\n\n    vec3 normalEC = normalize(v_normalEC);\n    materialInput.normalEC = u_normalDirection * normalEC;\n\n    czm_material material = czm_getMaterial(materialInput);\n    // czm_lightDirectionEC在cesium1.66开始加入的\n    return mix(czm_phong(normalize(positionToEyeEC), material, czm_lightDirectionEC), vec4(material.diffuse, material.alpha), 0.4);\n\n}\n\nbool isOnBoundary(float value, float epsilon)\n{\n    float width = getIntersectionWidth();\n    float tolerance = width * epsilon;\n\n#ifdef GL_OES_standard_derivatives\n    float delta = max(abs(dFdx(value)), abs(dFdy(value)));\n    float pixels = width * delta;\n    float temp = abs(value);\n    // There are a couple things going on here.\n    // First we test the value at the current fragment to see if it is within the tolerance.\n    // We also want to check if the value of an adjacent pixel is within the tolerance,\n    // but we don\'t want to admit points self are obviously not on the surface.\n    // For example, if we are looking for "value" to be close to 0, but value is 1 and the adjacent value is 2,\n    // then the delta would be 1 and "temp - delta" would be "1 - 1" which is zero even though neither of\n    // the points is close to zero.\n    return temp < tolerance && temp < pixels || (delta < 10.0 * tolerance && temp - delta < tolerance && temp < pixels);\n#else\n    return abs(value) < tolerance;\n#endif\n}\n\nvec4 shade(bool isOnBoundary)\n{\n    if (u_showIntersection && isOnBoundary)\n    {\n        return getIntersectionColor();\n    }\n    if(u_type == 1.0){\n        return getLineColor();\n    }\n    return getColor(u_radius, v_positionEC);\n}\n\nfloat ellipsoidSurfaceFunction(vec3 point)\n{\n    vec3 scaled = czm_ellipsoidInverseRadii * point;\n    return dot(scaled, scaled) - 1.0;\n}\n\nvoid main()\n{\n    vec3 sensorVertexWC = czm_model[3].xyz;      // (0.0, 0.0, 0.0) in model coordinates\n    vec3 sensorVertexEC = czm_modelView[3].xyz;  // (0.0, 0.0, 0.0) in model coordinates\n\n    //vec3 pixDir = normalize(v_position);\n    float positionX = v_position.x;\n    float positionY = v_position.y;\n    float positionZ = v_position.z;\n\n    vec3 zDir = vec3(0.0, 0.0, 1.0);\n    vec3 lineX = vec3(positionX, 0 ,positionZ);\n    vec3 lineY = vec3(0, positionY, positionZ);\n    float resX = dot(normalize(lineX), zDir);\n    if(resX < cos(u_xHalfAngle)-0.00001){\n        discard;\n    }\n    float resY = dot(normalize(lineY), zDir);\n    if(resY < cos(u_yHalfAngle)-0.00001){\n        discard;\n    }\n\n\n    float ellipsoidValue = ellipsoidSurfaceFunction(v_positionWC);\n\n    // Occluded by the ellipsoid?\n\tif (!u_showThroughEllipsoid)\n\t{\n\t    // Discard if in the ellipsoid\n\t    // PERFORMANCE_IDEA: A coarse check for ellipsoid intersection could be done on the CPU first.\n\t    if (ellipsoidValue < 0.0)\n\t    {\n            discard;\n\t    }\n\n\t    // Discard if in the sensor\'s shadow\n\t    if (inSensorShadow(sensorVertexWC, v_positionWC))\n\t    {\n\t        discard;\n\t    }\n    }\n\n    // Notes: Each surface functions should have an associated tolerance based on the floating point error.\n    bool isOnEllipsoid = isOnBoundary(ellipsoidValue, czm_epsilon3);\n    //isOnEllipsoid = false;\n    //if((resX >= 0.8 && resX <= 0.81)||(resY >= 0.8 && resY <= 0.81)){\n    /*if(false){\n        out_FragColor = vec4(1.0,0.0,0.0,1.0);\n    }else{\n        out_FragColor = shade(isOnEllipsoid);\n    }\n*/\nif(u_type==0.){\n    out_FragColor =u_color;\n}else{\n    out_FragColor =u_lineColor;\n} \n}',
    'diffuse',
    '_lightPositionEC',
    '_setMove',
    '_getDefaultShader',
    'hierarchy',
    'clockRange',
    'fromUniformScale',
    'multiplyByScalar',
    '拖拽整体平移',
    'videoProject',
    '_setMaterial',
    'defineProperties',
    '矩形文本标签',
    'updateMoveAllPosition',
    'Rectangle',
    '_getSiteTimes',
    'rgba(255,255,0,1)',
    'Test',
    '_matrix4',
    '_endFovV',
    '_swallowTailFactor',
    '漫游线',
    '着色器特效多彩点',
    'userSelect',
    '反射水面',
    'rgba(255, 60, 0, 0.43)',
    '纯文本点',
    '_getImage',
    'resourceToken',
    'shaderEffetNeonPoint',
    'CircleSpiralMaterialProperty',
    'getElementsByTagName',
    ';background: linear-gradient(',
    'fractionDigits',
    '\n            // "Type 2 Supernova" by Duke\n            // https://www.shadertoy.com/view/lsyXDK\n            // VR ready version of this shader is here https://www.shadertoy.com/view/XdKSDV (thanks RavenWorks for help)\n            //-------------------------------------------------------------------------------------\n            // Based on "Supernova remnant" (https://www.shadertoy.com/view/MdKXzc)\n            // "Volumetric explosion" (https://www.shadertoy.com/view/lsySzd)\n            // otaviogood\'s "Alien Beacon" (https://www.shadertoy.com/view/ld2SzK)\n            // Shane\'s "Cheap Cloud Flythrough" (https://www.shadertoy.com/view/Xsc3R4) shaders\n            // Some ideas came from other shaders from this wonderful site\n            // Press 1-2-3 to zoom in and zoom out.\n            // License: Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License\n            //-------------------------------------------------------------------------------------\n\n            #define DITHERING\n            #define BACKGROUND\n\n            //#define TONEMAPPING\n\n            //-------------------\n            #define pi 3.14159265\n            #define R(p, a) p=cos(a)*p+sin(a)*vec2(p.y, -p.x)\n\n            // iq\'s noise\n            float noise( in vec3 x )\n            {\n                vec3 p = floor(x);\n                vec3 f = fract(x);\n                f = f*f*(3.0-2.0*f);\n                vec2 uv = (p.xy+vec2(37.0,17.0)*p.z) + f.xy;\n                vec2 rg = texture( iChannel0, (uv+ 0.5)/256.0, 0.0 ).yx;\n                return 1. - 0.82*mix( rg.x, rg.y, f.z );\n            }\n\n            float fbm( vec3 p )\n            {\n            return noise(p*.06125)*.5 + noise(p*.125)*.25 + noise(p*.25)*.125 + noise(p*.4)*.2;\n            }\n\n            float Sphere( vec3 p, float r )\n            {\n                return length(p)-r;\n            }\n\n            //==============================================================\n            // otaviogood\'s noise from https://www.shadertoy.com/view/ld2SzK\n            //--------------------------------------------------------------\n            // This spiral noise works by successively adding and rotating sin waves while increasing frequency.\n            // It should work the same on all computers since it\'s not based on a hash function like some other noises.\n            // It can be much faster than other noise functions if you\'re ok with some repetition.\n            const float nudge = 4.;\t// size of perpendicular vector\n            float normalizer = 1.0 / sqrt(1.0 + nudge*nudge);\t// pythagorean theorem on that perpendicular to maintain scale\n            float SpiralNoiseC(vec3 p)\n            {\n                float n = 1.-mod(iTime * 0.1,-1.); // noise amount\n                float iter = 2.0;\n                for (int i = 0; i < 8; i++)\n                {\n                    // add sin and cos scaled inverse with the frequency\n                    n += -abs(sin(p.y*iter) + cos(p.x*iter)) / iter;\t// abs for a ridged look\n                    // rotate by adding perpendicular and scaling down\n                    p.xy += vec2(p.y, -p.x) * nudge;\n                    p.xy *= normalizer;\n                    // rotate on other axis\n                    p.xz += vec2(p.z, -p.x) * nudge;\n                    p.xz *= normalizer;\n                    // increase the frequency\n                    iter *= 1.733733;\n                }\n                return n;\n            }\n\n            float VolumetricExplosion(vec3 p)\n            {\n                float final = Sphere(p,4.);\n                final += noise(p*20.)*.4;\n                final += SpiralNoiseC(p.zxy*fbm(p*10.))*2.5; //1.25;\n\n                return final;\n            }\n\n            float map(vec3 p) \n            {\n                R(p.xz, iMouse.x*0.008*pi+iTime*0.1);\n\n                float VolExplosion = VolumetricExplosion(p/(1.+mod(iTime * 0.1,-1.)))*(1.+mod(iTime * 0.1,-1.)); // scale\n                \n                return abs(VolExplosion)+0.07;\n            }\n            //--------------------------------------------------------------\n\n            // assign color to the media\n            vec3 computeColor( float density, float radius )\n            {\n                // color based on density alone, gives impression of occlusion within\n                // the media\n                vec3 result = mix( vec3(1.0,0.9,0.8), vec3(0.4,0.15,0.1), density );\n                \n                // color added to the media\n                vec3 colCenter = 7.*vec3(0.8,1.0,1.0);\n                vec3 colEdge = 1.5*vec3(0.48,0.53,0.5);\n                result *= mix( colCenter, colEdge, min( (radius+.05)/.9, 1.15 ) );\n                \n                return result;\n            }\n\n            bool RaySphereIntersect(vec3 org, vec3 dir, out float near, out float far)\n            {\n                float b = dot(dir, org);\n                float c = dot(org, org) - 8.;\n                float delta = b*b - c;\n                if( delta < 0.0) \n                    return false;\n                float deltasqrt = sqrt(delta);\n                near = -b - deltasqrt;\n                far = -b + deltasqrt;\n                return far > 0.0;\n            }\n\n            // Applies the filmic curve from John Hable\'s presentation\n            // More details at : http://filmicgames.com/archives/75\n            vec3 ToneMapFilmicALU(vec3 _color)\n            {\n                _color = max(vec3(0), _color - vec3(0.004));\n                _color = (_color * (6.2*_color + vec3(0.5))) / (_color * (6.2 * _color + vec3(1.7)) + vec3(0.06));\n                return pow(_color, vec3(2.2));\n            }\n\n            void main()\n            {  \n                const float KEY_1 = 49.5/256.0;\n                const float KEY_2 = 50.5/256.0;\n                const float KEY_3 = 51.5/256.0;\n                float key = 0.0;\n                key += 0.7*texture(iChannel1, vec2(KEY_1,0.25)).x;\n                key += 0.7*texture(iChannel1, vec2(KEY_2,0.25)).x;\n                key += 0.7*texture(iChannel1, vec2(KEY_3,0.25)).x;\n\n                // vec2 uv = fragCoord/iResolution.xy;\n                vec2 uv = 2.*v_st.xy - vec2(1., 1.);\n                // ro: ray origin\n                // rd: direction of the ray\n                vec3 rd = normalize(vec3(uv.xy , 1.0));\n                vec3 ro = vec3(0., 0., -6.+key*1.6);\n                \n                // ld, td: local, total density \n                // w: weighting factor\n                float ld=0., td=0., w=0.;\n\n                // t: length of the ray\n                // d: distance function\n                float d=1., t=0.;\n                \n                const float h = 0.1;\n            \n                vec4 sum = vec4(0.0);\n            \n                float min_dist=0.0, max_dist=0.0;\n\n                if(RaySphereIntersect(ro, rd, min_dist, max_dist))\n                {\n                \n                t = min_dist*step(t,min_dist);\n            \n                // raymarch loop\n                for (int i=0; i<86; i++)\n                {\n                \n                    vec3 pos = ro + t*rd;\n            \n                    // Loop break conditions.\n                    if(td>0.9 || d<0.11*t || t>10. || sum.a > 0.99 || t>max_dist) break;\n                    \n                    // evaluate distance function\n                    float d = map(pos);\n                    \n                    // point light calculations\n                    vec3 ldst = vec3(0.0)-pos;\n                    float lDist = max(length(ldst), 0.001);\n\n                    // the color of light \n                    vec3 lightColor=vec3(1.0,0.5,0.25);\n                    \n                    sum.rgb+=(vec3(0.67,0.75,1.00)/(lDist*lDist*15.)/100.); // star itself\n                    sum.rgb+=(lightColor/exp(lDist*lDist*lDist*.08)/30.); // bloom\n                    \n                    if (d<h) \n                    {\n                        // compute local density \n                        ld = h - d;\n                        \n                        // compute weighting factor \n                        w = (1. - td) * ld;\n                \n                        // accumulate density\n                        td += w + 1./200.;\n                    \n                        vec4 col = vec4( computeColor(td,lDist), td );\n                        \n                        // emission\n                        sum += sum.a * vec4(sum.rgb, 0.0) * 0.2 / lDist;\t\n                        \n                        // uniform scale density\n                        col.a *= 0.2;\n                        // colour by alpha\n                        col.rgb *= col.a;\n                        // alpha blend in contribution\n                        sum = sum + col*(1.0 - sum.a);  \n                \n                    }\n                \n                    td += 1./70.;\n\n                    #ifdef DITHERING\n                    // idea from https://www.shadertoy.com/view/lsj3Dw\n                    vec2 uvd = uv;\n                    uvd.y*=120.;\n                    uvd.x*=280.;\n                    d=abs(d)*(.8+0.08*texture(iChannel2,vec2(uvd.y,-uvd.x+0.5*sin(4.*iTime+uvd.y*4.0))).r);\n                    #endif \n                    \n                    // trying to optimize step size\n                    t += max(d * 0.1 * max(min(length(ldst),length(ro)),0.1), 0.01);\n\n                }\n                \n                // simple scattering\n                sum *= 1. / exp( ld * 0.2 ) * 0.8;\n                    \n                sum = clamp( sum, 0.0, 1.0 );\n            \n                sum.xyz = sum.xyz*sum.xyz*(3.0-2.0*sum.xyz);\n                \n                }\n                \n                #ifdef BACKGROUND\n                // stars background\n                if (td<.8)\n                {\n                    vec3 stars = vec3(noise(rd*500.0)*0.5+0.5);\n                    vec3 starbg = vec3(0.0);\n                    starbg = mix(starbg, vec3(0.8,0.9,1.0), smoothstep(0.99, 1.0, stars)*clamp(dot(vec3(0.0),rd)+0.75,0.0,1.0));\n                    starbg = clamp(starbg, 0.0, 1.0);\n                    sum.xyz += starbg; \n                }\n                #endif\n            \n                #ifdef TONEMAPPING\n                fragColor = vec4(ToneMapFilmicALU(sum.xyz),1.0);\n                #else\n                fragColor = vec4(sum.xyz,1.0);\n                #endif\n            }\n\n        ',
    'uri',
    'rgba(255,255,255,0)',
    'abs',
    'createShadowMap',
    'repeat',
    '_reflectorViewMatrix',
    'isPointLight',
    'SingleTileImageryProvider',
    'duration',
    '\nvec2 rotate2D(vec2 _st, float _angle){\n    _st -= 0.5;\n    _st =  mat2(cos(_angle),-sin(_angle),\n                sin(_angle),cos(_angle)) * _st;\n    _st += 0.5;\n    return _st;\n} \nczm_material czm_getMaterial(czm_materialInput materialInput){\n  czm_material material = czm_getDefaultMaterial(materialInput);\n  vec2 st = materialInput.st  ;\n  float angle = czm_frameNumber * speed / 1000.0; \n  //angle=mod(angle,360.);如果要限定 旋转最大最小值 可通过entity实现\n  st = rotate2D(st,(angle));\n  vec4 colorImage = texture(image, fract(st));\n  material.alpha =colorImage.a ;\n  material.diffuse = color.rgb * 1.5;\n  return material;\n}',
    'asynchronous',
    'pointerEvents',
    'moveVertexEntity',
    'preUpdate',
    'margin:5px 1px;display:flex',
    'rectPyramid',
    '</div>\n            <div class=\'title-line\'></div>\n            <div class="label-content"> \n                ',
    'shaderEffetMagicBall',
    '\nconst float PI = 3.14159265359;\nconst float TWO_PI = 6.28318530718;\nconst int N = 3;\t\t\t\t// triangle polygons please\nconst float r0 = 0.01;\t\t\t// size of centre circle\nconst float r_blue = 0.025;\t\t// size of blue radar blips\nconst float r_red = 0.015;\t\t// size of red radar blips\nconst float edge = 0.95;\t\t// overall size\nconst float offset = 0.05;\n\nfloat plot(const vec2 st, const float pct, const float width)\n{\n        return smoothstep(pct - width, pct, st.y) -\n            smoothstep(pct, pct + width, st.y);\n    }\n\nfloat drawPolygon(const vec2 polygonCenter, const int N, const float radius, vec2 pos)\n{\n    pos = pos - polygonCenter;\n    float d = 0.0;\n    float a = atan(pos.x, pos.y);\n    float r = TWO_PI / float(N);\n    d = cos(floor(0.5 + a / r)*r - a)*length(pos);\n    return (1.0 - smoothstep(radius, radius + radius/10.0, d));\n}\n\nfloat gradations(const float a, const float gradNum, const float outRad, const float tickLen, const float tickWidth, const float r, const float move)\n{\n    float f = step(0.0, cos((a + move)*gradNum) - tickWidth)*tickLen + (outRad - tickLen);\n    return 1.0 - step(f, r) * 1.0 - step(r, outRad - tickLen);\n}\n\nczm_material czm_getMaterial(czm_materialInput materialInput){\n    czm_material material = czm_getDefaultMaterial(materialInput);\n    vec2 v_st = materialInput.st;\n    float time = czm_frameNumber / 60.0;\n    float iTime = u_speed * time;\n    vec2 pos = 2.0*v_st.xy - vec2(1., 1.) ; // center what being drawn\n    pos /=2.05;\n    vec4 grndSpd = vec4(0.0, iTime/5.0, 0.0, 0.0);\n    vec4 mapcol = vec4 (1.0, 1.0, 1.0, 1.0);\n    \n    vec3 color = vec3(0.0, 0.0, 0.0);\n\n    float r = length(pos) * 2.0;\n    float a = atan(pos.y, pos.x); // angle of pixel\n    float an = PI - mod(iTime / 1.0, TWO_PI); // angle of radar sweep\n        float blipSpd = 3.0; // Blip / Trace speed\n    vec2 translate1 = vec2(cos(iTime/ blipSpd), sin(iTime / blipSpd));\n    vec2 translate2 = vec2(sin(iTime / blipSpd), cos(iTime / blipSpd));\n    vec2 left1 = translate1 * 0.35;\n    vec2 right1 = -translate1 * 0.30;\n    vec2 left2 = translate2 * 0.15;\n    vec2 right2 = -translate2 * 0.25;\n        \n    // Radar Sweep\n        float sn = step(PI/2.0, an) * step(-PI/2.0, (a + an)) * step(r, edge) * (1.0 - 0.55 * (a + (TWO_PI) - an));\n    float sw = step(an, a) * step(r, edge);\n    float s_blade = sw * (1.0 - (a - an) * 20.0);\n    float s = sw * (1.0 - 0.55 * (a - an));\n    s = max(sn,s);\n    float se = step(r, edge - 0.05);\n    \n    // Center point\n    float s1 = smoothstep(edge - 0.00, edge + 0.01, r)* smoothstep(edge + 0.02, edge + 0.01, r);   \n    \n    // Circular concentric rings\n    float s0 = 1.0 - smoothstep(r0 / 2.0, r0, length(pos));\n        float smb = (1.0 - smoothstep(0.2, 0.2 + 0.01, length(pos))) * (1.0 - smoothstep(0.2 +0.01, 0.2, length(pos)));\n        float smr = (1.0 - smoothstep(0.3, 0.3 + 0.01, length(pos))) * (1.0 - smoothstep(0.3 +0.01, 0.3, length(pos)));\n        \n    // Circular concentric gradations\n    float gradNum = 120.0;\n    float tickWidth = 0.9;\n    const float tickLen = 0.04;\n    float outRad = edge;\n    float move = 0.0;\n    float sm = 0.75*gradations(a, gradNum, outRad, tickLen, tickWidth, r, move);   \n    \n    gradNum = 36.0;\n    tickWidth = 0.95;\n    outRad = 0.6;\n    move = sin(iTime/10.0);\n    smr += 0.5*gradations(a, gradNum, outRad, tickLen, tickWidth, r, move);\n\n    outRad = 0.4;\n    move = cos(iTime/10.0);\n    smb += 0.5*gradations(a, gradNum, outRad, tickLen, tickWidth, r, move);\n\n    // Radial spoke gradations \n    float sr = plot(pos, pos.x, 0.003) * step(r, edge - 0.06);\n    sr += plot(vec2(0.0, 0.0), pos.x, 0.002) * step(r, edge - 0.06);\n    sr += plot(vec2(0.0, 0.0), pos.y, 0.003) * step(r, edge - 0.06);\n    sr += plot(-pos, pos.x, 0.003) * step(r, edge - 0.06);\n        sr *= 0.75;\n\n    // Blue circular radar blip traces\n    vec2 st_trace1 = left2;\n    float s_trace1 = s * (1.0 - smoothstep(r_blue / 10.0, r_blue, length(pos - st_trace1)));\n    s_trace1 += s * (1.0 - smoothstep(r_blue / 10.0, r_blue, length(pos - st_trace1 + vec2(+offset, +offset))));\n    s_trace1 += s * (1.0 - smoothstep(r_blue / 10.0, r_blue, length(pos - st_trace1 + vec2(+2.0 *offset, +2.0 *offset))));\n\n    vec2 st_trace2 = right1;\n    float s_trace2 = s * (1.0 - smoothstep(r_blue / 10.0, r_blue, length(pos - st_trace2)));\n\n    // Red Trianglular radar flight blip trace \n    vec2 st_trace3 = left1;\n    float st1 = s * (drawPolygon(st_trace3, N, r_red , pos));\n    st1 += s * (drawPolygon(st_trace3 + vec2(-offset, -offset), N, r_red, pos));\n    st1 += s * (drawPolygon(st_trace3 + vec2(+offset, -offset), N, r_red, pos));\n\n    vec2 st_trace4 = right2;\n    float st2 = s * (drawPolygon(st_trace4, N, r_red, pos));  \n        \n    // Lets add all the bits together and send them to screen\n    float s_grn = max(s * mapcol.y, s_blade);\n    s_grn = max(s_grn, (s0 +  sr + sm));\n    s_grn += s1 / 1.5  + smb + smr;\n\n    float s_red = st1*2.0 + st2*2.0 + smr;\n        \n    float s_blue = max(s_trace1 + s_trace2, s_blade) + smb;\n\n    if (s_trace1 > 0.0 || s_trace2 > 0.0) { s_blue = max(s, s_blue); s_grn = max(s_grn, s_blue); }\n\n    color += vec3(s_red , s_grn, s_blue);   \n        \n    vec4 texColor = mapcol * s;\n    \n    // Output to screen   \n    vec4 fragColor = vec4(color, .8);//Set the screen pixel to that color\n    material.alpha = fragColor.r + fragColor.g + fragColor.b;\n    material.diffuse = fragColor.rgb * u_color.rgb  + color.rgb;\n    return material;\n}',
    'refractionWater',
    'onTick',
    'lineOpts',
    'uniformState',
    'CornerType',
    '规则墙体',
    'flat',
    'source',
    'hiddenColor',
    'target',
    'RadarScanFS',
    'sampler',
    'rectangleLabel',
    'baseHeight',
    'fromDegrees',
    'probeRadar',
    'rgba(227,108,9,0.5)',
    'innerConeCos',
    'FireworksFS',
    '_drawCommands',
    '_handleMouseMove',
    '_overObject',
    'separator',
    '_subSegmentV',
    'visibilitychange',
    '图片旋转',
    'rgba(28,25,125,0.99)',
    '_innerFovRadiusPairs',
    'flipY',
    'doubleArrow',
    'childNodes',
    '_remove',
    'eyeOffset',
    'off',
    'then',
    '_inverseViewProjectionDirty',
    'LESS',
    '_initProperty',
    'screenSpaceCameraController',
    '_lineInfo',
    'setLineDash',
    'CoolBallFS',
    'typeName',
    'scene',
    'heightReference',
    '_domeVA',
    '_videoExt',
    'PolyGradientMaterialProperty',
    'deactivate',
    'appendFromUrl',
    'projectionMatrix',
    'removeEventListener',
    'createObjectURL',
    '着色器特效发光盒子',
    '简单图标点',
    'worldToWindowCoordinates',
    'GeoJsonDataSource',
    'u_map',
    'pitch',
    'fromRotationZ',
    '_hiddenProperty',
    'filter',
    '_polygonEntity',
    '_graphics',
    '着色器特效光锥',
    'DoubleLoopFS',
    '_computeCommandList',
    'update',
    'createGuid',
    '_radiusCatch',
    'drawImage',
    'GlowBoxFS',
    'outlineColor',
    'waterColumn',
    'startFovV',
    'writeText',
    '\n    void main()\n    {\n       out_FragColor = vec4(1.0,.0,.0,.8);\n    }',
    'waveStrength',
    '_postRender',
    'translucencyMode',
    'CircleScanType_4',
    'PolyElevationContourType',
    'fontSize',
    'values',
    '_totalTime',
    'simpleWall',
    'fromCamera',
    'shaderEffect',
    'white',
    '_emitterModelMatrix',
    '_createTextureWebGL',
    'rgba(255, 255, 255, 0.0)',
    '_filter',
    'tabindex',
    'CircleLightRingType',
    'deg)',
    '_createLeftCrossSectionCommand',
    'raiseEvent',
    '10004upwMri',
    '_getTailPoints',
    'layerOpts',
    '.xt3d-divgraphic{position:absolute;left:0;bottom:0}.xt3d-divgraphic-show{opacity:1;transition:opacity 2s}.xt3d-divgraphic-hide{opacity:0;transition:opacity 2s}.xt3d-fillet-label-content{text-align:center;padding:5px 30px;margin:0;border-radius:5px;max-width:180px;max-height:134px;text-wrap:nowrap}.xt3d-fillet-label-line{content:"";position:absolute;bottom:-60px;left:calc(50% - 3px);display:block;width:3px;height:60px}.xt3d-divgraphic-animation,.xt3d-divgraphic-animation:after,.xt3d-divgraphic-animation:before,.xt3d-divgraphic-animation p,.xt3d-divgraphic-animation p:after,.xt3d-divgraphic-animation p:before{margin:0;padding:0;box-sizing:border-box}.xt3d-divgraphic-animation{width:10px;height:10px;border-radius:50%;cursor:pointer;color:#0ff;background:currentColor;z-index:3;left:50%;top:50%;transform:translate(-50%,-50%);box-shadow:0 0 2em currentColor,0 0 .5em currentColor;position:absolute}.xt3d-divgraphic-animation:after,.xt3d-divgraphic-animation:before,.xt3d-divgraphic-animation p:after,.xt3d-divgraphic-animation p:before{content:"";position:absolute;width:100%;height:100%;left:50%;top:50%;border-radius:50%;transform:translate(-50%,-50%)}.xt3d-divgraphic-animation:after,.xt3d-divgraphic-animation:before{border:1px solid;animation:xt3d-mapAni 1s ease infinite}.xt3d-divgraphic-animation p:before{border:1px solid}.xt3d-divgraphic-animation p{position:absolute;left:50%;top:50%;width:0;height:0;border-radius:50%;transform:translate(-50%,-50%);animation:xt3d-mapAni 2s ease infinite}@-webkit-keyframes xt3d-mapAni{0%{width:0;height:0;opacity:1;filter:alpha(opacity=1)}25%{width:12px;height:12px;opacity:.7;filter:alpha(opacity=70)}50%{width:20px;height:20px;opacity:.5;filter:alpha(opacity=50)}75%{width:30px;height:30px;opacity:.2;filter:alpha(opacity=20)}to{width:40px;height:40px;opacity:0;filter:alpha(opacity=0)}}@-moz-keyframes xt3d-mapAni{0%{width:0;height:0;opacity:1;filter:alpha(opacity=1)}25%{width:12px;height:12px;opacity:.7;filter:alpha(opacity=70)}50%{width:20px;height:20px;opacity:.5;filter:alpha(opacity=50)}75%{width:30px;height:30px;opacity:.2;filter:alpha(opacity=20)}to{width:40px;height:40px;opacity:0;filter:alpha(opacity=0)}}@-o-keyframes xt3d-mapAni{0%{width:0;height:0;opacity:1;filter:alpha(opacity=1)}25%{width:12px;height:12px;opacity:.7;filter:alpha(opacity=70)}50%{width:20px;height:20px;opacity:.5;filter:alpha(opacity=50)}75%{width:30px;height:30px;opacity:.2;filter:alpha(opacity=20)}to{width:40px;height:40px;opacity:0;filter:alpha(opacity=0)}}@-ms-keyframes xt3d-mapAni{0%{width:0;height:0;opacity:1;filter:alpha(opacity=1)}25%{width:12px;height:12px;opacity:.7;filter:alpha(opacity=70)}50%{width:20px;height:20px;opacity:.5;filter:alpha(opacity=50)}75%{width:30px;height:30px;opacity:.2;filter:alpha(opacity=20)}to{width:40px;height:40px;opacity:0;filter:alpha(opacity=0)}}@keyframes xt3d-mapAni{0%{width:0;height:0;opacity:1;filter:alpha(opacity=1)}25%{width:12px;height:12px;opacity:.7;filter:alpha(opacity=70)}50%{width:20px;height:20px;opacity:.5;filter:alpha(opacity=50)}75%{width:30px;height:30px;opacity:.2;filter:alpha(opacity=20)}to{width:40px;height:40px;opacity:0;filter:alpha(opacity=0)}}.xt3d-divgraphic-brighten{position:absolute;left:0;bottom:0;display:block}.divpoint-wrap{position:relative;padding:30px;overflow:hidden}.divpoint .arrow{position:absolute;bottom:0;left:0;width:45px;height:2px;transform:rotate(-45deg) translate(5px,-15px)}.divpoint .area{position:relative;min-width:180px;min-height:150px}.divpoint .b-t-l{position:absolute;top:0;left:0;width:1px;height:62px;transform:rotate(45deg) translate(52px,-22px);z-index:10}.divpoint .b-b-r{position:absolute;bottom:0;right:0;width:1px;height:62px;transform:rotate(45deg) translate(-52px,22px);z-index:10}.divpoint .b-t{position:absolute;top:0;left:44px;right:0;height:1px;z-index:10}.divpoint .b-r{position:absolute;top:0;right:0;bottom:44px;width:1px;z-index:10}.divpoint .b-b{position:absolute;left:0;right:44px;bottom:0;height:1px;z-index:10}.divpoint .b-l{position:absolute;top:44px;left:0;bottom:0;width:1px;z-index:10}.divpoint .label-wrap{padding-left:12px;color:#fff;font-size:16px;white-space:nowrap;overflow:hidden}.divpoint .title{margin-top:20px;padding:0 12px 0 30px;height:36px;line-height:36px;position:relative}.divpoint .title-line{width:calc(100% - 2px);height:2px;background:#4984ed;margin:2px 0}@keyframes dynamic-border-label-animate{0%,to{clip:rect(0px,var(--clip-width-1),2px,0px)}25%{clip:rect(0px,2px,var(--clip-height-1),0px)}50%{clip:rect(var(--clip-height-2),var(--clip-width-1),var(--clip-width-1),0px)}75%{clip:rect(0px,var(--clip-width-1),var(--clip-height-1),var(--clip-width-2))}}.dynamic-border-label-container{position:absolute;left:0;bottom:0;cursor:pointer;--clip-height-1: 40px;--clip-height-2: 38px;--clip-width-1: 165px;--clip-width-2: 163px;--text-left-position: -75px;--animation-name: dynamic-border-label-animate}.xt3d-component-animate-marker_boder{margin:auto;color:#15d1f2;box-shadow:inset 0 0 0 1px #15d1f2}.xt3d-component-animate-marker_text{color:#fff;font-size:14px;display:flex;width:100%;height:100%;align-items:center;justify-content:center;font-weight:bolder;-webkit-user-select:none;user-select:none;cursor:pointer}.xt3d-component-animate-marker_boder,.xt3d-component-animate-marker_boder:before,.xt3d-component-animate-marker_boder:after{position:absolute;top:0;bottom:0;left:0;right:0}.xt3d-component-animate-marker_boder:before,.xt3d-component-animate-marker_boder:after{content:"";margin:-5%;box-shadow:inset 0 0 0 2px;animation:var(--animation-name) 8s linear infinite}.xt3d-component-animate-marker_boder:before{animation-delay:-4s}.xt3d-divgraphic-multiline{min-width:96px;min-height:35px;background:linear-gradient(0deg,#1e202a 0%,#0d1013 100%);cursor:default;border:1px solid #14171c;box-shadow:0 2px 21px #2122278c;border-radius:4px;box-sizing:border-box}.xt3d-divgraphic-multiline-img_0{width:calc(100% + 22px);height:38px;position:absolute;bottom:-39px;left:-22px;background-position:0px 0px}.xt3d-divgraphic-multiline-img_1{width:100%;height:1px;position:absolute;bottom:0}.xt3d-divgraphic-multiline-img_1:after{content:"";display:block;position:absolute;width:0;height:0;left:50%;margin-left:-11px;bottom:-9px;border-left:11px solid transparent;border-right:11px solid transparent;border-top:11px solid #1e202a}.xt3d-divgraphic-multiline-text{width:100%;height:100%;text-align:center;padding:8px 20px;font-size:14px;font-family:MicrosoftYaHei;font-weight:400;-webkit-box-sizing:border-box;box-sizing:border-box;border-top:1px solid #303030;text-wrap:nowrap}.xt3d-rectangle-label-text{width:100%;height:100%;text-align:center;padding:5px 20px;font-size:14px;font-family:MicrosoftYaHei;color:#fff;box-sizing:border-box;text-wrap:nowrap}.upright-label-container{width:15px;text-align:center;background:transparent;font-size:9px;color:#fff;font-family:Microsoft YaHei,Helvetica Neue For Number,-apple-system,BlinkMacSystemFont,Segoe UI,Roboto,Hiragino Sans GB,PingFang SC,Helvetica Neue,Helvetica,Arial,sans-serif!important;box-sizing:border-box;font-weight:600;display:flex;flex-direction:column;align-items:center}.upright-label-item{writing-mode:vertical-lr;font-size:16px;letter-spacing:4px}.pre-topCard-list-item-line{display:block;height:100px;width:1px;margin-top:3px;background-color:#fff}.pre-topCard-list-item-circle{width:10px;height:10px;background-color:#fff;border-radius:50%;margin-top:-10px}.shining{color:#cce7f8;font-size:2.5rem;-webkit-animation:shining .5s alternate infinite;animation:shining .5s alternate infinite}@keyframes shining{0%{text-shadow:0 0 10px lightblue,0 0 20px lightblue,0 0 30px lightblue,0 0 40px skyblue,0 0 50px skyblue,0 0 60px skyblue}to{text-shadow:0 0 5px lightblue,0 0 10px lightblue,0 0 15px lightblue,0 0 20px skyblue,0 0 25px skyblue,0 0 30px skyblue}}.maske{letter-spacing:.2rem;font-size:32px;font-weight:700;background-image:-webkit-linear-gradient(left,#147B96,#E6D205 25%,#147B96 50%,#E6D205 75%,#147B96);-webkit-text-fill-color:transparent;-webkit-background-clip:text;background-clip:text;background-size:200% 100%;animation:maskedAnimation 4s infinite linear}@keyframes maskedAnimation{0%{background-position:0 0}to{background-position:-100% 0}}.flame{text-shadow:0 0 20px #fefcc9,10px -10px 30px #feec85,-20px -20px 40px #ffae34,20px -40px 50px #ec760c,0 -80px 70px #f38e1c;font-weight:700;font-size:32px;font-family:微软雅黑;color:#fff;text-align:center;animation:flame 2s infinite}@keyframes flame{0%{text-shadow:0 0 20px #fefcc9,10px -10px 30px #feec85,-20px -20px 40px #ffae34,20px -40px 50px #ec760c,0 -80px 70px #f38e1c}30%{text-shadow:0 0 20px #fefcc9,10px -10px 30px #feec85,-20px -20px 40px #ffae34,20px -40px 50px #ec760c,10px -90px 80px #f38e1c}60%{text-shadow:0 0 20px #fefcc9,10px -10px 30px #feec85,-20px -20px 40px #ffae34,20px -40px 50px #ec760c,0px -80px 100px #f38e1c}to{text-shadow:0 0 20px #fefcc9,10px -10px 30px #feec85,-20px -20px 40px #ffae34,20px -40px 50px #ec760c,0 -80px 70px #f38e1c}}.divIndicator-fixed{width:10px;height:10px;border:2px solid #fff;border-radius:50%;background:blue;position:absolute;z-index:100}.divIndicator-drag{position:absolute;cursor:pointer;margin:0}.divIndicator-line{position:absolute;width:0;top:6px;left:6px;z-index:99}.graphic-draw-tip-container{position:absolute;background:rgba(0,0,0,.637);padding:6px;color:#fff;pointer-events:none;white-space:nowrap;z-index:5}.graphic-draw-tip-container:before{position:absolute;content:"";top:calc(50% - 10px);left:-10px;border-bottom:10px solid transparent;border-top:10px solid transparent;border-right:10px solid rgba(0,0,0,.637)}.xt3d-point-picker-container{background:#e3eee59e;position:absolute;left:0;bottom:0;cursor:default;padding:0 2px 5px 5px;border:1px solid #9c9944e8}.xt3d-point-picker-container:before{position:absolute;content:"";left:50%;bottom:-50px;height:48px;border-left:2px dashed #ffffff70;pointer-events:none}.xt3d-point-picker-container .info-item-close{position:relative;height:15px;cursor:pointer}.xt3d-point-picker-container .info-close{display:block;height:15px;width:15px;position:absolute;top:1px;right:0;background:#dc9e9e4d;line-height:15px;text-align:center}.xt3d-point-picker-container .info-close:hover{background:#e01c1c6c}.xt3d-point-picker-container .info-input{height:18px;width:140px;line-height:18px;font-size:18px}.xt3d-point-picker-container .info-input1{height:18px;width:50px;line-height:18px;font-size:18px}.xt3d-point-picker-container .value-type{background:#dadae3;padding:0 3px;border:1px solid #5f5f66;cursor:pointer;margin-left:1px}.xt3d-point-picker-container .value-type-activate,.xt3d-point-picker-container .value-type:hover{background:#64ada0a1;color:#fff}.xt3d-point-picker-container textarea{resize:none;text-wrap:nowrap;overflow:hidden}.xt3d-point-picker-container button{outline:none;border:0;display:block;margin:2px 5px 0;cursor:pointer;background:#b7b1b196;box-shadow:inset 0 0 0 1px #575757,0 2px 21px #2122278c;opacity:.96;border-radius:2px}.xt3d-point-picker-container button:hover{background:#64ada0a1;color:#fff}.xt3d-camera-picker-container{background:#e3eee59e;position:absolute;left:0;bottom:0;cursor:default;padding:0 2px 5px 5px;border:1px solid #9c9944e8}.xt3d-camera-picker-container .info-item-close{position:relative;height:15px;cursor:pointer}.xt3d-camera-picker-container .info-close{display:block;height:15px;width:15px;position:absolute;top:1px;right:0;background:#dc9e9e4d;line-height:15px;text-align:center}.xt3d-camera-picker-container .info-close:hover{background:#e01c1c6c}.xt3d-camera-picker-container .info-input{height:18px;width:140px;line-height:18px;font-size:18px}.xt3d-camera-picker-container .info-input1{height:18px;width:50px;line-height:18px;font-size:18px}.xt3d-camera-picker-container textarea{resize:none;text-wrap:nowrap;overflow:hidden}.xt3d-camera-picker-container button{outline:none;border:0;display:block;margin:2px 5px 0;cursor:pointer;background:#b7b1b196;box-shadow:inset 0 0 0 1px #575757,0 2px 21px #2122278c;opacity:.96;border-radius:2px}.xt3d-camera-picker-container button:hover{background:#64ada0a1;color:#fff}.xt3d-popup-container{position:absolute;left:0;bottom:0}.xt3d-popup-show{opacity:1;transition:opacity 2s}.xt3d-popup-hide{opacity:0;transition:opacity 2s}.xt3d-popup-close-button{position:absolute;top:0;right:0;padding:4px 4px 0 0;text-align:center;width:18px;height:14px;font:16px/14px Tahoma,Verdana,sans-serif;color:#c3c3c3;text-decoration:none;font-weight:700;background:transparent;cursor:pointer;-webkit-user-select:none;user-select:none}.xt3d-popup-close-button:hover{color:#c3080893}.xt3d-popup-header{height:24px;border-bottom:1px solid #939393;padding:4px}.xt3d-popup-wrapper{text-align:center;min-height:50px;min-width:100px;overflow-y:auto;background-color:#3f4854e0;box-shadow:0 3px 14px #0006;text-align:left;border-radius:4px;padding:8px;color:#fff}.xt3d-popup-tip-container{margin:0 auto;width:40px;height:20px;position:relative;overflow:hidden}.xt3d-popup-tip{background-color:#3f4854e0;box-shadow:0 3px 14px #0006;width:17px;height:17px;padding:1px;margin:-10px auto 0;transform:rotate(45deg)}',
    'remarks',
    '.json',
    'rippleSize',
    '_off',
    '_graphic',
    'TimeInterval',
    'text/json',
    'shader',
    '_shaderGeometryOpts',
    'TextureUniform',
    'SampledPositionProperty',
    'NEAREST',
    '_mouseMove',
    '_createPrimitives',
    'entity',
    '_setSvgDom',
    'minHeight',
    'handleLeftDown',
    'lightShadowMapCube',
    'fromDegreesArrayHeights',
    'color:  ',
    'destroyObject',
    '_getOuterFovRadiusPairs',
    'hasEdit',
    '_cameraDistance',
    'translucencyByDistance',
    'verticalPoistion',
    'hasOwnProperty',
    'out_FragColor = vec4(out_FragColor.rgb + u_color, (out_FragColor.g + out_FragColor.r) * .6 );',
    'floatMarker',
    'VEC3',
    'deselectedGraphic',
    'OPAQUE',
    'sampler2D',
    'eType',
    '_updateLineStyle',
    'bgFS',
    '_createCircleInstance',
    'setAttribute',
    '_coordinates',
    '_nextStopsIndex',
    '_markerBg',
    'shaderEffetNeonLight',
    'PolylineFlowMaterialProperty',
    'PolylineFlowType',
    '  \nin vec3 v_position; \nin vec3 v_normal;  \nuniform vec4  defaultColor; \nuniform float specular; \nuniform float shininess; \nuniform vec3  emission; \nin vec2 v_st; \nuniform bool isLine; \nuniform float glowPower; \nuniform vec4 czm_pickColor;\nvoid main() {\nvec3 positionToEyeEC = -v_position; \nvec3 normalEC =normalize(v_normal);\nvec4 color=defaultColor;\n \n//if(v_st.x<0.5){\n//    color.a =0.75-v_st.x; \n//}\n//else  {\n//    color.a =v_st.x-0.25; \n//}\nczm_material material;\nmaterial.specular = specular;\nmaterial.shininess = shininess;\n// material.normal =  normalEC;\nmaterial.emission =emission;//vec3(0.2,0.2,0.2);\nmaterial.diffuse = color.rgb ;\nif(isLine){\nmaterial.alpha = 1.0; \n}\nelse{\nmaterial.alpha =  color.a; \n}\n//float glow = glowPower / abs(v_st.t  ) - (glowPower / 0.5); \n// \n//material.emission = max(vec3(glow - 1.0 + color.rgb), color.rgb); \n//if(isLine)\n//    material.alpha = clamp(0.0, 1.0, glow) * color.a; \n \n\nif(v_st.x==0.0){ \n  \n    out_FragColor = color ;\n\n}else { \n\n    out_FragColor =color ; \n} \n\n}\n',
    'ShinyRingFS',
    'upright-label-container',
    'addColorStop',
    'drivenDistance',
    'offsetLeft',
    'divIndicator-line',
    'customTemplate',
    'ColorMaterialProperty',
    'resize',
    'RGBA',
    'polyline',
    'range',
    'startPosition',
    'inverse',
    'pre-topCard-list-item-circle',
    '_writeLeftDow',
    '着色器特效基础对象',
    '_createSvgEditAnchors',
    'Matrix3',
    'heatStyle',
    'backgroundColor',
    'fontBold',
    'CylinderGlowGradientWallType',
    'MITERED',
    'asset',
    '_modelMatrixNeedsUpdate',
    '_videoCamera',
    '_onUpdate',
    '_shadowMapMatrix',
    '_poiIcon',
    '_initEvents',
    'texture(image, fract(materialInput.st)).a * alpha',
    'topWidth',
    '#ffffff6b',
    'liquidfill',
    '#00ffff',
    'function',
    'specular',
    'lineWidth',
    '_colorTexture',
    ';box-shadow: 0 0 10px 2px ',
    '当前绘制类型：',
    'box',
    '_lineDom',
    'custom',
    'getElementsByClassName',
    '当前测量类型：高度测量',
    '_bottomWidth',
    'specularIntensity',
    'value-type-activate',
    '光椎体',
    'fromVertices',
    'rgba(255, 255, 255, 0.14)',
    '_createLabelDom',
    '_hasEdit',
    '__addHook',
    '点击鼠标左键确定第1个点的位置',
    'GIF图标点',
    'TEXTURE_2D',
    '_boundingSphereWC',
    '_vertexEntity_1',
    'ShadowMode',
    'subSegmentV',
    '_createHLSVideoEle',
    '_addSelectedDom',
    'fromArray',
    '_hdr',
    'plane',
    '64px Microsoft YaHei',
    'shaderEffetWavePetals',
    'height:',
    'imageryLayers',
    '着色器特效爆炸',
    'readAsArrayBuffer',
    'negate',
    'brightness',
    'indexOf',
    'StencilFunction',
    'body',
    'u_radius',
    'UrlTemplateImageryProvider',
    'rgba(255, 255, 255, 1.0)',
    'download',
    'sin',
    'FrustumOutlineGeometry',
    'bold ',
    '着色器特效烟花',
    '_subSegmentH',
    'padding-left: 10px;',
    'GlowPyramidFS',
    'xt3d-component-animate-marker_boder',
    'proxy',
    'HeightReference',
    '_currentRadius',
    'raiseToTop',
    'createPickId',
    ' / ',
    'depthMask',
    'polyGradient',
    '__proto__',
    '着色器特效多彩星链',
    'get',
    '_htmlGraphics',
    'availability',
    'DrawCommand',
    'upright-label-item',
    '_getColor',
    '方位值：',
    'billboard',
    '376FZukjr',
    'intersectionWidth',
    '\n        uniform sampler2D colorTexture;\n        uniform sampler2D depthTexture;\n         \n        uniform vec4 lightPositionEC;\n        uniform float intensity;\n        uniform vec3 lightColor;\n        uniform vec3 direction;\n        uniform float outerConeCos;\n        uniform float innerConeCos;\n        uniform vec2 shadowMapDarknessType;\n        \n        in vec2 v_textureCoordinates;\n        \n        const float M_PI = 3.141592653589793;\n        \n        vec3 getEyeCoordinate3FromWindowCoordinateMars3D(vec2 fragCoord, float logDepthOrDepth) {\n          vec4 eyeCoordinate = czm_windowToEyeCoordinates(fragCoord, logDepthOrDepth);\n          return eyeCoordinate.xyz / eyeCoordinate.w;\n        }\n        \n        vec3 vectorFromOffsetMars3D(vec4 eyeCoordinate, vec2 positiveOffset) {\n          vec2 glFragCoordXY = v_textureCoordinates.xy * czm_viewport.zw;\n          float upOrRightLogDepth = czm_unpackDepth(texture(depthTexture, (glFragCoordXY + positiveOffset) / czm_viewport.zw));\n          float downOrLeftLogDepth = czm_unpackDepth(texture(depthTexture, (glFragCoordXY - positiveOffset) / czm_viewport.zw));\n        \n          bvec2 upOrRightInBounds = lessThan(glFragCoordXY + positiveOffset, czm_viewport.zw);\n          float useUpOrRight = float(upOrRightLogDepth > 0.0 && upOrRightInBounds.x && upOrRightInBounds.y);\n          float useDownOrLeft = float(useUpOrRight == 0.0);\n          vec3 upOrRightECMars3D = getEyeCoordinate3FromWindowCoordinateMars3D(glFragCoordXY + positiveOffset, upOrRightLogDepth);\n          vec3 downOrLeftEC = getEyeCoordinate3FromWindowCoordinateMars3D(glFragCoordXY - positiveOffset, downOrLeftLogDepth);\n          return (upOrRightECMars3D - (eyeCoordinate.xyz / eyeCoordinate.w)) * useUpOrRight + ((eyeCoordinate.xyz / eyeCoordinate.w) - downOrLeftEC) * useDownOrLeft;\n        }\n        \n        float getRangeAttenuationMars3D(float range, float d) {\n          if(range <= 0.0) {\n            return 1.0 / pow(d, 2.0);\n          }\n          return max(min(1.0 - pow(d / range, 4.0), 1.0), 0.0) / pow(d, 2.0);\n        }\n        \n        float getSpotAttenuationMars3D(vec3 pointToLight, vec3 direction, float outerConeCos, float innerConeCos) {\n          float actualCos = dot(normalize(direction), normalize(-pointToLight));\n          if(actualCos > outerConeCos) {\n            if(actualCos < innerConeCos) {\n              return smoothstep(outerConeCos, innerConeCos, actualCos);\n            }\n            return 1.0;\n          }\n          return 0.0;\n        }\n        \n        vec3 getLightIntensityMars3D(vec3 color, float intensity, float type, float range, vec3 pointToLight, vec3 direction, float outerConeCos, float innerConeCos) {\n          float rangeAttenuation = 1.0;\n          float spotAttenuation = 1.0;\n          rangeAttenuation = getRangeAttenuationMars3D(range, length(pointToLight));\n          if(type == 2.0) {\n            spotAttenuation = getSpotAttenuationMars3D(pointToLight, direction, outerConeCos, innerConeCos);\n          }\n          return rangeAttenuation * spotAttenuation * intensity * color;\n        }\n        \n        void main() {\n          vec4 color = texture(colorTexture, v_textureCoordinates);\n          float logDepthOrDepth = czm_unpackDepth(texture(depthTexture, v_textureCoordinates));\n          if(logDepthOrDepth >= 1.0) {\n            out_FragColor =  color;\n            return;\n          }  \n          vec4 eyeCoordinate = czm_windowToEyeCoordinates(v_textureCoordinates.xy * czm_viewport.zw, logDepthOrDepth);\n          vec3 downUp = vectorFromOffsetMars3D(eyeCoordinate, vec2(0.0, 1.0));\n          vec3 leftRight = vectorFromOffsetMars3D(eyeCoordinate, vec2(1.0, 0.0));\n          vec3 normalEC = normalize(cross(leftRight, downUp));\n          vec3 positionEC = eyeCoordinate.xyz / eyeCoordinate.w;\n        \n          vec3 totalColor = vec3(0.0);\n        \n          \n            vec4 lightPEC = lightPositionEC;\n            vec2 shadowMapDT = shadowMapDarknessType;\n        \n            vec3 pointToLightEC = positionEC - lightPEC.xyz;\n            float pointToLightECLength = length(pointToLightEC);\n            vec3 l = normalize(pointToLightEC);\n            float NdotL = clamp(dot(-normalEC, l), 0.0, 1.0);\n        \n            float type = shadowMapDT.y;\n            vec3 colorIntensity = getLightIntensityMars3D(lightColor, intensity, type, lightPEC.w, pointToLightEC, direction, outerConeCos, innerConeCos);\n            totalColor += NdotL * colorIntensity;\n          \n        \n          out_FragColor = vec4(color.xyz + totalColor, 1.0);\n        } ',
    'u_invert',
    'rayPlane',
    'framebuffer',
    '#8a1ce9',
    'TextureWrap',
    'fromTranslationRotationScale',
    'pixelOffset',
    '_radius',
    '点击鼠标右键结束绘制',
    'value',
    '_vertexEntity_',
    'userData',
    'orientation',
    'cartesianToCartographic',
    '_lineEntity',
    'shaderEffetStellarChain',
    'u_type',
    '_editAnchorIndex',
    '\nczm_material czm_getMaterial(czm_materialInput materialInput)\n{\n     czm_material material = czm_getDefaultMaterial(materialInput);\n     vec2 st = materialInput.st;\n     float time =czm_frameNumber * speed / 1000.0;\n     vec4 colorImage0 = texture(image0, vec2(fract(st.t), st.t));\n     vec4 colorImage1 = texture(image1, vec2(fract(st.t - time), st.t));\n     material.alpha = (colorImage0.a+colorImage1.a)*color.a;\n     material.diffuse =  1.5 * color.rgb  ;  \n     return material;\n } \n ',
    '_fragmentShader',
    'createGeometry',
    '_pickRS',
    '_offCenterFrustum',
    '墙体渐变材质',
    '_getInnerFovRadiusParire',
    'combine',
    'with',
    '_initDrawCommands',
    'attackArrow',
    'fromRotationTranslation',
    '_heatStyle',
    'lessThanOrEquals',
    'info-close',
    '\n            const int NUM_STEPS = 8;\n            const float PI     = 3.141592;\n            const float EPSILON  = 1e-3;\n            //#define EPSILON_NRM (0.1 / iResolution.x)\n            #define EPSILON_NRM (0.1 / 200.0)\n            // sea\n            const int ITER_GEOMETRY = 3;\n            const int ITER_FRAGMENT = 5;\n            const float SEA_HEIGHT = 0.6;\n            const float SEA_CHOPPY = 4.0;\n            const float SEA_SPEED = 1.8;\n            const float SEA_FREQ = 0.16;\n            const vec3 SEA_BASE = vec3(0.1,0.19,0.22);\n            const vec3 SEA_WATER_COLOR = vec3(0.8,0.9,0.6);\n            //#define SEA_TIME (1.0 + iTime * SEA_SPEED)\n            const mat2 octave_m = mat2(1.6,1.2,-1.2,1.6);\n            float iTime;\n            // math\n            mat3 fromEuler(vec3 ang) {\n              vec2 a1 = vec2(sin(ang.x),cos(ang.x));\n              vec2 a2 = vec2(sin(ang.y),cos(ang.y));\n              vec2 a3 = vec2(sin(ang.z),cos(ang.z));\n              mat3 m;\n              m[0] = vec3(a1.y*a3.y+a1.x*a2.x*a3.x,a1.y*a2.x*a3.x+a3.y*a1.x,-a2.y*a3.x);\n              m[1] = vec3(-a2.y*a1.x,a1.y*a2.y,a2.x);\n              m[2] = vec3(a3.y*a1.x*a2.x+a1.y*a3.x,a1.x*a3.x-a1.y*a3.y*a2.x,a2.y*a3.y);\n              return m;\n            }\n            float hash( vec2 p ) {\n              float h = dot(p,vec2(127.1,311.7));\n              return fract(sin(h)*43758.5453123);\n            }\n            float noise( in vec2 p ) {\n              vec2 i = floor( p );\n              vec2 f = fract( p );\n              vec2 u = f*f*(3.0-2.0*f);\n              return -1.0+2.0*mix( mix( hash( i + vec2(0.0,0.0) ),\n                       hash( i + vec2(1.0,0.0) ), u.x),\n                    mix( hash( i + vec2(0.0,1.0) ),\n                       hash( i + vec2(1.0,1.0) ), u.x), u.y);\n            }\n            // lighting\n            float diffuse(vec3 n,vec3 l,float p) {\n              return pow(dot(n,l) * 0.4 + 0.6,p);\n            }\n            float specular(vec3 n,vec3 l,vec3 e,float s) {\n              float nrm = (s + 8.0) / (PI * 8.0);\n              return pow(max(dot(reflect(e,n),l),0.0),s) * nrm;\n            }\n            // sky\n            vec3 getSkyColor(vec3 e) {\n              e.y = max(e.y,0.0);\n              return vec3(pow(1.0-e.y,2.0), 1.0-e.y, 0.6+(1.0-e.y)*0.4);\n            }\n            // sea\n            float sea_octave(vec2 uv, float choppy) {\n              uv += noise(uv);\n              vec2 wv = 1.0-abs(sin(uv));\n              vec2 swv = abs(cos(uv));\n              wv = mix(wv,swv,wv);\n              return pow(1.0-pow(wv.x * wv.y,0.65),choppy);\n            }\n            float map(vec3 p) {\n              float freq = SEA_FREQ;\n              float amp = SEA_HEIGHT;\n              float choppy = SEA_CHOPPY;\n              vec2 uv = p.xz; uv.x *= 0.75;\n              float d, h = 0.0;\n              float SEA_TIME = 1.0 + iTime * SEA_SPEED;\n              for(int i = 0; i < ITER_GEOMETRY; i++) {\n                d = sea_octave((uv+SEA_TIME)*freq,choppy);\n                d += sea_octave((uv-SEA_TIME)*freq,choppy);\n                h += d * amp;\n                uv *= octave_m; freq *= 1.9; amp *= 0.22;\n                choppy = mix(choppy,1.0,0.2);\n              }\n              return p.y - h;\n            }\n            float map_detailed(vec3 p) {\n              float freq = SEA_FREQ;\n              float amp = SEA_HEIGHT;\n              float choppy = SEA_CHOPPY;\n              vec2 uv = p.xz; uv.x *= 0.75;\n              float SEA_TIME = 1.0 + iTime * SEA_SPEED;\n              float d, h = 0.0;\n              for(int i = 0; i < ITER_FRAGMENT; i++) {\n                d = sea_octave((uv+SEA_TIME)*freq,choppy);\n                d += sea_octave((uv-SEA_TIME)*freq,choppy);\n                h += d * amp;\n                uv *= octave_m; freq *= 1.9; amp *= 0.22;\n                choppy = mix(choppy,1.0,0.2);\n              }\n              return p.y - h;\n            }\n            vec3 getSeaColor(vec3 p, vec3 n, vec3 l, vec3 eye, vec3 dist) {\n              float fresnel = clamp(1.0 - dot(n,-eye), 0.0, 1.0);\n              fresnel = pow(fresnel,3.0) * 0.65;\n              vec3 reflected = getSkyColor(reflect(eye,n));\n              vec3 refracted = SEA_BASE + diffuse(n,l,80.0) * SEA_WATER_COLOR * 0.12;\n              vec3 color = mix(refracted,reflected,fresnel);\n              float atten = max(1.0 - dot(dist,dist) * 0.001, 0.0);\n              color += SEA_WATER_COLOR * (p.y - SEA_HEIGHT) * 0.18 * atten;\n              color += vec3(specular(n,l,eye,60.0));\n              return color;\n            }\n            // tracing\n            vec3 getNormal(vec3 p, float eps) {\n              vec3 n;\n              n.y = map_detailed(p);\n              n.x = map_detailed(vec3(p.x+eps,p.y,p.z)) - n.y;\n              n.z = map_detailed(vec3(p.x,p.y,p.z+eps)) - n.y;\n              n.y = eps;\n              return normalize(n);\n            }\n            float heightMapTracing(vec3 ori, vec3 dir, out vec3 p) {\n              float tm = 0.0;\n              float tx = 1000.0;\n              float hx = map(ori + dir * tx);\n              if(hx > 0.0) return tx;\n              float hm = map(ori + dir * tm);\n              float tmid = 0.0;\n              for(int i = 0; i < NUM_STEPS; i++) {\n                tmid = mix(tm,tx, hm/(hm-hx));\n                p = ori + dir * tmid;\n                float hmid = map(p);\n                if(hmid < 0.0) {\n                  tx = tmid;\n                  hx = hmid;\n                } else {\n                  tm = tmid;\n                  hm = hmid;\n                }\n              }\n              return tmid;\n            }\n              czm_material czm_getMaterial(czm_materialInput materialInput)\n                 {\n                  czm_material material = czm_getDefaultMaterial(materialInput);\n                  vec2 uv = materialInput.st;\n                  uv = uv * 2.0 - 1.0;\n                  iTime = czm_frameNumber /500.0; \n                  // ray\n                  vec3 ang = vec3(0, 1.2, 0.0);\n                    vec3 ori = vec3(0.0,3.5,0);\n                  vec3 dir = normalize(vec3(uv.xy,-2.0)); dir.z += length(uv) * 0.15;\n                  dir = normalize(dir) * fromEuler(ang);\n                  // tracing\n                  vec3 p;\n                  heightMapTracing(ori,dir,p);\n                  vec3 dist = p - ori;\n                  vec3 n = getNormal(p, dot(dist,dist) * EPSILON_NRM);\n                  vec3 light = normalize(vec3(0.0,1.0,0.8)); \n                  // color\n                  vec3 color = mix(\n                    getSkyColor(dir),\n                    getSeaColor(p,n,czm_lightDirectionWC,dir,dist),\n                    pow(smoothstep(0.0,-0.05,dir.y),0.3)); \n                    material.diffuse = pow(color,vec3(0.75));\n                    material.alpha =0.99;\n                     return material;\n                 }\n              ',
    '_setCartesian3Array',
    'topShow',
    'loadFromGeoJson',
    'headingPitchRollQuaternion',
    'fromPoints',
    'WallFlowMaterialProperty',
    'VirusFS',
    'endPosition',
    '_clear',
    '点光源',
    'editEnd',
    'PLANE',
    'findIndex',
    '_mousePoint',
    'showIntersection',
    '_separator',
    'breathe',
    '_originalWorldPosition',
    'u_shadowMapCube',
    '_tileset',
    '_matrix',
    '_speed',
    'block',
    '_initAnchorPositions',
    'BLUE',
    'cartesian3Array',
    'drawingBufferToWgs84Coordinates',
    'onDomClick',
    'curvyLine',
    '_conetntDom',
    'xt3d-heatmap-container',
    'distanceDisplayCondition',
    '_lookAtState',
    'createRadialGradient',
    'appearanceOpts',
    'outerConeCos',
    'URL',
    ',\n                y: ',
    'MaterialAppearance',
    'trailLine',
    '_entity',
    'xt3d-component-animate-marker_text',
    'isConstant',
    'drivenLineWidth',
    'SphereElectricMaterialProperty',
    '_content',
    '_domeLineCommand',
    'CircleScanType_3',
    'vertexShader',
    '_createRingCanvas',
    '_reflectMatrix',
    '_renderState',
    'multiply',
    'segmentV',
    '\n#ifdef GL_ES\nprecision mediump float;\n#endif\n  in vec2 v_st;\n  vec2  resolution;\n  vec2  mouse;\n  float time;\n \n\nconst int num_x = 5;\nconst int num_y = 5;\n\nvec4 draw_ball(int i, int j) {\n\n    float w = resolution.x/10.0;\n    float h = resolution.y/5.0;\n \n\n\tfloat t = time*0.2;\n\tfloat x = w  * (1.0 + cos(1.5 * t + float(2*i+4*j))) ; \n\tfloat y = h  * (1.0 + sin(2.3 * t + float(3*i+4*j))); \n\tfloat size = 4.0 - 2.0 * sin(t);\n\tvec2 pos = vec2(x, y);\n\tfloat dist = length(v_st.xy * 600.0 - 100.0 - pos);\n\tfloat intensity = pow(size/(dist/2.0), 2.0);\n\tvec4 color = vec4(0.0);\n\tcolor.r = 0.5 + cos(t*float(i));\n\tcolor.g = 0.5 + sin(t*float(j));\n\tcolor.b = 0.5 + sin(float(j));\n\treturn color*intensity;\n}\n\nvoid main() {\n    mouse=vec2(1.0,1.0);\n    time = czm_frameNumber /100.0;\n    resolution = vec2(1920.0,1080.0);//czm_viewport.zw;\n\n\tvec4 color = vec4(0.0);\n\tfor (int i = 0; i < num_x; ++i) {\n\t\tfor (int j = 0; j < num_y; ++j) {\n\t\t\tcolor += draw_ball(i, j);\n\t\t}\n\t}\n\t \n    float alpha=(color.x+ color.y + color.z) ; \n    if(alpha < 0.2){ \n        discard; \n    } \n\tout_FragColor = color*2.0;// + shadow; \n\tout_FragColor.a =alpha; //1.0;\n}',
    'postProcess',
    'createEditAnchor',
    'properties',
    '87507bklKWL',
    'PolylineTrailMaterialProperty',
    '_ringImage',
    '_mouseMoveEvent',
    'velocity',
    '_headTailFactor',
    '_createInnerCurveCommand',
    'FILL_AND_OUTLINE',
    '图片标签',
    '_spaceDistance',
    '_topWidth',
    '_computeVAO',
    'routeChange',
    'lineColor',
    'getContext',
    '_materialAttachTarget',
    'CylinderGlowFlowWallSource',
    'fromQuaternion',
    'VerticalOrigin',
    '_boundingSphere',
    'ExplosionFS',
    '_computeVertexNormals',
    'reverse',
    'substring',
    'drivenLineColor',
    'createIfNeeded',
    'stopPropagation',
    '_screenSpaceEventHandler',
    '_superGif',
    '视频投影',
    '_materialOpts',
    'enable',
    '着色器特效辐射花瓣',
    'title',
    'wallGradient',
    'maxScaleHeight',
    'lastIndexOf',
    'WallGradientMaterialProperty',
    '拖拽旋转',
    'reflectorProjectionMatrix',
    '){\n                    material.diffuse =color.rgb;\n                    return; \n                }\n               \n                // 动态光环\n                float time = fract(czm_frameNumber / 360.0);\n                time = abs(time - 0.5) * 2.0;\n                float diff = step(0.005, abs( clamp(position.',
    '_computeMatrix',
    'DYNAMIC_DRAW',
    'angleT',
    'rgba(255, 255, 255,255)',
    '着色器特效火焰云',
    'PostProcessStageSampleMode',
    'colorTexture',
    'upright',
    'alignedAxis',
    '_trackedView',
    '_colorSubscription',
    'VertexArray',
    'infiniteProjectionMatrix',
    'graphicClassType',
    'provider',
    'fabric',
    'normalMap',
    'imageryProvider',
    'circleScan_1',
    'u_shadowMapDarkness',
    'PolylineTrailType',
    '_bottomHeight',
    'setInputAction',
    'fromRadians',
    'color',
    'click',
    'endColor',
    'MOUSE_MOVE',
    'sphere',
    '_appearance',
    '尾迹线',
    '_labelField',
    '_mergeOpts',
    'parseFromString',
    'imageSize',
    'span',
    'fromType',
    'EllipsoidOutlineGeometry',
    '_createParticle',
    '_event',
    'matrixMove',
    '基础底图图层',
    '_labelOpts',
    'maximumHeights',
    'topPsts',
    '_graphicDrawManager',
    '_loadMaterial',
    'uprightLabel',
    '_controlPositions',
    'LINES',
    '\n            #define HIGH_QUALITY\n            #ifdef HIGH_QUALITY\n            #define STEPS 130\n            #define ALPHA_WEIGHT 0.015\n            #define BASE_STEP 0.025\n            #else\n            #define STEPS 50\n            #define ALPHA_WEIGHT 0.05\n            #define BASE_STEP 0.1\n            #endif\n\n            #define time iTime\n            vec2 mo;\n            vec2 rot(in vec2 p, in float a){float c = cos(a), s = sin(a);return p*mat2(c,s,-s,c);}\n            float hash21(in vec2 n){ return fract(sin(dot(n, vec2(12.9898, 4.1414))) * 43758.5453); }\n            float noise(in vec3 p)\n            {\n            vec3 ip = floor(p), fp = fract(p);\n                fp = fp*fp*(3.0-2.0*fp);\n            vec2 tap = (ip.xy+vec2(37.0,17.0)*ip.z) + fp.xy;\n            vec2 cl = texture( iChannel0, (tap + 0.5)/256.0, 0.0 ).yx;\n            return mix(cl.x, cl.y, fp.z);\n            }\n\n            float fbm(in vec3 p, in float sr)\n            {\n                p *= 3.5;\n                float rz = 0., z = 1.;\n                for(int i=0;i<4;i++)\n                {\n                    float n = noise(p-time*.6);\n                    rz += (sin(n*4.4)-.45)*z;\n                    z *= .47;\n                    p *= 3.5;\n                }\n                return rz;\n            }\n\n            vec4 map(in vec3 p)\n            {\n                float dtp = dot(p,p);\n            p = .5*p/(dtp + .2);\n                p.xz = rot(p.xz, p.y*2.5);\n                p.xy = rot(p.xz, p.y*2.);\n                \n                float dtp2 = dot(p, p);\n                p = (mo.y + .6)*3.*p/(dtp2 - 5.);\n                float r = clamp(fbm(p, dtp*0.1)*1.5-dtp*(.35-sin(time*0.3)*0.15), 0. ,1.);\n                vec4 col = vec4(.5,1.7,.5,.96)*r;\n                \n                float grd = clamp((dtp+.7)*0.4,0.,1.);\n                col.b += grd*.6;\n                col.r -= grd*.5;    \n                vec3 lv = mix(p,vec3(0.3),2.);\n                grd = clamp((col.w - fbm(p+lv*.05,1.))*2., 0.01, 1.5 );\n                col.rgb *= vec3(.5, 0.4, .6)*grd + vec3(4.,0.,.4);\n                col.a *= clamp(dtp*2.-1.,0.,1.)*0.07+0.87;\n                \n                return col;\n            }\n\n            vec4 vmarch(in vec3 ro, in vec3 rd)\n            {\n            vec4 rz = vec4(0);\n            float t = 2.5;\n                t += 0.03*hash21(gl_FragCoord.xy);\n            for(int i=0; i<STEPS; i++)\n            {\n                if(rz.a > 0.99 || t > 6.)break;\n                vec3 pos = ro + t*rd;\n                    vec4 col = map(pos);\n                    float den = col.a;\n                    col.a *= ALPHA_WEIGHT;\n                col.rgb *= col.a*1.7;\n                rz += col*(1. - rz.a);\n                    t += BASE_STEP - den*(BASE_STEP-BASE_STEP*0.015);\n            }\n                return rz;\n            }\n\n            void main()\n            {\n            vec2 p = 2.*v_st.xy - vec2(1., 1.);\n            mo =vec2(.5,.5);\n            mo = (mo==vec2(.0))?mo=vec2(0.5,1.):mo;\n            \n            vec3 ro = 4.*normalize(vec3(cos(2.75-2.0*(mo.x+time*0.05)), sin(time*0.22)*0.2, sin(2.75-2.0*(mo.x+time*0.05))));\n            vec3 eye = normalize(vec3(0) - ro);\n            vec3 rgt = normalize(cross(vec3(0,1,0), eye));\n            vec3 up = cross(eye,rgt);\n            vec3 rd = normalize(p.x*rgt + p.y*up + (3.3-sin(time*0.3)*.7)*eye);\n            \n            vec4 col = clamp(vmarch(ro, rd),0.,1.);\n                col.rgb = pow(col.rgb, vec3(.9));\n                /*col.rb = rot(col.rg, 0.35);\n                col.gb = rot(col.gb, -0.1);*/\n                \n                fragColor = vec4(col.rgb * u_glow, 1.0);\n            }\n    ',
    'u_shadowMapTSDBANSS',
    'MipmapHint',
    'PrimitiveType',
    'ImageMaterialProperty',
    'shape',
    '_lockView',
    'flowLine',
    'fillText',
    'topline',
    'drivenPositions',
    'filletLabel',
    'isArray',
    'videoEle',
    'clipboard',
    'ellipse',
    '_setMatrix',
    'max',
    '_getProperties',
    '_graphicClassType',
    'lifetime',
    'initLeftDownEventHandler',
    '_availability',
    'svg格式错误！',
    'div',
    'lineVolume',
    'geometryInstances',
    'out_FragColor = vec4(out_FragColor.rgb , (out_FragColor.g + out_FragColor.r) *.3);if(out_FragColor.a<0.05)discard;',
    'onchange',
    'pickPosition',
    '相控阵雷达',
    'ColorGeometryInstanceAttribute',
    '\n                      void main(){\n                        iTime = czm_frameNumber / 60.0;\n                        iTime = u_speed * iTime;\n                        iResolution = czm_viewport.zw;\n                        customFsMain();  \n                       ',
    'offsetTop',
    'generateMipmap',
    'cartesian3',
    'CircleColorfulMaterialProperty',
    '_tickEventHandler',
    'parentElement',
    'uuid',
    '_seletedDom',
    '_arrowPrimitive',
    '_frameState',
    'reject',
    'shaderEffetRadarScan',
    'radii',
    'multiplyByVector',
    'roll',
    '_show',
    'VERTEX_FORMAT',
    'debugFrustum',
    'shaderEffetFlameCloud',
    'Ellipsoid',
    'toRadians',
    '_quaternion',
    '_postProcess',
    'replaceMain',
    'rgba(0,0,0,0)',
    'Color',
    'shaderEffetGlowStar',
    'xt3d-divgraphic-multiline-img_1',
    'WallScrollMaterialProperty',
    'graphicType',
    '\n#ifdef GL_ES\nprecision mediump float;\n#endif\n\n#extension GL_OES_standard_derivatives : enable \nuniform vec2 resolution;\nin vec2 v_st;\n\nvec3 hsv(float hue) {\n\treturn clamp(abs(fract(hue+vec3(0,2,1)/3.)*6.-3.)-1.,0.,1.);\n}\n\nvoid main(   ) {\n    float time = czm_frameNumber /60.0;\n    vec2  st= v_st * 2.0 - 1.0;\n\tvec2 uv = st;\n\tfloat d = length(uv) - 0.3;\n\tfloat bright = 0.001 / (d * d);\n    vec3 color = hsv(time * 001.2) + vec3(1, 1, 1);\n    float alpha=1.0; \n    vec3 c=color * bright; \n    float a=c.x+c.y+c.z;\n    alpha= a/4.0; \n\tout_FragColor = vec4(c, alpha);\n}',
    'sunDirectionWC',
    'HIGHLIGHT',
    '_circle',
    'fromPointNormal',
    'addOuterCylinder',
    '_name',
    'percent',
    '\n\n    uniform vec4 czm_pickColor; \n    void main() \n    { \n        czm_old_main(); \n        if (out_FragColor.a == 0.0) { \n        discard; \n        } \n        out_FragColor = czm_pickColor; \n    }',
    'waterAlpha',
    'Property',
    'circleScan_2',
    'shaderEffetMagicRing',
    'LEFT_DOWN',
    'shaderGeometryOpts',
    'executeInClosestFrustum',
    '374KtJkjN',
    '_setClusterEvent',
    'sources',
    '_initByImage',
    'xt3d-divgraphic-multiline-text',
    'mergeOpts',
    'fineArrow',
    '_measureResultEntity',
    'shaderEffetSwim',
    'toDataURL',
    '_outlineGeometry',
    '该对象已经被销毁，即调用了destroy()方法',
    'VaryingType',
    '_textureSize',
    '_createImage',
    ' 30px, ',
    '_videoTexture',
    'DISABLED',
    'divText',
    'REPLACE_MATERIAL',
    'u_bgColor',
    '_addRotateAnchor',
    '\nuniform vec4 color;\nuniform float speed;\nfloat circle(vec2 uv, float r, float blur) {\n  float d = length(uv) * 2.0;\n  float c = smoothstep(r+blur, r, d);\n  return c;\n}\nczm_material czm_getMaterial(czm_materialInput materialInput)\n{\n  czm_material material = czm_getDefaultMaterial(materialInput);\n  vec2 st = materialInput.st - 0.5;\n  material.diffuse =2.8 * color.rgb;\n  material.emission = vec3(0);\n  float t =fract(czm_frameNumber * speed / 1000.0);\n  float s = 0.3;\n  float radius1 = smoothstep(.0, s, t) * 0.5;\n  float alpha1 = circle(st, radius1, 0.01) * circle(st, radius1, -0.01);\n  float alpha2 = circle(st, radius1, 0.01 - radius1) * circle(st, radius1, 0.01);\n  float radius2 = 0.5 + smoothstep(s, 1.0, t) * 0.5;\n  float alpha3 = circle(st, radius1, radius2 + 0.01 - radius1) * circle(st, radius1, -0.01);\n  material.alpha = smoothstep(1.0, s, t) * (alpha1 + alpha2*0.1 + alpha3*0.1);\n  material.alpha *=color.a ;\n  return material;\n}  ',
    'heading',
    'px; width: ',
    '_initConsts',
    '_loadData',
    ',\n                roll: ',
    '_float',
    '事件类型无效！',
    'routeStopArrival',
    'setView',
    '_circleImage',
    'czm_pickColor',
    '_container_line',
    '直箭头',
    'maximumParticleLife',
    '平方米',
    'waterPrimitive',
    '_selectedEnable',
    'ceil',
    '_isPointLight',
    '/2F98C27F46D34E6D89E0A0E30B56F59E.png',
    'preRender',
    'getEvent',
    'WallLightType',
    '_segmentV',
    'expression',
    'GlowDiamondFS',
    'analyse',
    'fromGeometry',
    '_computeProperty',
    '_glow',
    '/C295257FB6934A52B51259670E5EB673.png',
    '_neckAngle',
    '_cache',
    '_imageUrl',
    '_farDepthFromNearPlusOne',
    '_time',
    'silhouetteColor',
    '_near',
    'positionWC',
    '); // 颜色\n                color *= vec4(vec3(position.',
    'processTexture',
    'vertexArray',
    'BOTTOM',
    'angle',
    'BoundingSphere',
    '_tailWidthFactor',
    '_entities',
    'fromHeadingPitchRoll',
    'UNIT_Z',
    '_shadowMapTexture',
    'distortionScale',
    '图标加载失败！',
    'count',
    'zoomToView',
    '_graphicOpts',
    'setPosition',
    'fromDimensions',
    '_hierarchy',
    'readFeatures',
    'px; background: #ff000042; border: 1px dashed red; position: absolute;  left: -5px; top: -5px;pointer-events:none;z-index:-1',
    'directionWC',
    '_provider',
    'clientX',
    'CircleScanMaterialProperty_4',
    '_originalWidth',
    'rgba(0,183,239,0.5)',
    '_generatePositions',
    '_proxy',
    '拖拽改变位置',
    'updatePosition',
    'toFixed',
    '按下鼠标左键确定位置',
    'rectSensor',
    'postRender',
    'cos',
    'updateEnvironment',
    'showThroughEllipsoid',
    'getSize',
    '\n#ifdef GL_ES\nprecision highp float;\n#endif\n\nin vec3 position;\nin vec2 st;\nin vec3 normal;\nuniform mat4 modelViewMatrix;\nuniform mat3 normalMatrix;\nuniform mat4 projectionMatrix;\nout vec3 v_position;\nout vec3 v_normal;\nout vec2 v_st;\n\nout vec3 v_light0Direction;\n\nvoid main(void)\n{\n    vec4 pos = czm_modelView * vec4(position, 1.0);\n    v_normal = normalMatrix * normal;\n    v_st = st;\n    v_position = pos.xyz;\n    v_light0Direction = mat3(modelViewMatrix) * vec3(1.0, 1.0, 1.0);\n    gl_Position = projectionMatrix * pos;\n} ',
    'background',
    'minimumSpeed',
    'CircleImageRotateMaterialProperty',
    'out_FragColor = vec4(out_FragColor.rgb + u_color, (out_FragColor.g + out_FragColor.r) *.4);',
    'modelView',
    'LEFT_CLICK',
    'Resource',
    '当前半径',
    'uniforms',
    'gravity',
    'textAlign',
    'maximumDistance',
    'glowEnable',
    'getMonth',
    'font',
    'near',
    'Lobster Two',
    'circleImageDiffuse',
    'zox',
    '_material',
    'customProjectionMatrix',
    'mouseup',
    'headLength',
    '\nuniform samplerCube u_shadowMapCube;\nuniform vec4 u_shadowMapTSDBANSS;\nuniform float u_shadowMapDarkness;\nuniform vec4 u_lightPositionEC;\n\nstruct shadowParameters {\n  vec3 texCoordsCube;\n  float depthBias;\n  float depth;\n  float nDotL;\n  vec2 texelStepSize;\n  float normalShadingSmooth;\n  float darkness;\n};\n\nfloat shadowVisibilityCube(samplerCube shadowMap, shadowParameters shadowParameters) {\n  float depthBias = shadowParameters.depthBias;\n  float depth = shadowParameters.depth;\n  float nDotL = shadowParameters.nDotL;\n  float normalShadingSmooth = shadowParameters.normalShadingSmooth;\n  float darkness = shadowParameters.darkness;\n  vec3 uvw = shadowParameters.texCoordsCube;\n\n  depth -= depthBias;\n  return czm_shadowDepthCompare(shadowMap, uvw, depth);\n}\n\nczm_material czm_getMaterial(czm_materialInput materialInput) {\n  czm_material material = czm_getDefaultMaterial(materialInput);\n\n  vec3 positionEC = materialInput.positionToEyeEC;\n  vec3 pointToLightEC = positionEC - u_lightPositionEC.xyz;\n  float pointToLightECLength = length(pointToLightEC);\n  vec3 l = normalize(pointToLightEC);\n\n  shadowParameters shadowParameters;\n  shadowParameters.texelStepSize = u_shadowMapTSDBANSS.xy;\n  shadowParameters.depthBias = u_shadowMapTSDBANSS.z;\n  shadowParameters.normalShadingSmooth = u_shadowMapTSDBANSS.w;\n  shadowParameters.darkness = u_shadowMapDarkness;\n  shadowParameters.depth = pointToLightECLength / u_radius;\n  shadowParameters.texCoordsCube = czm_inverseViewRotation * l;\n  float visibility = shadowVisibilityCube(u_shadowMapCube, shadowParameters);\n\n  if(visibility == 1.0) {\n    material.diffuse = u_visibleColor.rgb;\n    material.alpha = u_visibleColor.a;\n  } else {\n    material.diffuse = u_hiddenColor.rgb;\n    material.alpha = u_hiddenColor.a;\n  }\n  return material;\n}',
    'shaderEffetGradientRing',
    '#070acbbf',
    'fov',
    'className',
    'height: 10px; width: 10px; background: #01ff01; border-radius: 50%;border: 1px solid black;position:absolute;pointer-events: none;z-index:0',
    'moveDomLeft',
    '粗直箭头',
    'visibleColor',
    'RIGHT_CLICK',
    'MeasureAreaResult',
    '_projection',
    'mode',
    'LINE_STRIP',
    'TimeIntervalCollection',
    'DEFAULT',
    '螺旋材质',
    'horizontalPoistion',
    'fontColor',
    'minimumClock',
    'rgba(255, 255, 255, 0.5)',
    'xt3d-rectangle-label-text',
    'shaderEffetGlowDiamond',
    'mergeImage',
    'repeatDis',
    '_texture',
    '_maximumHeights',
    'graphic-draw-tip-container',
    'createPropertyDescriptor',
    '_createMp4VideoEle',
    '_isDestroyed',
    '_heatContainer',
    '_setContent',
    'floatingState',
    'upWC',
    '5zITefg',
    '1411377XUgKmg',
    'appearanceType',
    'LEFT_UP',
    'bold',
    '_line',
    '_createGraphic',
    'rect',
    'get_canvas',
    'needsUpdate',
    'sqrt',
    'primitives',
    '雷达四棱锥体',
    'ReflectionWater',
    'appendFromGeoJson',
    'options',
    'simpleLine',
    '当前测量类型：水平面积',
    'colorBlendAmount',
    'strokeStyle',
    'CylinderGeometry',
    'selectedGraphic',
    '_createCamera',
    'u_yHalfAngle',
    '_translation',
    '_domeFrontCommand',
    'green',
    'addGraphic',
    'fragmentShaderSource',
    '_shaderType',
    '_camera',
    'waterConservancy',
    'totalDistance',
    '燕尾进攻方向',
    'arrowAxis',
    'dimensions',
    'addMaterial',
    '_resetAnchorStyle',
    '_segmentH',
    'fourPposition',
    '_editAnchor',
    'wallScroll',
    '_createEntity',
    '_cartesian3Array',
    'normalShadingSmooth',
    'Entity',
    '_setSvgStyle',
    'UNSIGNED_BYTE',
    'PointGraphic',
    'divIndicator',
    '_imageryProvider',
    'pointLightRadius',
    '_layerOpts',
    '.arrow',
    'SphereScanMaterialProperty',
    'darkness',
    'camera',
    'container',
    'IDENTITY',
    'tailedAttackArrow',
    'CircleWaveMaterialProperty',
    'svg',
    '_currentHeight',
    '_geometry',
    'defineProperty',
    'setRotate',
    '_getCameraInfo',
    'FRONT',
    '_editType',
    'POSITION_NORMAL_AND_ST',
    '绘制图元名称：',
    '_sampler',
    '_createAppearence',
    '_createGraphicEditor',
    '_videoPlayer',
    '_getMaxHeight',
    'headRadius',
    'rectangle',
    '_computeShape',
    '_createMaterial',
    'fourPindices',
    'rgb(255,255,255)',
    '_bottomCircle',
    '_setComponentVisible',
    '<div class="divpoint divpoint-theme">\n      <div class="divpoint-wrap">\n        <div class="area">\n          <div class="arrow-lt"></div>\n          <div class="b-t"></div>\n          <div class="b-r"></div>\n          <div class="b-b"></div>\n          <div class="b-l"></div>\n          <div class="arrow-rb"></div>\n          <div class="label-wrap">\n            <div class="title">',
    'fromTranslation',
    'pointLight',
    'sphereScan',
    'emitterModelMatrix',
    'PolygonGraphic',
    '_removeSvgEditAnchors',
    '_getPositoinByAngle',
    'floor',
    'createVertexBuffer',
    'min',
    '_contentDom',
    'CircleImageRotateType',
    '发光面板',
    'shouldAnimate',
    'images',
    'BoxGeometry',
    '_type',
    'BufferUsage',
    '_removeGraphicEditor',
    'endFovV',
    '108552ZWEqMK',
    '_moveHeightPositions',
    '\n     in vec3 v_positionEC;\n     uniform float globalAlpha;\n     \n     void main() {\n       czm_materialInput materialInput;\n       materialInput.positionToEyeEC = v_positionEC;\n       czm_material material = czm_getMaterial(materialInput);\n       out_FragColor = vec4(material.diffuse + material.emission, material.alpha * globalAlpha);\n     }\n     ',
    '\n#extension GL_OES_standard_derivatives : enable\n\nprecision highp float;\n\nuniform float time;\nuniform vec2 mouse;\nuniform vec2 resolution; \nuniform vec3 color; \nin vec2 v_st;\nvoid main( void ) {\n\n    vec2 p = v_st;\n    vec2 q = p - vec2(0.5, 0.5);\n    float time = czm_frameNumber / 100.0;\n    // Time varying pixel color\n    vec3 col = 0.5 + 0.5 * sin(vec3(0,2,4) + time / 5.);\n    \n    float r = 0.4 + 0.1 * cos(atan(-p.x + 2., q.y / 2.) * 60.0 + 20.0 * -q.x );\n\n    float r1 = 0.2 + 0.1 * cos(atan(-p.x + 2., -p.y / 2.) * 5.0 + 20.0 * -q.x + time * 2.);\n    float r2 = 0.4 + 0.1 * cos(atan(-p.x + 4., -q.y / 2.) * 10.0 + 20.0 * -q.x + time * 4.);\n    float r3 = 0.6 + 0.1 * cos(atan(-p.x + 6., -q.y / 2.) * 20.0 + 20.0 * -q.x + time * 6.);\n    \n    col /= smoothstep(r / 100., r, length(q)) / .2;\n\n    col.x /= smoothstep(r1 / 100., r1, length(q)) / .2;\n    col.y /= smoothstep(r2 / 100., r2, length(q)) / .2;\n    col.z /= smoothstep(r3 / 100., r3, length(q)) / .2;\n\n    out_FragColor = vec4(col,0.3);\n    if((col.x+col.y+col.z)<0.8){\n        out_FragColor.a=0.01;\n        // discard\n    };\n    vec3 rgb= vec3(0.0,0.0,0.0);\n    out_FragColor = vec4(out_FragColor.rgb +rgb, out_FragColor.a);\n}\n',
    '_updateGraphic',
    '_view',
    'lerp',
    'fromDate',
    'rgba(255, 255, 255, 0.9)',
    '_centerWorldPosition',
    '#10eedc',
    '_createDrawCommands',
    '_defaultView',
    'extend2Earth',
    'endFovH',
    'zReverse',
    '当前绘制类型：面，至少需要3个点',
    'ImageryLayer',
    'glowCone',
    '_initUniforms',
    '_onVisibilityChange',
    'domEle',
    '_originalHeight',
    '_on',
    '_getMatrix',
    'PRE_MULTIPLIED_ALPHA_BLEND',
    'moveDomTop',
    '_tempEntity',
    'CustomShader',
    'Appearance',
    '图元图层',
    '个点的位置',
    'Quaternion',
    'depthTest',
    'getValue',
    '#006ab4',
    '_materialType',
    ' \nczm_material czm_getMaterial( czm_materialInput cmi )\n{\n    czm_material material = czm_getDefaultMaterial(cmi);\n    vec2 st = cmi.st;\n    float t = fract(czm_frameNumber / 1000.0) * 1.;\n    vec2 st1 = vec2(fract(st.s - t),st.t); \n    float alpha =  1.- st.t;\n     \n    vec4 color_ = vec4(color.rgb * color.a, alpha * 1.2);\n    material.diffuse =color_.rgb * 2.5;\n    material.alpha = color_.a;\n\n    if(breathe){\n       float powerRatio = fract(czm_frameNumber / 30.0) + 1.0; \n       material.alpha  = pow(1.0 -  st.t, powerRatio);\n    } \n    return material;\n} ',
    '_remarks',
    '\nuniform vec4 color;\nuniform float speed; \nczm_material czm_getMaterial(czm_materialInput materialInput){\n  czm_material material = czm_getDefaultMaterial(materialInput);\n  vec2 st = materialInput.st * 2.0 - 1.0;\n  float t = czm_frameNumber * speed / 1000.0 ;\n  vec3 col = vec3(0.0);\n  vec2 p = vec2(sin(t), cos(t));\n  float d = length(st - dot(p, st) * p);\n  if (dot(st, p) < 0.) {\n    d = length(st);\n  }\n\n  col = .006 / d * color.rgb;\n\n  if(distance(st,vec2(0)) >  0.99 ){\n    col =color.rgb;\n  }\n\n  material.alpha  = pow(length(col),2.0);\n  material.diffuse = col * 3.0 ;\n  return material;\n} ',
    '20px 微软雅黑',
    'clientWidth',
    'show',
    '_particleGraphics',
    'PerspectiveFrustum',
    '_drawOpts',
    '_setPosition',
    '_uri',
    'gamma',
    'attachMediaElement',
    '_modelViewProjectionDirty',
    'horizontalOrigin',
    '_overlayCommandList',
    'canvas',
    '_createFramebuffer',
    'px Microsoft YaHei',
    '_isSelected',
    '_clampToGround',
    '_gravityScratch',
    'moveForward',
    'concat',
    '_preUpdate',
    '_inputP',
    'getMinutes',
    '_createTexture',
    'bind',
    '_globeShow',
    'xt3d-fillet-label-line',
    'boardGraphic',
    'emissionRate',
    '_createEntities',
    'preventDefault',
    'runAnimations',
    'showBottomCircle',
    '_viewProjectionDirty',
    '_destroyResource',
    '_height',
    '，需要',
    'routeStart',
    'activate',
    'setContent',
    'strokeWidth',
    '开敞度分析图元',
    'glowPower',
    'VideoGraphic',
    'equals',
    '_shaderGeometryType',
    '\n#define NUM_RAYS 13.\n \n#define VOLUMETRIC_STEPS 19\n\n#define MAX_ITER 35\n#define FAR 6.\n\n#define time iTime*1.1\n\n\nmat2 mm2(in float a){float c = cos(a), s = sin(a);return mat2(c,-s,s,c);}\nfloat noise( in float x ){return texture(iChannel0, vec2(x*.01,1.),0.0).x;}\n\nfloat hash( float n ){return fract(sin(n)*43758.5453);}\n\nfloat noise(in vec3 p)\n{\n    vec3 ip = floor(p);\n    vec3 fp = fract(p);\n    fp = fp*fp*(3.0-2.0*fp);\n    \n    vec2 tap = (ip.xy+vec2(37.0,17.0)*ip.z) + fp.xy;\n    vec2 rg = texture( iChannel0, (tap + 0.5)/256.0, 0.0 ).yx;\n    return mix(rg.x, rg.y, fp.z);\n}\n\nmat3 m3 = mat3( 0.00,  0.80,  0.60,\n                -0.80,  0.36, -0.48,\n                -0.60, -0.48,  0.64 );\n\n\n//See: https://www.shadertoy.com/view/XdfXRj\nfloat flow(in vec3 p, in float t)\n{\n    float z=2.;\n    float rz = 0.;\n    vec3 bp = p;\n    for (float i= 1.;i < 5.;i++ )\n    {\n    p += time*.1;\n    rz+= (sin(noise(p+t*0.8)*6.)*0.5+0.5) /z;\n    p = mix(bp,p,0.6);\n    z *= 2.;\n    p *= 2.01;\n        p*= m3;\n    }\n    return rz;\t\n}\n\n//could be improved\nfloat sins(in float x)\n{\n    float rz = 0.;\n    float z = 2.;\n    for (float i= 0.;i < 3.;i++ )\n    {\n        rz += abs(fract(x*1.4)-0.5)/z;\n        x *= 1.3;\n        z *= 1.15;\n        x -= time*.65*z;\n    }\n    return rz;\n}\n\nfloat segm( vec3 p, vec3 a, vec3 b)\n{\n    vec3 pa = p - a;\n    vec3 ba = b - a;\n    float h = clamp( dot(pa,ba)/dot(ba,ba), 0.0, 1. );\t\n    return length( pa - ba*h )*.5;\n}\n\nvec3 path(in float i, in float d)\n{\n    vec3 en = vec3(0.,0.,1.);\n    float sns2 = sins(d+i*0.5)*0.22;\n    float sns = sins(d+i*.6)*0.21;\n    en.xz *= mm2((hash(i*10.569)-.5)*6.2+sns2);\n    en.xy *= mm2((hash(i*4.732)-.5)*6.2+sns);\n    return en;\n}\n\nvec2 map(vec3 p, float i)\n{\n    float lp = length(p);\n    vec3 bg = vec3(0.);   \n    vec3 en = path(i,lp);\n    \n    float ins = smoothstep(0.11,.46,lp);\n    float outs = .15+smoothstep(.0,.15,abs(lp-1.));\n    p *= ins*outs;\n    float id = ins*outs;\n    \n    float rz = segm(p, bg, en)-0.011;\n    return vec2(rz,id);\n}\n\nfloat march(in vec3 ro, in vec3 rd, in float startf, in float maxd, in float j)\n{\n    float precis = 0.001;\n    float h=0.5;\n    float d = startf;\n    for( int i=0; i<MAX_ITER; i++ )\n    {\n        if( abs(h)<precis||d>maxd ) break;\n        d += h*1.2;\n        float res = map(ro+rd*d, j).x;\n        h = res;\n    }\n    return d;\n}\n\n//volumetric marching\nvec3 vmarch(in vec3 ro, in vec3 rd, in float j, in vec3 orig)\n{   \n    vec3 p = ro;\n    vec2 r = vec2(0.);\n    vec3 sum = vec3(0);\n    float w = 0.;\n    for( int i=0; i<VOLUMETRIC_STEPS; i++ )\n    {\n        r = map(p,j);\n        p += rd*.03;\n        float lp = length(p);\n        \n        vec3 col = sin(vec3(1.05,2.5,1.52)*3.94+r.y)*.85+0.4;\n        col.rgb *= smoothstep(.0,.015,-r.x);\n        col *= smoothstep(0.04,.2,abs(lp-1.1));\n        col *= smoothstep(0.1,.34,lp);\n        sum += abs(col)*5. * (1.2-noise(lp*2.+j*13.+time*5.)*1.1) / (log(distance(p,orig)-2.)+.75);\n    }\n    return sum;\n}\n\n//returns both collision dists of unit sphere\nvec2 iSphere2(in vec3 ro, in vec3 rd)\n{\n    vec3 oc = ro;\n    float b = dot(oc, rd);\n    float c = dot(oc,oc) - 1.;\n    float h = b*b - c;\n    if(h <0.0) return vec2(-1.);\n    else return vec2((-b - sqrt(h)), (-b + sqrt(h)));\n}\n\nvoid main()\n{\t\n    // vec2 p = fragCoord.xy/iResolution.xy-0.5;\n    // p.x*=iResolution.x/iResolution.y;\n    vec2 p = 2.*v_st.xy - vec2(1., 1.);\n    vec2 um = vec2(.5,.5);\n    \n    //camera\n    vec3 ro = vec3(0.,0.,5.);\n    vec3 rd = normalize(vec3(p*.7,-1.5));\n    mat2 mx = mm2(time*.4+um.x*6.);\n    mat2 my = mm2(time*0.3+um.y*6.); \n    ro.xz *= mx;rd.xz *= mx;\n    ro.xy *= my;rd.xy *= my;\n    \n    vec3 bro = ro;\n    vec3 brd = rd;\n    \n    vec3 col = vec3(0.0125,0.,0.025);\n    #if 1\n    for (float j = 1.;j<NUM_RAYS+1.;j++)\n    {\n        ro = bro;\n        rd = brd;\n        mat2 mm = mm2((time*0.1+((j+1.)*5.1))*j*0.25);\n        ro.xy *= mm;rd.xy *= mm;\n        ro.xz *= mm;rd.xz *= mm;\n        float rz = march(ro,rd,2.5,FAR,j);\n    if ( rz >= FAR)continue;\n        vec3 pos = ro+rz*rd;\n        col = max(col,vmarch(pos,rd,j, bro));\n    }\n    #endif\n    \n    ro = bro;\n    rd = brd;\n    vec2 sph = iSphere2(ro,rd);\n    \n    if (sph.x > 0.)\n    {\n        vec3 pos = ro+rd*sph.x;\n        vec3 pos2 = ro+rd*sph.y;\n        vec3 rf = reflect( rd, pos );\n        vec3 rf2 = reflect( rd, pos2 );\n        float nz = (-log(abs(flow(rf*1.2,time)-.01)));\n        float nz2 = (-log(abs(flow(rf2*1.2,-time)-.01)));\n        col += (0.1*nz*nz* vec3(0.12,0.12,.5) + 0.05*nz2*nz2*vec3(0.55,0.2,.55))*0.8;\n    }\n    \n    fragColor = vec4(col*u_glow, 1.0);\n}',
    '点击鼠标左键结束，点击鼠标右键取消',
    'glow',
    '_createAppearance',
    'MaterialType',
    'xt3d-divgraphic-show',
    '_visibilityChangeEvent',
    'shaderEffetMagicBox',
    'WallFlowType',
    'routeGraphic',
    'CustomColorType',
    'tileset',
    '_createImageryProvider',
    'toStringTag',
    '_scanePlaneXHalfAngle',
    'point',
    '_initBySvg',
    'rgba(255, 0, 0, 0.4)',
    '_addHook',
    'createEvent',
    '墙体发光材质',
    'attributes',
    'readFeature',
    '_synToWindowPosition',
    'enableRotate',
    '_zReverse',
    'shaderEffetGlowCross',
    'poiMarker',
    'clear',
    'conditions',
    '_virtualCamera',
    'CLAMP_TO_EDGE',
    'VelocityOrientationProperty',
    '_shadowMapCamera',
    'ScreenSpaceEventHandler',
    '拖拽缩放',
    'REPEAT',
    'StellarChainFS',
    'ellipsoid',
    'getTime',
    'merge',
    '_overWrite',
    'getTooltip',
    '_id',
    'GlowCrossFS',
    'vertexFormat',
    'startColor',
    'minimumCone',
    'ConeEmitter',
    'WebGLConstants',
    'xt3d材质',
    '扫描材质_4',
    'defaultFS',
    '_plane',
    'shaderEffetRealFlame',
    'gifMarker',
    '_setChartOption',
    '_computeModelMatrix',
    '_computedScanPlaneModelMatrix',
    '椎体雷达',
    'czm_material czm_getMaterial(czm_materialInput materialInput)\n {\n      czm_material material = czm_getDefaultMaterial(materialInput);\n      vec2 st = materialInput.st;\n      vec4 colorImage = texture(image, vec2(fract(st.t - time), st.t));\n      material.alpha = colorImage.a * color.a;\n      material.diffuse =  1.9 * color.rgb  ;\n      return material;\n  }',
    'jpe',
    'geoJsonData',
    '\n            precision highp float;\n            precision highp int;\n            uniform float u_speed;\n            uniform vec3 u_color;\n            uniform float u_time;\n            uniform float u_glow;\n            in vec3 v_positionEC;\n            in vec3 v_normalEC;\n            in vec2 v_st;\n            mat4 mat  = mat4 ( vec4 ( 1.0 , 0.0 , 0.0 , 0.0 ),\n                vec4 ( 0.0 , 1.0 , 0.0 , 0.0 ),\n                vec4 ( 0.0 , 0.0 , 1.0 , 0.0 ),\n                vec4 ( 0.0 , 0.0 , 0.0 , 1.0 ) );\n\n            vec2 pos;\n\n            vec4 col = vec4 ( .3, .1, 0.1, 1000.0 );\n\n            void Line2 ( vec2 a, vec2 b );\n            void Line2 ( vec2 a, vec2 b ) {\n            float d = distance ( pos , a ) + distance ( pos , b ) - distance ( a , b ) + 1e-5;\n            col += max ( 1. - pow ( d * 14. , 0.1 ) , -0.01 );\n            }\n\n            void Line4 ( vec4 a, vec4 b );\n            void Line4 ( vec4 a, vec4 b ) {\n            a = mat * a;\n            a.xyz /= 1.5 + a.w * 2.;\n            b = mat * b;\n            b.xyz /= 1.5 + b.w * 2.;\n            Line2 ( a.xy , b.xy );\n            }\n\n            void Point ( vec4 p );\n            void Point ( vec4 p ) {\n            p = mat * p;\n            p.xyz /= 1.5 + p.w * 2.;\n            \n            float d = distance ( pos , p.xy );\n            \n            if ( d < .3 )\n            if ( p.z < col.a ) {\n                col.b += max ( 1.0 - pow ( d * 5.0 , .1 ) , 0.0 );\n            }\n            }\n\n            void Rotate ( float angle, float d1, float d2, float d3, float d4);\n            void Rotate ( float angle, float d1, float d2, float d3, float d4) {\n            float c = cos (angle), s = sin (angle);\n            mat *= mat4 ( vec4 (  c*d1+(1.-d1),  s * d2 * d1 , -s * d3 * d1 ,  s * d4 * d1 ),\n                    vec4 ( -s * d1 * d2 ,  c*d2+(1.-d2),  s * d3 * d2 , -s * d4 * d2 ),\n                    vec4 (  s * d1 * d3 , -s * d2 * d3 ,  c*d3+(1.-d3),  s * d4 * d3 ),\n                    vec4 ( -s * d1 * d4 ,  s * d2 * d4 , -s * d3 * d4 ,  c*d4+(1.-d4)) );\n            }\n\n            void main( void ) {\n            float time = czm_frameNumber / 60.0;\n            time = u_speed * time;\n            pos = v_st - 0.5;\n            \n            Rotate ( time,      0.0, 1.0, 1.0, 0.0 );\n            Rotate ( time * .5, 1.0, 0.0, 1.0, 0.0 );\n            Rotate ( time * .3, 1.0, 1.0, 0.0, 0.0 );\n            \n            Line4 ( vec4 ( .2, .2, .2, .2 ), vec4 (-.2, .2, .2, .2 ) );\n            Line4 ( vec4 ( .2, .2, .2, .2 ), vec4 ( .2,-.2, .2, .2 ) );\n            Line4 ( vec4 ( .2, .2, .2, .2 ), vec4 ( .2, .2,-.2, .2 ) );\n            Line4 ( vec4 ( .2, .2, .2, .2 ), vec4 ( .2, .2, .2,-.2 ) );\n            \n            Line4 ( vec4 ( .2, .2, .2,-.2 ), vec4 (-.2, .2, .2,-.2 ) );\n            Line4 ( vec4 ( .2, .2, .2,-.2 ), vec4 ( .2,-.2, .2,-.2 ) );\n            Line4 ( vec4 ( .2, .2, .2,-.2 ), vec4 ( .2, .2,-.2,-.2 ) );\n            \n            Line4 ( vec4 ( .2, .2,-.2, .2 ), vec4 (-.2, .2,-.2, .2 ) );\n            Line4 ( vec4 ( .2, .2,-.2, .2 ), vec4 ( .2,-.2,-.2, .2 ) );\n            \n            Line4 ( vec4 ( .2, .2,-.2,-.2 ), vec4 (-.2, .2,-.2,-.2 ) );\n            Line4 ( vec4 ( .2, .2,-.2,-.2 ), vec4 ( .2,-.2,-.2,-.2 ) );\n            Line4 ( vec4 ( .2, .2,-.2,-.2 ), vec4 ( .2, .2,-.2, .2 ) );\n            \n            Line4 ( vec4 ( .2,-.2, .2, .2 ), vec4 (-.2,-.2, .2, .2 ) );\n            Line4 ( vec4 ( .2,-.2, .2, .2 ), vec4 ( .2,-.2,-.2, .2 ) );\n            Line4 ( vec4 ( .2,-.2, .2, .2 ), vec4 ( .2,-.2, .2,-.2 ) );\n            \n            Line4 ( vec4 ( .2,-.2, .2,-.2 ), vec4 (-.2,-.2, .2,-.2 ) );\n            Line4 ( vec4 ( .2,-.2, .2,-.2 ), vec4 ( .2,-.2,-.2,-.2 ) );\n            \n            Line4 ( vec4 ( .2,-.2,-.2, .2 ), vec4 (-.2,-.2,-.2, .2 ) );\n            Line4 ( vec4 ( .2,-.2,-.2, .2 ), vec4 ( .2,-.2,-.2,-.2 ) );\n            \n            Line4 ( vec4 ( .2,-.2,-.2,-.2 ), vec4 (-.2,-.2,-.2,-.2 ) );\n            \n            \n            Line4 ( vec4 (-.2, .2, .2, .2 ), vec4 (-.2,-.2, .2, .2 ) );\n            Line4 ( vec4 (-.2, .2, .2, .2 ), vec4 (-.2, .2,-.2, .2 ) );\n            Line4 ( vec4 (-.2, .2, .2, .2 ), vec4 (-.2, .2, .2,-.2 ) );\n            \n            Line4 ( vec4 (-.2, .2, .2,-.2 ), vec4 (-.2,-.2, .2,-.2 ) );\n            Line4 ( vec4 (-.2, .2, .2,-.2 ), vec4 (-.2, .2,-.2,-.2 ) );\n            \n            Line4 ( vec4 (-.2, .2,-.2, .2 ), vec4 (-.2,-.2,-.2, .2 ) );\n            \n            Line4 ( vec4 (-.2, .2,-.2,-.2 ), vec4 (-.2,-.2,-.2,-.2 ) );\n            Line4 ( vec4 (-.2, .2,-.2,-.2 ), vec4 (-.2, .2,-.2, .2 ) );\n            \n            Line4 ( vec4 (-.2,-.2, .2, .2 ), vec4 (-.2,-.2,-.2, .2 ) );\n            Line4 ( vec4 (-.2,-.2, .2, .2 ), vec4 (-.2,-.2, .2,-.2 ) );\n            \n            Line4 ( vec4 (-.2,-.2, .2,-.2 ), vec4 (-.2,-.2,-.2,-.2 ) );\n            \n            Line4 ( vec4 (-.2,-.2,-.2, .2 ), vec4 (-.2,-.2,-.2,-.2 ) );\n            \n            Point ( vec4 ( .2, .2, .2, .2 ) );\n            Point ( vec4 ( .2, .2, .2,-.2 ) );\n            Point ( vec4 ( .2, .2,-.2, .2 ) );\n            Point ( vec4 ( .2, .2,-.2,-.2 ) );\n            Point ( vec4 ( .2,-.2, .2, .2 ) );\n            Point ( vec4 ( .2,-.2, .2,-.2 ) );\n            Point ( vec4 ( .2,-.2,-.2, .2 ) );\n            Point ( vec4 ( .2,-.2,-.2,-.2 ) );\n            \n            Point ( vec4 (-.2, .2, .2, .2 ) );\n            Point ( vec4 (-.2, .2, .2,-.2 ) );\n            Point ( vec4 (-.2, .2,-.2, .2 ) );\n            Point ( vec4 (-.2, .2,-.2,-.2 ) );\n            Point ( vec4 (-.2,-.2, .2, .2 ) );\n            Point ( vec4 (-.2,-.2, .2,-.2 ) );\n            Point ( vec4 (-.2,-.2,-.2, .2 ) );\n            Point ( vec4 (-.2,-.2,-.2,-.2 ) );\n            \n            out_FragColor = vec4( col.xyz * u_glow, 1.0 );\n            out_FragColor = vec4(out_FragColor.rgb + u_color, out_FragColor.r* out_FragColor.r);\n            }\n    ',
    'serializeToString',
    'axis',
    '#080707',
    '_createVideoPlane',
    '_createLabel',
    'uprightLine',
    'ClassificationType',
    'polygonOpts',
    '#00FFFF',
    'isSpotLight',
    'setTargetPosition',
    'forEach',
    '波束雷达',
    'u_intersectionColor',
    '_mouseDown',
    'NICEST',
    'PlaneGeometry',
    '_pickPosition',
    'circleSpiral',
    'lookAtTransform',
    'boundingSphere',
    '墙体图片流动材质',
    '_dimensions',
    '按下鼠标右键取消绘制',
    'FrustumGeometry',
    '_scanePlaneYHalfAngle',
    'maskImage',
    'GlowPointFS',
    '_sectorBackCommand',
    'MeasureDistanceVertex',
    'appendChild',
    'layer',
    '.b-l',
    'createLinearGradient',
    'rangeRadar',
    'red',
    'defaultTexture',
    'version',
    'shaderEffetDoubleLoop',
    '_param',
    'czm_material czm_getMaterial(czm_materialInput materialInput){\n    czm_material material = czm_getDefaultMaterial(materialInput);\n    vec4 tColor = u_color;\n    vec2 st = materialInput.st;\n    vec2 center = st - vec2(0.5,0.5);\n    float length = length(center)/0.5;\n    float time = 1. - abs(czm_frameNumber / 360. - 0.5);\n\n    float param = 1. - step(length, 0.6);//大于0.6模糊，rate = 0.6\n    float scale = param * length;// 0.6< length 返回0，反之返回1.\n    float alpha = param * (1.0 - abs(scale - 0.8) / 0.2);// 0.8 < length 返回0，反之返回1.\n\n    float param1 = step(length, 0.7);//小于0.5模糊\n    float scale1 = param1 * length;// 0.6< length 返回0，反之返回1.\n    alpha += param1 * (1.0 - abs(scale1 - 0.35) / 0.35);// 0.8 < length 返回0，反之返回1.\n\n    material.diffuse = u_color.rgb * vec3(u_color.a);\n    material.alpha = pow(alpha, 4.0);\n    return material;\n}',
    'Polygon',
    '点击鼠标左键新增点，点击鼠标右键取消',
    '_labelEntity',
    'KeyboardEventModifier',
    '_modelMatrix',
    'Material',
    'removeInputAction',
    'customShader',
    'type',
    'Texture',
    '_removeComponent',
    'buildModuleUrl',
    'flyToBoundingSphere',
    'updateAndExecuteCommands',
    'pick',
    'round',
    'eastNorthUpToFixedFrame',
    '_getArrowBodyPoints',
    'cylinderRadar',
    '\n#ifdef GL_ES\nprecision mediump float;\n#endif\n\n#extension GL_OES_standard_derivatives : enable\nin vec2 v_st;\n  float time;\n  vec2 mouse;\n  vec2 resolution;\n\nconst int num_balls = 2000; // change this, but be careful \nconst float coordinate_scale = 1000.;\nconst float structure_size = 1.20; // Size from 0 to 1, max size = 1\nconst float glow_decay = 2.00; \nconst float trail_len = 151.0;\n#define pi 3.14159265358\nconst float speed = 0.189;\nconst float rot_speed = speed*10.;\nconst float starting_pt = 2.0; // This is a good number\n\nvec4 draw_ball(float i, float j, float size) {\n\tfloat balls = float(num_balls);\n\tfloat dt = starting_pt + time * speed;\n\t// Map coordinates to 0-1\n\tvec2 coord = v_st *  1.2 - .1;//gl_FragCoord.xy/resolution.xy;\n\t//map coordinates to -coord_scale -> +coord_scale\n\tcoord = coord *coordinate_scale-coordinate_scale/2.;\n\tcoord -= vec2(coord.x/2.,coord.y/2.);\n\t\n\t//Controls motion of balls\n\tfloat spacing = (2.*pi)/balls;\n\tfloat xi = sin(dt + i * spacing) * structure_size;\n\tfloat yi = cos(dt + i *spacing) * structure_size;\n\tfloat x =  (sin(dt*i*spacing)*100. - cos(dt*j*spacing)) * structure_size;\n\tfloat y =  (cos(dt*j*spacing)*100. - sin(dt*i*spacing)) * structure_size;\n\ty *= ((xi- dt)/-dt) + sin(i*spacing + dt*sin(dt/100.));\n\tx *= ((yi - dt)/dt)  + sin(i*spacing - dt*cos(dt/100.));\n\t//Correct aspect ratio\n\tcoord.x *= resolution.x/resolution.y; \n\tvec2 pos = vec2(x,y);\n\tmat2 rot = mat2(cos(dt*rot_speed), -sin(dt*rot_speed), sin(dt*rot_speed), cos(dt*rot_speed));\n\tpos *= rot;\n\tfloat dist = length(coord - pos);\n\t\n\t//Controls how quickly brightness falls off\n\tfloat intensity = pow(size/dist, glow_decay);\n\t\n\tvec4 color = vec4(vec3(1.0) * abs(sin(vec3(time*0.25,time/2.,time/3.))), 1.0);\n\treturn color * intensity;\n}\n\n\nvoid main( void ) {\n\n    mouse=vec2(1.0,1.0);\n    time = czm_frameNumber /100.0;\n    resolution = czm_viewport.zw;\n\n\tvec4 col = vec4(0.0);\n\tfor (int i = 0; i < num_balls; i++) {\n\t\tvec2 pt = vec2(float(i),float(i));\n\t\tcol += draw_ball(float(i),float(i), 2.1-distance(pt,vec2(0.))/coordinate_scale);\n\t}\n   \n  float alpha=(col.x+ col.y + col.z) ; \n  if(alpha < 0.2){ \n      discard; \n  } \n  out_FragColor = col;\n}',
    'dataSources',
    'default',
    '_outlineColor',
    'rgba(255,0,0,1)',
    '_measureType',
    'direction',
    'depthBias',
    'getValueOrClonedDefault',
    'texture',
    'RED',
    'image/octet-stream',
    'maximumCone',
    'MAX_VALUE',
    'JulianDate',
    'querySelector',
    'rgba(255, 242, 0, 0.34)',
    'tailedSquadCombat',
    'call',
    '\n\n#ifdef GL_ES\nprecision mediump float;\n#endif\n\nuniform float time; \nin vec2 v_st;\n\n#define PI 3.14159265358979\n#define N 12\nvoid main( void ) {\n    float time = czm_frameNumber /60.0;\n\tfloat size = 0.30;\n\tfloat dist = 0.120;\n\tfloat ang = 0.0;\n\tvec2 pos = vec2(0.0,0.0);\n\tvec3 color = vec3(0.1);\n    vec2 st=v_st * 2.0 - 1.0;\n    \n\tfor(int i=0; i<N; i++){\n\t\tfloat r = 0.3;\n\t\tang += PI / (float(N)*0.1)+time/120.;\n\t\tpos = vec2(cos(ang),sin(ang))*r*sin(time+ang/.6);\t\t\t\t  \n\t\tdist += size / distance(pos,st);\n\t\tvec3 c = vec3(0.03, 0.05, .1);\n\t\tcolor = c*dist;\n    }\n    \n    float alpha=1.0; \n    vec3 c=color;\n    float a=c.x+c.y+c.z;\n    alpha= c.r/3. ; //a; \n    if(a<1.5)discard;\n\tout_FragColor = vec4(c*vec3(1,1,1), alpha);\n}',
    'GeometryPipeline',
    'resolve',
    'Assets/Textures/waterNormals.jpg',
    'initMouseEvent',
    '_removeEditAnchor',
    'subtract',
    'wallFlow',
    '#ff4500',
    'style',
    'RectPyramid',
    'out_FragColor = vec4(out_FragColor.rgb + u_color, out_FragColor.r + out_FragColor.g + out_FragColor.b);',
    'editAnchorType',
    'createAttributeLocations',
    '笛卡尔',
    '_uniforms',
    '_addHWAnchors',
    'xt3d-divgraphic-brighten',
    'timeInfo',
    'renderState',
    '_convertPositions',
    'CircleScanType_2',
    '_chartDom',
    'fromRotationY',
    'shadowMap',
    '_scanPlaneVA',
    'imageSubRegion',
    '_createDrawTool',
    '_leftClickEvent',
    '_debugFrustum',
    'series',
    'onload',
    '光柱体',
    'scissorTest',
    'udpateHeightPosition',
    '.area',
    '_createPolygon',
    '\n                        uniform vec4 color;\n                        uniform float repeat;\n                        uniform float offset;\n                        uniform float thickness;\n                        czm_material czm_getMaterial(czm_materialInput materialInput)\n                        { \n                            czm_material material = czm_getDefaultMaterial(materialInput);\n                            material.diffuse = 1.5 * color.rgb;\n                            vec2 st = materialInput.st;\n                            float dis = distance(st,vec2(0.5,0.5));\n                            float per = fract(czm_frameNumber/100.);\n                            if(dis > per * 0.5){\n                                discard;\n                            }else {\n                                material.alpha = color.a * dis / per / 2.0;\n                            }\n                            return material;\n                        }\n                    ',
    'sign',
    'left',
    'now',
    'input',
    '_beforDestroy',
    'materialOpts',
    'info-input',
    'out_FragColor = vec4(out_FragColor.rgb, out_FragColor.r + out_FragColor.r);',
    'CircleColorfulType',
    'accept',
    'SphereGeometry',
    'init',
    'uniformMap',
    '_renderStateOptions',
    'BOTH',
    'emission',
    'addInnereCylinder',
    'mainFS',
    '_inverseProjectionDirty',
    'prototype',
    '箭头轴图元',
    '_createAnimateEntity',
    '_stopsTimes',
    'freeze',
    'CircleImageDiffuseMaterialProperty',
    '_reAdd',
    'DivGraphic',
    'beginPath',
    '_dom',
    '_setCamera',
    'ZERO',
    'lineCommand',
    '_fetchSvg',
    '燕尾分队战斗',
    '\nfloat noise(vec3 p) //Thx to Las^Mercury\n{\nvec3 i = floor(p);\nvec4 a = dot(i, vec3(1., 57., 21.)) + vec4(0., 57., 21., 78.);\nvec3 f = cos((p-i)*acos(-1.))*(-.5)+.5;\na = mix(sin(cos(a)*a),sin(cos(1.+a)*(1.+a)), f.x);\na.xy = mix(a.xz, a.yw, f.y);\nreturn mix(a.x, a.y, f.z);\n}\n\nfloat sphere(vec3 p, vec4 spr)\n{\nreturn length(spr.xyz-p) - spr.w;\n}\n\nfloat flame(vec3 p)\n{\nfloat d = sphere(p*vec3(1.,.5,1.), vec4(.0,-1.,.0,1.));\nreturn d + (noise(p+vec3(.0,iTime*2.,.0)) + noise(p*3.)*.5)*.25*(p.y) ;\n}\n\nfloat scene(vec3 p)\n{\nreturn min(100.-length(p) , abs(flame(p)) );\n}\n\nvec4 raymarch(vec3 org, vec3 dir)\n{\nfloat d = 0.0, glow = 0.0, eps = 0.02;\nvec3  p = org;\nbool glowed = false;\n\nfor(int i=0; i<64; i++)\n{\n    d = scene(p) + eps;\n    p += d * dir;\n    if( d>eps )\n    {\n    if(flame(p) < .0)\n        glowed=true;\n    if(glowed)\n            glow = float(i)/64.;\n    }\n}\nreturn vec4(p,glow);\n}\n\nvoid main( )\n{\n// vec2 v = -1.0 + 2.0 * fragCoord.xy / iResolution.xy;\n// v.x *= iResolution.x/iResolution.y;\n\nvec2 v = 2.*v_st.xy - vec2(1., 1.);\nv.x /= 0.6;\nv.y /= 0.4;\nvec3 org = vec3(0., -2., 4.);\nvec3 dir = normalize(vec3(v.x*1.6, -v.y, -1.5));\n\nvec4 p = raymarch(org, dir);\nfloat glow = p.w;\n\nvec4 col = mix(vec4(1.,.5,.1,1.), vec4(0.1,.5,1.,1.), p.y*.02+.4);\n\nfragColor = mix(vec4(0.), col*u_glow, pow(glow*2.,4.));  \n}\n',
    'u_visibleColor',
    'src',
    'color("',
    'Geometry',
    'shadowMaps',
    'url(',
    'CullFace',
    'disableDepthTestDistance',
    'colorStops',
    'calc(60% - ',
    'Cesium3DTileStyle',
    'move',
    '_translate',
    '6enuJjc',
    '\nuniform vec4 color; \nczm_material czm_getMaterial(czm_materialInput materialInput)\n{\n    czm_material material = czm_getDefaultMaterial(materialInput); \n    vec2 st = materialInput.st;          \n    float alpha = distance(st,center); \n    material.alpha = color.a  * alpha  * 1.5; \n    material.diffuse = color.rgb * 1.3;                            \n    return material;\n} ',
    'mozHidden',
    'imagery',
    '_startFovH',
    '_createLineEntity',
    '_setStyleHook',
    'CylinderGlowFlowWall',
    'getFragmentShader',
    'PolygonGeometry',
    '_geoJsonData',
    'polylineOpts',
    '_add',
    'animationSpeed',
    'textMarker',
    '_createCircleEntitiy',
    'none',
    'color.a',
    '_setRotate',
    'water_0',
    'czm_material czm_getMaterial(czm_materialInput materialInput)\n{\n     czm_material material = czm_getDefaultMaterial(materialInput);\n     vec2 st = materialInput.st;\n     vec4 colorImage = texture(image, vec2(fract(st.t - time), st.t));\n     material.alpha = colorImage.a * color.a;\n     material.diffuse =  1.9 * color.rgb  ;\n     return material;\n }',
    'u_showThroughEllipsoid',
    '\n            const float PI = 3.14159265359;\n            const float TWO_PI = 6.28318530718;\n            const int N = 3;\t\t\t\t// triangle polygons please\n            const float r0 = 0.01;\t\t\t// size of centre circle\n            const float r_blue = 0.025;\t\t// size of blue radar blips\n            const float r_red = 0.015;\t\t// size of red radar blips\n            const float edge = 0.95;\t\t// overall size\n            const float offset = 0.05;\n\n            float plot(const vec2 st, const float pct, const float width)\n              {\n                    return smoothstep(pct - width, pct, st.y) -\n                          smoothstep(pct, pct + width, st.y);\n                }\n\n            float drawPolygon(const vec2 polygonCenter, const int N, const float radius, vec2 pos)\n              {\n                pos = pos - polygonCenter;\n                float d = 0.0;\n                float a = atan(pos.x, pos.y);\n                float r = TWO_PI / float(N);\n                d = cos(floor(0.5 + a / r)*r - a)*length(pos);\n                return (1.0 - smoothstep(radius, radius + radius/10.0, d));\n              }\n\n            float gradations(const float a, const float gradNum, const float outRad, const float tickLen, const float tickWidth, const float r, const float move)\n              {\n                float f = step(0.0, cos((a + move)*gradNum) - tickWidth)*tickLen + (outRad - tickLen);\n                  return 1.0 - step(f, r) * 1.0 - step(r, outRad - tickLen);\n              }\n\n            void main( )\n            {\n              vec2 pos = 2.0*v_st.xy - vec2(1., 1.) ; // center what being drawn\n              pos /=2.05;\n              vec4 grndSpd = vec4(0.0, iTime/5.0, 0.0, 0.0);\n              vec4 mapcol = vec4 (1.0, 1.0, 1.0, 1.0);\n              \n              vec3 color = vec3(0.0, 0.0, 0.0);\n\n              float r = length(pos) * 2.0;\n              float a = atan(pos.y, pos.x); // angle of pixel\n              float an = PI - mod(iTime / 1.0, TWO_PI); // angle of radar sweep\n                float blipSpd = 3.0; // Blip / Trace speed\n              vec2 translate1 = vec2(cos(iTime/ blipSpd), sin(iTime / blipSpd));\n              vec2 translate2 = vec2(sin(iTime / blipSpd), cos(iTime / blipSpd));\n              vec2 left1 = translate1 * 0.35;\n              vec2 right1 = -translate1 * 0.30;\n              vec2 left2 = translate2 * 0.15;\n              vec2 right2 = -translate2 * 0.25;\n                \n            // Radar Sweep\n                float sn = step(PI/2.0, an) * step(-PI/2.0, (a + an)) * step(r, edge) * (1.0 - 0.55 * (a + (TWO_PI) - an));\n              float sw = step(an, a) * step(r, edge);\n              float s_blade = sw * (1.0 - (a - an) * 20.0);\n              float s = sw * (1.0 - 0.55 * (a - an));\n              s = max(sn,s);\n              float se = step(r, edge - 0.05);\n              \n            // Center point\n              float s1 = smoothstep(edge - 0.00, edge + 0.01, r)* smoothstep(edge + 0.02, edge + 0.01, r);   \n              \n            // Circular concentric rings\n              float s0 = 1.0 - smoothstep(r0 / 2.0, r0, length(pos));\n                float smb = (1.0 - smoothstep(0.2, 0.2 + 0.01, length(pos))) * (1.0 - smoothstep(0.2 +0.01, 0.2, length(pos)));\n                float smr = (1.0 - smoothstep(0.3, 0.3 + 0.01, length(pos))) * (1.0 - smoothstep(0.3 +0.01, 0.3, length(pos)));\n                \n            // Circular concentric gradations\n              float gradNum = 120.0;\n              float tickWidth = 0.9;\n              const float tickLen = 0.04;\n              float outRad = edge;\n              float move = 0.0;\n              float sm = 0.75*gradations(a, gradNum, outRad, tickLen, tickWidth, r, move);   \n              \n              gradNum = 36.0;\n              tickWidth = 0.95;\n              outRad = 0.6;\n              move = sin(iTime/10.0);\n              smr += 0.5*gradations(a, gradNum, outRad, tickLen, tickWidth, r, move);\n\n              outRad = 0.4;\n              move = cos(iTime/10.0);\n              smb += 0.5*gradations(a, gradNum, outRad, tickLen, tickWidth, r, move);\n\n            // Radial spoke gradations \n              float sr = plot(pos, pos.x, 0.003) * step(r, edge - 0.06);\n              sr += plot(vec2(0.0, 0.0), pos.x, 0.002) * step(r, edge - 0.06);\n              sr += plot(vec2(0.0, 0.0), pos.y, 0.003) * step(r, edge - 0.06);\n              sr += plot(-pos, pos.x, 0.003) * step(r, edge - 0.06);\n                sr *= 0.75;\n\n            // Blue circular radar blip traces\n              vec2 st_trace1 = left2;\n              float s_trace1 = s * (1.0 - smoothstep(r_blue / 10.0, r_blue, length(pos - st_trace1)));\n              s_trace1 += s * (1.0 - smoothstep(r_blue / 10.0, r_blue, length(pos - st_trace1 + vec2(+offset, +offset))));\n              s_trace1 += s * (1.0 - smoothstep(r_blue / 10.0, r_blue, length(pos - st_trace1 + vec2(+2.0 *offset, +2.0 *offset))));\n\n              vec2 st_trace2 = right1;\n              float s_trace2 = s * (1.0 - smoothstep(r_blue / 10.0, r_blue, length(pos - st_trace2)));\n\n            // Red Trianglular radar flight blip trace \n              vec2 st_trace3 = left1;\n              float st1 = s * (drawPolygon(st_trace3, N, r_red , pos));\n              st1 += s * (drawPolygon(st_trace3 + vec2(-offset, -offset), N, r_red, pos));\n              st1 += s * (drawPolygon(st_trace3 + vec2(+offset, -offset), N, r_red, pos));\n\n              vec2 st_trace4 = right2;\n              float st2 = s * (drawPolygon(st_trace4, N, r_red, pos));  \n                \n            // Lets add all the bits together and send them to screen\n              float s_grn = max(s * mapcol.y, s_blade);\n              s_grn = max(s_grn, (s0 +  sr + sm));\n              s_grn += s1 / 1.5  + smb + smr;\n\n              float s_red = st1*2.0 + st2*2.0 + smr;\n                \n              float s_blue = max(s_trace1 + s_trace2, s_blade) + smb;\n\n              if (s_trace1 > 0.0 || s_trace2 > 0.0) { s_blue = max(s, s_blue); s_grn = max(s_grn, s_blue); }\n\n              color += vec3(s_red , s_grn, s_blue);   \n                \n                vec4 texColor = mapcol * s;\n                \n                // Output to screen   \n                fragColor = vec4(color * u_glow, 0.8);//Set the screen pixel to that color\n\n            }\n    ',
    '_preRender',
    'loadLocalFile',
    '_cameraFrustum',
    '_drivenDistance',
    'dataset',
    'object',
    'ALPHA_BLEND',
    '当前绘制类型：线，至少需要2个点',
    'czm_material czm_getMaterial( czm_materialInput cmi )\n{\n    czm_material material = czm_getDefaultMaterial(cmi);\n    vec2 st = cmi.st;\n    float t = fract(speed *  czm_frameNumber / 1000.0);\n    vec2 st1 = vec2(fract(st.s * repeat - t),st.t); \n    float alpha = 1.-st.t;\n    float value = fract(st1.s/1.);\n    alpha *= sin(value * 3.1415926);\n    vec4 color_ = vec4(color.rgb * color.a, alpha * 1.2);\n    material.diffuse =color_.rgb*2.5;\n    material.alpha = color_.a;\n    return material;\n} ',
    'VertexFormat',
    'toArray',
    'POI文本点',
    '_domeLineVA',
    'error',
    'moveAll',
    'blendColor',
    '_measureResult',
    '_imageryLayer',
    'xt3d-fillet-label-content',
    '着色器特效魔法盒',
    'normalTexture',
    '按住alt拖拽改变整体高度',
    'topSteps',
    'south',
    'onclick',
    'editMode',
    'CylinderGlowGradientWall',
    'udpateControlPosition',
    'px; background: #ff000042; border: 1px dashed red; position: absolute;  left: -7px; top: 4px;pointer-events:none;z-index:-1',
    '_createEndEntity',
    'passes',
    'dynamicBorder',
    'cssText',
    '进攻方向',
    'modelViewMatrix',
    'fromFramebuffer',
    '\n#ifdef GL_ES\nprecision mediump float;\n#endif\n\n#extension GL_OES_standard_derivatives : enable\n\nfloat time;\nvec2 resolution;\nin vec2 v_st;\nvoid main( void ) {\n      time = czm_frameNumber /60.0;\n      resolution=czm_viewport.zw;\n     \n    vec2 p =v_st.xy* 2.0 - 1.0;// (v_st.xy * 2.0 - resolution) / min(resolution.x, resolution.y);\n    float rotTime = sin(time) + 1.0;\n    vec3 destColor = vec3(3.0 * rotTime, 1.0 + rotTime, 4.0 * rotTime);\n    float f = 0.0;\n    for(float i = 0.0; i < 64.0; i++){\n        float s = sin((time / 4.0) + i * -1.54321) * 1.5;\n        float c = cos((time / 4.0) + i * 1.54321) * 1.5;\n        f += .00025 / abs(length(p / vec2(c, s)) - 0.5);\n    }\n\n    vec3 color=vec3(destColor * f);\n    float alpha=(color.x+ color.y + color.z) / 3.; \n    if(alpha < 0.5){ \n        discard; \n    } \n    out_FragColor = vec4(color*1.2,alpha);\n\n}\n',
    'pickEllipsoid',
    '_addScaleAnchors',
    'debugFrustumColor',
    'perPositionHeight',
    'shaderEffetCoolBall',
    '_vertexEntity',
    '_container_drag',
    'circleImageRotate',
    'tan',
    '_createStartEntity',
    'outline',
    'ScreenSpaceEventType',
    '_totalDistance',
    'Feature',
    '_lastTime',
    '\n#ifdef GL_ES\n    precision highp float;\n        #endif \n    in vec3 position; \n    in vec2 st; \n    in vec3 normal; \n    uniform mat4 modelViewMatrix; \n    uniform mat3 normalMatrix; \n    uniform mat4 projectionMatrix; \n    out vec3 v_position; \n    out vec3 v_normal; \n    out vec2 v_st; \n    out vec3 v_light0Direction;  \n    void main(void)  {\n    vec4 pos =  modelViewMatrix * vec4( position,1.0);\n    v_normal =  normalMatrix *  normal;\n    v_st = st;\n    v_position = pos.xyz;\n    v_light0Direction = mat3( modelViewMatrix) * vec3(1.0,1.0,1.0);\n    gl_Position =  projectionMatrix * pos;\n}',
    'glowRange',
    'MagicBoxFS',
    "\nuniform vec4 u_intersectionColor;\nuniform float u_intersectionWidth;\nuniform vec4 u_lineColor;  \nbool inSensorShadow(vec3 coneVertexWC, vec3 pointWC)\n{\n    // Diagonal matrix from the unscaled ellipsoid space to the scaled space.    \n    vec3 D = czm_ellipsoidInverseRadii;\n\n    // Sensor vertex in the scaled ellipsoid space\n    vec3 q = D * coneVertexWC;\n    float qMagnitudeSquared = dot(q, q);\n    float test = qMagnitudeSquared - 1.0;\n    \n    // Sensor vertex to fragment vector in the ellipsoid's scaled space\n    vec3 temp = D * pointWC - q;\n    float d = dot(temp, q);\n    \n    // Behind silhouette plane and inside silhouette cone\n    return (d < -test) && (d / length(temp) < -sqrt(test));\n}\n \nvec4 getLineColor()\n{\n    return u_lineColor;\n}\n\nvec4 getIntersectionColor()\n{\n    return u_intersectionColor;\n}\n\nfloat getIntersectionWidth()\n{\n    return u_intersectionWidth;\n}\n\nvec2 sensor2dTextureCoordinates(float sensorRadius, vec3 pointMC)\n{\n    // (s, t) both in the range [0, 1]\n    float t = pointMC.z / sensorRadius;\n    float s = 1.0 + (atan(pointMC.y, pointMC.x) / czm_twoPi);\n    s = s - floor(s);\n    \n    return vec2(s, t);\n}",
    '着色器特效双环',
    'zIndex',
    'rotate',
    '着色器特效魔幻星星',
    '_mouseMoveALTHook',
    'CylinderGlowCircleType',
    '_initTipContent',
    'itemStyle',
    '获取材质失败！',
    '_url',
    'sphereElectric',
    'ClockRange',
    'PolylineGlowMaterialProperty',
    '_svgEditAnchors',
    'distance',
    'frustum',
    'border-right:3px solid ',
    'CircleImageDiffuseType',
    'SAMPLER_2D',
    'maximumSpeed',
    'arc',
    '_topOutlineGeometry',
    'SceneMode',
    'border-left: 2px dashed ',
    'DEPTH_STENCIL',
    'initMouseMoveEventHandler',
    'add',
    'getMoveByMatrix',
    'transform',
    'Pass',
    'Point',
    '_topHeight',
    'getPoint',
    'multiplyByUniformScale',
    '_getFixedFrameToEastNorthUpTransformFromWorldMatrix',
    'BlendingState',
    '着色器特效星链',
    'dot',
    '_startFovV',
    '上下拖拽改变高度',
    'showBackground',
    'removeAll',
    '.b-b',
    'blending',
    'replaceCache',
    'rgba(238,103,98,0.5)',
    'ALWAYS',
    '按住alt拖拽改变高度',
    'right',
    '3.0',
    '_buildingStyle',
    'toCssColorString',
    '_sectorLineVA',
    'event',
    'trackedEntity',
    'WallImageFlowType',
    '_createComponent',
    'value-type',
    'origin',
    '2262452ERbZDA',
    '_createMarkerBg',
    'remove',
    '_needUpdate',
    '个点，点击鼠标左键确定第',
    'ellipseOpts',
    'border',
    'toJson',
    'clone',
    'info-item',
    '_createImage_3',
    '未命名图元',
    'HALF_FLOAT',
    'czm_material czm_getMaterial(czm_materialInput materialInput){\n    czm_material material = czm_getDefaultMaterial(materialInput);\n    vec2 st = materialInput.st;\n    vec2 center = st - vec2(0.5,0.5);\n    float time = -czm_frameNumber * 3.1415926 / 180.;//扫描速度1度\n    float sin_t = sin(time);\n    float cos_t = cos(time);\n    vec2 center_rotate = vec2(center.s*cos_t-center.t*sin_t+0.5,center.s*sin_t+center.t*cos_t+0.5);\n    vec4 color = texture(image,center_rotate);\n    vec3 tColor = color.rgb * u_color.rgb;\n    tColor *= u_color.a;\n    material.diffuse = tColor;\n    float length = 2. - length(center)/0.5;\n    material.alpha = color.a * pow(length, 0.5);//color.r = 0 或1\n    return material;\n}',
    'czm_old_main',
    'textStyle',
    'shaderEffetColorfulStellarChain',
    'polyElevationCountour',
    '_position',
    'position:absolute;top:0px;left:0px;color:white;pointer-events: none;',
    'render',
    'glowLine',
    'name',
    'px;background: linear-gradient(to right,',
    'dispatchEvent',
    'DepthFunction',
    '_reAddHook',
    'regularWall',
    'content',
    'simpleMarker',
    'atan',
    'initLeftUpEventHandler',
    'siteTimes',
    'fog',
    'viewer是必须的！',
    'outFragColorBody',
    'toGeoJson',
    'setBuilding',
    '_measureTool',
    '/6910AEE9470643969E7A489A5C3CAD8C.png',
    '_primitive',
    'addScaleAnchor',
    'POSITION_ONLY',
    'rgba(140,0,58,0.99)',
    '按下鼠标右键结束绘制',
    '_floatingState',
    '点击鼠标左键开始，点击鼠标右键取消',
    '_graphicDraw',
    'image/png',
    '_pickId',
    'multilineLabel',
    'createElement',
    'flame',
    'muted',
    'GeometryAttribute',
    'isLine',
    'm3u8',
    'create',
    'createDrawCommand',
    '_swallowTailPnt',
    'isEditing',
    '_primitives',
    'shaderEffetFireworks',
    '_origin',
    'clientY',
    '_videoDom',
    'defaultValue',
    'cull',
    '_ellipseOpts',
    '\nczm_material czm_getMaterial(czm_materialInput materialInput)  { \n    czm_material material = czm_getDefaultMaterial(materialInput); \n    vec2 st = materialInput.st; \n    float t =fract(czm_frameNumber * speed / 1000.0);\n    vec4 colorImage =texture(image, vec2(fract( repeat * st.s - t),fract(st.t))); \n    material.alpha =  colorImage.a * color.a; \n    material.diffuse =(colorImage.rgb + color.rgb)* 1.5; \n    if(sampling==1.){\n        material.diffuse =color.rgb * 1.5; \n    } \n    if(sampling==2.){\n        material.diffuse =colorImage.rgb * 1.5; \n    } \n    return material;\n }',
    '_addOriginPoint',
    '_handlePickObject',
    'SphereElectricType',
    '_drawEntity',
    'showOrigin',
    '_computeHierarchy',
    '_outerCylinder',
    'lightingModel',
    'activateHook',
    'addEventListener',
    'union',
    '_currentFrustum',
    'FlameRingFS',
    '\n            void vertexMain(VertexInput vsInput, inout czm_modelVertexOutput vsOutput){\n                v_normalMC = vsInput.attributes.normalMC;\n              }',
    'onerror',
    'graphics',
    'Axis',
    'radius',
    'text',
    'transparent',
    '_updateAnchor',
    '_createHeatContainer',
    '_vertexEntity_0',
    'projection',
    '_computedModelMatrix',
    'spotLight',
    'verticalOrigin',
    'alertMarker',
    '着色器特效星座链',
    'globeShow',
    'appearance',
    '_polylineOpts',
    'mouseVal',
    'subSegmentH',
    'setCursor',
    '闭合曲线面',
    'scanPlaneColor',
    '_length',
    '.0), 1.0); // 渐变\n                if(!',
    'atan2',
    '_frustum',
    'customFsMain',
    '_leftUp',
    'Ray',
    'fromCssColorString',
    '当前绘制类型：点，需要一个点',
    'createViewportQuadCommand',
    'join',
    'CircleScanMaterialProperty_2',
    '_counter',
    'shaderEffetGlowBox',
    'change',
    'deactivateHook',
    'shaderGeometryType',
    '_scene',
    '.b-t-l',
    'urlTemplate',
    'updateFrameState',
    'scaleWidth',
    '_getDOMContainer',
    '抛物线',
    '_layerFilter',
    '_singleEnd',
    '_getArrowPoints',
    '释放鼠标完成修改',
    '_connPoint',
    'drawingBufferWidth',
    '_activate',
    '\n            </div>\n          </div>\n        </div>\n        <div class="b-t-l"></div>\n        <div class="b-b-r"></div>\n      </div>\n      <div class="arrow"></div> \n      </div>',
    'getPickRay',
    'blue',
    '_scale',
    ' \nprecision highp float;\n \n//uniform float time; \nfloat time;\nconst float PI = acos(-1.);\nconst float TAU = PI * 2.;\n\n#define saturate(x) clamp(x,0.,1.)\n#define _tail2x(p,n) (mod(p,2.)-1.)\n\nfloat Hash( vec2 p, in float s ){\n    return fract(sin(dot(vec3(p.xy,10.0 * abs(sin(s))),vec3(27.1,61.7, 12.4)))*273758.5453123);\n}\n\nfloat noise(in vec2 p, in float s){\n  vec2 i = floor(p);\n  vec2 f = fract(p);\n  return mix(\n    mix(Hash(i + vec2(0.,0.), s), Hash(i + vec2(1.,0.), s),f.x),\n    mix(Hash(i + vec2(0.,1.), s), Hash(i + vec2(1.,1.), s),f.x),f.y) * s;\n}\n\nfloat fbm(vec2 p){\n  float v = 0.0;\n  v += noise(p*34., .1);\n  v += noise(p*20., .04);\n  return v;\n}\n\nvec2 mPolar(vec2 p){\n  float a = atan(p.y, p.x);\n  float r = length(p);\n  return vec2(a, r);\n}\n\nvec2 tailY2x(vec2 p,float n){p*=n;return vec2(p.x,_tail2x(p.y,n));}\nmat2 rot(float a){float c=cos(a),s=sin(a);return mat2(c,-s,s,c);}\n\nhighp float rand(vec2 p){\n  highp float a = 12.9898;\n  highp float b = 78.233;\n  highp float c = 43758.5453;\n  highp float dt= dot(p ,vec2(a,b));\n  highp float sn= mod(dt,3.14);\n  return fract(sin(sn) * c);\n}\n\n// signed distance\nfloat sd(float d,float r){return r-d;} \nfloat sd(float d){return 1.-d;} \n// glow + fill\nfloat gf(float d,float r){return r/d;} \nfloat gf(float d){return 1./d;} \n\nfloat fill_na(float d){return step(0.,d);}\nfloat fill(float d){return smoothstep(0.,0.01,d);}\nfloat stroke(float d,float w){return 1.-smoothstep(w,w+0.01,abs(d));}\nfloat strokeInner(float d,float w){return stroke(d-w,w);}\nfloat strokeOuter(float d,float w){return stroke(d+w,w);}\n\nfloat lSquare(vec2 p){p = abs(p);return max(p.x,p.y);}     \n\nfloat lPoly(vec2 p,float n){\n  float a = atan(p.x,p.y)+PI;\n  float r = TAU/n;\n  return cos(floor(.5+a/r)*r-a)*length(p)/cos(r*.5);\n}\n\nfloat strokeStar(vec2 p,float n,float w){\n  float l =strokeInner(sd(lPoly(p,n*.5)),w);\n  l+=strokeInner(sd(lPoly(mod(n,2.)!=0.?vec2(-p.x,p.y):p*rot(TAU/n),n*.5)),w);\n  return l;\n}\n\nvec2 mPoly(vec2 p,float n,float s){\n  float r = TAU / n;\n  float a = floor(atan(p.y,p.x)/r)*r+r*.5;\n  return (vec2(cos(a),sin(a))*s-p)*rot(-a-PI*.5);\n}\n\nfloat wsaw(float x){return fract(x*.5+.5)*2.-1.;}\nfloat wtri(float x){return abs(2.*fract(x*.5-.25)-1.)*2.-1.;}\nfloat utri(float x){return abs(2.*fract(x*.5-.5)-1.);}\nfloat wtrz(float x,float w){return clamp(wtri(x*2.)*w,-1.,1.);} // 台形波 trapezoidal wave\n\n// ease\nfloat o2(float t){t=1.-t;return 1.-t*t;}\nfloat oN(float t,float n){return 1.-pow(1.-t,n);}\n\nfloat dot2(vec2 p){return dot(p,p);}\n\nvec2 mSimplePerspective(vec2 p){p.y+=.2;p.y*=3.;return p;}\n\nfloat ring(vec2 p,float t){\n  float alpha =    fract(-t);\n  float l = 0.;\n  vec2 p3=mPoly(p*rot(PI*.5),10.,1.);\n  l+=saturate(gf(abs(p3.x),.03)*fill(sd(length(p),1.1+fract(t)))*(1.-fill(sd(length(p),.9+fract(t))))*alpha);\n  \n  l+=saturate(.02/abs(sd(length(p),1.1+fract(t)))*alpha);\n  vec2 p4=mPolar(p*(.57-oN(t,1.3)*.28)).yx;\n  p4.x-=.65;\n  l+= saturate(abs(1./((p4.x + fbm( p4 + vec2(sin(t*.2),t*0.1))) * 50.0))*sd(dot2(tailY2x(p4+vec2(.1,0.),12.)),.9)*alpha);\n  return l;\n}\n\nfloat summoningCircle(vec2 p){\n  float l=0.;\n  l+=fill(sd(lSquare(p*rot(PI/3.*1.5)*vec2(100.,1.)),1.));\n  l+=fill(sd(lSquare(p*rot(PI/3.*2.5)*vec2(100.,1.)),1.));\n  l+=fill(sd(lSquare(p*rot(PI/3.*3.5)*vec2(100.,1.)),1.));\n  l=saturate(l);\n  l-=fill(sd(lPoly(p,3.)));\n  l=saturate(l);\n  float r = atan(p.y,p.x);\n  l+=strokeOuter(sd(length(p),.98),.008+wtrz(r/TAU*3.,12.)*.005);\n  l+=strokeInner(sd(length(p),.95),.005);\n  l+=strokeInner(sd(lPoly(p,3.)),.01);\n  l+=strokeInner(sd(lPoly(p,3.),.88),.02);\n  l+=strokeInner(sd(lPoly(p,6.),.53),.01);\n  vec2 q=mPoly(p*rot(PI*.5),3.,.5);\n  l+=fill(sd(lPoly(q,3.),.3));\n  vec2 q2=mPoly(p*rot(PI/3.+PI*.5),3.,.7);\n  l+=fill(sd(lPoly(q2,3.),.1));\n  l+=strokeInner(sd(lPoly(p*rot(PI),3.),.5),.02);\n  l+=fill(sd(length(p),.05));\n  vec2 q3=mPoly(p*rot(PI*.5),3.,1.);\n  l=saturate(l);\n  l-=fill(sd(length(q3),.2));\n  l=saturate(l);\n  l+=strokeInner(sd(length(q3),.18),.005);\n  l+=strokeInner(sd(length(q3),.15),.005);\n  l+=strokeStar(q3*rot(PI)*7.,6.,.1);\n  return l;\n}\n\nfloat render2(vec2 p){\n\tfloat time = time;\n\tvec2 _p = p;\n\t\n  //p=mSimplePerspective(p);\n  p*=rot(time);\n  p*=2.;\n  float tt = time*.75;\n  float l2 = ring(p,o2(fract(tt)));\n  l2+=ring(p*rot(PI/3.),o2(fract(tt+.5)));\n  float l=0.;\n  l = summoningCircle(p*=rot(floor(time*12.)/3.));\n  float return1 = l2;\n\t\n\t\n\ttime *= atan(-1.)/10.;\n\t_p *= 1.8;\n\t\n\ttime += pow(length(_p), .1)*15.;\n\t\n\t\n  //p=mSimplePerspective(p);\n\tp = _p;\n  p*=rot(time);\n  p*=2.;\n  tt = time*.75;\n  l2 = ring(p,o2(fract(tt)));\n  l2+=ring(p*rot(PI/3.),o2(fract(tt+.5)));\n  l=0.;\n  l = summoningCircle(p*=rot(floor(time*12.)/3.));\n  return max(return1, l2);\n}\nfloat render(vec2 p){\n\treturn render2(p);\n  //p=mSimplePerspective(p);\n  p*=rot(time);\n  p*=2.;\n  float tt = time*.75;\n  float l2 = ring(p,o2(fract(tt)));\n  l2+=ring(p*rot(PI/3.),o2(fract(tt+.5)));\n  float l=0.;\n  l = summoningCircle(p*=rot(floor(time*12.)/3.));\n  return l2;\n}\n\nin vec2 v_st;\n\nvoid main( ) {\n    time = czm_frameNumber /100.0;\n    vec2 resolution = czm_viewport.zw;\n    vec2  st =v_st  * 2.5 - 1.25;\n\n  vec2 p = st;//(gl_FragCoord.xy * 2.0 - resolution) / 200.;\n  float l=0.;\n  l = (render(p)+render(p+vec2(0.,1./min(resolution.x, resolution.y))))*.5;\n\n  vec3 color=l*vec3(0.75, 0.5, .05) * 2.;\n  float alpha=(color.x+ color.y + color.z) / 3.0; \n//   if(alpha < 0.05){ \n//       discard; \n//   } \n  out_FragColor = vec4(color,alpha);\n}',
    'polyWater',
    'buildingStyle',
    'destroy',
    'viewDome',
    'labelOpts',
    '#15d1f2',
    'rgba(149,0,235,0.99)',
    'px ',
    '_title',
    '_aspectRatio',
    'Buffer',
    'planeHeat',
    '_createVideoCamera',
    'ALL',
    'removeGraphic',
    'animateLabel',
    '拖拽改变半径',
    'splice',
    '_drawing',
    'MagicBallFS',
    '_lastUpdate',
    'Transforms',
    '弹跳点',
    '_setTipContent',
    'UNSIGNED_SHORT',
    'satellite',
    '_createDrawCommand',
    'CENTER',
    'shaderEffetExplosion',
    'Matrix4',
    'viewport',
    '_editGraphic',
    '_clusterEnable',
    'measureEnd',
    'showScanPlane',
    'stopTime',
    'play',
    '无效的背景色格式',
    'CircleSpiralType',
    'absolute',
    '_setInputContent',
    'RangeRadar',
    'attributeLocations',
    'clampToGround',
    'getElementById',
    '_getMoveHeight',
    'xt3d-camera-picker-container',
    'positionChange',
    'u_showIntersection',
    '_appearanceOpts',
    'polygon',
    'lookAtState',
    'topHeight',
    '追加失败！类型错误，请确保是通过toGeoJson（）方法导出的数据',
    '_outerFovRadiusPairs',
    '\nvec2 rotate2D(vec2 _st, float _angle){\n  _st -= 0.5;\n  _st =  mat2(cos(_angle),-sin(_angle),\n              sin(_angle),cos(_angle)) * _st;\n  _st += 0.5;\n  return _st;\n} \nczm_material czm_getMaterial(czm_materialInput materialInput){\nczm_material material = czm_getDefaultMaterial(materialInput);\nvec2 st = materialInput.st;\nfloat angle = czm_frameNumber * speed / 1000.0; \n//angle=mod(angle,360.);//如果要限定 旋转最大最小值 可通过entity实现\nst = rotate2D(st,(angle));\nvec4 colorImage = texture(image, st);\nmaterial.alpha =colorImage.a * color.a ;\nmaterial.diffuse = color.rgb;\nreturn material;\n} ',
    'particleGraphic',
    '\n                        uniform vec4 color;\n                        uniform float repeat;\n                        uniform float offset;\n                        uniform float thickness;\n                        czm_material czm_getMaterial(czm_materialInput materialInput)\n                        {\n                            czm_material material = czm_getDefaultMaterial(materialInput);\n                            float sp = 1.0/repeat;\n                            vec2 st = materialInput.st;\n                            float dis = distance(st, vec2(0.5));\n                            float m = mod(dis - czm_frameNumber/1000., sp);\n                            float a = step(sp*(1.0-thickness), m);\n                            material.diffuse = color.rgb;\n                            material.alpha = a * color.a;\n                            return material;\n                        } \n                    ',
    'baseWaterColor',
    'waterColor',
    '_createTopGeometry',
    'videoFusion',
    'SphereScanType',
    '_bIcon',
    '90%',
    'arcHeat',
    'backgroundPadding',
    'minimumParticleLife',
    'bounce',
    'longitude',
    'RectSensor',
    'shininess',
    '_moveAllPosition',
    '_minRadius',
    'BoundingRectangle',
    'margin-left: 19px;',
    'lineRadius',
    'windowPosition',
    '着色器特效发光棱锥',
    'clusterEnable',
    'intersectionColor',
    '_scanePlaneSP',
    '_label',
    '_getOrientation',
    'PositionPick',
    'layerId',
    'trackedView',
    'flyToById',
    '_slices',
    '_chart',
    'offsetHeight',
    '_createDebugFrustum',
    '_textureAtlasGUID',
    'colorMap',
    '_reflectorProjectionMatrix',
    'dom',
    'getValueOrUndefined',
    'createPlaneGeometry',
    'canvasHeight',
    'shaderEffetVirus',
    '_index',
    'flyTo',
    'shaderEffetGlowPyramid',
    'lineTo',
    'Cartesian4',
    '_heatmap',
    '_value',
    'height: ',
    'ParticleSystem',
    'CylinderGlowGradientWallSource',
    'AQUA',
    'PI_OVER_TWO',
    'north',
    '集结地',
    'createIndexBuffer',
    '_loadGif',
    '_editAnchors',
    'positionProperty',
    'GeometryInstance',
    'vjs-default-skin',
    '\nprecision highp float;                      \nczm_material czm_getMaterial(czm_materialInput materialInput)\n{ \n    czm_material material = czm_getDefaultMaterial(materialInput);  \n    float logDepthOrDepth = czm_unpackDepth(texture(czm_globeDepthTexture, gl_FragCoord.xy / czm_viewport.zw));\n\n    if(logDepthOrDepth>0.0){\n        vec4 eyeCoordinate = czm_windowToEyeCoordinates(gl_FragCoord.xy, logDepthOrDepth);\n        vec4 worldCoordinate4 = czm_inverseView * eyeCoordinate;\n        vec3 worldCoordinate = worldCoordinate4.xyz / worldCoordinate4.w;\n        //当前高度\n        float l=length(worldCoordinate.xyz);\n        l-= 6378137.;\n        materialInput.height=l;  \n   float distanceToContour = mod(materialInput.height, spacing);\n\n  #if (__VERSION__ == 300 || defined(GL_OES_standard_derivatives))\n      float dxc = abs(dFdx(materialInput.height));\n      float dyc = abs(dFdy(materialInput.height));\n      float dF = max(dxc, dyc) * czm_pixelRatio * width;\n      float alpha = (distanceToContour < dF) ? 1.0 : 0.0;\n  #else\n      // If no derivatives available (IE 10?), use pixel ratio\n      float alpha = (distanceToContour < (czm_pixelRatio * width)) ? 1.0 : 0.0;\n  #endif\n\n       vec4 outColor = czm_gammaCorrect(vec4(color.rgb, alpha * color.a));\n       material.diffuse = outColor.rgb;\n       material.alpha = outColor.a; \n    } \n    return material;\n}\n',
    'selectedObjectChanged',
    'yellow',
    'image2',
    '_initContent',
    'Camera',
    'matrixState',
    'button',
    'postProcessStages',
    'Constant',
    'shaderEffetShield',
    'REPLACE',
    '_frontFaceRS',
    '_circleEntity',
    'url',
    'boundingVolume',
    'fromRotationMatrix',
    'pixelMax',
    'shaderEffetFlameRing',
    'SceneTransforms',
    'sizeInMeters',
    '_applyGravity',
    '_drawCommand',
    'rgba(0, 255, 255, 0.31)',
    '电弧穹顶材质',
    '_udpateRectangleGeometry',
    'bounceHeight',
    'NONE',
    '_addPostRender',
    'Primitive',
    '_html5_api',
    'shaderEffetShinyRing',
    '\nuniform vec4 color;\nuniform float speed;\n#define pi 3.1415926535\n#define PI2RAD 0.01745329252\n#define TWO_PI (2. * PI)\nfloat rands(float p){\n  return fract(sin(p) * 10000.0);\n}\nfloat noise(vec2 p){\n  float time = fract( czm_frameNumber * speed / 1000.0);\n  float t = time / 20000.0;\n  if(t > 1.0) t -= floor(t);\n  return rands(p.x * 14. + p.y * sin(t) * 0.5);\n}\nvec2 sw(vec2 p){\n  return vec2(floor(p.x), floor(p.y));\n}\nvec2 se(vec2 p){\n  return vec2(ceil(p.x), floor(p.y));\n}\nvec2 nw(vec2 p){\n  return vec2(floor(p.x), ceil(p.y));\n}\nvec2 ne(vec2 p){\n  return vec2(ceil(p.x), ceil(p.y));\n}\nfloat smoothNoise(vec2 p){\n  vec2 inter = smoothstep(0.0, 1.0, fract(p));\n  float s = mix(noise(sw(p)), noise(se(p)), inter.x);\n  float n = mix(noise(nw(p)), noise(ne(p)), inter.x);\n  return mix(s, n, inter.y);\n}\nfloat fbm(vec2 p){\n  float z = 2.0;\n  float rz = 0.0;\n  vec2 bp = p;\n  for(float i = 1.0; i < 6.0; i++){\n    rz += abs((smoothNoise(p) - 0.5)* 2.0) / z;\n    z *= 2.0;\n    p *= 2.0;\n  }\n  return rz;\n}\nczm_material czm_getMaterial(czm_materialInput materialInput)\n{\n  czm_material material = czm_getDefaultMaterial(materialInput);\n  vec2 st = materialInput.st;\n  vec2 st2 = materialInput.st;\n  float time = fract( czm_frameNumber * speed / 1000.0);\n  if (st.t < 0.5) {\n    discard;\n  }\n  st *= 4.;\n  float rz = fbm(st);\n  st /= exp(mod( time * 2.0, pi));\n  rz *= pow(15., 0.9);\n  vec4 temp = vec4(0);\n  temp = mix( color  / rz, vec4(color.rgb, 0.1), 0.2);\n  if (st2.s < 0.05) {\n    temp = mix(vec4(color.rgb, 0.1), temp, st2.s / 0.05);\n  }\n  if (st2.s > 0.95){\n    temp = mix(temp, vec4(color.rgb, 0.1), (st2.s - 0.95) / 0.05);\n  }\n  material.diffuse = temp.rgb;\n  material.alpha = temp.a * 2.0;\n  return material;\n}',
    'sampling',
    '\n     in vec3 position3DHigh;\n     in vec3 position3DLow;\n     in float batchId;\n     \n     out vec3 v_positionEC;\n     \n     void main() {\n       vec4 p = czm_computePosition();\n       v_positionEC = (czm_modelViewRelativeToEye * p).xyz;\n       gl_Position = czm_modelViewProjectionRelativeToEye * p;\n     }',
    '_svgDom',
    '#4984ed',
    'singleTile',
    '\n                 uniform sampler2D colorTexture; // 反射贴图 \n                 uniform sampler2D normalTexture; // 法线贴图\n                 float time;\n                 \n                 uniform mat4 fixedFrameToEastNorthUpTransform; // 水面的东北天矩阵的逆矩阵\n                 \n                 // 从顶点着色器传来的\n                 in vec4 v_worldPosition; // 当前像素的世界坐标\n                 in vec4 v_uv; // 原本的纹理坐标乘以贴图矩阵\n                 \n                 // 可配置的参数\n                 uniform float rippleSize; // 波纹大小（数值越大波纹越密集）\n                 uniform vec4 waterColor; // 水面颜色\n                 uniform float waterAlpha; // 水面透明度\n                 uniform float reflectivity; // 水面反射率\n                 uniform vec3 lightDirection; // 光照方向\n                 uniform float sunShiny; // 光照强度\n                 uniform float distortionScale; // 倒影的扭曲程度\n                 \n                 vec3 sunDirection ;//= normalize( lightDirection );\n                 vec3 sunColor;// = vec3( 1.0 );\n                 \n                 \n                 // 获取噪声\n                 // vec4 czm_getWaterNoise(sampler2D normalMap, vec2 uv, float time, float angleInRadians)\n                 vec4 getNoise( sampler2D normalMap, vec2 uv ) {\n                     vec2 uv0 = ( uv / 103.0 ) + vec2( time / 17.0, time / 29.0 );\n                     vec2 uv1 = uv / 107.0 - vec2( time / -19.0, time / 31.0 );\n                     vec2 uv2 = uv / vec2( 8907.0, 9803.0 ) + vec2( time / 101.0, time / 97.0 );\n                     vec2 uv3 = uv / vec2( 1091.0, 1027.0 ) - vec2( time / 109.0, time / -113.0 );\n                     vec4 noise = texture( normalMap, uv0 ) +\n                         texture( normalMap, uv1 ) +\n                         texture( normalMap, uv2 ) +\n                         texture( normalMap, uv3 );\n                     return noise * 0.5 - 1.0;\n                 }\n                 \n                 void sunLight( const vec3 surfaceNormal, const vec3 eyeDirection, float shiny, float spec, float diffuse, inout vec3 diffuseColor, inout vec3 specularColor ) {\n                     vec3 reflection = normalize( reflect( -sunDirection, surfaceNormal ) );  // 获得太阳对表面法线的反射向量\n                     float direction = max( 0.0, dot( eyeDirection, reflection ) );  // 当太阳反射方向和眼睛的方向一致时，direction 最大，为 1，当角度大于 90度时最小，最小为 0\n                     specularColor += pow( direction, shiny ) * sunColor * spec;\n                     diffuseColor += max( dot( sunDirection, surfaceNormal ), 0.0 ) * sunColor * diffuse;\n                 }\n                 \n                 czm_material czm_getMaterial(czm_materialInput materialInput) {\n                     czm_material material = czm_getDefaultMaterial(materialInput);\n                     time= czm_frameNumber / 100.0;\n                     // 通过法线贴图计算新的表面法线\n                     vec2 transformedSt = materialInput.st * 2.0 - 1.0;  // [0, 1] => [-1, 1]\n                     vec4 noise = getNoise( normalTexture, transformedSt * rippleSize );\n                     vec3 surfaceNormal = normalize( noise.xzy );  // [0, +1]，Y up\n                 \n                     // 漫反射光\n                     vec3 diffuseLight = vec3( 0.0 );\n                     // 高光\n                     vec3 specularLight = vec3( 0.0 );\n                 \n                     // 获取视线方向（世界坐标）\n                     vec3 eye = ( czm_inverseView * vec4( vec3(0.0), 1.0 ) ).xyz;\n                     // 获取视线方向（水面的本地坐标）\n                     eye = ( fixedFrameToEastNorthUpTransform * vec4( eye, 1.0) ).xyz;\n                     // 当前像素的本地坐标\n                     vec3 world = ( fixedFrameToEastNorthUpTransform * vec4( v_worldPosition.xyz, 1.0) ).xyz;\n                 \n                     vec3 worldToEye = eye - world;  // east, north, up\n                     worldToEye = vec3( worldToEye.x, worldToEye.z, -worldToEye.y );  // Y up\n                     vec3 eyeDirection = normalize( worldToEye );\n                 \n                     float shiny = sunShiny;\n                     float spec = 2.0;\n                     float diffuse = 0.5;\n                     sunLight( surfaceNormal, eyeDirection, shiny, spec, diffuse, diffuseLight, specularLight );\n                 \n                     float distance = length( worldToEye );\n                     float distortionScale = distortionScale;\n                     vec2 distortion = surfaceNormal.xz * ( 0.001 + 1.0 / distance ) * distortionScale;\n                     vec3 reflectionSample = vec3( texture( colorTexture, (v_uv.xy / v_uv.w) * 0.5 + 0.5 + distortion ) );\n                 \n                     float theta = max( dot( eyeDirection, surfaceNormal ), 0.0 );\n                     float reflectivity = reflectivity;\n                     float reflectance = mix( reflectivity, 1.0, pow( 1.0 - theta, 5.0 ) );\n                 \n                     vec3 waterColor = waterColor.rgb;\n                 \n                     // surfaceNormal 是以反射平面为 X-Y 平面的，\n                     // 所以 eyeDirection 也得是以反射平面为 X-Y 平面。\n                     vec3 scatter = max( 0.0, dot( surfaceNormal, eyeDirection ) ) * waterColor;\n                     vec3 albedo = mix(\n                         sunColor * diffuseLight * 0.3 + scatter,\n                         vec3( 0.1 ) + reflectionSample * 0.9 + reflectionSample * specularLight,\n                         reflectance\n                     );\n                     material.diffuse = albedo.rgb;\n                     material.alpha = waterAlpha;\n                 \n                     return material;\n                 }',
    'cross',
    'squadCombat',
    'borderColor',
    '多行文本标签',
    'STATIC_DRAW',
    'xt3d',
    '_sectorFrontCommand',
    'Sampler',
    'bgMap',
    '_addWindowEvent',
    'CylinderGlowCircleSource',
    '分队战斗',
    'cesiumColor',
    '_frustumPlanes',
    '20241110',
    'features',
    'passState',
    '_editAnchorType',
    '_createGeometry',
    'pixelOffsetScaleByDistance',
    'IndexDatatype',
    'rgba(255,255,255,1)',
    'data:image/svg+xml;base64,',
    'amd',
    '视频格式错误！',
    '_measureEnd',
    'handleLeftUp',
    'top',
    '_handleOverObject',
    'replace',
    '_model',
    'increment',
    '值类型：',
    'PolygonHierarchy',
    'distortTexture',
    '_setStyle',
    '\n    \n    ',
    '_highlightDom',
    "\n#define MODEL_ROTATION vec2(.3, .25)\n#define CAMERA_ROTATION vec2(.5, .5)\n\n// 0: Defaults\n// 1: Model\n// 2: Camera\n#define MOUSE_CONTROL 1\n\n//#define DEBUG\n\n// 1, 2, or 3\n//#define LOOP 1\n\n\n// --------------------------------------------------------\n// HG_SDF\n// https://www.shadertoy.com/view/Xs3GRB\n// --------------------------------------------------------\n\nvoid pR(inout vec2 p, float a) {\n    p = cos(a)*p + sin(a)*vec2(p.y, -p.x);\n}\n\nfloat pReflect(inout vec3 p, vec3 planeNormal, float offset) {\n    float t = dot(p, planeNormal)+offset;\n    if (t < 0.) {\n        p = p - (2.*t)*planeNormal;\n    }\n    return sign(t);\n}\n\nfloat smax(float a, float b, float r) {\n    float m = max(a, b);\n    if ((-a < r) && (-b < r)) {\n        return max(m, -(r - sqrt((r+a)*(r+a) + (r+b)*(r+b))));\n    } else {\n        return m;\n    }\n}\n\n\n// --------------------------------------------------------\n// Icosahedron domain mirroring\n// Adapted from knighty https://www.shadertoy.com/view/MsKGzw\n// --------------------------------------------------------\n\n#define PI 3.14159265359\n\nvec3 facePlane;\nvec3 uPlane;\nvec3 vPlane;\n\nint Type=5;\nvec3 nc;\nvec3 pab;\nvec3 pbc;\nvec3 pca;\n\nvoid initIcosahedron() {//setup folding planes and vertex\n    float cospin=cos(PI/float(Type)), scospin=sqrt(0.75-cospin*cospin);\n    nc=vec3(-0.5,-cospin,scospin);//3rd folding plane. The two others are xz and yz planes\n    pbc=vec3(scospin,0.,0.5);//No normalization in order to have 'barycentric' coordinates work evenly\n    pca=vec3(0.,scospin,cospin);\n    pbc=normalize(pbc); pca=normalize(pca);//for slightly better DE. In reality it's not necesary to apply normalization :)\n\tpab=vec3(0,0,1);\n    \n    facePlane = pca;\n    uPlane = cross(vec3(1,0,0), facePlane);\n    vPlane = vec3(1,0,0);\n}\n\nvoid pModIcosahedron(inout vec3 p) {\n    p = abs(p);\n    pReflect(p, nc, 0.);\n    p.xy = abs(p.xy);\n    pReflect(p, nc, 0.);\n    p.xy = abs(p.xy);\n    pReflect(p, nc, 0.);\n}\n\n\n// --------------------------------------------------------\n// Triangle tiling\n// Adapted from mattz https://www.shadertoy.com/view/4d2GzV\n// --------------------------------------------------------\n\nconst float sqrt3 = 1.7320508075688772;\nconst float i3 = 0.5773502691896258;\n\nconst mat2 cart2hex = mat2(1, 0, i3, 2. * i3);\nconst mat2 hex2cart = mat2(1, 0, -.5, .5 * sqrt3);\n\n#define PHI (1.618033988749895)\n#define TAU 6.283185307179586\n\nstruct TriPoints {\n\tvec2 a;\n    vec2 b;\n    vec2 c;\n    vec2 center;\n    vec2 ab;\n    vec2 bc;\n    vec2 ca;\n};\n\nTriPoints closestTriPoints(vec2 p) {    \n    vec2 pTri = cart2hex * p;\n    vec2 pi = floor(pTri);\n    vec2 pf = fract(pTri);\n    \n    float split1 = step(pf.y, pf.x);\n    float split2 = step(pf.x, pf.y);\n    \n    vec2 a = vec2(split1, 1);\n    vec2 b = vec2(1, split2);\n    vec2 c = vec2(0, 0);\n\n    a += pi;\n    b += pi;\n    c += pi;\n\n    a = hex2cart * a;\n    b = hex2cart * b;\n    c = hex2cart * c;\n    \n    vec2 center = (a + b + c) / 3.;\n    \n\tvec2 ab = (a + b) / 2.;\n    vec2 bc = (b + c) / 2.;\n    vec2 ca = (c + a) / 2.;\n\n    return TriPoints(a, b, c, center, ab, bc, ca);\n}\n\n\n// --------------------------------------------------------\n// Geodesic tiling\n// --------------------------------------------------------\n\nstruct TriPoints3D {\n\tvec3 a;\n    vec3 b;\n    vec3 c;\n\tvec3 center;\n    vec3 ab;\n    vec3 bc;\n    vec3 ca;\n};\n\nvec3 intersection(vec3 n, vec3 planeNormal, float planeOffset) {\n    float denominator = dot(planeNormal, n);\n    float t = (dot(vec3(0), planeNormal ) + planeOffset) / -denominator;\n    return n * t;\n}\n\n//// Edge length of an icosahedron with an inscribed sphere of radius of 1\n//float edgeLength = 1. / ((sqrt(3.) / 12.) * (3. + sqrt(5.)));\n//// Inner radius of the icosahedron's face\n//float faceRadius = (1./6.) * sqrt(3.) * edgeLength;\nfloat faceRadius = 0.3819660112501051;\n\n// 2D coordinates on the icosahedron face\nvec2 icosahedronFaceCoordinates(vec3 p) {\n    vec3 pn = normalize(p);\n    vec3 i = intersection(pn, facePlane, -1.);\n    return vec2(dot(i, uPlane), dot(i, vPlane));\n}\n\n// Project 2D icosahedron face coordinates onto a sphere\nvec3 faceToSphere(vec2 facePoint) {\n\treturn normalize(facePlane + (uPlane * facePoint.x) + (vPlane * facePoint.y));\n}\n\nTriPoints3D geodesicTriPoints(vec3 p, float subdivisions) {\n    // Get 2D cartesian coordiantes on that face\n    vec2 uv = icosahedronFaceCoordinates(p);\n    \n    // Get points on the nearest triangle tile\n\tfloat uvScale = subdivisions / faceRadius / 2.;\n    TriPoints points = closestTriPoints(uv * uvScale);\n    \n    // Project 2D triangle coordinates onto a sphere \n    vec3 a = faceToSphere(points.a / uvScale);\n    vec3 b = faceToSphere(points.b / uvScale);\n    vec3 c = faceToSphere(points.c / uvScale);\n    vec3 center = faceToSphere(points.center / uvScale);\n    vec3 ab = faceToSphere(points.ab / uvScale);\n    vec3 bc = faceToSphere(points.bc / uvScale);\n    vec3 ca = faceToSphere(points.ca / uvScale);\n    \n    return TriPoints3D(a, b, c, center, ab, bc, ca);\n}\n\n\n// --------------------------------------------------------\n// Spectrum colour palette\n// IQ https://www.shadertoy.com/view/ll2GD3\n// --------------------------------------------------------\n\nvec3 pal( in float t, in vec3 a, in vec3 b, in vec3 c, in vec3 d ) {\n    return a + b*cos( 6.28318*(c*t+d) );\n}\n\nvec3 spectrum(float n) {\n    return pal( n, vec3(0.5,0.5,0.5),vec3(0.5,0.5,0.5),vec3(1.0,1.0,1.0),vec3(0.0,0.33,0.67) );\n}\n\n\n// --------------------------------------------------------\n// Model/Camera Rotation\n// --------------------------------------------------------\n\nmat3 sphericalMatrix(float theta, float phi) {\n    float cx = cos(theta);\n    float cy = cos(phi);\n    float sx = sin(theta);\n    float sy = sin(phi);\n    return mat3(\n        cy, -sy * -sx, -sy * cx,\n        0, cx, sx,\n        sy, cy * -sx, cy * cx\n    );\n}\n\nmat3 mouseRotation(bool enable, vec2 xy) {\n    if (enable) {\n        vec2 mouse = iMouse.xy / iResolution.xy;\n\n        if (mouse.x != 0. && mouse.y != 0.) {\n            xy.x = mouse.x;\n            xy.y = mouse.y;\n        }\n    }\n    float rx, ry;\n    \n    rx = (xy.y + .5) * PI;\n    ry = (-xy.x) * 2. * PI;\n    \n    return sphericalMatrix(rx, ry);\n}\n\nmat3 modelRotation() {\n    mat3 m = mouseRotation(MOUSE_CONTROL==1, MODEL_ROTATION);\n    return m;\n}\n\nmat3 cameraRotation() {\n    mat3 m = mouseRotation(MOUSE_CONTROL==2, CAMERA_ROTATION);\n    return m;\n}\n\n\n// --------------------------------------------------------\n// Animation \n// --------------------------------------------------------\n\nconst float SCENE_DURATION = 6.;\nconst float CROSSFADE_DURATION = 2.;\n\nfloat time;\n\nstruct HexSpec {\n    float roundTop;\n    float roundCorner;\n\tfloat height;\n    float thickness;\n    float gap;    \n};\n    \nHexSpec newHexSpec(float subdivisions) {\n\treturn HexSpec(\n        .05 / subdivisions,\n        .1 / subdivisions,\n        2.,\n        2.,\n        .005\n    );\n}\n    \n// Animation 1\n    \nfloat animSubdivisions1() {\n\treturn mix(2.4, 3.4, cos(time * PI) * .5 + .5);\n}\n\nHexSpec animHex1(vec3 hexCenter, float subdivisions) {\n    HexSpec spec = newHexSpec(subdivisions);\n    \n    float offset = time * 3. * PI;\n    offset -= subdivisions;\n    float blend = dot(hexCenter, pca);\n    blend = cos(blend * 30. + offset) * .5 + .5;\n    spec.height = mix(1.75, 2., blend);\n\n    spec.thickness = spec.height;\n\n    return spec;\n}\n\n// Animation 2\n\nfloat animSubdivisions2() {\n    return mix(1., 2.3, sin(time * PI/2.) * .5 + .5);\n}\n\nHexSpec animHex2(vec3 hexCenter, float subdivisions) {\n    HexSpec spec = newHexSpec(subdivisions);\n    \n    float blend = hexCenter.y;\n    spec.height = mix(1.6, 2., sin(blend * 10. + time * PI) * .5 + .5);\n    \n    spec.roundTop = .02 / subdivisions;\n    spec.roundCorner = .09 / subdivisions;\n    spec.thickness = spec.roundTop * 4.;\n    spec.gap = .01;\n\n    return spec;\n}\n\n// Animation 3\n\nfloat animSubdivisions3() {\n\treturn 5.;\n}\n\nHexSpec animHex3(vec3 hexCenter, float subdivisions) {\n    HexSpec spec = newHexSpec(subdivisions);\n    \n    float blend = acos(dot(hexCenter, pab)) * 10.;\n    blend = cos(blend + time * PI) * .5 + .5;\n    spec.gap = mix(.01, .4, blend) / subdivisions;\n\n    spec.thickness = spec.roundTop * 2.;\n\n\treturn spec;\n}\n\n// Transition between animations\n\nfloat sineInOut(float t) {\n  return -0.5 * (cos(PI * t) - 1.0);\n}\n\nfloat transitionValues(float a, float b, float c) {\n    #ifdef LOOP\n        #if LOOP == 1\n            return a;\n        #endif\n        #if LOOP == 2\n            return b;\n        #endif\n        #if LOOP == 3\n            return c;\n        #endif\n    #endif\n    float t = time / SCENE_DURATION;\n    float scene = floor(mod(t, 3.));\n    float blend = fract(t);\n    float delay = (SCENE_DURATION - CROSSFADE_DURATION) / SCENE_DURATION;\n    blend = max(blend - delay, 0.) / (1. - delay);\n    blend = sineInOut(blend);\n    float ab = mix(a, b, blend);\n    float bc = mix(b, c, blend);\n    float cd = mix(c, a, blend);\n    float result = mix(ab, bc, min(scene, 1.));\n    result = mix(result, cd, max(scene - 1., 0.));\n    return result;\n}\n \nHexSpec transitionHexSpecs(HexSpec a, HexSpec b, HexSpec c) {\n    float roundTop = transitionValues(a.roundTop, b.roundTop, c.roundTop);\n    float roundCorner = transitionValues(a.roundCorner, b.roundCorner, c.roundCorner);\n\tfloat height = transitionValues(a.height, b.height, c.height);\n    float thickness = transitionValues(a.thickness, b.thickness, c.thickness);\n    float gap = transitionValues(a.gap, b.gap, c.gap);\n\treturn HexSpec(roundTop, roundCorner, height, thickness, gap);\n}\n\n\n// --------------------------------------------------------\n// Modelling \n// --------------------------------------------------------\n\nconst vec3 FACE_COLOR = vec3(.9,.9,1.);\nconst vec3 BACK_COLOR = vec3(.1,.1,.15);\nconst vec3 BACKGROUND_COLOR = vec3(.0, .005, .03);\n\nstruct Model {\n    float dist;\n    vec3 albedo;\n    float glow;\n};\n\nModel hexModel(\n    vec3 p,\n    vec3 hexCenter,\n    vec3 edgeA,\n    vec3 edgeB,\n    HexSpec spec\n) {\n    float d;\n\n    float edgeADist = dot(p, edgeA) + spec.gap;\n    float edgeBDist = dot(p, edgeB) - spec.gap;\n    float edgeDist = smax(edgeADist, -edgeBDist, spec.roundCorner);\n\n    float outerDist = length(p) - spec.height;\n    d = smax(edgeDist, outerDist, spec.roundTop);\n\n    float innerDist = length(p) - spec.height + spec.thickness;\n    d = smax(d, -innerDist, spec.roundTop);\n    \n    vec3 color;\n\n    float faceBlend = (spec.height - length(p)) / spec.thickness;\n    faceBlend = clamp(faceBlend, 0., 1.);\n    color = mix(FACE_COLOR, BACK_COLOR, step(.5, faceBlend));\n    \n    vec3 edgeColor = spectrum(dot(hexCenter, pca) * 5. + length(p) + .8);    \n\tfloat edgeBlend = smoothstep(-.04, -.005, edgeDist);\n    color = mix(color, edgeColor, edgeBlend); \n\n    return Model(d, color, edgeBlend);\n}\n\n// checks to see which intersection is closer\nModel opU( Model m1, Model m2 ){\n    if (m1.dist < m2.dist) {\n        return m1;\n    } else {\n        return m2;\n    }\n}\n\nModel geodesicModel(vec3 p) {\n\n    pModIcosahedron(p);\n    \n    float subdivisions = transitionValues(\n        animSubdivisions1(),\n        animSubdivisions2(),\n        animSubdivisions3()\n   \t);\n\tTriPoints3D points = geodesicTriPoints(p, subdivisions);\n        \n\tvec3 edgeAB = normalize(cross(points.center, points.ab));\n\tvec3 edgeBC = normalize(cross(points.center, points.bc));\n    vec3 edgeCA = normalize(cross(points.center, points.ca));\n    \n    Model model, part;\n    HexSpec spec;\n\n\tspec = transitionHexSpecs(\n        animHex1(points.b, subdivisions),\n        animHex2(points.b, subdivisions),\n        animHex3(points.b, subdivisions)\n    );\n    part = hexModel(p, points.b, edgeAB, edgeBC, spec);\n    model = part;\n\n\tspec = transitionHexSpecs(\n        animHex1(points.c, subdivisions),\n        animHex2(points.c, subdivisions),\n        animHex3(points.c, subdivisions)\n    );\n    part = hexModel(p, points.c, edgeBC, edgeCA, spec);\n    model = opU(model, part);\n    \n\tspec = transitionHexSpecs(\n        animHex1(points.a, subdivisions),\n        animHex2(points.a, subdivisions),\n        animHex3(points.a, subdivisions)\n    );\n    part = hexModel(p, points.a, edgeCA, edgeAB, spec);\n    model = opU(model, part);\n    \n\treturn model;\n}\n\nModel map( vec3 p ){\n    mat3 m = modelRotation();\n    p *= m;  \n    #ifndef LOOP\n    \tpR(p.xz, time * PI/16.);\n    #endif\n    Model model = geodesicModel(p);\n    return model;\n}\n\n// --------------------------------------------------------\n// LIGHTING\n// Adapted from IQ https://www.shadertoy.com/view/Xds3zN\n// --------------------------------------------------------\n\nvec3 doLighting(Model model, vec3 pos, vec3 nor, vec3 ref, vec3 rd) {\n    vec3 lightPos = normalize(vec3(.5,.5,-1.));\n    vec3 backLightPos = normalize(vec3(-.5,-.3,1));\n    vec3 ambientPos = vec3(0,1,0);\n    \n    vec3  lig = lightPos;\n    float amb = clamp((dot(nor, ambientPos) + 1.) / 2., 0., 1.);\n    float dif = clamp( dot( nor, lig ), 0.0, 1.0 );\n    float bac = pow(clamp(dot(nor, backLightPos), 0., 1.), 1.5);\n    float fre = pow( clamp(1.0+dot(nor,rd),0.0,1.0), 2.0 );\n    \n    vec3 lin = vec3(0.0);\n    lin += 1.20 * dif * vec3(.9);\n    lin += 0.80 * amb * vec3(.5, .7, .8);\n    lin += 0.30 * bac * vec3(.25);\n    lin += 0.20 * fre * vec3(1);\n    \n    vec3 albedo = model.albedo;\n    vec3 col = mix(albedo * lin, albedo, model.glow);    \n\n    return col;\n}   \n\n\n// --------------------------------------------------------\n// Ray Marching\n// Adapted from cabbibo https://www.shadertoy.com/view/Xl2XWt\n// --------------------------------------------------------\n\nconst float MAX_TRACE_DISTANCE = 8.; // max trace distance\nconst float INTERSECTION_PRECISION = .001; // precision of the intersection\nconst int NUM_OF_TRACE_STEPS = 100;\nconst float FUDGE_FACTOR = .9; // Default is 1, reduce to fix overshoots\n\nstruct CastRay {\n    vec3 origin;\n    vec3 direction;\n};\n\nstruct Ray {\n    vec3 origin;\n    vec3 direction;\n    float len;\n};\n\nstruct Hit {\n    Ray ray;\n    Model model;\n    vec3 pos;\n    bool isBackground;\n    vec3 normal;\n    vec3 color;\n};\n\nvec3 calcNormal( in vec3 pos ){\n    vec3 eps = vec3( 0.001, 0.0, 0.0 );\n    vec3 nor = vec3(\n        map(pos+eps.xyy).dist - map(pos-eps.xyy).dist,\n        map(pos+eps.yxy).dist - map(pos-eps.yxy).dist,\n        map(pos+eps.yyx).dist - map(pos-eps.yyx).dist );\n    return normalize(nor);\n}\n    \nHit raymarch(CastRay castRay){\n\n    float currentDist = INTERSECTION_PRECISION * 2.0;\n    Model model;\n    \n    Ray ray = Ray(castRay.origin, castRay.direction, 0.);\n\n    for( int i=0; i< NUM_OF_TRACE_STEPS ; i++ ){\n        if (currentDist < INTERSECTION_PRECISION || ray.len > MAX_TRACE_DISTANCE) {\n            break;\n        }\n        model = map(ray.origin + ray.direction * ray.len);\n        currentDist = model.dist;\n        ray.len += currentDist * FUDGE_FACTOR;\n    }\n    \n    bool isBackground = false;\n    vec3 pos = vec3(0);\n    vec3 normal = vec3(0);\n    vec3 color = vec3(0);\n    \n    if (ray.len > MAX_TRACE_DISTANCE) {\n        isBackground = true;\n    } else {\n        pos = ray.origin + ray.direction * ray.len;\n        normal = calcNormal(pos);\n    }\n\n    return Hit(ray, model, pos, isBackground, normal, color);\n}\n\n\n// --------------------------------------------------------\n// Rendering\n// --------------------------------------------------------\n\nvoid shadeSurface(inout Hit hit){\n    \n    vec3 color = BACKGROUND_COLOR;\n    \n    if (hit.isBackground) {\n        hit.color = color;\n        return;\n    }\n\n    vec3 ref = reflect(hit.ray.direction, hit.normal);\n\n    #ifdef DEBUG\n        color = hit.normal * 0.5 + 0.5;\n    #else \n        color = doLighting(\n            hit.model,\n            hit.pos,\n            hit.normal,\n            ref,\n            hit.ray.direction\n        );\n    #endif\n\n    hit.color = color;\n}\n\nvec3 render(Hit hit){\n    shadeSurface(hit);\n\treturn hit.color;\n}\n\n\n// --------------------------------------------------------\n// Camera\n// https://www.shadertoy.com/view/Xl2XWt\n// --------------------------------------------------------\n\nmat3 calcLookAtMatrix( in vec3 ro, in vec3 ta, in float roll )\n{\n    vec3 ww = normalize( ta - ro );\n    vec3 uu = normalize( cross(ww,vec3(sin(roll),cos(roll),0.0) ) );\n    vec3 vv = normalize( cross(uu,ww));\n    return mat3( uu, vv, ww );\n}\n\nvoid doCamera(out vec3 camPos, out vec3 camTar, out float camRoll, in float time, in vec2 mouse) {\n    float dist = 5.5;\n    camRoll = 0.;\n    camTar = vec3(0,0,0);\n    camPos = vec3(0,0,-dist);\n    camPos *= cameraRotation();\n    camPos += camTar;\n}\n\n\n// --------------------------------------------------------\n// Gamma\n// https://www.shadertoy.com/view/Xds3zN\n// --------------------------------------------------------\n\nconst float GAMMA = 2.2;\n\nvec3 gamma(vec3 color, float g) {\n    return pow(color, vec3(g));\n}\n\nvec3 linearToScreen(vec3 linearRGB) {\n    return gamma(linearRGB, 1.0 / GAMMA);\n}\n\nvoid main()\n{\n    time = iTime;\n\n    #ifdef LOOP\n        #if LOOP == 1\n            time = mod(time, 2.);   \n        #endif\n        #if LOOP == 2\n            time = mod(time, 4.);   \n        #endif\n        #if LOOP == 3\n            time = mod(time, 2.);\n    \t#endif\n    #endif\n    \n    initIcosahedron();\n    \n    vec2 p =  2.*v_st.xy - vec2(1., 1.);//(-iResolution.xy + 2.0*fragCoord.xy)/iResolution.y;\n    vec2 m = iMouse.xy / iResolution.xy;\n\n    vec3 camPos = vec3( 0., 0., 2.);\n    vec3 camTar = vec3( 0. , 0. , 0. );\n    float camRoll = 0.;\n    \n    // camera movement\n    doCamera(camPos, camTar, camRoll, iTime, m);\n    \n    // camera matrix\n    mat3 camMat = calcLookAtMatrix( camPos, camTar, camRoll );  // 0.0 is the camera roll\n    \n    // create view ray\n    vec3 rd = normalize( camMat * vec3(p.xy,2.0) ); // 2.0 is the lens length\n    \n    Hit hit = raymarch(CastRay(camPos, rd));\n\n    vec3 color = render(hit);\n    \n    #ifndef DEBUG\n        color = linearToScreen(color);\n    #endif\n\n    fragColor = vec4(color,1.0);\n}\n",
    'rgba(255, 0, 0, 1)',
    '_mouseVal'
  ]
  return (a0_0x270b = function () {
    return t
  })()
}
;(function (t, e) {
  const i = a0_0x3b79,
    s = a0_0x270b()
  for (;;)
    try {
      if (
        818706 ===
        -parseInt(i(1336)) / 1 +
          (-parseInt(i(751)) / 2) * (parseInt(i(1756)) / 3) +
          -parseInt(i(1894)) / 4 +
          (parseInt(i(1335)) / 5) * (parseInt(i(2329)) / 6) +
          (-parseInt(i(1005)) / 7) * (-parseInt(i(910)) / 8) +
          (parseInt(i(2585)) / 9) * (parseInt(i(511)) / 10) +
          (-parseInt(i(1175)) / 11) * (-parseInt(i(1440)) / 12)
      )
        break
      s.push(s.shift())
    } catch (t) {
      s.push(s.shift())
    }
})(),
  (function () {
    'use strict'
    const t = a0_0x3b79
    try {
      if ('undefined' != typeof document) {
        var e = document.createElement(t(1679))
        e[t(1621)](document.createTextNode(t(754))), document[t(131)][t(1621)](e)
      }
    } catch (e) {
      console[t(1792)](t(2511), e)
    }
  })(),
  (function (t, e) {
    const i = a0_0x3b79
    i(1784) == typeof exports && 'undefined' != typeof module
      ? (module[i(408)] = e())
      : 'function' == typeof define && define[i(2219)]
        ? define(e)
        : ((t = 'undefined' != typeof globalThis ? globalThis : t || self)[i(2201)] = e())
  })(this, function () {
    'use strict'
    const t = a0_0x3b79,
      e = {}
    e[t(922)] = t(256)
    const i = 100,
      s = Math.PI / 2,
      n = 3e-6,
      r = 2 * Math.PI,
      h = Object[t(1731)](
        Object[t(1399)](
          {
            __proto__: null,
            ERADIUS: 6378137,
            FITTING_COUNT: i,
            HALF_PI: s,
            PRADIUS: 6356725,
            TWO_PI: r,
            ZERO_TOLERANCE: n
          },
          Symbol[t(1540)],
          e
        )
      ),
      l = {}
    l[t(922)] = 'Module'
    const c = {}
    ;(c[t(1362)] = t(1362)),
      (c[t(2051)] = t(2051)),
      (c.selectedGraphic = t(1356)),
      (c[t(786)] = t(786)),
      (c[t(2279)] = t(2279)),
      (c[t(182)] = t(182)),
      (c[t(289)] = 'drawCancel'),
      (c[t(308)] = 'editStart'),
      (c[t(957)] = 'editEnd'),
      (c[t(355)] = t(355)),
      (c[t(1518)] = 'routeStart'),
      (c[t(2414)] = t(2414)),
      (c[t(1205)] = t(1205)),
      (c[t(1017)] = t(1017)),
      (c[t(2084)] = t(2084)),
      (c.matrixMove = t(1086)),
      (c[t(244)] = t(244)),
      (c[t(2070)] = t(2070)),
      (c.selectedObjectChanged = t(2158))
    const m = Object[t(1731)](
        Object[t(1399)](
          {
            __proto__: null,
            fromAssetId: function (t) {
              return new Promise((e, i) => {
                const s = a0_0x3b79
                let n = s(439) + t + s(315) + xt3d[s(627)]
                fetch(n)
                  [s(687)]((t) => t.json())
                  [s(687)]((t) => {
                    const i = s,
                      n = {}
                    ;(n[i(2524)] = xt3d.resourceToken),
                      200 == t[i(2421)] &&
                        e(new Cesium[i(1283)]({ url: t[i(825)][i(2171)], queryParameters: n }))
                  })
                  [s(472)]((t) => {
                    i(t)
                  })
              })
            }
          },
          Symbol[t(1540)],
          l
        )
      ),
      p = c
    class d {
      constructor() {
        this._cache = {}
      }
      _on(e, i, s) {
        const n = t
        if (!p[e]) throw new Cesium.DeveloperError(n(1204))
        let o,
          r = this[n(1219)](e)
        return (
          r || ((r = new Cesium[n(525)]()), (this[n(1230)][e] = r)),
          r && i && (o = r[n(1973)](i, s || this)),
          o
        )
      }
      [t(758)](e, i, s) {
        const n = t
        let o = this.getEvent(e),
          r = !1
        return o && i && (r = o[n(704)](i, s || this)), r
      }
      _fire(e, i) {
        const s = t
        let n = this.getEvent(e)
        n && n[s(750)](i)
      }
      on(e, i, s) {
        return this[t(1463)](e, i, s)
      }
      [t(534)](e, i, s) {
        let n = this[t(1463)](
          e,
          (t) => {
            i(t), n && n()
          },
          s
        )
      }
      [t(686)](e, i, s) {
        return this[t(758)](e, i, s)
      }
      [t(450)](t, e) {
        this._fire(t, e)
      }
      [t(1219)](e) {
        return this[t(1230)][e] || void 0
      }
    }
    const f = {}
    ;(f.color = t(1070)),
      (f[t(259)] = t(259)),
      (f[t(899)] = t(899)),
      (f.polyWater = 'polyWater'),
      (f[t(1911)] = 'polyElevationCountour'),
      (f[t(1103)] = t(1103)),
      (f[t(986)] = t(986)),
      (f[t(1915)] = t(1915)),
      (f[t(326)] = t(326)),
      (f.circleImageDiffuse = t(1294)),
      (f.circleImageRotate = 'circleImageRotate'),
      (f.circleLightRing = 'circleLightRing'),
      (f[t(1064)] = 'circleScan_1'),
      (f[t(1170)] = 'circleScan_2'),
      (f[t(479)] = t(479)),
      (f.circleScan_4 = t(150)),
      (f[t(1609)] = 'circleSpiral'),
      (f.circleWave = 'circleWave'),
      (f[t(1845)] = t(1845)),
      (f[t(1422)] = 'sphereScan'),
      (f[t(1677)] = t(1677)),
      (f[t(1039)] = t(1039)),
      (f[t(2344)] = t(2344)),
      (f[t(533)] = 'wallImageFlow'),
      (f[t(1376)] = t(1376)),
      (f.wallCool = t(497))
    const C = f
    class v {
      constructor(e) {
        const i = t
        ;(this[i(568)] = e), (this[i(1085)] = new d()), (this[i(1330)] = !1)
      }
      isDestroyed() {
        return this[t(1330)]
      }
      get [t(1888)]() {
        return this._event
      }
      on(t, e, i) {
        return this._event.on(t, e, i)
      }
      [t(534)](t, e, i) {
        this._event.once(t, e, i)
      }
      [t(450)](e, i, s) {
        const n = t
        return this._event[n(450)](e, i, s)
      }
      [t(686)](e, i, s) {
        const n = t
        return this[n(1085)][n(686)](e, i, s)
      }
      [t(2039)]() {
        const e = t
        this[e(1712)]()
        let i = ((t) => {
          const i = e
          let s = new Set(),
            n = t
          do {
            Object.getOwnPropertyNames(n)[i(2512)]((t) => s[i(1861)](t))
          } while ((n = Object[i(385)](n)))
          return [...s.keys()].filter((e) => 'function' == typeof t[e])
        })(this)
        const s = () => {
          const t = e
          throw new Cesium[t(2360)](t(1186))
        }
        for (let t = 0; t < i[e(277)]; t++) this[i[t]] = s
        for (const t in this) e(837) != typeof this[t] && (this[t] = null)
        this[e(1330)] = !0
      }
      [t(1712)]() {}
    }
    function _(t, e) {
      return new Promise((e, i) => {
        const s = a0_0x3b79
        let n = document[s(1945)](s(1711))
        ;(n[s(1640)] = s(2554)), (n[s(1717)] = '.' + t)
        let o = document[s(1546)](s(156))
        o[s(1674)](s(1071), !0, !1, window, 0, 0, 0, 0, 0, !1, !1, !1, !1, 0, null),
          n.dispatchEvent(o),
          n[s(1973)](s(2015), (i) => {
            const n = s
            let o = i[n(662)][n(2465)][0],
              r = new FileReader()
            ;(r[n(1701)] = (t) => {
              const i = n
              e(t[i(662)][i(551)])
            }),
              n(2241) == t ? r[n(233)](o) : r[n(874)](o)
          })
      })
    }
    function g(e, i) {
      const s = t
      e = JSON.stringify(e)
      const n = {}
      n[s(1640)] = s(761)
      let o = new Blob([e], n),
        r = document[s(1546)](s(156)),
        a = document[s(1945)]('a')
      ;(a[s(883)] = i + s(756)),
        (a[s(297)] = window[s(983)][s(705)](o)),
        (a[s(1783)][s(409)] = [s(761), a.download, a[s(297)]][s(2011)](':')),
        r[s(1674)](s(1071), !0, !1, window, 0, 0, 0, 0, 0, !1, !1, !1, !1, 0, null),
        a[s(1918)](r)
    }
    const y = {}
    ;(y[t(900)] = null), (y.downloadToLocalFile = g), (y[t(1780)] = _)
    const w = {}
    w[t(922)] = t(256)
    const x = Object[t(1731)](Object[t(1399)](y, Symbol[t(1540)], w))
    class b extends v {
      constructor(e = {}) {
        const i = t
        super(e),
          (this[i(1570)] = Cesium[i(1960)](e.id, Cesium[i(721)]())),
          (this[i(528)] = ''),
          (this._show = Cesium[i(1960)](e[i(1482)], !0)),
          (this._name = Cesium[i(1960)](e[i(1916)], i(1905))),
          (this._remarks = Cesium[i(1960)](e[i(755)], '')),
          (this._style = Cesium[i(1960)](e.style, {})),
          (this._drawEnd = !0),
          (this[i(2520)] = !1),
          (this[i(1496)] = !1),
          (this[i(324)] = 0),
          (this[i(350)] = 0),
          (this[i(314)] = null),
          (this._coordinates = []),
          (this[i(1912)] = null),
          (this[i(2238)] = null),
          (this[i(2469)] = []),
          (this._cartesian3Array = []),
          (this[i(216)] = Cesium[i(1960)](e[i(924)], {})),
          (this._boundingSphere = new Cesium[i(1242)](Cesium[i(310)][i(1738)], 1))
      }
      get id() {
        return this[t(1570)]
      }
      get [t(1059)]() {
        return this[t(1115)]
      }
      get [t(1158)]() {
        return this._graphicType
      }
      get typeName() {
        return this[t(528)]
      }
      get [t(1916)]() {
        return this._name
      }
      set [t(1916)](e) {
        this[t(1165)] = e
      }
      get remarks() {
        return this[t(1478)]
      }
      set remarks(e) {
        this[t(1478)] = e
      }
      get [t(1482)]() {
        return this[t(1144)]
      }
      set [t(1482)](e) {
        const i = t
        ;(this[i(1144)] = e), this[i(415)](e)
      }
      get [t(468)]() {
        return this[t(1496)]
      }
      set [t(468)](e) {
        const i = t
        ;(this[i(1496)] = e), this[i(2231)]()
      }
      get [t(1804)]() {
        return this[t(2520)]
      }
      set [t(1804)](t) {
        ;(this._editMode = t), this._setEditMode()
      }
      get userData() {
        return this[t(216)]
      }
      set [t(924)](e) {
        this[t(216)] = e
      }
      get [t(2251)]() {
        return this._position
      }
      get [t(1131)]() {
        return this[t(2238)]
      }
      get positions() {
        return this[t(2469)]
      }
      get [t(972)]() {
        return this[t(1378)]
      }
      get [t(153)]() {
        return this._geometryType
      }
      get [t(2365)]() {
        return this[t(794)]
      }
      get [t(293)]() {
        return this[t(324)]
      }
      get fixPointCount() {
        return this[t(350)]
      }
      get [t(1679)]() {
        return this[t(343)]
      }
      set style(e) {
        const i = t
        ;(this[i(343)] = e), this[i(2231)]()
      }
      get [t(1611)]() {
        return this[t(1024)]
      }
      [t(200)]() {}
      [t(415)](t) {}
      setPosition(t) {}
      [t(2515)](t) {}
      [t(2231)]() {}
      [t(1930)]() {
        const e = t
        let i = this[e(1114)]()
        const s = {}
        return (
          (s.type = this[e(153)]),
          (s[e(2365)] = this.coordinates),
          {
            type: e(1829),
            geometry: s,
            properties: {
              id: this[e(1570)],
              name: this[e(1165)],
              remarks: this[e(1478)],
              graphicClassType: this[e(1115)],
              graphicType: this[e(2336)],
              style: this._getStyle(),
              userData: this[e(216)],
              ...i
            }
          }
        )
      }
      [t(393)]() {
        return this[t(343)]
      }
      downloadToLocalFile(e) {
        g(this[t(1930)](), e)
      }
      [t(1114)]() {}
      [t(1545)](t) {}
      [t(2389)](t) {}
    }
    const S = {}
    ;(S[t(1107)] = 'filletLabel'),
      (S[t(2052)] = t(2052)),
      (S.brightenDiv = 'brightenDiv'),
      (S[t(1810)] = t(1810)),
      (S[t(2352)] = t(2352)),
      (S[t(835)] = 'liquidfill'),
      (S[t(1944)] = t(1944)),
      (S[t(665)] = t(665)),
      (S[t(1093)] = t(1093)),
      (S[t(807)] = t(807)),
      (S.divText = 'divText'),
      (S[t(1384)] = t(1384)),
      (S.simpleMarker = t(1923)),
      (S[t(143)] = t(143)),
      (S[t(784)] = t(784)),
      (S.gifMarker = t(1582)),
      (S.alertMarker = t(1991)),
      (S[t(1554)] = 'poiMarker'),
      (S.textMarker = t(1770)),
      (S[t(1396)] = t(1396)),
      (S[t(1351)] = t(1351)),
      (S[t(975)] = 'curvyLine'),
      (S[t(1596)] = t(1596)),
      (S[t(1121)] = 'lineVolume'),
      (S.polygon = t(2087)),
      (S[t(2584)] = t(2584)),
      (S[t(1074)] = 'sphere'),
      (S[t(2048)] = 'planeHeat'),
      (S[t(2102)] = t(2102)),
      (S.gatheringPlace = t(181)),
      (S[t(2286)] = t(2286)),
      (S[t(2568)] = t(2568)),
      (S.attackArrow = t(941)),
      (S[t(1394)] = 'tailedAttackArrow'),
      (S.squadCombat = t(2197)),
      (S[t(1668)] = t(1668)),
      (S[t(2433)] = t(2433)),
      (S[t(1181)] = 'fineArrow'),
      (S[t(682)] = t(682)),
      (S[t(2098)] = 'videoFusion'),
      (S[t(608)] = 'videoProject'),
      (S[t(738)] = t(738)),
      (S[t(1921)] = 'regularWall'),
      (S[t(549)] = 'diffuseWall'),
      (S[t(1421)] = t(1421)),
      (S.spotLight = t(1989)),
      (S[t(2316)] = t(2316)),
      (S.model = t(441)),
      (S[t(265)] = t(265)),
      (S.arrowAxis = 'arrowAxis'),
      (S[t(2500)] = t(2500)),
      (S[t(1946)] = t(1946)),
      (S[t(726)] = t(726)),
      (S[t(649)] = t(649)),
      (S[t(1625)] = t(1625)),
      (S[t(1270)] = t(1270)),
      (S[t(427)] = 'beamRadar'),
      (S[t(668)] = t(668)),
      (S.cylinderRadar = t(1650)),
      (S[t(1458)] = t(1458)),
      (S.glowCylinder = t(203)),
      (S[t(868)] = t(868)),
      (S[t(2513)] = 'shaderEffetBase'),
      (S.shaderEffetRealFlame = t(1581)),
      (S[t(2139)] = t(2139)),
      (S.shaderEffetMagicBall = t(651)),
      (S.shaderEffetMagicBox = t(1534)),
      (S[t(1147)] = 'shaderEffetFlameCloud'),
      (S.shaderEffetCoolBall = t(1820)),
      (S.shaderEffetStellarChain = t(928)),
      (S[t(870)] = t(870)),
      (S[t(1956)] = t(1956)),
      (S[t(1140)] = t(1140)),
      (S[t(2175)] = t(2175)),
      (S.shaderEffetNeonLight = t(797)),
      (S.shaderEffetGlowBox = t(2014)),
      (S[t(1322)] = 'shaderEffetGlowDiamond'),
      (S[t(2065)] = 'shaderEffetExplosion'),
      (S[t(1301)] = 'shaderEffetGradientRing'),
      (S.shaderEffetColorfulPoint = t(2580)),
      (S[t(2188)] = t(2188)),
      (S[t(1183)] = t(1183)),
      (S[t(1629)] = t(1629)),
      (S[t(2136)] = t(2136)),
      (S[t(1171)] = t(1171)),
      (S[t(243)] = 'shaderEffetPetal'),
      (S[t(1553)] = t(1553)),
      (S[t(1910)] = t(1910)),
      (S.shaderEffetConstellationChain = t(482)),
      (S[t(628)] = 'shaderEffetNeonPoint'),
      (S[t(2498)] = t(2498)),
      (S[t(1155)] = t(1155)),
      (S[t(2167)] = 'shaderEffetShield'),
      (S[t(1213)] = 'waterPrimitive'),
      (S[t(2417)] = t(2417)),
      (S.refractionWater = t(653)),
      (S[t(2040)] = t(2040))
    const P = t(1383),
      M = t(209),
      A = t(1121),
      T = t(1424),
      E = t(2425),
      z = t(465),
      D = t(1734),
      I = t(325),
      k = t(165),
      F = t(596),
      R = t(1524),
      L = 'SphereGraphic',
      O = t(591),
      B = t(1508),
      V = t(572),
      N = t(1536),
      H = t(536),
      G = t(2093),
      W = t(2062),
      U = t(1529),
      j = t(958),
      q = t(740),
      Y = t(1366),
      X = t(1224),
      Q = S
    function Z(e, i, s, n) {
      const o = t
      i = Cesium[o(2066)][o(1902)](i)
      let r = Cesium[o(2066)].clone(i, new Cesium.Matrix4()),
        a = e[o(1391)][o(1236)],
        h = Cesium[o(2066)][o(2483)](r, a, new Cesium[o(310)]())
      Cesium[o(310)].normalize(h, h)
      let l = Cesium[o(310)][o(1676)](s, a, new Cesium[o(310)]()),
        c = Cesium[o(310)][o(1676)](n, a, new Cesium[o(310)]()),
        u = new Cesium[o(2007)](a, l),
        m = new Cesium[o(2007)](a, c),
        p = new Cesium.Cartesian3()
      const d = new Cesium.Cartesian3(),
        f = Cesium[o(310)][o(1246)]
      Cesium[o(310)][o(2196)](h, f, d),
        Cesium[o(310)][o(2196)](f, d, p),
        Cesium[o(310)][o(379)](p, p)
      const C = new Cesium[o(2314)](p, 0)
      Cesium[o(2314)][o(1863)](C, i, C)
      const v = Cesium[o(2525)][o(914)](u, C),
        _ = Cesium[o(2525)].rayPlane(m, C)
      if (!Cesium[o(2330)](v) || !Cesium.defined(_)) return
      const g = Cesium[o(2066)][o(2483)](r, v, new Cesium[o(310)]()),
        y = Cesium[o(2066)][o(2483)](r, _, new Cesium[o(310)]()),
        w = new Cesium[o(310)]()
      return (
        Cesium[o(310)][o(1676)](y, g, w),
        (w.x = w.y = 0),
        Cesium[o(2066)].fromTranslation(w, new Cesium[o(2066)]())
      )
    }
    function K(e, i) {
      const s = t,
        n = e.scene
      let o = n[s(1125)](i)
      if (!o) {
        const t = new Cesium.Cartesian2()
        ;(t.x = i.x),
          (t.y = n[s(1493)][s(281)] - i.y),
          (o = Cesium.SceneTransforms[s(516)](n, t, new Cesium[s(194)]())),
          (o = Cesium[s(2176)][s(973)](n, t, 0.1, new Cesium[s(310)]()))
      }
      return o
    }
    function J(e, i, s) {
      const n = t
      let o = Cesium[n(310)].magnitude(s),
        r = new Cesium[n(1148)](o, o, o)
      return e[n(696)].camera[n(1816)](i, r)
    }
    function $(e, i) {
      const s = t
      let n = new Cesium[s(310)](),
        o = Cesium[s(2058)][s(1648)](e)
      return (
        Cesium.Matrix4[s(814)](o, o),
        Cesium[s(2066)][s(2483)](o, i, n),
        Cesium[s(310)][s(379)](n, n),
        Cesium[s(475)][s(363)](Math[s(2003)](n.x, n.y))
      )
    }
    class tt extends b {
      constructor(e) {
        const i = t
        super(e),
          (this[i(528)] = i(1347)),
          (this[i(314)] = i(1865)),
          (this._graphicClassType = W),
          (this[i(2336)] = Q[i(649)]),
          (this[i(350)] = 1),
          (this._style[i(1981)] = Cesium[i(1960)](this[i(343)][i(1981)], 100)),
          (this[i(343)].hAngle = Cesium[i(1960)](this[i(343)][i(2519)], 22)),
          (this[i(343)][i(2493)] = Cesium[i(1960)](this[i(343)][i(2493)], 22)),
          (this._style.heading = Cesium[i(1960)](this._style.heading, 0)),
          (this._style[i(711)] = Cesium.defaultValue(this[i(343)][i(711)], 0)),
          (this[i(343)][i(1143)] = Cesium[i(1960)](this[i(343)].roll, 0)),
          (this[i(343)].color = Cesium[i(1960)](this._style[i(1070)], i(581))),
          (this._style[i(506)] = Cesium[i(1960)](this[i(343)].selectedColor, i(1544))),
          (this[i(343)][i(2248)] = Cesium.defaultValue(this[i(343)][i(2248)], !0)),
          (this[i(343)][i(1018)] = Cesium[i(1960)](this[i(343)][i(1018)], i(1320))),
          (this._style.topShow = Cesium[i(1960)](this[i(343)][i(948)], !0)),
          (this[i(343)][i(1105)] = Cesium[i(1960)](this[i(343)][i(1105)], !0)),
          (this._style[i(1801)] = Cesium.defaultValue(this[i(343)][i(1801)], 2)),
          (this._drawCommands = [])
        let s = Cesium[i(1960)](e.position, [111, 28, 0])
        ;(this[i(2238)] = Cesium[i(310)].fromDegrees(s[0], s[1], s[2])),
          (this[i(1912)] = s),
          (this[i(1636)] = Cesium[i(2066)][i(1902)](Cesium[i(2066)].IDENTITY)),
          (this[i(1150)] = new Cesium[i(1472)]()),
          (this[i(1359)] = new Cesium[i(310)]()),
          (this._scale = new Cesium[i(310)](1, 1, 1)),
          (this[i(967)] = new Cesium[i(2066)]()),
          (this[i(1897)] = !0),
          (this[i(1024)] = new Cesium[i(1242)](this[i(2238)], 1.5 * this._style[i(1981)]))
      }
      [t(1253)](e) {
        const i = t
        ;(this[i(1912)] = e),
          (this[i(2238)] = Cesium[i(310)][i(667)](e[0], e[1], e[2])),
          (this[i(794)] = this[i(1912)]),
          (this[i(1897)] = !0)
      }
      [t(1601)](e, i) {
        const s = t
        if (!e) return
        let n = e,
          o = Cesium[s(310)][s(667)](n[0], n[1], n[2]),
          r = o[s(1902)](),
          a = Cesium[s(310)].subtract(this[s(2238)], o, o)
        a = Cesium[s(310)].normalize(a, a)
        let h = Cesium[s(310)][s(379)](this[s(2238)], new Cesium[s(310)]()),
          l = $(this._cartesian3, r),
          c = this._style
        i && (c[s(1981)] = Cesium.Cartesian3[s(1849)](this[s(2238)], r)),
          (c[s(1198)] = l - 90),
          (c[s(1143)] = 0),
          (c[s(711)] = Cesium.Math[s(363)](Cesium[s(310)][s(354)](a, h)) - 180),
          (this[s(1679)] = c)
      }
      _setStyle() {
        const e = t
        ;(this[e(1897)] = !0),
          (this[e(1024)] = new Cesium[e(1242)](this[e(2238)], this[e(343)][e(1981)]))
      }
      _computeVAO(e) {
        const i = t
        let s = this._style[i(1981)],
          n = e[i(1015)] / 2,
          o = e[i(1866)] / 2,
          r = [],
          a = [],
          h = [],
          l = [],
          c = [],
          u = [],
          m = [],
          p = [],
          d = new Cesium.Cartesian3(-n, -o, s),
          f = new Cesium[i(310)](n, -o, s),
          C = new Cesium.Cartesian3(-n, o, s),
          v = new Cesium[i(310)](n, o, s)
        m[i(2553)](0, 0, 0),
          m.push(d.x, d.y, d.z),
          m[i(2553)](C.x, C.y, C.z),
          m[i(2553)](v.x, v.y, v.z),
          m[i(2553)](f.x, f.y, f.z),
          c[i(2553)](0, 1, 2),
          c[i(2553)](0, 2, 3),
          c[i(2553)](0, 3, 4),
          c[i(2553)](0, 4, 1),
          u[i(2553)](0, 1),
          u[i(2553)](0, 2),
          u.push(0, 3),
          u[i(2553)](0, 4),
          u[i(2553)](1, 2),
          u[i(2553)](2, 3),
          u[i(2553)](3, 4),
          u[i(2553)](4, 1)
        for (let t = this[i(343)][i(1801)], e = 0, s = 0; s <= t; s++) {
          let n = Cesium.Cartesian3[i(1446)](d, C, s / t, new Cesium[i(310)]()),
            o = Cesium.Cartesian3[i(1446)](f, v, s / t, new Cesium[i(310)]()),
            h = []
          for (let s = 0; s <= t; s++) {
            let l = Cesium[i(310)][i(1446)](n, o, s / t, new Cesium[i(310)]())
            r[i(2553)](l.x, l.y, l.z), a.push(1, 1), h[i(2553)](e++)
          }
          p.push(h)
        }
        for (let t = 1; t < p[i(277)]; t++)
          for (let e = 1; e < p[t].length; e++) {
            let s = p[t - 1][e - 1],
              n = p[t][e - 1],
              o = p[t][e],
              r = p[t - 1][e]
            h.push(s, n, o), h[i(2553)](s, o, r)
          }
        for (let t = 0; t < p[i(277)]; t++) l.push(p[t][0]), l[i(2553)](p[t][p[t][i(277)] - 1])
        let _ = p[i(277)]
        for (let t = 0; t < p[0][i(277)]; t++) l[i(2553)](p[0][t]), l[i(2553)](p[_ - 1][t])
        return {
          topPpositions: new Float32Array(r),
          topPindices: new Int32Array(h),
          topPsts: new Float32Array(a),
          topOindices: new Int32Array(l),
          fourPposition: new Float32Array(m),
          fourPindices: new Int32Array(c),
          fourOindices: new Int32Array(u)
        }
      }
      _createGeometry(e, i, s, n, o) {
        const r = t
        let a = {
            position: new Cesium.GeometryAttribute({
              componentDatatype: Cesium[r(2319)].DOUBLE,
              componentsPerAttribute: 3,
              values: i
            }),
            st: new Cesium[r(1948)]({
              componentDatatype: Cesium[r(2319)][r(2441)],
              componentsPerAttribute: 2,
              values: s
            })
          },
          h = Cesium[r(1242)][r(852)](i),
          l = new Cesium.Geometry({
            attributes: a,
            indices: e,
            primitiveType: n,
            boundingSphere: h
          })
        return (
          (function (t) {
            const e = r
            let i = t[e(495)],
              s = t[e(1548)],
              n = i.length
            if (s[e(2251)]) {
              let o = s[e(2251)][e(736)]
              if (s[e(226)]) {
                let o = s[e(226)][e(736)]
                for (let t = 0; t < n; t++) o[t] = 0
                let r,
                  a,
                  h,
                  l = s[e(226)].values,
                  c = new Cesium[e(310)](),
                  u = new Cesium[e(310)](),
                  m = new Cesium[e(310)](),
                  p = new Cesium[e(310)](),
                  d = new Cesium[e(310)]()
                for (let t = 0; t < n; t += 3)
                  (r = 3 * i[t]),
                    (a = 3 * i[t + 1]),
                    (h = 3 * i[t + 2]),
                    Cesium[e(310)][e(866)](o, r, c),
                    Cesium[e(310)].fromArray(o, a, u),
                    Cesium.Cartesian3[e(866)](o, h, m),
                    Cesium[e(310)][e(1676)](m, u, p),
                    Cesium[e(310)][e(1676)](c, u, d),
                    Cesium.Cartesian3[e(2196)](p, d, p),
                    (l[r] += p.x),
                    (l[r + 1] += p.y),
                    (l[r + 2] += p.z),
                    (l[a] += p.x),
                    (l[a + 1] += p.y),
                    (l[a + 2] += p.z),
                    (l[h] += p.x),
                    (l[h + 1] += p.y),
                    (l[h + 2] += p.z)
                !(function (i) {
                  const s = e
                  let n,
                    o,
                    r,
                    a,
                    h = t.attributes[s(226)].values
                  for (let t = 0; t < h[s(277)]; t += 3)
                    (n = h[t]),
                      (o = h[t + 1]),
                      (r = h[t + 2]),
                      (a = 1 / Math.sqrt(n * n + o * o + r * r)),
                      (h[t] = n * a),
                      (h[t + 1] = o * a),
                      (h[t + 2] = r * a)
                })(),
                  (s.normal[e(1344)] = !0)
              } else
                s[e(226)] = new Cesium[e(1948)]({
                  componentDatatype: Cesium[e(2319)][e(2441)],
                  componentsPerAttribute: 3,
                  values: new Float32Array(o.length)
                })
            }
          })(l),
          l
        )
      }
      [t(2270)]() {
        const e = t
        let i = (function (t, e, i, s) {
            const n = a0_0x3b79,
              o = {}
            ;(o[n(277)] = i),
              (o.zReverse = !0),
              (o[n(508)] = i),
              (o.bottomWidth = i),
              (o[n(2089)] = i),
              (o.topWidth = i)
            let r = o
            return (
              (t = Cesium[n(475)][n(1149)](t)),
              (e = Cesium.Math[n(1149)](e)),
              (r[n(508)] = 0),
              (r[n(298)] = 0),
              (r[n(2089)] = i * Math[n(1824)](t)),
              (r[n(833)] = i * Math[n(1824)](e)),
              new et(r)
            )
          })(this[e(343)].hAngle, this[e(343)].vAngle, this[e(343)][e(1981)]),
          s = this[e(1016)](i)
        ;(this[e(1398)] = this._createGeometry(
          s[e(1415)],
          s[e(1374)],
          s[e(1090)],
          Cesium[e(1099)].TRIANGLES
        )),
          (this[e(222)] = this[e(2214)](
            s.topPindices,
            s.topPpositions,
            s[e(1090)],
            Cesium[e(1099)][e(2367)]
          )),
          (this[e(1856)] = this[e(2214)](
            s[e(377)],
            s[e(2318)],
            s[e(1090)],
            Cesium[e(1099)][e(1095)]
          )),
          (this._outlineGeometry = this[e(2214)](
            s[e(2292)],
            s[e(1374)],
            s[e(1090)],
            Cesium.PrimitiveType[e(1095)]
          ))
        let n = new Float32Array(this[e(1398)][e(1548)][e(2251)][e(736)][e(277)])
        for (let t = 0; t < n[e(277)]; t++) n[t] = this[e(1398)][e(1548)].position[e(736)][t]
        this[e(672)] &&
          this._drawCommands[e(277)] &&
          this[e(672)][e(1602)](function (t) {
            const i = e
            t[i(1239)] = t[i(1239)] && t[i(1239)][i(2039)]()
          })
      }
      _computeMatrix() {
        const e = t
        return (
          (this[e(223)] = Cesium[e(2058)][e(1648)](this._cartesian3)),
          Cesium[e(1472)][e(1245)](
            Cesium[e(183)][e(667)](this[e(343)][e(1198)], this._style[e(711)], this._style.roll),
            this._quaternion
          ),
          (this[e(967)] = Cesium[e(2066)][e(235)](
            this[e(1359)],
            this[e(1150)],
            new Cesium[e(310)](1, 1, 1),
            this[e(967)]
          )),
          Cesium[e(2066)][e(332)](this[e(223)], this[e(967)], this[e(967)]),
          this[e(967)]
        )
      }
      _createDrawCommand(e, i) {
        const s = t
        let n,
          o = this[s(1943)]
        const r = {}
        ;(r[s(2416)] = this),
          (r.id = this[s(1570)]),
          o || ((o = i[s(2325)].createPickId(r)), (this[s(1943)] = o)),
          (n = o[s(1070)])
        let a = i[s(2325)],
          h = new Cesium[s(310)]()
        Cesium[s(2066)][s(2483)](this[s(967)], e[s(1611)][s(365)], h)
        let l = new Cesium.BoundingSphere(h, e[s(1611)][s(1981)]),
          c = new Cesium[s(905)]({
            modelMatrix: this[s(967)],
            owner: this,
            primitiveType: e[s(2283)],
            pass: Cesium.Pass[s(2403)],
            boundingVolume: l,
            pickId: 'czm_pickColor'
          }),
          u = this,
          m = Cesium[s(1671)][s(1683)](e)
        return (
          (c.vertexArray = Cesium.VertexArray.fromGeometry({
            context: a,
            geometry: e,
            attributeLocations: m,
            bufferUsage: Cesium[s(1437)].STATIC_DRAW
          })),
          (c[s(1239)][s(577)] = m),
          (c[s(187)] = Cesium[s(178)][s(1879)]({
            context: a,
            vertexShaderSource: s(1276),
            fragmentShaderSource:
              '\nin vec3 v_position;\nin vec3 v_normal; \nuniform vec4  color; \nuniform float specular;\nuniform float shininess;\nuniform vec3  emission;\nin vec2 v_st;\nuniform bool isLine;\nuniform float glowPower;\nuniform vec4 czm_pickColor;\nvoid main() {\n    vec3 positionToEyeEC = -v_position;\n    vec3 normalEC = normalize(v_normal); \n    \n    czm_material material;\n    material.specular = specular;\n    material.shininess = shininess;\n    material.normal = normalEC;\n    material.emission = emission;//vec3(0.2,0.2,0.2);\n    material.diffuse = color.rgb;\n    if (isLine) {\n        material.alpha = 1.0;\n    }\n    else {\n        material.alpha = color.a;\n    } \n    if (v_st.x == 0.0) {\n        out_FragColor = color;\n    } else {\n        out_FragColor = czm_phong(normalize(positionToEyeEC), material, czm_lightDirectionEC);\n    }\n    out_FragColor = color;\n} ',
            attributeLocations: m
          })),
          (c[s(1689)] = Cesium.RenderState[s(2407)]({
            blending: Cesium[s(1870)][s(1785)],
            depthTest: { enabled: !0, func: Cesium[s(1919)][s(689)] },
            cull: { enabled: !1, face: Cesium.CullFace[s(274)] }
          })),
          (c.uniformMap = {}),
          (c[s(1720)][s(1208)] = function () {
            return n
          }),
          (c[s(1720)].projectionMatrix = function () {
            const t = s
            return i[t(1391)][t(1850)].projectionMatrix
          }),
          (c[s(1720)][s(1813)] = function () {
            const t = s
            return i.context[t(656)][t(1281)]
          }),
          (c[s(1720)].shininess = function () {
            const t = s
            return u[t(2108)] || (u.shininess = 0), u[t(2108)]
          }),
          (c[s(1720)][s(1723)] = function () {
            const t = s
            return u[t(1723)] || (u[t(1723)] = new Cesium[t(310)](0.2, 0.2, 0.2)), u[t(1723)]
          }),
          (c[s(1720)].specular = function () {
            const t = s
            return u[t(838)] || (u[t(838)] = 0), u[t(838)]
          }),
          (c[s(1720)][s(1949)] = function () {
            const t = s
            return (
              e.primitiveType == Cesium[t(1099)][t(1095)] ||
              e[t(2283)] == Cesium[t(1099)].LINE_STRIP
            )
          }),
          (c[s(1720)][s(1070)] = function () {
            const t = s
            return e[t(2283)] == Cesium[t(1099)].LINES || e[t(2283)] == Cesium[t(1099)][t(1313)]
              ? u[t(1654)]
              : u._color
          }),
          (c[s(1720)][s(2347)] = function () {
            const t = s
            return i.context[t(656)][t(226)]
          }),
          (c[s(1720)][s(1523)] = function () {
            return 0.25
          }),
          c
        )
      }
      [t(1451)](e) {
        const i = t
        this[i(1046)](), this[i(2270)]()
        const s = this._isSelected ? this[i(343)][i(506)] : this[i(343)][i(1070)]
        if (
          ((this._color = Cesium[i(1154)].fromCssColorString(s)),
          (this[i(1654)] = Cesium[i(1154)][i(2008)](this[i(343)][i(1018)])),
          (this[i(1398)][i(1611)] = Cesium.BoundingSphere[i(852)](
            this[i(1398)][i(1548)].position.values
          )),
          (this[i(672)] = []),
          this[i(672)].push(this[i(2063)](this[i(1398)], e)),
          this[i(343)].lineShow && this[i(672)][i(2553)](this[i(2063)](this[i(1185)], e)),
          this._style[i(948)])
        ) {
          let t = this[i(2063)](this[i(222)], e)
          if ((this[i(672)].push(t), this[i(343)][i(1105)])) {
            let t = this[i(2063)](this[i(1856)], e)
            this._drawCommands[i(2553)](t)
          }
        }
      }
      update(e) {
        const i = t
        this[i(1144)] &&
          (this._needUpdate && (this[i(1451)](e), (this[i(1897)] = !1)),
          e[i(1312)] == Cesium[i(1857)].SCENE3D &&
            this[i(672)].forEach(function (t) {
              e[i(331)].push(t)
            }))
      }
      [t(2039)]() {}
      _addHook(e) {
        const i = t
        ;(this[i(2462)] = e),
          (this[i(2122)] = e.id),
          (this[i(2366)] = this[i(1570)]),
          e[i(245)][i(696)][i(1346)][i(1861)](this)
      }
      _removeHook(e) {
        const i = t
        e._viewer[i(696)][i(1346)][i(1896)](this),
          this[i(1943)] && this[i(1943)][i(2039)](),
          this[i(672)].forEach(function (t) {
            const e = i
            t[e(1239)] = t.vertexArray && t[e(1239)].destroy()
          }),
          (this[i(672)] = [])
      }
    }
    class et {
      constructor(e) {
        const i = t
        ;(this[i(848)] = e[i(298)]),
          (this[i(1067)] = e[i(508)]),
          (this[i(1015)] = e[i(833)]),
          (this[i(1866)] = e.topHeight),
          (this[i(2001)] = e.length),
          (this[i(1552)] = e[i(1455)]),
          (this[i(2125)] = Cesium[i(1960)](e[i(212)], 8))
      }
    }
    function it(e) {
      const i = t
      ;(this[i(277)] = e[i(277)]),
        (this.topRadius = e[i(523)]),
        (this[i(2395)] = e[i(2395)]),
        (this[i(212)] = e.slices ? e[i(212)] : 64),
        (this.zReverse = e[i(1455)])
    }
    let st = new Cesium[t(2007)]()
    ;(it[t(2214)] = function (e) {
      const i = t
      var s = e.length,
        n = e[i(523)],
        o = e[i(2395)],
        r = e[i(212)],
        a = (2 * Math.PI) / (r - 1),
        h = [],
        l = [],
        c = [],
        u = [],
        m = [o, n],
        p = [0, e[i(1455)] ? -s : s],
        d = 0,
        f = Math.atan2(o - n, s),
        C = new Cesium[i(194)]()
      C.z = Math[i(884)](f)
      for (var v = Math[i(1272)](f), _ = 0; _ < p.length; _++) {
        u[_] = []
        for (var g = m[_], y = 0; y < r; y++) {
          u[_][i(2553)](d++)
          var w = a * y,
            x = g * Math[i(1272)](w),
            b = g * Math.sin(w)
          h.push(x, b, p[_]),
            (x = v * Math[i(1272)](w)),
            (b = v * Math[i(884)](w)),
            l[i(2553)](x, b, C.z),
            c.push(_ / (p[i(277)] - 1), 0)
        }
      }
      var S = []
      for (_ = 1; _ < p[i(277)]; _++)
        for (y = 1; y < r; y++) {
          var P = u[_ - 1][y - 1],
            M = u[_][y - 1],
            A = u[_][y],
            T = u[_ - 1][y]
          S[i(2553)](A),
            S[i(2553)](T),
            S[i(2553)](P),
            S[i(2553)](A),
            S.push(P),
            S[i(2553)](M),
            y == u[_][i(277)] - 1 &&
              ((P = u[_ - 1][y]),
              (M = u[_][y]),
              (A = u[_][0]),
              (T = u[_ - 1][0]),
              S[i(2553)](A),
              S[i(2553)](T),
              S[i(2553)](P),
              S[i(2553)](A),
              S[i(2553)](P),
              S.push(M))
        }
      ;(S = new Int16Array(S)),
        (h = new Float32Array(h)),
        (l = new Float32Array(l)),
        (c = new Float32Array(c))
      var E = {
          position: new Cesium[i(1948)]({
            componentDatatype: Cesium[i(2319)].DOUBLE,
            componentsPerAttribute: 3,
            values: h
          }),
          normal: new Cesium[i(1948)]({
            componentDatatype: Cesium.ComponentDatatype[i(2441)],
            componentsPerAttribute: 3,
            values: l
          }),
          st: new Cesium[i(1948)]({
            componentDatatype: Cesium[i(2319)].FLOAT,
            componentsPerAttribute: 2,
            values: c
          })
        },
        z = Cesium[i(1242)][i(852)](h),
        D = new Cesium.Geometry({
          attributes: E,
          indices: S,
          primitiveType: Cesium[i(1099)][i(2367)],
          boundingSphere: z
        })
      return (h = []), (S = []), (c = []), D
    }),
      (it[t(933)] = function (e, i) {
        const s = t
        if (!i) return this._createGeometry(e)
        Cesium[s(2066)][s(2483)](i, Cesium.Cartesian3[s(1738)], Cesium[s(310)]),
          Cesium[s(310)][s(1902)](st[s(1893)])
        var n = e[s(277)],
          o = e.topRadius,
          r = (e.bottomRadius, e[s(212)]),
          h = (2 * Math.PI) / (r - 1),
          l = [],
          c = [],
          m = [],
          p = [],
          d = [0, e[s(1455)] ? -n : n],
          f = 0
        ;(f = 0), l[s(2553)](0, 0, 0), c[s(2553)](1, 1), f++
        for (var C = new Cesium.Cartesian3(), v = o / 15, _ = 0; _ < 16; _++) {
          for (var g = v * _, y = [], w = 0; w < r; w++) {
            var x = h * w,
              b = g * Math[s(1272)](x),
              S = g * Math[s(884)](x)
            ;(C.x = b), (C.y = S), (C.z = d[1])
            var P = (0, a[s(1453)])(C, i, st)
            P
              ? (y.push(f), l[s(2553)](b, S, d[1]), c.push(_ / 15, 1), f++)
              : ((P = u), y[s(2553)](-1))
          }
          p[s(2553)](y)
        }
        for (var M, A, T = [0, p.length - 1], E = 0; E < T[s(277)]; E++)
          for (_ = T[E], w = 1; w < p[_].length; w++)
            (M = p[_][w - 1]), (A = p[_][w]), M >= 0 && A >= 0 && m[s(2553)](0, M, A)
        ;(l = new Float32Array(l)), (m = new Int32Array(m)), (c = new Float32Array(c))
        var z = {
            position: new Cesium[s(1948)]({
              componentDatatype: Cesium[s(2319)][s(340)],
              componentsPerAttribute: 3,
              values: l
            }),
            st: new Cesium[s(1948)]({
              componentDatatype: Cesium.ComponentDatatype.FLOAT,
              componentsPerAttribute: 2,
              values: c
            })
          },
          D = Cesium[s(1242)][s(852)](l),
          I = new Cesium[s(1746)]({
            attributes: z,
            indices: m,
            primitiveType: Cesium[s(1099)][s(2367)],
            boundingSphere: D
          })
        return _computeVertexNormals(I), (l = []), (m = []), I
      }),
      (it[t(306)] = function (e) {
        const i = t
        var s = e[i(277)],
          n = e[i(523)],
          o = e.bottomRadius,
          r = e[i(212)],
          a = (2 * Math.PI) / (r - 1),
          h = [],
          l = [],
          c = [],
          u = [],
          m = [o, n],
          p = [0, e[i(1455)] ? -s : s],
          d = 0,
          f = Math[i(2003)](o - n, s),
          C = new Cesium[i(194)]()
        C.z = Math[i(884)](f)
        for (var v = Math.cos(f), _ = 0; _ < p.length; _++) {
          u[_] = []
          for (var g = m[_], y = 0; y < r; y++) {
            u[_][i(2553)](d++)
            var w = a * y,
              x = g * Math.cos(w),
              b = g * Math.sin(w)
            h[i(2553)](x, b, p[_]),
              (x = v * Math.cos(w)),
              (b = v * Math.sin(w)),
              l[i(2553)](x, b, C.z),
              c[i(2553)](_ / (p[i(277)] - 1), 0)
          }
        }
        var S = []
        for (_ = 1; _ < p.length; _++)
          for (y = 1; y < r; y += 1) {
            var P = u[_ - 1][y - 1],
              M = u[_][y - 1]
            u[_][y], u[_ - 1][y], y % 8 == 1 && S[i(2553)](P, M)
          }
        ;(S = new Int16Array(S)),
          (h = new Float32Array(h)),
          (l = new Float32Array(l)),
          (c = new Float32Array(c))
        var A = {
            position: new Cesium.GeometryAttribute({
              componentDatatype: Cesium[i(2319)][i(340)],
              componentsPerAttribute: 3,
              values: h
            }),
            normal: new Cesium[i(1948)]({
              componentDatatype: Cesium.ComponentDatatype.FLOAT,
              componentsPerAttribute: 3,
              values: l
            }),
            st: new Cesium[i(1948)]({
              componentDatatype: Cesium[i(2319)][i(2441)],
              componentsPerAttribute: 2,
              values: c
            })
          },
          T = Cesium[i(1242)][i(852)](h),
          E = new Cesium[i(1746)]({
            attributes: A,
            indices: S,
            primitiveType: Cesium[i(1099)][i(1095)],
            boundingSphere: T
          })
        return (h = []), (S = []), (c = []), E
      })
    class nt extends b {
      constructor(e) {
        const i = t
        super(e),
          (this[i(350)] = 1),
          (this[i(528)] = i(1586)),
          (this._geometryType = i(1865)),
          (this[i(2336)] = Q.cylinderRadar),
          (this._graphicClassType = W),
          (this[i(343)][i(1981)] = Cesium.defaultValue(this[i(343)][i(1981)], 100)),
          (this[i(343)][i(1241)] = Cesium[i(1960)](this._style.angle, 22)),
          (this[i(343)][i(1048)] = 90 - this[i(343)][i(1241)]),
          (this._style.heading = Cesium.defaultValue(this._style.heading, 0)),
          (this._style[i(711)] = Cesium.defaultValue(this[i(343)].pitch, 0)),
          (this[i(343)][i(1143)] = Cesium[i(1960)](this[i(343)][i(1143)], 0)),
          (this[i(343)].color = Cesium[i(1960)](this[i(343)].color, i(400))),
          (this[i(343)][i(506)] = Cesium[i(1960)](this[i(343)][i(506)], i(624))),
          (this[i(343)][i(2248)] = Cesium.defaultValue(this._style[i(2248)], !0)),
          (this[i(343)][i(1018)] = Cesium[i(1960)](this[i(343)][i(1018)], i(149))),
          (this[i(343)][i(948)] = Cesium.defaultValue(this[i(343)][i(948)], !0)),
          (this[i(343)].topline = Cesium.defaultValue(this[i(343)][i(1105)], !0)),
          (this._drawCommands = [])
        let s = Cesium[i(1960)](e[i(2251)], [120, 40, 0])
        ;(this[i(2238)] = Cesium[i(310)][i(667)](s[0], s[1], s[2])),
          (this[i(1912)] = s),
          (this[i(1636)] = Cesium[i(2066)][i(1902)](Cesium.Matrix4.IDENTITY)),
          (this[i(1150)] = new Cesium[i(1472)]()),
          (this[i(1359)] = new Cesium[i(310)]()),
          (this[i(2035)] = new Cesium.Cartesian3(1, 1, 1)),
          (this[i(967)] = new Cesium[i(2066)]()),
          (this[i(1897)] = !0),
          (this[i(1024)] = new Cesium[i(1242)](this[i(2238)], 1.2 * this[i(343)][i(1981)]))
      }
      [t(1601)](e, i) {
        const s = t
        if (!e) return
        let n = e,
          o = Cesium.Cartesian3[s(667)](n[0], n[1], n[2]),
          r = o[s(1902)](),
          a = Cesium[s(310)][s(1676)](this[s(2238)], o, o)
        a = Cesium[s(310)][s(379)](a, a)
        let h = Cesium[s(310)][s(379)](this[s(2238)], new Cesium.Cartesian3()),
          l = $(this._cartesian3, r),
          c = this._style
        i && (c[s(1981)] = Cesium[s(310)][s(1849)](this._cartesian3, r)),
          (c[s(1198)] = l - 90),
          (c.roll = 0),
          (c[s(711)] = Cesium[s(475)][s(363)](Cesium.Cartesian3[s(354)](a, h)) - 180),
          (this[s(1679)] = c)
      }
      setPosition(e) {
        const i = t
        ;(this[i(1912)] = e),
          (this._cartesian3 = Cesium[i(310)].fromDegrees(e[0], e[1], e[2])),
          (this[i(794)] = this[i(1912)]),
          (this._needUpdate = !0)
      }
      [t(2231)]() {
        const e = t
        ;(this[e(343)][e(1048)] = 90 - this[e(343)][e(1241)]),
          (this[e(1897)] = !0),
          (this._boundingSphere = new Cesium.BoundingSphere(
            this[e(2238)],
            1.2 * this[e(343)].radius
          ))
      }
      [t(1046)](e, i) {
        const s = t
        return (
          (this[s(223)] = Cesium[s(2058)].eastNorthUpToFixedFrame(this._cartesian3)),
          Cesium[s(1472)].fromHeadingPitchRoll(
            Cesium[s(183)][s(667)](
              this[s(343)].heading,
              this[s(343)][s(711)],
              this[s(343)][s(1143)]
            ),
            this[s(1150)]
          ),
          (this[s(967)] = Cesium[s(2066)][s(235)](
            this._translation,
            this._quaternion,
            new Cesium[s(310)](1, 1, 1),
            this[s(967)]
          )),
          Cesium.Matrix4[s(332)](this[s(223)], this[s(967)], this[s(967)]),
          this[s(967)]
        )
      }
      [t(2097)]() {
        const e = t
        let i = [],
          s = [],
          n = [],
          o = []
        for (
          var r = this._style.radius,
            a = 90 - parseInt(this._style.angleT),
            h = a < 1 ? a / 8 : 1,
            l = (2 * Math.PI) / 127,
            c = 0,
            u = this[e(343)][e(1048)];
          u < 91;
          u += h
        ) {
          var m = Cesium[e(475)][e(1149)](u < 90 ? u : 90)
          m = Math.cos(m) * r
          for (var p = [], d = 0; d < 128; d++) {
            let t = l * d,
              n = m * Math[e(1272)](t),
              o = m * Math.sin(t),
              a = Math[e(1345)](r * r - n * n - o * o)
            i[e(2553)](n, o, a), s[e(2553)](1, 1), p[e(2553)](c++)
          }
          o[e(2553)](p)
        }
        for (u = 1; u < o[e(277)]; u++)
          for (d = 1; d < o[u][e(277)]; d++) {
            let t = o[u - 1][d - 1],
              i = o[u][d - 1],
              s = o[u][d],
              r = o[u - 1][d]
            n[e(2553)](t, i, s), n.push(t, s, r)
          }
        ;(i = new Float32Array(i)), (n = new Int32Array(n)), (s = new Float32Array(s))
        let f = {
            position: new Cesium[e(1948)]({
              componentDatatype: Cesium[e(2319)].DOUBLE,
              componentsPerAttribute: 3,
              values: i
            }),
            st: new Cesium[e(1948)]({
              componentDatatype: Cesium.ComponentDatatype[e(2441)],
              componentsPerAttribute: 2,
              values: s
            })
          },
          C = Cesium[e(1242)].fromVertices(i),
          v = new Cesium[e(1746)]({
            attributes: f,
            indices: n,
            primitiveType: Cesium.PrimitiveType.TRIANGLES,
            boundingSphere: C
          })
        return this._computeVertexNormals(v), v
      }
      _createTopOutlineGeometry() {
        const e = t
        for (
          var i = this[e(343)][e(1981)],
            s = [],
            n = [],
            o = [],
            r = [],
            a = 90 - parseInt(this[e(343)][e(1048)]),
            h = a < 1 ? a / 8 : 1,
            l = (2 * Math.PI) / 127,
            c = 0,
            u = this._style.angleT;
          u < 91;
          u += h
        ) {
          var m = Cesium[e(475)][e(1149)](u < 90 ? u : 90)
          m = Math[e(1272)](m) * i
          for (var p = [], d = 0; d < 128; d++) {
            let t = l * d,
              o = m * Math[e(1272)](t),
              r = m * Math[e(884)](t),
              a = Math.sqrt(i * i - o * o - r * r)
            s[e(2553)](o, r, a), n[e(2553)](1, 1), p[e(2553)](c++)
          }
          r[e(2553)](p)
        }
        for (u = 1; u < r.length; u++)
          for (d = 1; d < r[u][e(277)]; d++) {
            let t = r[u - 1][d - 1],
              e = r[u][d - 1],
              i = r[u][d]
            r[u - 1][d], d % 8 == 1 && o.push(t, e), u % 8 == 1 && o.push(e, i)
          }
        ;(s = new Float32Array(s)), (o = new Int32Array(o)), (n = new Float32Array(n))
        let f = {
            position: new Cesium[e(1948)]({
              componentDatatype: Cesium[e(2319)][e(340)],
              componentsPerAttribute: 3,
              values: s
            }),
            st: new Cesium[e(1948)]({
              componentDatatype: Cesium[e(2319)].FLOAT,
              componentsPerAttribute: 2,
              values: n
            })
          },
          C = Cesium[e(1242)][e(852)](s),
          v = new Cesium[e(1746)]({
            attributes: f,
            indices: o,
            primitiveType: Cesium[e(1099)].LINES,
            boundingSphere: C
          })
        return this[e(1026)](v), v
      }
      [t(2270)]() {
        const e = t
        ;(this[e(1398)] = it[e(933)](
          new it({
            topRadius:
              this[e(343)].radius * Math[e(1272)](Cesium.Math[e(1149)](this[e(343)][e(1048)])),
            bottomRadius: 0,
            length:
              this[e(343)][e(1981)] * Math[e(884)](Cesium[e(475)].toRadians(this._style[e(1048)]))
          })
        )),
          (this[e(222)] = this[e(2097)]()),
          (this[e(1856)] = this[e(207)]()),
          (this[e(1185)] = it[e(306)](
            new it({
              topRadius:
                this[e(343)][e(1981)] * Math[e(1272)](Cesium[e(475)][e(1149)](this[e(343)].angleT)),
              bottomRadius: 0,
              slices: 128,
              length:
                this._style[e(1981)] * Math[e(884)](Cesium[e(475)][e(1149)](this[e(343)][e(1048)]))
            })
          ))
        let i = new Float32Array(this._geometry[e(1548)][e(2251)].values[e(277)])
        for (var s = 0; s < i[e(277)]; s++) i[s] = this._geometry[e(1548)][e(2251)][e(736)][s]
        this[e(672)][e(1602)]((t) => {
          const i = e
          t.vertexArray = t[i(1239)] && t.vertexArray[i(2039)]()
        }),
          this[e(672)][e(2054)](0, this[e(672)][e(277)])
      }
      [t(2063)](e, i) {
        const s = t
        let n,
          o = this[s(1943)]
        const r = {}
        ;(r.primitive = this),
          (r.id = this._id),
          o || ((o = i[s(2325)][s(896)](r)), (this._pickId = o)),
          (n = o[s(1070)])
        let a = this
        var h = i[s(2325)]
        let l = new Cesium[s(310)]()
        Cesium[s(2066)][s(2483)](this[s(967)], e[s(1611)][s(365)], l)
        var c = new Cesium[s(1242)](l, e.boundingSphere.radius)
        let u = new Cesium.DrawCommand({
            modelMatrix: this._matrix,
            owner: this,
            primitiveType: e.primitiveType,
            pass: Cesium.Pass[s(2403)],
            boundingVolume: c,
            pickId: s(1208)
          }),
          m = this,
          p = Cesium[s(1671)][s(1683)](e)
        return (
          (u[s(1239)] = Cesium[s(1057)].fromGeometry({
            context: h,
            geometry: e,
            attributeLocations: p,
            bufferUsage: Cesium[s(1437)][s(2200)]
          })),
          (u[s(1239)][s(577)] = p),
          (u[s(187)] = Cesium[s(178)].replaceCache({
            context: h,
            vertexShaderSource: s(1831),
            fragmentShaderSource: s(800),
            attributeLocations: p
          })),
          (u[s(1689)] = Cesium.RenderState[s(2407)]({
            blending: Cesium[s(1870)][s(1785)],
            depthTest: { enabled: !0, func: Cesium[s(1919)][s(689)] },
            cull: { enabled: !1, face: Cesium[s(1749)][s(274)] }
          })),
          (u.uniformMap = {}),
          (u[s(1720)][s(1208)] = function () {
            return n
          }),
          (u.uniformMap[s(703)] = function () {
            const t = s
            return i[t(2325)][t(656)][t(1987)]
          }),
          (u[s(1720)][s(1813)] = function () {
            const t = s
            return i.context[t(656)].modelView
          }),
          (u[s(1720)][s(2108)] = function () {
            const t = s
            return m[t(2108)] || (m[t(2108)] = 0), m[t(2108)]
          }),
          (u[s(1720)][s(1723)] = function () {
            const t = s
            return m[t(1723)] || (m.emission = new Cesium.Cartesian3(0.2, 0.2, 0.2)), m[t(1723)]
          }),
          (u[s(1720)][s(838)] = function () {
            const t = s
            return m[t(838)] || (m[t(838)] = 0), m[t(838)]
          }),
          (u.uniformMap[s(1949)] = function () {
            const t = s
            return (
              e[t(2283)] == Cesium.PrimitiveType[t(1095)] || e[t(2283)] == Cesium[t(1099)][t(1313)]
            )
          }),
          (u[s(1720)][s(262)] = function () {
            const t = s
            let i = Cesium[t(1154)][t(971)]
            return (i = e.primitiveType == Cesium[t(1099)][t(1095)] ? a[t(1654)] : a[t(2390)])
          }),
          (u.uniformMap[s(2347)] = function () {
            const t = s
            return i[t(2325)][t(656)][t(226)]
          }),
          (u[s(1720)][s(1523)] = function () {
            return 0.25
          }),
          u
        )
      }
      [t(1451)](e) {
        const i = t
        this[i(1046)](), this[i(2270)]()
        let s = this[i(1496)] ? this._style[i(506)] : this[i(343)][i(1070)]
        if (
          ((this[i(2390)] = Cesium[i(1154)].fromCssColorString(s)),
          (this[i(1654)] = Cesium[i(1154)][i(2008)](this[i(343)][i(1018)])),
          (this[i(1398)][i(1611)] = Cesium[i(1242)].fromVertices(
            this[i(1398)][i(1548)][i(2251)][i(736)]
          )),
          (this[i(672)] = []),
          this[i(672)].push(this[i(2063)](this[i(1398)], e)),
          this[i(343)][i(2248)] && this[i(672)][i(2553)](this[i(2063)](this._outlineGeometry, e)),
          this[i(343)][i(948)])
        ) {
          let t = this._createDrawCommand(this._topGeometry, e)
          if ((this[i(672)][i(2553)](t), this._style.topline)) {
            let t = this[i(2063)](this._topOutlineGeometry, e)
            this[i(672)].push(t)
          }
        }
      }
      [t(720)](e) {
        const i = t
        this._show &&
          (this[i(1897)] && (this[i(1451)](e), (this[i(1897)] = !1)),
          e.mode == Cesium[i(1857)][i(339)] &&
            this._drawCommands[i(1602)](function (t) {
              e[i(331)].push(t)
            }))
      }
      [t(2039)](e) {
        const i = t
        this[i(672)][i(1602)](function (t) {
          const e = i
          t[e(1239)] = t[e(1239)] && t[e(1239)][e(2039)]()
        }),
          (this[i(672)] = [])
      }
      [t(1026)](e) {
        const i = t
        var s = e[i(495)],
          n = e[i(1548)],
          o = s[i(277)]
        if (n[i(2251)]) {
          var r = n[i(2251)].values
          if (n[i(226)]) {
            r = n[i(226)][i(736)]
            for (var a = 0; a < o; a++) r[a] = 0
            for (
              var h,
                l,
                c,
                u = n[i(226)][i(736)],
                m = new Cesium[i(310)](),
                p = new Cesium[i(310)](),
                d = new Cesium[i(310)](),
                f = new Cesium.Cartesian3(),
                C = new Cesium[i(310)](),
                v = 0;
              v < o;
              v += 3
            )
              (h = 3 * s[v]),
                (l = 3 * s[v + 1]),
                (c = 3 * s[v + 2]),
                Cesium[i(310)].fromArray(r, h, m),
                Cesium.Cartesian3[i(866)](r, l, p),
                Cesium[i(310)][i(866)](r, c, d),
                Cesium[i(310)][i(1676)](d, p, f),
                Cesium[i(310)][i(1676)](m, p, C),
                Cesium.Cartesian3[i(2196)](f, C, f),
                (u[h] += f.x),
                (u[h + 1] += f.y),
                (u[h + 2] += f.z),
                (u[l] += f.x),
                (u[l + 1] += f.y),
                (u[l + 2] += f.z),
                (u[c] += f.x),
                (u[c + 1] += f.y),
                (u[c + 2] += f.z)
            !(function (t) {
              const s = i
              for (var n, o, r, a, h = e.attributes[s(226)][s(736)], l = 0; l < h[s(277)]; l += 3)
                (n = h[l]),
                  (o = h[l + 1]),
                  (r = h[l + 2]),
                  (a = 1 / Math[s(1345)](n * n + o * o + r * r)),
                  (h[l] = n * a),
                  (h[l + 1] = o * a),
                  (h[l + 2] = r * a)
            })(),
              (n.normal[i(1344)] = !0)
          } else
            n[i(226)] = new Cesium[i(1948)]({
              componentDatatype: Cesium[i(2319)].FLOAT,
              componentsPerAttribute: 3,
              values: new Float32Array(r[i(277)])
            })
        }
        return e
      }
      _addHook(e) {
        const i = t
        ;(this[i(2462)] = e),
          (this[i(2122)] = e.id),
          (this[i(2366)] = this[i(1570)]),
          e[i(245)].scene[i(1346)].add(this)
      }
      [t(2389)](e) {
        const i = t
        e[i(245)].scene[i(1346)][i(1896)](this),
          this[i(1943)] && this._pickId[i(2039)](),
          this[i(672)][i(1602)](function (t) {
            const e = i
            t[e(1239)] = t.vertexArray && t[e(1239)][e(2039)]()
          }),
          (this[i(672)] = [])
      }
    }
    class ot extends b {
      constructor(e = {}) {
        const i = t
        super(e),
          (this._typeName = i(1603)),
          (this[i(314)] = 'Point'),
          (this[i(2336)] = Q[i(427)]),
          (this[i(1115)] = W),
          (this[i(350)] = 1),
          (this._style.color = Cesium[i(1960)](this[i(343)][i(1070)], i(1655))),
          (this[i(343)][i(506)] = Cesium[i(1960)](this[i(343)][i(506)], i(615))),
          (this[i(343)][i(277)] = Cesium[i(1960)](this[i(343)][i(277)], 2e3)),
          (this[i(343)][i(1981)] = Cesium[i(1960)](this._style[i(1981)], 100)),
          (this[i(343)][i(1198)] = Cesium[i(1960)](this[i(343)][i(1198)], 0)),
          (this[i(343)][i(711)] = Cesium[i(1960)](this._style[i(711)], 0)),
          (this[i(343)][i(1143)] = Cesium[i(1960)](this._style[i(1143)], 0))
        let s = Cesium[i(1960)](e[i(2251)], [111, 28, 0])
        ;(this[i(2238)] = Cesium.Cartesian3[i(667)](s[0], s[1], s[2])),
          (this[i(794)] = this[i(1912)] = s),
          (this[i(1934)] = this._createPrimitive()),
          (this[i(2549)] = !1),
          (this[i(1024)] = new Cesium[i(1242)](this[i(2238)], 1.2 * this[i(343)][i(277)]))
      }
      get primitive() {
        return this[t(1934)]
      }
      _setVisible(e) {
        const i = t
        this[i(1934)][i(1482)] = e
      }
      _setStyle() {
        const e = t
        this[e(2549)] &&
          (this[e(684)](this[e(2462)]),
          (this._primitive = this[e(438)]()),
          this[e(1768)](this[e(2462)]))
        const i = this[e(1496)] ? this[e(343)][e(506)] : this[e(343)][e(1070)]
        ;(this[e(1934)][e(1994)].material[e(1285)][e(1070)] = Cesium[e(1154)][e(2008)](i)),
          (this._boundingSphere = new Cesium[e(1242)](this._cartesian3, 1.2 * this[e(343)][e(277)]))
      }
      [t(1253)](e) {
        const i = t
        ;(this[i(1912)] = e),
          (this[i(2238)] = Cesium[i(310)][i(667)](e[0], e[1], e[2])),
          (this._coordinates = this._position),
          this[i(2231)]()
      }
      [t(1601)](e, i) {
        const s = t
        if (!e) return
        let n = e,
          o = Cesium.Cartesian3[s(667)](n[0], n[1], n[2]),
          r = o[s(1902)](),
          a = Cesium.Cartesian3[s(1676)](this[s(2238)], o, o)
        a = Cesium[s(310)][s(379)](a, a)
        let h = Cesium[s(310)][s(379)](this[s(2238)], new Cesium.Cartesian3()),
          l = $(this[s(2238)], r),
          c = this._style
        i && (c[s(277)] = Cesium[s(310)][s(1849)](this[s(2238)], r)),
          (c[s(1198)] = l - 90),
          (c[s(1143)] = 0),
          (c.pitch = Cesium.Math.toDegrees(Cesium[s(310)].angleBetween(a, h))),
          (this[s(1679)] = c)
      }
      [t(438)]() {
        const e = t
        let i = (function (t, e, i, s, n, o, r) {
            const a = a0_0x3b79
            ;(e = Cesium.Math[a(1149)](e)),
              (i = Cesium[a(475)][a(1149)](i)),
              (s = Cesium[a(475)][a(1149)](s)),
              (n = n || 0),
              (o = o || 0),
              (r = r || 0)
            var h = Cesium[a(1472)][a(1245)](new Cesium[a(183)](e, i, s)),
              l = Cesium.Matrix4[a(235)](
                new Cesium[a(310)](0, 0, 0),
                h,
                new Cesium.Cartesian3(1, 1, 1)
              ),
              c = Cesium[a(2058)][a(1648)](t, void 0, new Cesium[a(2066)]()),
              u = Cesium.Matrix4[a(999)](c, l, new Cesium[a(2066)]()),
              m = new Cesium[a(2066)]()
            return Cesium[a(2066)][a(229)](u, new Cesium[a(310)](n, o, r), m), m
          })(
            this[e(2238)],
            this[e(343)].heading,
            this[e(343)][e(711)],
            this[e(343)][e(1143)],
            0,
            0,
            -this[e(343)][e(277)] / 2
          ),
          s = new Cesium[e(1355)]({
            length: this[e(343)][e(277)],
            topRadius: 0,
            bottomRadius: this[e(343)][e(1981)],
            slices: 128,
            shadows: Cesium[e(862)][e(1192)],
            vertexFormat: Cesium.MaterialAppearance[e(2272)].TEXTURED[e(1572)]
          }),
          n = new Cesium[e(2155)]({ geometry: s, modelMatrix: i })
        ;(this[e(1636)] = i), (this[e(2457)] = n)
        let o = this[e(1530)]()
        const r = {}
        return (
          (r[e(1122)] = [n]),
          (r.appearance = o),
          (r.releaseGeometryInstances = !1),
          (r[e(644)] = !1),
          new Cesium[e(2186)](r)
        )
      }
      [t(1530)]() {
        const e = t
        return new Cesium[e(985)]({
          material: new Cesium[e(1637)]({
            fabric: {
              uniforms: {
                color: Cesium[e(1154)][e(2008)](this[e(343)].color),
                repeat: 1,
                offset: 0,
                thickness: 0.8
              },
              source: e(1707)
            },
            translucent: !0
          }),
          faceForward: !1,
          closed: !0
        })
      }
      _addHook(t) {
        this._add(t)
      }
      [t(2389)](e) {
        this[t(684)](e)
      }
      [t(1768)](e) {
        const i = t
        ;(this[i(1934)][i(2122)] = e.id),
          (this[i(1934)][i(2366)] = this[i(1570)]),
          e[i(245)].scene[i(1346)][i(1861)](this[i(1934)]),
          (this._isAdd = !0),
          (this[i(2462)] = e)
      }
      [t(684)](e) {
        const i = t
        this._isAdd && (e._viewer[i(696)][i(1346)][i(1896)](this._primitive), (this[i(2549)] = !1))
      }
    }
    class rt extends ot {
      constructor(e = {}) {
        const i = t
        ;(e.style = Cesium[i(1960)](e[i(1679)], {})),
          (e[i(1679)][i(638)] = Cesium[i(1960)](e[i(1679)][i(638)], 20)),
          (e[i(1679)].thickness = Cesium[i(1960)](e[i(1679)][i(171)], 0.3)),
          super(e),
          (this._typeName = '探测雷达'),
          (this[i(2336)] = Q.probeRadar)
      }
      [t(1530)]() {
        const e = t
        return new Cesium[e(985)]({
          material: new Cesium[e(1637)]({
            fabric: {
              uniforms: {
                color: Cesium[e(1154)][e(2008)](this[e(343)][e(1070)]),
                repeat: this[e(343)][e(638)],
                thickness: this[e(343)][e(171)]
              },
              source: e(2094)
            },
            translucent: !0
          }),
          faceForward: !1,
          closed: !0
        })
      }
    }
    const at = {}
    ;(at[t(2251)] = 0), (at[t(226)] = 1)
    let ht = at,
      lt = t(1834),
      ct =
        '\nin vec4 position;\nin vec3 normal; \nout vec3 v_position;\nout vec3 v_positionWC;\nout vec3 v_positionEC;\nout vec3 v_normalEC;\n\nvoid main()\n{\n    gl_Position = czm_modelViewProjection * position;\n    v_position = vec3(position);\n    v_positionWC = (czm_model * position).xyz;\n    v_positionEC = (czm_modelView * position).xyz;\n    v_normalEC = czm_normal * normal;\n}',
      ut = t(2266),
      mt = t(598)
    class pt extends b {
      constructor(e) {
        const i = t
        super(e),
          (this[i(528)] = i(1126)),
          (this[i(314)] = i(1865)),
          (this._graphicClassType = W),
          (this[i(2336)] = Q.rectSensor),
          (this[i(350)] = 1),
          (this._style.slice = Cesium[i(1960)](this[i(343)][i(123)], 32)),
          (this._modelMatrix = new Cesium.Matrix4()),
          (this[i(1988)] = new Cesium[i(2066)]()),
          (this._computedScanPlaneModelMatrix = new Cesium[i(2066)]()),
          (this._style.radius = Cesium[i(1960)](this[i(343)][i(1981)], 20)),
          (this[i(343)].hAngle = Cesium[i(1960)](this[i(343)][i(2519)], 70)),
          (this[i(343)][i(2493)] = Cesium.defaultValue(this[i(343)][i(2493)], 70)),
          (this[i(343)][i(1198)] = Cesium[i(1960)](this[i(343)][i(1198)], 0)),
          (this[i(343)].pitch = Cesium[i(1960)](this[i(343)][i(711)], 0)),
          (this[i(343)].roll = Cesium[i(1960)](this[i(343)][i(1143)], 0)),
          (this._style[i(1070)] = Cesium[i(1960)](this[i(343)].color, i(2180))),
          (this[i(343)][i(506)] = Cesium[i(1960)](this[i(343)][i(506)], i(456))),
          (this[i(343)].lineShow = Cesium[i(1960)](this[i(343)][i(2248)], !0)),
          (this[i(343)][i(1018)] = Cesium.defaultValue(this._style.lineColor, i(1320))),
          (this._style[i(948)] = Cesium[i(1960)](this._style[i(948)], !0)),
          (this._style[i(1105)] = Cesium[i(1960)](this[i(343)][i(1105)], !0)),
          (this[i(961)] = Cesium[i(1960)](this._style[i(961)], !0)),
          (this[i(2117)] = Cesium[i(1960)](this[i(343)][i(2117)], Cesium[i(1154)][i(588)])),
          (this[i(911)] = Cesium[i(1960)](this._style[i(911)], 5)),
          (this[i(1274)] = Cesium.defaultValue(this[i(343)][i(1274)], !1)),
          (this[i(343)][i(2071)] = Cesium[i(1960)](this[i(343)][i(2071)], !0)),
          (this[i(343)][i(2e3)] = Cesium[i(1960)](this._style[i(2e3)], i(582))),
          (this._style[i(2288)] = Cesium[i(1960)](this[i(343)].scanPlaneRate, 10)),
          (this[i(1296)] = Cesium[i(1637)][i(1082)](i(1154))),
          (this[i(1296)]._uniforms[i(1070)] = Cesium.Color[i(971)][i(329)](0.4)),
          (this[i(1541)] = 0),
          (this[i(1616)] = 0),
          (this[i(1233)] = Cesium[i(1665)][i(1710)]()),
          (this[i(860)] = new Cesium.BoundingSphere()),
          (this._scanRadialCommand = void 0),
          (this[i(1897)] = !0),
          (this[i(672)] = []),
          (this[i(2169)] = void 0),
          (this._backFaceRS = void 0),
          (this[i(2536)] = void 0)
        let s = Cesium[i(1960)](e[i(2251)], [120, 40, 0])
        this.setPosition(s)
      }
      [t(1253)](e) {
        const i = t
        ;(this[i(1912)] = e),
          (this[i(2238)] = Cesium.Cartesian3[i(667)](e[0], e[1], e[2])),
          (this[i(794)] = this[i(1912)]),
          (this[i(1897)] = !0),
          (this[i(1024)] = new Cesium.BoundingSphere(this._cartesian3, 1.2 * this[i(343)][i(1981)]))
      }
      [t(2231)]() {
        const e = t
        this._isAdd &&
          ((this._needUpdate = !0), this[e(684)](this[e(2462)]), this[e(1768)](this[e(2462)])),
          (this._boundingSphere = new Cesium[e(1242)](this[e(2238)], 1.2 * this._style[e(1981)]))
      }
      [t(940)](e) {
        const i = t,
          s = {}
        ;(s[i(2416)] = this), (s.id = this[i(1570)])
        let n = e[i(2325)][i(896)](s)
        ;(this[i(1943)] = n),
          (this[i(1208)] = this[i(1943)][i(1070)]),
          (this[i(2202)] = new Cesium[i(905)]({
            owner: this,
            primitiveType: Cesium[i(1099)].TRIANGLES,
            boundingVolume: this[i(860)],
            pickId: i(1208)
          })),
          (this[i(1619)] = new Cesium[i(905)]({
            owner: this,
            primitiveType: Cesium[i(1099)][i(2367)],
            boundingVolume: this[i(860)]
          })),
          (this._sectorVA = void 0),
          (this._sectorLineCommand = new Cesium[i(905)]({
            owner: this,
            primitiveType: Cesium[i(1099)][i(1095)],
            boundingVolume: this[i(860)]
          })),
          (this._sectorLineVA = void 0),
          (this[i(2521)] = new Cesium[i(905)]({
            owner: this,
            primitiveType: Cesium.PrimitiveType[i(1095)],
            boundingVolume: this[i(860)]
          })),
          (this[i(2566)] = void 0),
          (this[i(1360)] = new Cesium[i(905)]({
            owner: this,
            primitiveType: Cesium[i(1099)][i(2367)],
            boundingVolume: this._boundingSphereWC,
            pickId: 'czm_pickColor'
          })),
          (this[i(521)] = new Cesium[i(905)]({
            owner: this,
            primitiveType: Cesium.PrimitiveType.TRIANGLES,
            boundingVolume: this._boundingSphereWC
          })),
          (this[i(698)] = void 0),
          (this[i(993)] = new Cesium[i(905)]({
            owner: this,
            primitiveType: Cesium[i(1099)][i(1095)],
            boundingVolume: this[i(860)]
          })),
          (this._domeLineVA = void 0),
          (this[i(2573)] = new Cesium[i(905)]({
            owner: this,
            primitiveType: Cesium[i(1099)][i(2367)],
            boundingVolume: this[i(860)]
          })),
          (this._scanPlaneBackCommand = new Cesium[i(905)]({
            owner: this,
            primitiveType: Cesium[i(1099)][i(2367)],
            boundingVolume: this[i(860)]
          }))
      }
      [t(1459)](e) {
        const i = t
        let s = this
        const n = {}
        ;(n[i(929)] = function () {
          return 0
        }),
          (n.u_xHalfAngle = function () {
            return s[i(447)]
          }),
          (n[i(1358)] = function () {
            return s.yHalfAngle
          }),
          (n[i(880)] = function () {
            const t = i
            return s._style[t(1981)]
          }),
          (n[i(1777)] = function () {
            return s[i(1274)]
          }),
          (n[i(2085)] = function () {
            return s[i(961)]
          }),
          (n[i(1604)] = function () {
            return s[i(2117)]
          }),
          (n[i(2282)] = function () {
            return s[i(911)]
          }),
          (n[i(2551)] = function () {
            return 1
          }),
          (n[i(2339)] = function () {
            return s.lineColor
          }),
          (n[i(202)] = function () {
            const t = i
            return s[t(468)] ? s[t(506)] : s[t(1070)]
          }),
          (n.czm_pickColor = function () {
            return s.czm_pickColor
          })
        const o = {}
        ;(o[i(1208)] = function () {
          return s[i(1208)]
        }),
          (o.u_xHalfAngle = function () {
            return s[i(1541)]
          }),
          (o[i(1358)] = function () {
            return s._scanePlaneYHalfAngle
          }),
          (o[i(880)] = function () {
            const t = i
            return s[t(343)][t(1981)]
          }),
          (o[i(202)] = function () {
            return s[i(2e3)]
          }),
          (o[i(1777)] = function () {
            return s[i(1274)]
          }),
          (o.u_showIntersection = function () {
            return s[i(961)]
          }),
          (o.u_intersectionColor = function () {
            return s[i(2117)]
          }),
          (o[i(2282)] = function () {
            return s[i(911)]
          }),
          (o.u_normalDirection = function () {
            return 1
          }),
          (o[i(2339)] = function () {
            return s.lineColor
          }),
          (this._uniforms = n),
          (this._scanUniforms = o)
      }
      [t(1451)](e) {
        const i = t
        ;(this[i(447)] = Cesium[i(475)].toRadians(this._style.hAngle / 2)),
          (this.yHalfAngle = Cesium[i(475)].toRadians(this[i(343)].vAngle / 2)),
          (this[i(2071)] = this[i(343)].showScanPlane),
          (this[i(2e3)] = Cesium.Color[i(2008)](this[i(343)][i(2e3)])),
          (this[i(1070)] = Cesium[i(1154)].fromCssColorString(this[i(343)][i(1070)])),
          (this.selectedColor = Cesium[i(1154)][i(2008)](this._style[i(506)])),
          (this[i(2248)] = this[i(343)][i(2248)]),
          (this[i(1018)] = Cesium[i(1154)][i(2008)](this[i(343)][i(1018)])),
          (this[i(948)] = this[i(343)][i(948)]),
          (this[i(1105)] = this[i(343)][i(1105)]),
          this[i(1459)](e)
        let s = new Cesium[i(183)](
            Cesium[i(475)][i(1149)](this[i(343)].heading),
            Cesium.Math[i(1149)](this._style[i(711)]),
            Cesium.Math[i(1149)](this[i(343)].roll)
          ),
          n = Cesium.Transforms[i(950)](this[i(2238)], s),
          o = Cesium[i(2066)][i(942)](Cesium[i(819)][i(1022)](n), this[i(2238)]),
          r = new Cesium[i(1242)](Cesium[i(310)][i(1738)], this[i(343)][i(1981)])
        ;(this[i(1636)] = o),
          Cesium[i(2066)][i(1868)](this[i(1636)], this._style[i(1981)], this[i(1988)]),
          Cesium[i(1242)][i(1863)](r, this[i(1636)], this[i(860)])
        let a = this[i(1274)]
        !(function (t, s) {
          const n = i
          let o = e[n(2325)],
            r = (function (t, e) {
              const i = n
              let s = t[i(447)],
                o = t[i(2559)],
                r = e[i(2400)],
                a = e[i(1295)],
                h = [],
                l = Cesium[i(819)][i(1693)](s, dt)
              return (
                h[i(2553)](
                  r[i(2512)](function (t) {
                    const e = i
                    return Cesium.Matrix3[e(1142)](l, t, new Cesium[e(310)]())
                  })
                ),
                (l = Cesium[i(819)][i(2377)](-o, dt)),
                h.push(
                  a[i(2512)](function (t) {
                    const e = i
                    return Cesium[e(819)].multiplyByVector(l, t, new Cesium[e(310)]())
                  })[i(1027)]()
                ),
                (l = Cesium.Matrix3[i(1693)](-s, dt)),
                h[i(2553)](
                  r[i(2512)](function (t) {
                    const e = i
                    return Cesium[e(819)][e(1142)](l, t, new Cesium[e(310)]())
                  })[i(1027)]()
                ),
                (l = Cesium[i(819)][i(2377)](o, dt)),
                h[i(2553)](
                  a[i(2512)](function (t) {
                    const e = i
                    return Cesium[e(819)].multiplyByVector(l, t, new Cesium[e(310)]())
                  })
                ),
                h
              )
            })(t, Ct(t, t.xHalfAngle, t[n(2559)]))
          if (
            ((t._sectorVA = (function (t, e) {
              const i = n
              let s = Array[i(1727)][i(1500)][i(560)]([], e).length - e[i(277)],
                o = new Float32Array(18 * s),
                r = 0
              for (let t = 0, s = e[i(277)]; t < s; t++) {
                let s = e[t],
                  n = Cesium[i(310)][i(379)](Cesium[i(310)][i(2196)](s[0], s[s.length - 1], ft), ft)
                for (let t = 0, e = s[i(277)] - 1; t < e; t++)
                  (o[r++] = 0),
                    (o[r++] = 0),
                    (o[r++] = 0),
                    (o[r++] = -n.x),
                    (o[r++] = -n.y),
                    (o[r++] = -n.z),
                    (o[r++] = s[t].x),
                    (o[r++] = s[t].y),
                    (o[r++] = s[t].z),
                    (o[r++] = -n.x),
                    (o[r++] = -n.y),
                    (o[r++] = -n.z),
                    (o[r++] = s[t + 1].x),
                    (o[r++] = s[t + 1].y),
                    (o[r++] = s[t + 1].z),
                    (o[r++] = -n.x),
                    (o[r++] = -n.y),
                    (o[r++] = -n.z)
              }
              let a = Cesium[i(2047)][i(1428)]({
                  context: t,
                  typedArray: o,
                  usage: Cesium[i(1437)].STATIC_DRAW
                }),
                h = 6 * Float32Array[i(369)],
                l = [
                  {
                    index: ht[i(2251)],
                    vertexBuffer: a,
                    componentsPerAttribute: 3,
                    componentDatatype: Cesium[i(2319)].FLOAT,
                    offsetInBytes: 0,
                    strideInBytes: h
                  },
                  {
                    index: ht.normal,
                    vertexBuffer: a,
                    componentsPerAttribute: 3,
                    componentDatatype: Cesium.ComponentDatatype[i(2441)],
                    offsetInBytes: 3 * Float32Array[i(369)],
                    strideInBytes: h
                  }
                ]
              const c = {}
              return (c.context = t), (c.attributes = l), new Cesium[i(1057)](c)
            })(o, r)),
            t[n(2248)] &&
              (t[n(1887)] = (function (t, e) {
                const i = n
                let s = e[i(277)],
                  o = new Float32Array(9 * s),
                  r = 0
                for (let t = 0, s = e[i(277)]; t < s; t++) {
                  let i = e[t]
                  ;(o[r++] = 0),
                    (o[r++] = 0),
                    (o[r++] = 0),
                    (o[r++] = i[0].x),
                    (o[r++] = i[0].y),
                    (o[r++] = i[0].z)
                }
                let a = Cesium[i(2047)].createVertexBuffer({
                    context: t,
                    typedArray: o,
                    usage: Cesium.BufferUsage[i(2200)]
                  }),
                  h = 3 * Float32Array[i(369)],
                  l = [
                    {
                      index: ht[i(2251)],
                      vertexBuffer: a,
                      componentsPerAttribute: 3,
                      componentDatatype: Cesium[i(2319)].FLOAT,
                      offsetInBytes: 0,
                      strideInBytes: h
                    }
                  ]
                const c = {}
                return (c[i(2325)] = t), (c[i(1548)] = l), new Cesium[i(1057)](c)
              })(o, r)),
            t[n(2248)] &&
              (t._sectorSegmentLineVA = (function (t, e) {
                const i = n
                let s = Array[i(1727)][i(1500)][i(560)]([], e).length - e.length,
                  o = new Float32Array(9 * s),
                  r = 0
                for (let t = 0, s = e[i(277)]; t < s; t++) {
                  let s = e[t]
                  for (let t = 0, e = s[i(277)] - 1; t < e; t++)
                    (o[r++] = s[t].x),
                      (o[r++] = s[t].y),
                      (o[r++] = s[t].z),
                      (o[r++] = s[t + 1].x),
                      (o[r++] = s[t + 1].y),
                      (o[r++] = s[t + 1].z)
                }
                let a = Cesium[i(2047)][i(1428)]({
                    context: t,
                    typedArray: o,
                    usage: Cesium.BufferUsage.STATIC_DRAW
                  }),
                  h = 3 * Float32Array[i(369)],
                  l = [
                    {
                      index: ht[i(2251)],
                      vertexBuffer: a,
                      componentsPerAttribute: 3,
                      componentDatatype: Cesium[i(2319)][i(2441)],
                      offsetInBytes: 0,
                      strideInBytes: h
                    }
                  ]
                const c = {}
                return (c[i(2325)] = t), (c[i(1548)] = l), new Cesium[i(1057)](c)
              })(o, r)),
            t.topShow &&
              ((t[n(698)] = (function (t) {
                const e = n
                let i = Cesium[e(2487)][e(933)](
                  new Cesium.EllipsoidGeometry({
                    vertexFormat: Cesium[e(1788)][e(1936)],
                    stackPartitions: 32,
                    slicePartitions: 32
                  })
                )
                return Cesium[e(1057)][e(1225)]({
                  context: t,
                  geometry: i,
                  attributeLocations: ht,
                  bufferUsage: Cesium[e(1437)][e(2200)],
                  interleave: !1
                })
              })(o)),
              t.topline &&
                (t._domeLineVA = (function (t) {
                  const e = n
                  let i = Cesium[e(1083)][e(933)](
                    new Cesium[e(1083)]({
                      vertexFormat: Cesium.VertexFormat[e(1936)],
                      stackPartitions: 32,
                      slicePartitions: 32
                    })
                  )
                  return Cesium[e(1057)][e(1225)]({
                    context: t,
                    geometry: i,
                    attributeLocations: ht,
                    bufferUsage: Cesium[e(1437)][e(2200)],
                    interleave: !1
                  })
                })(o))),
            t.showScanPlane)
          )
            if (n(333) == t.scanPlaneMode) {
              let e = Ct(t, Cesium[n(475)][n(2148)], 0)
              t._scanPlaneVA = vt(o, e.zox)
            } else {
              let e = Ct(t, 0, Cesium[n(475)][n(2148)])
              t[n(1695)] = vt(o, e[n(2400)])
            }
        })(this),
          (function (t, e, s) {
            const n = i,
              o = {}
            o[n(299)] = !e
            const r = {}
            r.enabled = !e
            const a = {}
            a[n(299)] = !e
            const h = {}
            h.enabled = !e
            const l = {}
            ;(l[n(1473)] = h), (l[n(898)] = !0)
            const c = {}
            c[n(299)] = !0
            const u = {}
            ;(u[n(1473)] = c),
              (u[n(898)] = !0),
              []
                ? ((t._frontFaceRS = Cesium[n(515)][n(2407)]({
                    depthTest: o,
                    depthMask: !1,
                    blending: Cesium.BlendingState[n(1192)],
                    cull: { enabled: !0, face: Cesium[n(1749)][n(274)] }
                  })),
                  (t[n(2499)] = Cesium[n(515)][n(2407)]({
                    depthTest: r,
                    depthMask: !1,
                    blending: Cesium[n(1870)].ALPHA_BLEND,
                    cull: { enabled: !0, face: Cesium[n(1749)][n(1402)] }
                  })),
                  (t[n(934)] = Cesium[n(515)][n(2407)]({
                    depthTest: a,
                    depthMask: !1,
                    blending: Cesium.BlendingState[n(1785)]
                  })))
                : ((t[n(2169)] = Cesium[n(515)][n(2407)](l)),
                  (t[n(934)] = Cesium[n(515)].fromCache(u)))
          })(this, a),
          (function (t, e, s) {
            const n = i
            ;(function (t, e, i) {
              const s = a0_0x3b79,
                n = {}
              n.sources = [lt, i[s(239)], mt]
              let o = e[s(2325)],
                r = ct,
                a = new Cesium[s(231)](n)
              t[s(2536)] = Cesium[s(178)][s(1879)]({
                context: o,
                shaderProgram: t._sp,
                vertexShaderSource: r,
                fragmentShaderSource: a,
                attributeLocations: ht
              })
            })(t, e, s),
              t[n(2071)] &&
                (function (t, e, i) {
                  const s = n,
                    o = {}
                  o[s(1177)] = [lt, i[s(239)], ut]
                  let r = e.context,
                    a = ct,
                    h = new Cesium[s(231)](o)
                  t[s(2118)] = Cesium.ShaderProgram[s(1879)]({
                    context: r,
                    shaderProgram: t[s(2118)],
                    vertexShaderSource: a,
                    fragmentShaderSource: h,
                    attributeLocations: ht
                  })
                })(t, e, s)
          })(this, e, this[i(1296)]),
          (function (t, e) {
            const s = i
            t._drawCommands[s(277)] = 0
            let n = e ? Cesium.Pass[s(2403)] : Cesium[s(1864)][s(787)]
            _t(
              t,
              t[s(2202)],
              t._sectorBackCommand,
              t[s(2169)],
              t[s(2499)],
              t._sp,
              t[s(267)],
              t[s(1685)],
              t._computedModelMatrix,
              e,
              n
            ),
              t[s(2248)] &&
                _t(
                  t,
                  t[s(411)],
                  void 0,
                  t[s(2169)],
                  t._backFaceRS,
                  t._sp,
                  t[s(1887)],
                  t._uniforms,
                  t._computedModelMatrix,
                  e,
                  n,
                  !0
                ),
              t[s(2248)] &&
                _t(
                  t,
                  t[s(2521)],
                  void 0,
                  t._frontFaceRS,
                  t[s(2499)],
                  t[s(2536)],
                  t[s(2566)],
                  t[s(1685)],
                  t[s(1988)],
                  e,
                  n,
                  !0
                ),
              t[s(948)] &&
                (_t(
                  t,
                  t._domeFrontCommand,
                  t[s(521)],
                  t[s(2169)],
                  t[s(2499)],
                  t[s(2536)],
                  t._domeVA,
                  t._uniforms,
                  t[s(1988)],
                  e,
                  n
                ),
                t[s(1105)] &&
                  _t(
                    t,
                    t[s(993)],
                    void 0,
                    t[s(2169)],
                    t._backFaceRS,
                    t[s(2536)],
                    t[s(1791)],
                    t[s(1685)],
                    t[s(1988)],
                    e,
                    n,
                    !0
                  )),
              t[s(2071)] &&
                _t(
                  t,
                  t[s(2573)],
                  t[s(2478)],
                  t[s(2169)],
                  t[s(2499)],
                  t[s(2118)],
                  t[s(1695)],
                  t[s(2372)],
                  t[s(1585)],
                  e,
                  n
                )
          })(this, !0)
      }
      [t(720)](e) {
        const i = t
        if (!this[i(1144)] || e[i(1312)] !== Cesium[i(1857)][i(339)]) return
        if (
          (this[i(1897)] && (this[i(940)](e), this._createDrawCommands(e), (this._needUpdate = !1)),
          this[i(343)][i(2071)])
        ) {
          let t = e.time,
            s = Cesium[i(1665)][i(2479)](t, this._time)
          s < 0 && (this._time = Cesium[i(1665)][i(1902)](t, this[i(1233)]))
          let n,
            o = Math[i(1113)]((s % this[i(343)].scanPlaneRate) / this[i(343)][i(2288)], 0),
            r = this[i(447)],
            a = this[i(2559)]
          n = 2 * r * o - r
          let h = Math[i(1824)](a),
            l = Math[i(1272)](n),
            c = Math[i(1924)](l * h)
          ;(this[i(1541)] = n),
            (this[i(1616)] = c),
            Cesium[i(819)][i(1693)](this[i(1541)], dt),
            Cesium[i(2066)][i(383)](this[i(1636)], dt, this._computedScanPlaneModelMatrix),
            Cesium[i(2066)].multiplyByUniformScale(
              this[i(1585)],
              this[i(343)][i(1981)],
              this[i(1585)]
            )
        }
        let s = e.commandList
        e.passes
        let n = this._drawCommands
        for (let t = 0, e = n.length; t < e; t++) {
          let e = n[t]
          s[i(2553)](e)
        }
      }
      _addHook(t) {
        this._add(t)
      }
      [t(1768)](e) {
        const i = t
        ;(this._layer = e),
          (this[i(2122)] = e.id),
          (this[i(2366)] = this._id),
          e[i(245)][i(696)][i(1346)][i(1861)](this),
          (this[i(2549)] = !0)
      }
      [t(2389)](e) {
        this[t(684)](e)
      }
      _remove(e) {
        const i = t
        this[i(2549)] &&
          (e[i(245)][i(696)][i(1346)][i(1896)](this),
          (this[i(2549)] = !1),
          this[i(1943)] && this._pickId[i(2039)]())
      }
      [t(2039)]() {
        const e = t
        this[e(672)].forEach(function (t) {
          const i = e
          t[i(1239)] = t[i(1239)] && t[i(1239)][i(2039)]()
        }),
          (this._drawCommands = [])
      }
    }
    let dt = new Cesium.Matrix3(),
      ft = new Cesium[t(310)]()
    function Ct(e, i, s) {
      const n = t
      let o = e[n(1679)].slice,
        r = Math[n(1272)](s),
        a = Math[n(1824)](s),
        h = Math[n(1272)](i),
        l = Math[n(1824)](i),
        c = Math.atan(h * a),
        u = Math[n(1924)](r * l),
        m = []
      for (let t = 0; t < o; t++) {
        let e = (2 * c * t) / (o - 1) - c
        m.push(new Cesium.Cartesian3(0, Math[n(884)](e), Math[n(1272)](e)))
      }
      let p = []
      for (let t = 0; t < o; t++) {
        let e = (2 * u * t) / (o - 1) - u
        p[n(2553)](new Cesium[n(310)](Math[n(884)](e), 0, Math[n(1272)](e)))
      }
      const d = {}
      return (d[n(2400)] = m), (d.zox = p), d
    }
    function vt(e, i) {
      const s = t
      let n = i[s(277)] - 1,
        o = new Float32Array(9 * n),
        r = 0
      for (let t = 0; t < n; t++)
        (o[r++] = 0),
          (o[r++] = 0),
          (o[r++] = 0),
          (o[r++] = i[t].x),
          (o[r++] = i[t].y),
          (o[r++] = i[t].z),
          (o[r++] = i[t + 1].x),
          (o[r++] = i[t + 1].y),
          (o[r++] = i[t + 1].z)
      let a = Cesium[s(2047)][s(1428)]({
          context: e,
          typedArray: o,
          usage: Cesium[s(1437)][s(2200)]
        }),
        h = 3 * Float32Array[s(369)],
        l = [
          {
            index: ht[s(2251)],
            vertexBuffer: a,
            componentsPerAttribute: 3,
            componentDatatype: Cesium[s(2319)].FLOAT,
            offsetInBytes: 0,
            strideInBytes: h
          }
        ]
      const c = {}
      return (c[s(2325)] = e), (c.attributes = l), new Cesium[s(1057)](c)
    }
    function _t(e, i, s, n, o, r, a, h, l, c, u, m) {
      const p = t
      c &&
        s &&
        ((s.vertexArray = a),
        (s[p(1689)] = o),
        (s.shaderProgram = r),
        (s.uniformMap = Cesium[p(938)](h, e[p(1296)]._uniforms)),
        (s[p(1720)][p(2551)] = function () {
          return -1
        }),
        (s[p(134)] = u),
        (s.modelMatrix = l),
        e[p(672)].push(s)),
        (i[p(1239)] = a),
        (i[p(1689)] = n),
        (i[p(187)] = r),
        (i[p(1720)] = Cesium[p(938)](h, e[p(1296)][p(1685)])),
        m &&
          (i[p(1720)][p(929)] = function () {
            return 1
          }),
        (i[p(134)] = u),
        (i.modelMatrix = l),
        e[p(672)].push(i)
    }
    class gt extends b {
      constructor(e) {
        const i = t
        super(e),
          (this[i(528)] = '范围雷达'),
          (this[i(314)] = 'Point'),
          (this._graphicClassType = W),
          (this[i(2336)] = Q[i(1625)]),
          (this._fixPointCount = 1),
          (this[i(343)].outerRadius = Cesium[i(1960)](this._style[i(2334)], 100)),
          (this._style.innerRadius = Cesium[i(1960)](this[i(343)][i(384)], 10)),
          (this[i(343)][i(145)] = Cesium[i(1960)](this[i(343)].startFovH, 0)),
          (this[i(343)].endFovH = Cesium[i(1960)](this._style[i(1454)], 360)),
          (this[i(343)][i(727)] = Cesium.defaultValue(this[i(343)].startFovV, 0)),
          (this[i(343)][i(1439)] = Cesium[i(1960)](this[i(343)][i(1439)], 90)),
          (this[i(343)][i(2349)] = Cesium.defaultValue(this[i(343)][i(2349)], 20)),
          (this._style[i(1e3)] = Cesium[i(1960)](this._style.segmentV, 20)),
          (this[i(343)][i(1997)] = Cesium[i(1960)](this[i(343)][i(1997)], 3)),
          (this[i(343)][i(863)] = Cesium[i(1960)](this._style[i(863)], 3)),
          (this[i(1760)] = 0),
          (this._endFovH = 0),
          (this[i(1873)] = 0),
          (this[i(618)] = 0),
          (this._segmentH = 1),
          (this[i(1221)] = 1),
          (this[i(888)] = 1),
          (this[i(676)] = 1),
          (this[i(343)][i(1198)] = Cesium[i(1960)](this[i(343)].heading, 0)),
          (this._style[i(711)] = Cesium.defaultValue(this[i(343)].pitch, 0)),
          (this[i(343)][i(1143)] = Cesium[i(1960)](this._style[i(1143)], 0)),
          (this[i(343)][i(1070)] = Cesium.defaultValue(this[i(343)][i(1070)], i(1667))),
          (this[i(343)].selectedColor = Cesium.defaultValue(this._style[i(506)], i(1544))),
          (this[i(343)][i(2248)] = Cesium[i(1960)](this._style[i(2248)], !0)),
          (this[i(343)][i(1018)] = Cesium[i(1960)](
            this[i(343)][i(1018)],
            'rgba(255, 255, 255, 1)'
          )),
          (this._drawCommands = [])
        let s = Cesium[i(1960)](e[i(2251)], [120, 40, 0])
        this[i(1253)](s),
          (this[i(1636)] = Cesium[i(2066)][i(1902)](Cesium[i(2066)][i(1393)])),
          (this[i(1150)] = new Cesium.Quaternion()),
          (this[i(1359)] = new Cesium.Cartesian3()),
          (this[i(2035)] = new Cesium.Cartesian3(1, 1, 1)),
          (this._matrix = new Cesium.Matrix4()),
          (this._needUpdate = !0)
      }
      [t(1253)](e) {
        const i = t
        ;(this[i(1912)] = e),
          (this._cartesian3 = Cesium[i(310)][i(667)](e[0], e[1], e[2])),
          (this[i(794)] = this[i(1912)]),
          (this._needUpdate = !0),
          (this[i(1024)] = new Cesium[i(1242)](this._cartesian3, this[i(343)].outerRadius))
      }
      [t(2231)]() {
        const e = t
        ;(this._needUpdate = !0),
          (this[e(1024)] = new Cesium[e(1242)](this[e(2238)], this[e(343)][e(2334)]))
      }
      [t(1046)]() {
        const e = t
        ;(this._modelMatrix = Cesium[e(2058)][e(1648)](this[e(2238)])),
          Cesium[e(1472)][e(1245)](
            Cesium[e(183)][e(667)](
              this._style[e(1198)],
              this[e(343)][e(711)],
              this[e(343)][e(1143)]
            ),
            this._quaternion
          ),
          (this[e(967)] = Cesium[e(2066)][e(235)](
            this[e(1359)],
            this._quaternion,
            new Cesium[e(310)](1, 1, 1),
            this[e(967)]
          )),
          Cesium.Matrix4[e(332)](this[e(1636)], this[e(967)], this[e(967)]),
          (this._modelMatrix = this[e(967)])
      }
      [t(2265)](e) {
        const i = t
        var s = this[i(888)] * this[i(1373)],
          n = this[i(676)] * this[i(1221)],
          o = xt(this[i(1760)], this[i(440)], this[i(1873)], this[i(618)], s, n, this[i(2091)]),
          r = xt(this[i(1760)], this[i(440)], this._startFovV, this._endFovV, s, n, this[i(2091)]),
          a = bt(s, n),
          h = yt(this._segmentH, this[i(1221)], this[i(888)], this[i(676)])
        return this[i(420)](e, o, r, a, h)
      }
      [t(1011)](e) {
        const i = t
        var s = this[i(888)] * this._segmentH,
          n = this[i(676)] * this[i(1221)],
          o = xt(this[i(1760)], this[i(440)], this[i(1873)], this[i(618)], s, n, this[i(680)]),
          r = xt(
            this[i(1760)],
            this[i(440)],
            this[i(1873)],
            this[i(618)],
            s,
            n,
            this._innerFovRadiusPairs
          ),
          a = bt(s, n),
          h = yt(this[i(1373)], this[i(1221)], this._subSegmentH, this[i(676)])
        return this[i(420)](e, o, r, a, h)
      }
      [t(749)](e) {
        const i = t
        var s = this[i(676)] * this[i(1221)],
          n = St(this._startFovH, this[i(1873)], this[i(618)], 10, s, this[i(680)], this[i(2091)]),
          o = St(this._startFovH, this[i(1873)], this._endFovV, 10, s, this[i(680)], this[i(2091)]),
          r = bt(10, s),
          a = yt(10, this[i(1221)], 1, this[i(676)])
        return this._createRawCommand(e, n, o, r, a)
      }
      _createRightCrossSectionCommand(e) {
        const i = t
        var s = this[i(676)] * this[i(1221)],
          n = St(this._endFovH, this[i(1873)], this[i(618)], 10, s, this[i(680)], this[i(2091)]),
          o = St(
            this[i(440)],
            this[i(1873)],
            this[i(618)],
            10,
            s,
            this[i(680)],
            this._outerFovRadiusPairs
          ),
          r = bt(10, s),
          a = yt(10, this._segmentV, 1, this[i(676)])
        return this._createRawCommand(e, n, o, r, a)
      }
      [t(420)](e, i, s, n, o) {
        const r = t
        ;(this.faceColor = this[r(468)]
          ? Cesium[r(1154)][r(2008)](this._style[r(506)])
          : Cesium[r(1154)].fromCssColorString(this._style.color)),
          (this[r(1018)] = Cesium.Color[r(2008)](this[r(343)][r(1018)]))
        let a,
          h = this[r(1943)]
        const l = {}
        ;(l[r(2416)] = this),
          (l.id = this[r(1570)]),
          h || ((h = e[r(896)](l)), (this._pickId = h)),
          (a = h[r(1070)])
        const c = {}
        c[r(1177)] = [r(2335)]
        const u = {
            sources: [
              '\n            in vec3 v_positionEC;\n            in vec3 v_normalEC;\n            // varying vec2 v_st;\n            \n            // uniform sampler2D myImage;\n            uniform vec4 u_color;\n            uniform vec4 czm_pickColor;\n            void main()\n            {\n                vec3 positionToEyeEC = -v_positionEC;\n            \n                vec3 normalEC = normalize(v_normalEC);\n                #ifdef FACE_FORWARD\n                    normalEC = faceforward(normalEC, vec3(0.0, 0.0, 1.0), -normalEC);\n                #endif\n            \n                czm_materialInput materialInput;\n                materialInput.normalEC = normalEC;\n                materialInput.positionToEyeEC = positionToEyeEC;\n                // materialInput.st = v_st;\n            \n                //czm_material material = czm_getMaterial(materialInput);\n                czm_material material = czm_getDefaultMaterial(materialInput);\n                // material.diffuse = texture2D(myImage, materialInput.st).rgb;\n                material.diffuse = u_color.rgb;\n                material.alpha = u_color.a;\n            \n                #ifdef FLAT\n                    out_FragColor = vec4(material.diffuse + material.emission, material.alpha);\n                #else\n                   out_FragColor = czm_phong(normalize(positionToEyeEC), material, czm_lightDirectionEC);\n                #endif\n            }'
            ]
          },
          m = {}
        ;(m[r(226)] = 1), (m[r(2251)] = 0)
        var p = this,
          d = Cesium[r(1469)].getDefaultRenderState(!0, !1, void 0),
          f = Cesium[r(515)].fromCache(d),
          C = new Cesium[r(231)](c),
          v = new Cesium[r(231)](u),
          _ = {
            u_color: function () {
              return p[r(213)]
            },
            czm_pickColor: function () {
              return a
            }
          },
          g = {
            u_color: function () {
              return p.lineColor
            },
            czm_pickColor: function () {
              return a
            }
          },
          y = Cesium[r(178)].fromCache({
            context: e,
            vertexShaderSource: C,
            fragmentShaderSource: v,
            attributeLocations: m
          }),
          w = Cesium[r(2047)][r(1428)]({
            context: e,
            typedArray: i,
            usage: Cesium[r(1437)].STATIC_DRAW
          }),
          x = Cesium[r(2047)][r(1428)]({
            context: e,
            typedArray: s,
            usage: Cesium[r(1437)][r(2200)]
          }),
          b = Cesium[r(2047)][r(2151)]({
            context: e,
            typedArray: n,
            usage: Cesium[r(1437)][r(2200)],
            indexDatatype: Cesium.IndexDatatype[r(2061)]
          }),
          S = Cesium[r(2047)][r(2151)]({
            context: e,
            typedArray: o,
            usage: Cesium[r(1437)][r(2200)],
            indexDatatype: Cesium[r(2216)][r(2061)]
          }),
          P = new Cesium[r(1057)]({
            context: e,
            attributes: [
              {
                index: 0,
                vertexBuffer: w,
                componentsPerAttribute: 3,
                componentDatatype: Cesium.ComponentDatatype[r(2441)]
              },
              {
                index: 1,
                vertexBuffer: x,
                componentsPerAttribute: 3,
                componentDatatype: Cesium[r(2319)][r(2441)]
              }
            ],
            indexBuffer: b
          }),
          M = new Cesium[r(1057)]({
            context: e,
            attributes: [
              {
                index: 0,
                vertexBuffer: w,
                componentsPerAttribute: 3,
                componentDatatype: Cesium[r(2319)][r(2441)]
              },
              {
                index: 1,
                vertexBuffer: x,
                componentsPerAttribute: 3,
                componentDatatype: Cesium[r(2319)].FLOAT
              }
            ],
            indexBuffer: S
          }),
          A = Cesium[r(1242)][r(852)](i)
        return {
          command: new Cesium.DrawCommand({
            vertexArray: P,
            primitiveType: Cesium.PrimitiveType[r(2367)],
            renderState: f,
            shaderProgram: y,
            uniformMap: _,
            owner: this,
            pass: Cesium[r(1864)][r(2403)],
            modelMatrix: new Cesium.Matrix4(),
            boundingVolume: new Cesium[r(1242)](),
            cull: !0,
            pickId: r(1208)
          }),
          lineCommand: new Cesium[r(905)]({
            vertexArray: M,
            primitiveType: Cesium.PrimitiveType[r(1095)],
            renderState: f,
            shaderProgram: y,
            uniformMap: g,
            owner: this,
            pass: Cesium[r(1864)].TRANSLUCENT,
            modelMatrix: new Cesium[r(2066)](),
            boundingVolume: new Cesium[r(1242)](),
            cull: !0
          }),
          initBoundingSphere: A
        }
      }
      [t(937)]() {
        const e = t
        let i = this._style[e(384)]
        return [
          { fov: Cesium.Math.toRadians(0), radius: i },
          { fov: Cesium[e(475)][e(1149)](10), radius: 0.9 * i },
          { fov: Cesium.Math[e(1149)](20), radius: 0.8 * i },
          { fov: Cesium[e(475)][e(1149)](30), radius: 0.7 * i },
          { fov: Cesium[e(475)][e(1149)](40), radius: 0.6 * i },
          { fov: Cesium[e(475)][e(1149)](50), radius: 0.5 * i },
          { fov: Cesium[e(475)][e(1149)](60), radius: 0.4 * i },
          { fov: Cesium.Math[e(1149)](70), radius: 0.3 * i },
          { fov: Cesium[e(475)][e(1149)](80), radius: 0.1 * i },
          { fov: Cesium.Math[e(1149)](90), radius: 0.01 * i }
        ]
      }
      [t(777)]() {
        const e = t
        let i = this[e(343)][e(2334)]
        return [
          { fov: Cesium[e(475)][e(1149)](0), radius: i },
          { fov: Cesium.Math.toRadians(10), radius: 0.9 * i },
          { fov: Cesium[e(475)].toRadians(20), radius: 0.8 * i },
          { fov: Cesium[e(475)][e(1149)](30), radius: 0.7 * i },
          { fov: Cesium.Math[e(1149)](40), radius: 0.6 * i },
          { fov: Cesium[e(475)][e(1149)](50), radius: 0.5 * i },
          { fov: Cesium[e(475)][e(1149)](60), radius: 0.4 * i },
          { fov: Cesium[e(475)][e(1149)](70), radius: 0.3 * i },
          { fov: Cesium[e(475)][e(1149)](80), radius: 0.1 * i },
          { fov: Cesium[e(475)].toRadians(90), radius: 0.01 * i }
        ]
      }
      [t(1451)](e) {
        const i = t
        ;(this._drawCommands = []),
          (this[i(680)] = this._getInnerFovRadiusParire()),
          (this[i(2091)] = this._getOuterFovRadiusPairs()),
          (this[i(1760)] = Cesium.Math[i(1149)](this[i(343)][i(145)])),
          (this[i(440)] = Cesium[i(475)][i(1149)](this[i(343)][i(1454)])),
          (this[i(1873)] = Cesium[i(475)][i(1149)](this[i(343)][i(727)])),
          (this[i(618)] = Cesium.Math.toRadians(this[i(343)][i(1439)])),
          (this._segmentH = this[i(343)][i(2349)]),
          (this[i(1221)] = this[i(343)][i(1e3)]),
          (this[i(888)] = this._style.subSegmentH),
          (this[i(676)] = this._style[i(863)]),
          this[i(1046)](),
          this[i(672)][i(2553)](this[i(2265)](e[i(2325)])),
          this._drawCommands[i(2553)](this[i(749)](e[i(2325)])),
          this[i(672)][i(2553)](this[i(2544)](e[i(2325)])),
          this._drawCommands[i(2553)](this[i(1011)](e.context)),
          this[i(672)][i(1602)]((t) => {
            const e = i
            ;(t[e(597)][e(223)] = this[e(1636)]),
              (t[e(597)][e(2172)] = new Cesium[e(1242)](this._cartesian3, this[e(343)][e(2334)])),
              (t.lineCommand.modelMatrix = Cesium[e(2066)][e(1393)]),
              (t.lineCommand[e(223)] = this[e(1636)]),
              (t.lineCommand.boundingVolume = new Cesium[e(1242)](
                this[e(2238)],
                this[e(343)][e(2334)]
              ))
          })
      }
      update(e) {
        const i = t
        this[i(1144)] &&
          (this._needUpdate && (this._createDrawCommands(e), (this[i(1897)] = !1)),
          e[i(1312)] == Cesium[i(1857)][i(339)] &&
            this._drawCommands[i(1602)]((t) => {
              const s = i
              t[s(597)] && e[s(331)][s(2553)](t[s(597)]),
                this._style[s(2248)] && t[s(1739)] && e[s(331)][s(2553)](t[s(1739)])
            }))
      }
      [t(2039)]() {}
      [t(1545)](e) {
        const i = t
        ;(this[i(2462)] = e),
          (this[i(2122)] = e.id),
          (this.graphicId = this[i(1570)]),
          e[i(245)].scene[i(1346)][i(1861)](this)
      }
      [t(2389)](e) {
        const i = t
        e[i(245)][i(696)][i(1346)].remove(this),
          this[i(1943)] && this._pickId[i(2039)](),
          this._drawCommands[i(1602)](function (t) {
            const e = i
            t.vertexArray = t[e(1239)] && t.vertexArray[e(2039)]()
          }),
          (this[i(672)] = [])
      }
    }
    function yt(t, e, i, s) {
      var n = t * i,
        o = e * s,
        r = new Uint16Array((t + 1) * (2 * o) + (e + 1) * (2 * n) + 8)
      for (c = 0; c < t + 1; ++c)
        for (var a = 0; a < o; ++a) {
          var h = c * i
          ;(r[2 * (c * o + a)] = a * (n + 1) + h), (r[2 * (c * o + a) + 1] = (a + 1) * (n + 1) + h)
        }
      for (var l = 2 * (t + 1) * o, c = 0; c < e + 1; ++c)
        for (a = 0; a < n; ++a) {
          var u = c * s
          ;(r[l + 2 * (a + c * n)] = u * (n + 1) + a),
            (r[l + 2 * (a + c * n) + 1] = u * (n + 1) + a + 1)
        }
      return r
    }
    function wt(e, i) {
      const s = t
      var n = i[s(959)](function (t) {
        return t[s(1303)] > e
      })
      if (n > 0) {
        var o = i[n - 1],
          r = i[n],
          a = (e - o[s(1303)]) / (r.fov - o.fov)
        return o.radius * (1 - a) + r[s(1981)] * a
      }
    }
    function xt(t, e, i, s, n, o, r) {
      for (var a = new Float32Array((n + 1) * (o + 1) * 3), h = 0; h < n + 1; ++h)
        for (var l = 0; l < o + 1; ++l) {
          var c = Pt(i, s, o, l),
            u = Mt(Pt(t, e, n, h), c),
            m = r ? wt(c, r) : 1
          ;(a[3 * (l * (n + 1) + h)] = u[0] * m),
            (a[3 * (l * (n + 1) + h) + 1] = u[1] * m),
            (a[3 * (l * (n + 1) + h) + 2] = u[2] * m)
        }
      return a
    }
    function bt(t, e) {
      for (var i = new Uint16Array(t * e * 6), s = 0; s < t; ++s)
        for (var n = 0; n < e; ++n) {
          var o = n * (t + 1) + s,
            r = n * (t + 1) + s + 1,
            a = (n + 1) * (t + 1) + s,
            h = (n + 1) * (t + 1) + s + 1,
            l = 6 * (n * t + s)
          ;(i[l + 0] = o),
            (i[l + 1] = r),
            (i[l + 2] = h),
            (i[l + 3] = o),
            (i[l + 4] = h),
            (i[l + 5] = a)
        }
      return i
    }
    function St(t, e, i, s, n, o, r) {
      for (var a = new Float32Array((s + 1) * (n + 1) * 3), h = 0; h < s + 1; ++h)
        for (var l = 0; l < n + 1; ++l) {
          var c = Pt(e, i, n, l),
            u = Mt(t, c),
            m = Pt(o ? wt(c, o) : 1, r ? wt(c, r) : 1, s, h)
          ;(a[3 * (l * (s + 1) + h)] = u[0] * m),
            (a[3 * (l * (s + 1) + h) + 1] = u[1] * m),
            (a[3 * (l * (s + 1) + h) + 2] = u[2] * m)
        }
      return a
    }
    function Pt(t, e, i, s) {
      return t + (s / i) * (e - t)
    }
    function Mt(e, i) {
      const s = t
      var n = e,
        o = i
      return [Math.cos(-n) * Math[s(1272)](o), Math[s(884)](-n) * Math[s(1272)](o), Math[s(884)](o)]
    }
    const At = {}
    ;(At[t(900)] = null),
      (At.BeamRadar = ot),
      (At[t(473)] = nt),
      (At.ProbeRadar = rt),
      (At[t(2078)] = gt),
      (At[t(1680)] = tt),
      (At[t(2107)] = pt)
    const Tt = {}
    Tt[t(922)] = t(256)
    const Et = Object[t(1731)](Object.defineProperty(At, Symbol[t(1540)], Tt))
    class zt extends b {
      constructor(e = {}) {
        const i = t
        super(e), (this[i(1115)] = P), (this[i(314)] = i(1865))
        let s = Cesium[i(1960)](e.position, [111, 28, 0])
        ;(this._entity = this[i(1377)]()), this[i(1253)](s), (this._fixPointCount = 1)
      }
      get [t(769)]() {
        return this[t(987)]
      }
      [t(1253)](e) {
        const i = t
        e &&
          ((this._position = e),
          (this[i(794)] = this[i(2251)]),
          (this._cartesian3 = Cesium[i(310)][i(667)](
            this[i(1912)][0],
            this[i(1912)][1],
            this[i(1912)][2]
          )),
          (this._entity[i(2251)] = this[i(2238)]),
          (this[i(1024)] = new Cesium.BoundingSphere(this[i(2238)], 5)))
      }
      [t(2515)](e) {
        const i = t
        Array[i(1108)](e) && this[i(1253)](e[0])
      }
      _createEntity() {
        const e = t,
          i = {}
        i[e(1482)] = !1
        const s = {}
        return (
          (s[e(2366)] = this[e(1570)]),
          (s[e(1542)] = i),
          (s[e(2240)] = {}),
          (s[e(909)] = {}),
          new Cesium[e(1380)](s)
        )
      }
      _setVisible(e) {
        const i = t
        this._entity[i(1482)] = e
      }
      [t(1545)](e) {
        const i = t
        ;(this._layer = e), (this[i(987)][i(2122)] = e.id), e._viewer[i(118)].add(this[i(987)])
      }
      [t(2389)](e) {
        const i = t
        e[i(245)].entities[i(1896)](this._entity)
      }
    }
    class Dt {
      constructor(e = {}) {
        const i = t
        ;(this[i(1482)] = Cesium[i(1960)](e[i(1482)], !0)),
          (this.text = Cesium[i(1960)](e.text, '')),
          (this.fontSize = Cesium[i(1960)](e[i(735)], 64)),
          (this[i(2527)] = Cesium[i(1960)](e[i(2527)], i(510))),
          (this[i(1679)] = Cesium[i(1960)](e[i(1679)], Cesium[i(268)][i(1012)])),
          (this[i(2341)] = Cesium[i(1960)](e[i(2341)], 0.5)),
          (this[i(1875)] = Cesium.defaultValue(e[i(1875)], !0)),
          (this.backgroundColor = Cesium[i(1960)](e.backgroundColor, i(853))),
          (this[i(2103)] = Cesium.defaultValue(e[i(2103)], new Cesium[i(194)](7, 5))),
          (this[i(919)] = Cesium.defaultValue(e[i(919)], new Cesium[i(194)](0, -40))),
          (this[i(685)] = Cesium.defaultValue(e[i(685)], Cesium.Cartesian3[i(1738)])),
          (this[i(1491)] = Cesium[i(1960)](e[i(1491)], Cesium.HorizontalOrigin.CENTER)),
          (this[i(1990)] = Cesium.defaultValue(e[i(1990)], Cesium[i(1023)][i(2064)])),
          (this[i(697)] = Cesium[i(1960)](e[i(697)], Cesium[i(893)][i(2184)])),
          (this[i(2262)] = Cesium[i(1960)](e.fillColor, i(287))),
          (this[i(725)] = Cesium.defaultValue(e.outlineColor, '#080707')),
          (this[i(376)] = Cesium[i(1960)](e.outlineWidth, 4)),
          (this[i(780)] = Cesium[i(1960)](e.translucencyByDistance, null)),
          (this[i(2215)] = Cesium.defaultValue(e[i(2215)], null)),
          (this[i(435)] = Cesium.defaultValue(e[i(435)], null)),
          (this[i(978)] = Cesium.defaultValue(e.distanceDisplayCondition, null)),
          (this[i(1750)] = Cesium.defaultValue(e.disableDepthTestDistance, null))
      }
      [t(1567)](e) {
        const i = t
        for (const t in this) {
          const i = this[t]
          null != i && null != i && (e[t] = i)
        }
        ;(e[i(1291)] = this[i(735)] + i(2044) + this[i(2527)]),
          (e.backgroundColor = Cesium.Color[i(2008)](this[i(821)])),
          (e[i(2262)] = Cesium[i(1154)].fromCssColorString(this[i(2262)])),
          (e[i(725)] = Cesium.Color[i(2008)](this[i(725)]))
      }
    }
    class It {
      constructor(e = {}) {
        const i = t
        ;(this.show = Cesium[i(1960)](e[i(1482)], !0)),
          (this.image = Cesium[i(1960)](e.image, '')),
          (this[i(1323)] = Cesium[i(1960)](e.mergeImage, !0)),
          (this.scale = Cesium[i(1960)](e[i(2341)], 1)),
          (this.pixelOffset = Cesium[i(1960)](e[i(919)], new Cesium[i(194)](0, 0))),
          (this[i(685)] = Cesium.defaultValue(e[i(685)], Cesium[i(310)].ZERO)),
          (this[i(1491)] = Cesium[i(1960)](e.horizontalOrigin, Cesium[i(2324)].CENTER)),
          (this[i(1990)] = Cesium[i(1960)](e[i(1990)], Cesium[i(1023)].BOTTOM)),
          (this[i(697)] = Cesium[i(1960)](e.heightReference, Cesium[i(893)].NONE)),
          (this[i(1070)] = Cesium[i(1960)](e[i(1070)], i(287))),
          (this[i(2307)] = Cesium[i(1960)](e[i(2307)], 0)),
          (this[i(1054)] = Cesium.defaultValue(e[i(1054)], Cesium[i(310)][i(1738)])),
          (this[i(2177)] = Cesium[i(1960)](e[i(2177)], !1)),
          (this[i(575)] = Cesium[i(1960)](e[i(575)], null)),
          (this[i(2306)] = Cesium[i(1960)](e[i(2306)], null)),
          (this.scaleByDistance = Cesium.defaultValue(e[i(435)], null)),
          (this[i(2215)] = Cesium.defaultValue(e[i(2215)], null)),
          (this[i(780)] = Cesium.defaultValue(e[i(780)], null)),
          (this[i(978)] = Cesium[i(1960)](e.distanceDisplayCondition, null)),
          (this[i(1750)] = Cesium.defaultValue(e[i(1750)], null)),
          (this.imageSubRegion = Cesium.defaultValue(e[i(1696)], null))
      }
      [t(1567)](e) {
        const i = t
        let s = { ...this }
        delete s[i(259)]
        for (const t in s) {
          const i = s[t]
          null != i && null != i && (e[t] = i)
        }
        this[i(1323)] && this[i(259)] && (e.image = this[i(259)]),
          (e[i(1070)] = Cesium[i(1154)][i(2008)](this.color)),
          (e[i(2307)] = Cesium[i(475)].toRadians(this[i(2307)])),
          (e[i(1482)] = this[i(1482)] && this[i(259)])
      }
    }
    class kt extends zt {
      constructor(e = {}) {
        const i = t
        super(e),
          (this._typeName = i(707)),
          (this._graphicType = Q[i(1923)]),
          (this[i(343)] = Cesium[i(1960)](e[i(1679)], {})),
          (this[i(343)].selectedColor = Cesium[i(1960)](this._style[i(506)], i(179))),
          (this[i(2119)] = new Dt(e.label)),
          this[i(2119)][i(1567)](this._entity[i(2240)]),
          (this[i(556)] = new It(e[i(909)])),
          this[i(556)][i(1567)](this[i(987)][i(909)])
      }
      get [t(2240)]() {
        return this[t(2119)]
      }
      set label(e) {
        const i = t
        ;(this._label = e), this[i(2119)][i(1567)](this._entity[i(2240)]), this[i(1444)]()
      }
      get [t(909)]() {
        return this[t(556)]
      }
      set [t(909)](e) {
        const i = t
        ;(this[i(556)] = e), this[i(556)][i(1567)](this[i(987)][i(909)]), this[i(1444)]()
      }
      [t(1444)]() {}
      [t(2231)]() {
        const e = t
        let i = this._isSelected ? this[e(343)].selectedColor : this._label[e(2262)]
        this[e(987)].label[e(2262)] = Cesium[e(1154)][e(2008)](i)
      }
      [t(1114)]() {
        const e = t,
          i = {}
        return (i[e(2240)] = this[e(2119)]), (i.billboard = this[e(556)]), i
      }
    }
    let Ft = {
      readFeature(e) {
        const i = t
        let s = e[i(138)][i(2365)],
          n = e.properties
        return (n[i(2251)] = s), this[i(1951)](n)
      },
      create(e) {
        const i = t
        switch (e[i(1158)]) {
          case Q[i(1923)]:
            return new kt(e)
          case Q[i(1991)]:
            return new (class extends kt {
              constructor(e = {}) {
                const i = t
                super(e),
                  (this[i(528)] = i(2369)),
                  (this[i(2336)] = Q[i(1991)]),
                  (this._style.color = Cesium.defaultValue(this[i(343)][i(1070)], '#FF0000')),
                  (this[i(343)].pixelSize = Cesium.defaultValue(this[i(343)][i(225)], 20)),
                  (this[i(343)][i(2174)] = Cesium[i(1960)](this[i(343)][i(2174)], 50)),
                  (this[i(343)][i(164)] = Cesium.defaultValue(this[i(343)][i(164)], 20)),
                  (this._billboard[i(1750)] = Cesium[i(1960)](this[i(556)][i(1750)], 1e3)),
                  this[i(2231)]()
              }
              [t(2231)]() {
                const e = t
                let i = 1,
                  s = !0,
                  n = this[e(343)][e(225)],
                  o = !0,
                  r = 0.7,
                  a = !0,
                  h = this[e(343)][e(2174)]
                const l = this[e(987)]
                let c = this[e(1496)] ? this[e(343)][e(506)] : this[e(343)][e(1070)]
                ;(c = Cesium[e(1154)][e(2008)](c)),
                  (l[e(1542)] = {
                    show: !0,
                    color: new Cesium[e(2569)](
                      () => (
                        s ? ((i -= 0.03), i <= 0 && (s = !1)) : ((i = 1), (s = !0)), c[e(329)](i)
                      ),
                      !1
                    ),
                    pixelSize: new Cesium[e(2569)](
                      (t, e) => (o ? ((n += 2), n >= h && (o = !1)) : ((n = 10), (o = !0)), n),
                      !1
                    ),
                    outlineColor: new Cesium[e(2569)](
                      () => (
                        a ? ((r -= 0.035), r <= 0 && (a = !1)) : ((r = 0.7), (a = !0)), c[e(329)](r)
                      ),
                      !1
                    ),
                    outlineWidth: this._style[e(164)]
                  }),
                  this._billboard[e(978)] && (l[e(1542)][e(978)] = this[e(556)][e(978)])
              }
            })(e)
          case Q.poiMarker:
            return new (class extends kt {
              constructor(e) {
                const i = t
                super(e), (this[i(528)] = i(1790)), (this[i(2336)] = Q.poiMarker)
                let s = Cesium[i(1960)](e.label, {})
                ;(this[i(2119)][i(2262)] = Cesium[i(1960)](s[i(1070)], i(287))),
                  (this[i(2119)][i(821)] = Cesium.defaultValue(s.backgroundColor, i(916))),
                  (this[i(2119)][i(1982)] = Cesium[i(1960)](s[i(1982)], i(1790))),
                  (this[i(2119)][i(735)] = Cesium[i(1960)](s[i(735)], 20)),
                  (this[i(2119)][i(376)] = Cesium[i(1960)](s[i(376)], 0)),
                  (this._billboard[i(1323)] = !1)
                const n = this._billboard.image,
                  o = n[i(1028)](n.indexOf('.') + 1)
                if (n && ['png', i(1588), 'jpeg', 'svg'][i(877)](o) < 0)
                  throw new Cesium[i(2360)](i(2436))
                this[i(1444)]()
              }
              _updateGraphic() {
                const e = t,
                  i = this[e(556)].image
                ;(this[e(987)].label[e(1982)] = ''),
                  (this._entity[e(2240)][e(1482)] = !1),
                  (this[e(987)][e(909)][e(1482)] = i),
                  (this._entity[e(909)][e(1491)] = 1),
                  (this[e(987)][e(909)][e(919)] = new Cesium.Cartesian2(-30, -10)),
                  i && (i.indexOf('.svg') > 0 ? this[e(1543)](i) : this._initByImage(i))
              }
              [t(2231)]() {
                this._updateGraphic()
              }
              _initBySvg(e) {
                const i = t
                fetch(e)
                  [i(687)]((t) => t.text())
                  .then((t) => {
                    const e = i
                    let s = new DOMParser()[e(1079)](t, 'image/svg+xml')[e(630)](e(1396))[0]
                    if (!(s instanceof SVGElement)) throw new Cesium[e(2360)](e(1119))
                    {
                      this[e(1381)](s)
                      let t = new XMLSerializer().serializeToString(s),
                        i = e(2218) + btoa(unescape(encodeURIComponent(t)))
                      this[e(1178)](i)
                    }
                  })
                  [i(472)]((t) => {
                    const e = i
                    throw new Cesium[e(2360)](e(453))
                  })
              }
              [t(1381)](e) {
                const i = t
                for (let t = 0; t < e.childNodes[i(277)]; t++)
                  if (e[i(683)][t][i(1679)]) {
                    let s = this._isSelected ? this[i(343)][i(506)] : this[i(2119)].backgroundColor
                    ;(e[i(683)][t][i(1679)][i(2327)] = s), (e[i(683)][t].style.stroke = s)
                  }
              }
              [t(1178)](e) {
                const i = t
                let s = document[i(1945)](i(2575))
                ;(s[i(1744)] = e),
                  (s.onload = (t) => {
                    const e = i
                    this[e(830)] = s
                    let n = this._createImage()
                    n[e(1019)]('2d')[e(723)](s, 14, 14, 32, 32), (this[e(987)][e(909)][e(259)] = n)
                  }),
                  (s[i(1978)] = (t) => {
                    const e = i
                    throw (
                      ((this[e(987)][e(909)][e(259)] = this[e(2507)]()),
                      new Cesium[e(2360)](e(1249)))
                    )
                  })
              }
              [t(1189)]() {
                const e = t
                let i = this[e(1895)]()
                this[e(796)] = i
                let s = this[e(2507)](),
                  n = document[e(1945)](e(1493))
                ;(n[e(575)] = i[e(575)] + s[e(575)] + 10),
                  (n[e(2306)] = Math[e(1113)](i[e(2306)], s[e(2306)]))
                let o = n.getContext('2d')
                return o[e(723)](i, 0, 0), o[e(723)](s, i[e(575)], 5), n
              }
              _createMarkerBg() {
                const e = t
                let i = document.createElement(e(1493)),
                  s = i.getContext('2d')
                ;(i.width = 60),
                  (i[e(2306)] = 80),
                  s[e(1735)](),
                  s[e(1855)](30, 30, 28, 0, 2 * Math.PI, !0),
                  s[e(418)]()
                let n = this[e(1496)] ? this[e(343)][e(506)] : this[e(2119)].backgroundColor
                return (
                  (s[e(537)] = n),
                  s[e(2327)](),
                  s[e(1735)](),
                  s.moveTo(6, 45),
                  s[e(2140)](30, 80),
                  s.lineTo(54, 45),
                  (s[e(537)] = n),
                  s[e(2327)](),
                  s[e(1735)](),
                  s.arc(30, 30, 23, 0, 2 * Math.PI, !0),
                  s[e(418)](),
                  (s[e(537)] = e(741)),
                  s[e(2327)](),
                  i
                )
              }
              [t(2507)]() {
                const e = t,
                  i = this._label[e(1982)]
                let s = this[e(2119)][e(735)],
                  n = Math.ceil(i[e(277)] / 15),
                  o = []
                if (n > 1)
                  for (let t = 1; t <= n; t++) {
                    const s = i[e(1028)](15 * (t - 1), 15 * t)
                    o[e(2553)](s)
                  }
                else o[e(2553)](i)
                let r = Math[e(1429)](15, i[e(277)]) * s
                r < 120 && (r = 120)
                let a = s * n
                a < 20 && (a = 20)
                let h = document[e(1945)](e(1493))
                ;(h[e(2306)] = a + 8 + 10), (h[e(575)] = r + 8)
                let l = h[e(1019)]('2d'),
                  c = [h[e(575)] / 2, h[e(2306)] / 2],
                  u = l.createRadialGradient(
                    c[0],
                    c[1],
                    h[e(2306)] / 2,
                    c[0],
                    c[1],
                    0.5 * h[e(575)]
                  ),
                  m = this[e(1496)] ? this._style.selectedColor : this._label[e(821)]
                u[e(803)](0, m),
                  u.addColorStop(1, e(834)),
                  (l[e(537)] = u),
                  l[e(412)](0, 0, h.width, h[e(2306)]),
                  (l.font = e(886) + s + e(1495))
                let p = this[e(2119)][e(2262)]
                ;(l[e(537)] = p),
                  (l[e(1287)] = 'center'),
                  (l[e(1354)] = '#d3d3d3'),
                  l[e(2382)](0, 0, h.width, h[e(2306)]),
                  (l[e(1354)] = this[e(2119)][e(725)]),
                  (l.lineWidth = this[e(2119)][e(376)])
                for (let t = 1; t <= n; t++) {
                  let i = t * s + 4 * (t - 1) + 4
                  l.strokeText(o[t - 1], c[0], i), l[e(1104)](o[t - 1], c[0], i)
                }
                return h
              }
            })(e)
          case Q.textMarker:
            return new (class extends kt {
              constructor(e) {
                const i = t
                super(e), (this._typeName = '纯文本点'), (this[i(2336)] = Q[i(1770)])
                let s = Cesium[i(1960)](e[i(2240)], {})
                ;(this[i(2119)][i(2262)] = Cesium[i(1960)](s[i(2262)], i(287))),
                  (this[i(2119)][i(821)] = Cesium[i(1960)](s.backgroundColor, '#0008D9')),
                  (this[i(2119)].text = Cesium[i(1960)](s[i(1982)], i(625))),
                  (this[i(2119)][i(735)] = Cesium[i(1960)](s[i(735)], 20)),
                  (this[i(2119)].outlineWidth = Cesium[i(1960)](s.outlineWidth, 0)),
                  (this[i(556)][i(1323)] = !1),
                  this[i(1444)]()
              }
              [t(1444)]() {
                const e = t
                ;(this[e(987)][e(2240)].text = ''),
                  (this[e(987)][e(2240)].show = !1),
                  (this[e(987)][e(909)].image = this[e(2507)]()),
                  (this[e(987)].billboard[e(1482)] = !0)
              }
              _setStyle() {
                const e = t
                this[e(987)].billboard[e(259)] = this._createText()
              }
              [t(2507)]() {
                const e = t,
                  i = this._label[e(1982)]
                let s = this._label[e(735)],
                  n = Math[e(1215)](i[e(277)] / 15),
                  o = []
                if (n > 1)
                  for (let t = 1; t <= n; t++) {
                    const s = i.substring(15 * (t - 1), 15 * t)
                    o[e(2553)](s)
                  }
                else o[e(2553)](i)
                let r = Math.min(15, i.length) * s
                r < 120 && (r = 120)
                let a = s * n
                a < 20 && (a = 20)
                let h = document[e(1945)]('canvas')
                ;(h[e(2306)] = a + 8 + 10), (h.width = r + 8)
                let l = h.getContext('2d'),
                  c = [h[e(575)] / 2, h[e(2306)] / 2],
                  u = l[e(980)](c[0], c[1], h.height / 2, c[0], c[1], 0.5 * h[e(575)])
                u[e(803)](0, this[e(2119)][e(821)]),
                  u[e(803)](1, e(834)),
                  this[e(2119)][e(1875)] ? (l[e(537)] = u) : (l[e(537)] = 'transparent'),
                  l[e(412)](0, 0, h[e(575)], h.height),
                  (l.font = e(886) + s + e(1495)),
                  (l[e(537)] = this._isSelected ? this._style[e(506)] : this[e(2119)][e(2262)]),
                  (l.textAlign = e(365)),
                  this[e(2119)][e(1875)] &&
                    ((l[e(1354)] = '#d3d3d3'), l[e(2382)](0, 0, h.width, h[e(2306)])),
                  (l[e(1354)] = this._label[e(725)]),
                  (l.lineWidth = this[e(2119)].outlineWidth)
                for (let t = 1; t <= n; t++) {
                  let i = t * s + 4 * (t - 1) + 4
                  l.strokeText(o[t - 1], c[0], i), l[e(1104)](o[t - 1], c[0], i)
                }
                return h
              }
            })(e)
          case Q.bounceMarker:
            return new (class extends kt {
              constructor(e) {
                const i = t
                super(e),
                  (this[i(528)] = i(2059)),
                  (this[i(2336)] = Q[i(143)]),
                  (this._style[i(2183)] = Cesium.defaultValue(this[i(343)].bounceHeight, 100)),
                  (this[i(343)].increment = Cesium[i(1960)](this[i(343)][i(2227)], 0.1)),
                  this[i(2105)]()
              }
              [t(415)](e) {
                const i = t
                ;(this[i(987)][i(1482)] = e), e && this.bounce()
              }
              [t(2105)]() {
                const e = t
                let i = this[e(343)][e(2183)]
                if (!i) return
                const s = Cesium[e(2285)].fromCartesian(this[e(2238)])
                let n = s[e(2306)] + i,
                  o = 0,
                  r = 0,
                  a = this[e(343)].increment,
                  h = 0
                this[e(987)][e(2251)] = new Cesium[e(2569)]((t) => {
                  const l = e
                  ;(o += r += a) > i && ((o = i), (r *= -1), (r *= 0.55)),
                    o == i && h++,
                    h > 10 && (this[l(987)][l(2251)] = this[l(2238)])
                  const c = n - o
                  return Cesium.Cartesian3.fromRadians(s[l(2106)], s[l(199)], c)
                })
              }
            })(e)
          case Q[i(784)]:
            return new (class extends kt {
              constructor(e) {
                const i = t
                super(e),
                  (this[i(528)] = i(260)),
                  (this[i(2336)] = Q[i(784)]),
                  (this._style.floatHeight = Cesium.defaultValue(this._style[i(390)], 10)),
                  (this[i(343)][i(2227)] = Cesium[i(1960)](this[i(343)][i(2227)], 0.02)),
                  this[i(2231)]()
              }
              [t(200)]() {
                const e = t
                this[e(2520)] || this[e(1203)]()
              }
              [t(1203)]() {
                const e = t
                let i = Cesium.Cartographic.fromCartesian(this._cartesian3),
                  s = i[e(2306)],
                  n = s,
                  o = !0,
                  r = i[e(2306)] + this[e(343)][e(390)]
                this._entity[e(2251)] = new Cesium[e(2569)](
                  (t) => (
                    o
                      ? ((n += this[e(343)].increment), n > r && (o = !1))
                      : ((n -= this._style[e(2227)]), n < s && (o = !0)),
                    Cesium[e(310)][e(1069)](i[e(2106)], i[e(199)], n)
                  )
                )
              }
              [t(1444)]() {
                this._float()
              }
              [t(2231)]() {
                const e = t
                this[e(1203)]()
                let i = this[e(1496)] ? this[e(343)][e(506)] : this[e(2119)][e(2262)]
                this._entity.label.fillColor = Cesium[e(1154)][e(2008)](i)
              }
            })(e)
          case Q[i(1582)]:
            return new (class extends kt {
              constructor(e) {
                const i = t
                super(e),
                  (this._typeName = i(858)),
                  (this._graphicType = Q[i(1582)]),
                  this[i(2231)]()
              }
              [t(2231)]() {
                const e = t
                let i = this[e(556)][e(259)]
                if (!i || i[e(877)](e(285)) < 0)
                  throw new Cesium[e(2360)]('image是必须的！，且仅支持gif')
                let s = this[e(1496)] ? this[e(343)].selectedColor : this._label.fillColor
                if (
                  ((this[e(987)].label[e(2262)] = Cesium[e(1154)][e(2008)](s)),
                  this[e(2100)] == i && this[e(1033)])
                )
                  return void (this[e(987)][e(909)][e(259)] = new Cesium[e(2569)](
                    () => this[e(1033)].get_canvas()[e(1184)](e(1942)),
                    !1
                  ))
                this[e(2100)] = i
                let n = document[e(1945)](e(2575))
                ;(n.src = i),
                  (n[e(1701)] = () => {
                    this[e(2152)](n)
                  })
              }
              [t(2152)](e) {
                const i = t,
                  s = {}
                ;(s[i(285)] = e),
                  (this[i(1033)] = new SuperGif(s)),
                  this[i(1033)][i(2338)](() => {
                    const t = i
                    this[t(987)][t(909)][t(259)] = new Cesium[t(2569)](
                      () => this[t(1033)][t(1343)]()[t(1184)](t(1942)),
                      !1
                    )
                  })
              }
            })(e)
        }
      }
    }
    const Rt = {}
    ;(Rt.diffuse = t(2518)),
      (Rt[t(2315)] = t(1773)),
      (Cesium[t(1637)].CustomColorType = t(1537)),
      Cesium[t(1637)][t(2257)][t(1371)](Cesium[t(1637)][t(1537)], {
        fabric: {
          type: Cesium[t(1637)][t(1537)],
          uniforms: { color: new Cesium[t(1154)](1, 0, 0, 0.7) },
          components: Rt
        },
        translucent: function (t) {
          return !0
        }
      }),
      (Cesium[t(1637)][t(1716)] = t(1716)),
      Cesium[t(1637)]._materialCache[t(1371)](Cesium[t(1637)][t(1716)], {
        fabric: {
          type: Cesium[t(1637)].CircleColorfulType,
          uniforms: { color: new Cesium[t(1154)](1, 0, 0, 0.7), speed: 3 },
          source:
            '\nuniform vec4 color;\nuniform float speed; \nczm_material czm_getMaterial(czm_materialInput materialInput){\nczm_material material = czm_getDefaultMaterial(materialInput);\nvec2 st = materialInput.st  * 2.0 - 1.0;\nfloat time =czm_frameNumber * speed / 1000.0;\nfloat radius = length(st);\nfloat angle = atan(st.y/st.x);\nfloat radius1 = sin(time * 2.0) + sin(40.0*angle+time)*0.01;\nfloat radius2 = cos(time * 3.0);\nvec3 fragColor = 0.2 + 0.5 * cos( time + color.rgb );\nfloat inten1 = 1.0 - sqrt(abs(radius1 - radius));\nfloat inten2 = 1.0 - sqrt(abs(radius2 - radius));\n\nmaterial.alpha = pow(inten1 + inten2 , 5.0)*.5 ;\nmaterial.diffuse = fragColor * (inten1 + inten2); \nreturn material;\n}'
        },
        translucent: function (t) {
          return !0
        }
      }),
      (Cesium[t(1637)].CircleWaveType = t(477)),
      Cesium[t(1637)][t(2257)][t(1371)](Cesium.Material[t(477)], {
        fabric: {
          type: Cesium[t(1637)][t(477)],
          uniforms: { color: new Cesium.Color(1, 0, 0, 0.7), speed: 3, count: 1, gradient: 0.1 },
          source: t(193)
        },
        translucent: function (t) {
          return !0
        }
      }),
      (Cesium[t(1637)][t(747)] = 'CircleLightRingType'),
      Cesium[t(1637)][t(2257)].addMaterial(Cesium.Material[t(747)], {
        fabric: {
          type: Cesium[t(1637)][t(747)],
          uniforms: { color: new Cesium[t(1154)](1, 0, 0, 0.7), speed: 3 },
          source: t(1197)
        },
        translucent: function (t) {
          return !0
        }
      }),
      (Cesium[t(1637)][t(2075)] = t(2075)),
      Cesium[t(1637)][t(2257)].addMaterial(Cesium[t(1637)][t(2075)], {
        fabric: {
          type: Cesium[t(1637)][t(2075)],
          uniforms: { color: new Cesium[t(1154)](1, 0, 0, 0.7), speed: 3 },
          source: t(423)
        },
        translucent: function (t) {
          return !0
        }
      }),
      (Cesium[t(1637)].CircleImageDiffuseType = t(1852)),
      Cesium[t(1637)]._materialCache[t(1371)](Cesium.Material[t(1852)], {
        fabric: {
          type: Cesium.Material[t(1852)],
          uniforms: { color: new Cesium[t(1154)](1, 0, 0, 0.7), speed: 3, image: '' },
          source: t(253)
        },
        translucent: function (t) {
          return !0
        }
      }),
      (Cesium[t(1637)][t(1431)] = t(1431)),
      Cesium[t(1637)][t(2257)][t(1371)](Cesium.Material[t(1431)], {
        fabric: {
          type: Cesium[t(1637)].CircleImageRotateType,
          uniforms: { color: new Cesium.Color(1, 0, 0, 0.7), speed: 3, image: '' },
          source: t(2092)
        },
        translucent: function (t) {
          return !0
        }
      }),
      (Cesium.Material[t(166)] = 'CircleScanType_1'),
      Cesium[t(1637)]._materialCache[t(1371)](Cesium.Material[t(166)], {
        fabric: {
          type: Cesium[t(1637)].CircleScanType_1,
          uniforms: { color: new Cesium.Color(1, 0, 0, 0.7), speed: 10 },
          source:
            '\nuniform vec4 color;\n uniform float speed; \n #define PI 3.14159265359 ;\n czm_material czm_getMaterial(czm_materialInput materialInput){\n czm_material material = czm_getDefaultMaterial(materialInput);\n vec2 st = materialInput.st;\n vec2 scrPt = st * 2.0 - 1.0;\n float time = czm_frameNumber * speed / 1000.0 ;\n vec3 col = vec3(0.0);\n mat2 rot;\n float theta = -time * 1.0 * PI - 2.2;\n float cosTheta, sinTheta;\n cosTheta = cos(theta);\n sinTheta = sin(theta);\n rot[0][0] = cosTheta;\n rot[0][1] = -sinTheta;\n rot[1][0] = sinTheta;\n rot[1][1] = cosTheta;\n vec2 scrPtRot = rot * scrPt;\n float angle = 1.0 - (atan(scrPtRot.y, scrPtRot.x) / 6.2831 + 0.5);\n float falloff = length(scrPtRot);\n material.alpha = pow(length(col + vec3(.5)),5.0);\n material.diffuse =  (0.5 +  pow(angle, 2.0) * falloff ) *color.rgb;\n return material;\n } '
        },
        translucent: function (t) {
          return !0
        }
      }),
      (Cesium[t(1637)][t(1691)] = 'CircleScanType_2'),
      Cesium[t(1637)]._materialCache[t(1371)](Cesium.Material[t(1691)], {
        fabric: {
          type: Cesium.Material.CircleScanType_2,
          uniforms: { color: new Cesium[t(1154)](1, 0, 0, 0.7), speed: 10 },
          source:
            '\nuniform vec4 color;\nuniform float speed;\n\n#define PI 3.14159265359\n\nfloat rand(vec2 co){\nreturn fract(sin(dot(co.xy ,vec2(12.9898,78.233))) * 43758.5453);\n}\n\nczm_material czm_getMaterial(czm_materialInput materialInput){\nczm_material material = czm_getDefaultMaterial(materialInput);\nvec2 st = materialInput.st;\nvec2 pos = st - vec2(0.5);\nfloat time = czm_frameNumber * speed / 1000.0 ;\nfloat r = length(pos);\nfloat t = atan(pos.y, pos.x) - time * 2.5;\nfloat a = (atan(sin(t), cos(t)) + PI)/(2.0*PI);\nfloat ta = 0.5;\nfloat v = smoothstep(ta-0.05,ta+0.05,a) * smoothstep(ta+0.05,ta-0.05,a);\nvec3 flagColor = color.rgb * v;\nfloat blink = pow(sin(time*1.5)*0.5+0.5, 0.8);\nflagColor = color.rgb *  pow(a, 8.0*(.2+blink))*(sin(r*500.0)*.5+.5) ;\nflagColor = flagColor * pow(r, 0.4);\nmaterial.alpha = length(flagColor) * 1.3;\nmaterial.diffuse = flagColor * 3.0;\nreturn material;\n} '
        },
        translucent: function (t) {
          return !0
        }
      }),
      (Cesium[t(1637)][t(994)] = t(994)),
      Cesium[t(1637)][t(2257)][t(1371)](Cesium[t(1637)][t(994)], {
        fabric: {
          type: Cesium.Material[t(994)],
          uniforms: { color: new Cesium[t(1154)](1, 0, 0, 0.7), speed: 10 },
          source: t(1479)
        },
        translucent: function (t) {
          return !0
        }
      }),
      (Cesium.Material[t(733)] = t(733)),
      Cesium[t(1637)][t(2257)].addMaterial(Cesium[t(1637)][t(733)], {
        fabric: {
          type: Cesium.Material[t(733)],
          uniforms: { u_color: new Cesium[t(1154)](1, 0, 0, 0.7), u_speed: 10 },
          source: t(652)
        },
        translucent: function (t) {
          return !0
        }
      }),
      (Cesium[t(1637)][t(1066)] = t(323)),
      Cesium[t(1637)]._materialCache.addMaterial(Cesium[t(1637)][t(1066)], {
        fabric: {
          type: Cesium[t(1637)].PolylineTrailType,
          uniforms: {
            color: new Cesium.Color(1, 0, 0, 0.5),
            speed: 0,
            gradient: 0.01,
            percent: 0.1
          },
          source: t(368)
        },
        translucent: function (t) {
          return !0
        }
      }),
      (Cesium[t(1637)][t(799)] = t(2388)),
      Cesium[t(1637)][t(2257)][t(1371)](Cesium.Material.PolylineFlowType, {
        fabric: {
          type: Cesium.Material[t(799)],
          uniforms: {
            color: new Cesium[t(1154)](1, 0, 0, 0.5),
            speed: 0,
            image: '',
            repeat: 3,
            sampling: 1
          },
          source: t(1963)
        },
        translucent: function (t) {
          return !0
        }
      })
    const Lt = {}
    ;(Lt[t(259)] = ''),
      (Cesium[t(1637)][t(1535)] = t(1535)),
      Cesium[t(1637)][t(2257)][t(1371)](Cesium[t(1637)][t(1535)], {
        fabric: {
          type: Cesium[t(1637)].WallFlowType,
          uniforms: {
            color: new Cesium.Color(1, 0, 0, 0.5),
            image0:
              'data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAADgAAAA4CAYAAACohjseAAAACXBIWXMAAAsTAAALEwEAmpwYAAAKTWlDQ1BQaG90b3Nob3AgSUNDIHByb2ZpbGUAAHjanVN3WJP3Fj7f92UPVkLY8LGXbIEAIiOsCMgQWaIQkgBhhBASQMWFiApWFBURnEhVxILVCkidiOKgKLhnQYqIWotVXDjuH9yntX167+3t+9f7vOec5/zOec8PgBESJpHmomoAOVKFPDrYH49PSMTJvYACFUjgBCAQ5svCZwXFAADwA3l4fnSwP/wBr28AAgBw1S4kEsfh/4O6UCZXACCRAOAiEucLAZBSAMguVMgUAMgYALBTs2QKAJQAAGx5fEIiAKoNAOz0ST4FANipk9wXANiiHKkIAI0BAJkoRyQCQLsAYFWBUiwCwMIAoKxAIi4EwK4BgFm2MkcCgL0FAHaOWJAPQGAAgJlCLMwAIDgCAEMeE80DIEwDoDDSv+CpX3CFuEgBAMDLlc2XS9IzFLiV0Bp38vDg4iHiwmyxQmEXKRBmCeQinJebIxNI5wNMzgwAABr50cH+OD+Q5+bk4eZm52zv9MWi/mvwbyI+IfHf/ryMAgQAEE7P79pf5eXWA3DHAbB1v2upWwDaVgBo3/ldM9sJoFoK0Hr5i3k4/EAenqFQyDwdHAoLC+0lYqG9MOOLPv8z4W/gi372/EAe/tt68ABxmkCZrcCjg/1xYW52rlKO58sEQjFu9+cj/seFf/2OKdHiNLFcLBWK8ViJuFAiTcd5uVKRRCHJleIS6X8y8R+W/QmTdw0ArIZPwE62B7XLbMB+7gECiw5Y0nYAQH7zLYwaC5EAEGc0Mnn3AACTv/mPQCsBAM2XpOMAALzoGFyolBdMxggAAESggSqwQQcMwRSswA6cwR28wBcCYQZEQAwkwDwQQgbkgBwKoRiWQRlUwDrYBLWwAxqgEZrhELTBMTgN5+ASXIHrcBcGYBiewhi8hgkEQcgIE2EhOogRYo7YIs4IF5mOBCJhSDSSgKQg6YgUUSLFyHKkAqlCapFdSCPyLXIUOY1cQPqQ28ggMor8irxHMZSBslED1AJ1QLmoHxqKxqBz0XQ0D12AlqJr0Rq0Hj2AtqKn0UvodXQAfYqOY4DRMQ5mjNlhXIyHRWCJWBomxxZj5Vg1Vo81Yx1YN3YVG8CeYe8IJAKLgBPsCF6EEMJsgpCQR1hMWEOoJewjtBK6CFcJg4Qxwicik6hPtCV6EvnEeGI6sZBYRqwm7iEeIZ4lXicOE1+TSCQOyZLkTgohJZAySQtJa0jbSC2kU6Q+0hBpnEwm65Btyd7kCLKArCCXkbeQD5BPkvvJw+S3FDrFiOJMCaIkUqSUEko1ZT/lBKWfMkKZoKpRzame1AiqiDqfWkltoHZQL1OHqRM0dZolzZsWQ8ukLaPV0JppZ2n3aC/pdLoJ3YMeRZfQl9Jr6Afp5+mD9HcMDYYNg8dIYigZaxl7GacYtxkvmUymBdOXmchUMNcyG5lnmA+Yb1VYKvYqfBWRyhKVOpVWlX6V56pUVXNVP9V5qgtUq1UPq15WfaZGVbNQ46kJ1Bar1akdVbupNq7OUndSj1DPUV+jvl/9gvpjDbKGhUaghkijVGO3xhmNIRbGMmXxWELWclYD6yxrmE1iW7L57Ex2Bfsbdi97TFNDc6pmrGaRZp3mcc0BDsax4PA52ZxKziHODc57LQMtPy2x1mqtZq1+rTfaetq+2mLtcu0W7eva73VwnUCdLJ31Om0693UJuja6UbqFutt1z+o+02PreekJ9cr1Dund0Uf1bfSj9Rfq79bv0R83MDQINpAZbDE4Y/DMkGPoa5hpuNHwhOGoEctoupHEaKPRSaMnuCbuh2fjNXgXPmasbxxirDTeZdxrPGFiaTLbpMSkxeS+Kc2Ua5pmutG003TMzMgs3KzYrMnsjjnVnGueYb7ZvNv8jYWlRZzFSos2i8eW2pZ8ywWWTZb3rJhWPlZ5VvVW16xJ1lzrLOtt1ldsUBtXmwybOpvLtqitm63Edptt3xTiFI8p0in1U27aMez87ArsmuwG7Tn2YfYl9m32zx3MHBId1jt0O3xydHXMdmxwvOuk4TTDqcSpw+lXZxtnoXOd8zUXpkuQyxKXdpcXU22niqdun3rLleUa7rrStdP1o5u7m9yt2W3U3cw9xX2r+00umxvJXcM970H08PdY4nHM452nm6fC85DnL152Xlle+70eT7OcJp7WMG3I28Rb4L3Le2A6Pj1l+s7pAz7GPgKfep+Hvqa+It89viN+1n6Zfgf8nvs7+sv9j/i/4XnyFvFOBWABwQHlAb2BGoGzA2sDHwSZBKUHNQWNBbsGLww+FUIMCQ1ZH3KTb8AX8hv5YzPcZyya0RXKCJ0VWhv6MMwmTB7WEY6GzwjfEH5vpvlM6cy2CIjgR2yIuB9pGZkX+X0UKSoyqi7qUbRTdHF09yzWrORZ+2e9jvGPqYy5O9tqtnJ2Z6xqbFJsY+ybuIC4qriBeIf4RfGXEnQTJAntieTE2MQ9ieNzAudsmjOc5JpUlnRjruXcorkX5unOy553PFk1WZB8OIWYEpeyP+WDIEJQLxhP5aduTR0T8oSbhU9FvqKNolGxt7hKPJLmnVaV9jjdO31D+miGT0Z1xjMJT1IreZEZkrkj801WRNberM/ZcdktOZSclJyjUg1plrQr1zC3KLdPZisrkw3keeZtyhuTh8r35CP5c/PbFWyFTNGjtFKuUA4WTC+oK3hbGFt4uEi9SFrUM99m/ur5IwuCFny9kLBQuLCz2Lh4WfHgIr9FuxYji1MXdy4xXVK6ZHhp8NJ9y2jLspb9UOJYUlXyannc8o5Sg9KlpUMrglc0lamUycturvRauWMVYZVkVe9ql9VbVn8qF5VfrHCsqK74sEa45uJXTl/VfPV5bdra3kq3yu3rSOuk626s91m/r0q9akHV0IbwDa0b8Y3lG19tSt50oXpq9Y7NtM3KzQM1YTXtW8y2rNvyoTaj9nqdf13LVv2tq7e+2Sba1r/dd3vzDoMdFTve75TsvLUreFdrvUV99W7S7oLdjxpiG7q/5n7duEd3T8Wej3ulewf2Re/ranRvbNyvv7+yCW1SNo0eSDpw5ZuAb9qb7Zp3tXBaKg7CQeXBJ9+mfHvjUOihzsPcw83fmX+39QjrSHkr0jq/dawto22gPaG97+iMo50dXh1Hvrf/fu8x42N1xzWPV56gnSg98fnkgpPjp2Snnp1OPz3Umdx590z8mWtdUV29Z0PPnj8XdO5Mt1/3yfPe549d8Lxw9CL3Ytslt0utPa49R35w/eFIr1tv62X3y+1XPK509E3rO9Hv03/6asDVc9f41y5dn3m978bsG7duJt0cuCW69fh29u0XdwruTNxdeo94r/y+2v3qB/oP6n+0/rFlwG3g+GDAYM/DWQ/vDgmHnv6U/9OH4dJHzEfVI0YjjY+dHx8bDRq98mTOk+GnsqcTz8p+Vv9563Or59/94vtLz1j82PAL+YvPv655qfNy76uprzrHI8cfvM55PfGm/K3O233vuO+638e9H5ko/ED+UPPR+mPHp9BP9z7nfP78L/eE8/sl0p8zAAAAIGNIUk0AAHolAACAgwAA+f8AAIDpAAB1MAAA6mAAADqYAAAXb5JfxUYAAAH0SURBVHja7JnNcoMgEIA/IlFMz23f/9EyvfUBOm2TSC9syzjRSIVEnWXGHDKy7Mf+gsZ7z5bHjo0PBVRABVRABVRABVRABVRABVRABVRABVRABVTA7MMej8dqZTon3XNaYA+Y1IlLARgZBvA2QO4WBJhrM34Ba6ACuokCzUwlTGZLjQF2qYBrctEdcBFAOwI4ZIFHKT5VdgWcLdCERHPuKe9HTP9IQD/R3S3wHQPagclTFTWFoMY8aEyXCtgJoANO/1SoRNLIIcsKoAtPtaCYyqHHHjBzAJcObuXHAQfgawOloW/BX8B2oPH2iXFYKrFMSWr9OY10MgJoZpQEX8A6ZqbMRjqZNrioSdid3C6WCjNlvVYAHfBUsI49Ki5b6WQO4ekLzt3J3ILxA8c2MwPwJC7aRr2oLwhyy1I+47oHadWkTHSZF0hxvxLrHoDPuNCf75gkcp/cr8lugTa24GWG3/tM75qMyagFXFwmLgvInj7DfBO56EcMeNpAe9aPQRd3Mt8rb679kIvKefCr0ELmARvjY0Cx4Ocda19pwC4wNXGZaK6cKJZ4V+onvuPEgnUEmXLDbe5oqRQZJgKs4xh0K82WNy0ogHViDC75qr8LXI18fKn5+wizplgbm9sAewu8hj+f2d54tcAb8BLKxJaGA95/BgD0bHrmTa2aKAAAAABJRU5ErkJggg==',
            image1: t(530),
            time: 10,
            speed: 30
          },
          source: t(931)
        },
        translucent: function (t) {
          return !0
        }
      }),
      (Cesium[t(1637)][t(2558)] = t(2558)),
      Cesium.Material[t(2257)][t(1371)](Cesium[t(1637)][t(2558)], {
        fabric: { type: Cesium[t(1637)].WallGradientType, uniforms: Lt, source: t(170) },
        translucent: function (t) {
          return !0
        }
      }),
      (Cesium[t(1637)][t(1890)] = t(1890)),
      Cesium.Material[t(2257)][t(1371)](Cesium.Material.WallImageFlowType, {
        fabric: {
          type: Cesium[t(1637)][t(1890)],
          uniforms: {
            color: new Cesium[t(1154)](1, 0, 0, 0.5),
            image: '',
            time: 0,
            speed: 0,
            repeat: 2,
            direction: 1
          },
          source:
            ' \nczm_material czm_getMaterial(czm_materialInput materialInput)\n{\n    czm_material material = czm_getDefaultMaterial(materialInput);\n    vec2 st = materialInput.st; \n\n    float time =czm_frameNumber * speed / 1000.0 ;\n    vec4 colorImage = texture(image, vec2(fract(repeat * st.s - time), st.t)); \n    if(direction==2.){ \n        colorImage = texture(image, vec2(fract(repeat * st.t - time), st.s)); \n    }\n    material.alpha = colorImage.a * color.a;\n    material.diffuse =  2.5 * color.rgb  ;\n    return material;\n}'
        },
        translucent: function (t) {
          return !0
        }
      }),
      (Cesium[t(1637)][t(1220)] = t(1220)),
      Cesium[t(1637)][t(2257)].addMaterial(Cesium[t(1637)][t(1220)], {
        fabric: {
          type: Cesium.Material[t(1220)],
          uniforms: { color: new Cesium[t(1154)](1, 0, 0, 0.5), breathe: !1 },
          source: t(1477)
        },
        translucent: function (t) {
          return !0
        }
      }),
      (Cesium[t(1637)][t(480)] = t(480)),
      Cesium[t(1637)]._materialCache.addMaterial(Cesium.Material.WallScrollType, {
        fabric: {
          type: Cesium[t(1637)][t(480)],
          uniforms: { color: Cesium[t(1154)][t(1661)], speed: 1, repeat: 5 },
          source: t(1787)
        },
        translucent: function (t) {
          return !0
        }
      }),
      (Cesium[t(1637)][t(2326)] = 'WallCoolType'),
      Cesium[t(1637)]._materialCache[t(1371)](t(2326), {
        fabric: {
          type: t(2326),
          uniforms: {
            color: new Cesium[t(1154)](1, 1, 1, 0.8),
            image: '',
            image2: '',
            image3: '',
            speed: 1,
            repeat: 3,
            glow: 1
          },
          source:
            '\nczm_material czm_getMaterial(czm_materialInput materialInput)\n{\n    czm_material material = czm_getDefaultMaterial(materialInput);\n    vec2 vUv = materialInput.st;\n    float time2 = czm_frameNumber / 100.0;\n    time2 = speed * time2;\n    vec4 colorb =texture(image2,vec2(vUv.s,vUv.t));\n    vec4 colora = texture(image,vec2(fract(vUv.s * repeat ),fract(vUv.t-time2)));\n    vec4 colorc= texture(image3,vec2(vUv.s,vUv.t));\n    material.alpha = colorb.a + colora.a + colorc.a;\n    material.diffuse = colorb.rgb + colora.rgb * colorc.rgb * color.rgb * glow;\n    return material;\n}'
        },
        translucent: function (t) {
          return !0
        }
      }),
      (Cesium.Material.PolyGradientType = t(382)),
      Cesium[t(1637)][t(2257)][t(1371)](Cesium[t(1637)][t(382)], {
        fabric: {
          type: Cesium.Material[t(382)],
          uniforms: {
            color: new Cesium[t(1154)](1, 0, 0, 0.5),
            center: new Cesium[t(194)](0.5, 0.5)
          },
          source: t(1757)
        },
        translucent: function (t) {
          return !0
        }
      }),
      (Cesium[t(1637)][t(734)] = t(734)),
      Cesium[t(1637)][t(2257)].addMaterial(Cesium[t(1637)][t(734)], {
        fabric: {
          type: Cesium[t(1637)][t(734)],
          uniforms: { color: new Cesium[t(1154)](1, 0, 0, 0.5), spacing: 200, width: 1 },
          source: t(2157)
        },
        translucent: function (t) {
          return !0
        }
      }),
      (Cesium.Material.SphereElectricType = t(1966)),
      Cesium[t(1637)][t(2257)][t(1371)](Cesium[t(1637)][t(1966)], {
        fabric: {
          type: Cesium[t(1637)][t(1966)],
          uniforms: { color: new Cesium[t(1154)](1, 0, 0, 0.7), speed: 10 },
          source: t(2189)
        },
        translucent: function (t) {
          return !0
        }
      }),
      (Cesium[t(1637)][t(2099)] = 'SphereScanType'),
      Cesium[t(1637)][t(2257)][t(1371)](Cesium.Material[t(2099)], {
        fabric: {
          type: Cesium[t(1637)][t(2099)],
          uniforms: { color: new Cesium.Color(1, 0, 0, 0.7), speed: 10 },
          source:
            '\nuniform vec4 color; \nczm_material czm_getMaterial(czm_materialInput materialInput)\n{\n    czm_material material = czm_getDefaultMaterial(materialInput);\n    vec2 st = materialInput.st;  \n    float t =  fract( czm_frameNumber * speed / 1000.0);; \n    t *= 1.03;\n    float alpha = smoothstep(t- 0.03, t, st.s) * step(-t, -st.s); \n    alpha += 0.1;\n    alpha *= step(-0.5, -abs(0.5-st.t));                             \n    material.diffuse = color.rgb;\n    material.alpha = alpha;\n    return material;\n}\n'
        },
        translucent: function (t) {
          return !0
        }
      })
    class Ot {
      constructor(e = {}) {
        const i = t
        ;(this[i(1035)] = e), (this[i(2273)] = new Cesium[i(525)]()), (this[i(1165)] = i(1577))
      }
      get [t(1713)]() {
        return this._materialOpts
      }
      get [t(1916)]() {
        return this[t(1165)]
      }
      set [t(1916)](t) {
        this._name = t
      }
      get [t(989)]() {
        return !1
      }
      get definitionChanged() {
        return this[t(2273)]
      }
      [t(2578)](t) {
        return null
      }
      [t(1474)](t, e) {
        return Cesium.defaultValue(e, {})
      }
      [t(1525)](t) {
        return this === t
      }
      [t(1180)]() {
        const e = t
        for (const t in this[e(1035)])
          if (Object[e(782)][e(1669)](this[e(1035)], t)) {
            const i = this[e(1035)][t]
            this[t] = i
          }
        ;(this[e(1035)][e(1916)] = this._name), this[e(1078)]()
      }
      [t(1078)]() {}
      [t(589)](t) {}
    }
    class Bt extends Ot {
      constructor(e = {}) {
        const i = t
        super(e), (this[i(1165)] = '多彩圆材质'), this[i(589)](e)
      }
      [t(589)](e) {
        const i = t
        ;(e[i(1070)] = Cesium[i(1960)](e[i(1070)], i(456))),
          (e.speed = Cesium[i(1960)](e[i(2535)], 3)),
          (this[i(1035)] = e),
          this.mergeOpts()
      }
      [t(1078)]() {
        const e = t
        this[e(1070)] = Cesium[e(1154)][e(2008)](this[e(1035)].color)
      }
      [t(2578)](e) {
        const i = t
        return Cesium[i(1637)][i(1716)]
      }
      [t(1474)](e, i) {
        const s = t
        return (
          ((i = Cesium.defaultValue(i, {}))[s(1070)] = Cesium[s(1169)][s(2133)](this[s(1070)], e)),
          (i[s(2535)] = this[s(2535)]),
          i
        )
      }
      [t(1525)](e) {
        const i = t
        return (
          this === e ||
          (e instanceof Bt &&
            Cesium[i(1169)][i(1525)](this[i(1070)], e[i(1070)]) &&
            Cesium[i(1169)][i(1525)](this[i(2535)], e[i(2535)]))
        )
      }
    }
    Object[t(610)](Bt[t(1727)], {
      color: Cesium[t(1328)](t(1070)),
      speed: Cesium[t(1328)]('speed')
    })
    class Vt extends Ot {
      constructor(e = {}) {
        const i = t
        super(e), (this[i(1165)] = i(513)), this[i(589)](e)
      }
      [t(589)](e) {
        const i = t
        ;(e.color = Cesium[i(1960)](e[i(1070)], i(456))),
          (e[i(2535)] = Cesium[i(1960)](e.speed, 1)),
          (e[i(1250)] = Cesium.defaultValue(e[i(1250)], 1)),
          (e[i(159)] = Cesium.defaultValue(e[i(159)], 0.1)),
          (this._materialOpts = e),
          this[i(1180)]()
      }
      getType(e) {
        const i = t
        return Cesium.Material[i(477)]
      }
      _mergeOpts() {
        const e = t
        this[e(1070)] = Cesium[e(1154)].fromCssColorString(this[e(1035)].color)
      }
      getValue(e, i) {
        const s = t
        return (
          ((i = Cesium[s(1960)](i, {}))[s(1070)] = Cesium[s(1169)][s(2133)](this.color, e)),
          (i[s(2535)] = this[s(2535)]),
          (i[s(1250)] = this[s(1250)]),
          (i[s(159)] = 1 + 10 * (1 - this[s(159)])),
          i
        )
      }
      [t(1525)](e) {
        const i = t
        return (
          this === e ||
          (e instanceof Vt &&
            Cesium.Property.equals(this.color, e.color) &&
            this.speed == e[i(2535)] &&
            this.count == e[i(1250)] &&
            this[i(159)] == e.gradient)
        )
      }
    }
    Object[t(610)](Vt[t(1727)], { color: Cesium.createPropertyDescriptor(t(1070)) })
    class Nt extends Ot {
      constructor(e = {}) {
        const i = t
        super(e), (this[i(1165)] = '光圈材质'), this[i(589)](e)
      }
      [t(589)](e) {
        const i = t
        ;(e[i(1070)] = Cesium[i(1960)](e[i(1070)], i(456))),
          (e.speed = Cesium[i(1960)](e[i(2535)], 1)),
          (this[i(1035)] = e),
          this[i(1180)]()
      }
      [t(2578)](e) {
        return Cesium[t(1637)].CircleLightRingType
      }
      [t(1078)]() {
        const e = t
        this.color = Cesium.Color.fromCssColorString(this[e(1035)][e(1070)])
      }
      [t(1474)](e, i) {
        const s = t
        return (
          ((i = Cesium.defaultValue(i, {}))[s(1070)] = Cesium[s(1169)][s(2133)](this[s(1070)], e)),
          (i[s(2535)] = this[s(2535)]),
          i
        )
      }
      equals(e) {
        const i = t
        return (
          this === e ||
          (e instanceof Nt &&
            Cesium[i(1169)].equals(this[i(1070)], e[i(1070)]) &&
            this[i(2535)] == e[i(2535)])
        )
      }
    }
    Object[t(610)](Nt[t(1727)], { color: Cesium.createPropertyDescriptor(t(1070)) })
    class Ht extends Ot {
      constructor(e = {}) {
        const i = t
        super(e), (this[i(1165)] = i(1316)), this[i(589)](e)
      }
      setOpts(e) {
        const i = t
        ;(e.color = Cesium[i(1960)](e[i(1070)], i(456))),
          (e[i(2535)] = Cesium[i(1960)](e[i(2535)], 1)),
          (this._materialOpts = e),
          this[i(1180)]()
      }
      [t(2578)](e) {
        const i = t
        return Cesium[i(1637)][i(2075)]
      }
      [t(1078)]() {
        const e = t
        this[e(1070)] = Cesium[e(1154)][e(2008)](this[e(1035)][e(1070)])
      }
      [t(1474)](e, i) {
        const s = t
        return (
          ((i = Cesium[s(1960)](i, {}))[s(1070)] = Cesium[s(1169)].getValueOrUndefined(
            this[s(1070)],
            e
          )),
          (i[s(2535)] = this.speed),
          i
        )
      }
      equals(e) {
        const i = t
        return (
          this === e ||
          (e instanceof Ht &&
            Cesium[i(1169)][i(1525)](this.color, e[i(1070)]) &&
            this[i(2535)] == e[i(2535)])
        )
      }
    }
    Object[t(610)](Ht.prototype, { color: Cesium.createPropertyDescriptor(t(1070)) })
    class Gt extends Ot {
      constructor(t = {}) {
        super(t), (this._name = '图片扩散'), this.setOpts(t)
      }
      [t(589)](e) {
        const i = t
        ;(e[i(1070)] = Cesium.defaultValue(e[i(1070)], i(456))),
          (e.speed = Cesium.defaultValue(e[i(2535)], 3)),
          (this[i(1035)] = e),
          this.mergeOpts()
      }
      [t(2578)](e) {
        const i = t
        return Cesium[i(1637)][i(1852)]
      }
      [t(1078)]() {
        const e = t
        this[e(1070)] = Cesium.Color.fromCssColorString(this[e(1035)][e(1070)])
      }
      [t(1474)](e, i) {
        const s = t
        return (
          ((i = Cesium[s(1960)](i, {})).color = Cesium[s(1169)][s(2133)](this[s(1070)], e)),
          (i.speed = this.speed),
          (i[s(259)] = this[s(259)]),
          i
        )
      }
      [t(1525)](e) {
        const i = t
        return (
          this === e ||
          (e instanceof Gt &&
            Cesium[i(1169)][i(1525)](this[i(1070)], e[i(1070)]) &&
            this.speed == e[i(2535)] &&
            this[i(259)] == e[i(259)])
        )
      }
    }
    Object[t(610)](Gt[t(1727)], { color: Cesium[t(1328)](t(1070)) })
    class Wt extends Ot {
      constructor(e = {}) {
        const i = t
        super(e), (this[i(1165)] = i(678)), this[i(589)](e)
      }
      [t(589)](e) {
        const i = t
        ;(e[i(1070)] = Cesium[i(1960)](e.color, 'rgba(255, 0, 0, 0.5)')),
          (e[i(2535)] = Cesium[i(1960)](e[i(2535)], 3)),
          (this[i(1035)] = e),
          this[i(1180)]()
      }
      [t(2578)](e) {
        const i = t
        return Cesium[i(1637)][i(1431)]
      }
      [t(1078)]() {
        const e = t
        this[e(1070)] = Cesium[e(1154)][e(2008)](this[e(1035)].color)
      }
      [t(1474)](e, i) {
        const s = t
        return (
          ((i = Cesium.defaultValue(i, {})).color = Cesium.Property[s(2133)](this[s(1070)], e)),
          (i[s(2535)] = this[s(2535)]),
          (i[s(259)] = this[s(259)]),
          i
        )
      }
      [t(1525)](e) {
        const i = t
        return (
          this === e ||
          (e instanceof Wt &&
            Cesium[i(1169)][i(1525)](this[i(1070)], e[i(1070)]) &&
            this.speed == e[i(2535)] &&
            this[i(259)] == e[i(259)])
        )
      }
    }
    Object[t(610)](Wt.prototype, { color: Cesium[t(1328)](t(1070)) })
    class Ut extends Ot {
      constructor(e = {}) {
        const i = t
        super(e), (this[i(1165)] = i(2415)), this[i(589)](e)
      }
      setOpts(e) {
        const i = t
        ;(e[i(1070)] = Cesium[i(1960)](e.color, i(456))),
          (e[i(2535)] = Cesium[i(1960)](e[i(2535)], 10)),
          (this._materialOpts = e),
          this[i(1180)]()
      }
      [t(2578)](e) {
        const i = t
        return Cesium.Material[i(166)]
      }
      [t(1078)]() {
        const e = t
        this[e(1070)] = Cesium[e(1154)][e(2008)](this._materialOpts[e(1070)])
      }
      [t(1474)](e, i) {
        const s = t
        return (
          ((i = Cesium[s(1960)](i, {})).color = Cesium[s(1169)][s(2133)](this.color, e)),
          (i[s(2535)] = this.speed),
          i
        )
      }
      [t(1525)](e) {
        const i = t
        return (
          this === e ||
          (e instanceof Ut &&
            Cesium[i(1169)][i(1525)](this[i(1070)], e[i(1070)]) &&
            this[i(2535)] == e[i(2535)])
        )
      }
    }
    Object[t(610)](Ut[t(1727)], { color: Cesium[t(1328)](t(1070)) })
    class jt extends Ot {
      constructor(e = {}) {
        const i = t
        super(e), (this[i(1165)] = i(580)), this[i(589)](e)
      }
      setOpts(e) {
        const i = t
        ;(e[i(1070)] = Cesium[i(1960)](e[i(1070)], i(456))),
          (e[i(2535)] = Cesium[i(1960)](e[i(2535)], 10)),
          (this[i(1035)] = e),
          this.mergeOpts()
      }
      [t(2578)](e) {
        const i = t
        return Cesium[i(1637)][i(1691)]
      }
      _mergeOpts() {
        const e = t
        this[e(1070)] = Cesium.Color[e(2008)](this[e(1035)][e(1070)])
      }
      [t(1474)](e, i) {
        const s = t
        return (
          ((i = Cesium[s(1960)](i, {}))[s(1070)] = Cesium.Property[s(2133)](this[s(1070)], e)),
          (i[s(2535)] = this.speed),
          i
        )
      }
      equals(e) {
        const i = t
        return (
          this === e ||
          (e instanceof jt &&
            Cesium[i(1169)][i(1525)](this[i(1070)], e.color) &&
            this[i(2535)] == e[i(2535)])
        )
      }
    }
    Object[t(610)](jt[t(1727)], { color: Cesium[t(1328)]('color') })
    class qt extends Ot {
      constructor(e = {}) {
        const i = t
        super(e), (this[i(1165)] = '扫描材质_3'), this[i(589)](e)
      }
      [t(589)](e) {
        const i = t
        ;(e[i(1070)] = Cesium.defaultValue(e[i(1070)], i(456))),
          (e[i(2535)] = Cesium[i(1960)](e[i(2535)], 10)),
          (this[i(1035)] = e),
          this.mergeOpts()
      }
      [t(2578)](e) {
        const i = t
        return Cesium[i(1637)][i(994)]
      }
      [t(1078)]() {
        const e = t
        this[e(1070)] = Cesium.Color.fromCssColorString(this[e(1035)][e(1070)])
      }
      getValue(e, i) {
        const s = t
        return (
          ((i = Cesium[s(1960)](i, {})).color = Cesium[s(1169)][s(2133)](this.color, e)),
          (i[s(2535)] = this[s(2535)]),
          i
        )
      }
      [t(1525)](e) {
        const i = t
        return (
          this === e ||
          (e instanceof qt &&
            Cesium[i(1169)][i(1525)](this[i(1070)], e[i(1070)]) &&
            this[i(2535)] == e[i(2535)])
        )
      }
    }
    Object[t(610)](qt[t(1727)], { color: Cesium.createPropertyDescriptor('color') })
    class Yt extends Ot {
      constructor(e = {}) {
        const i = t
        super(e), (this[i(1165)] = i(1578)), this[i(589)](e)
      }
      [t(589)](e) {
        const i = t
        ;(e[i(1070)] = Cesium[i(1960)](e[i(1070)], i(456))),
          (e[i(2535)] = Cesium[i(1960)](e[i(2535)], 10)),
          (this._materialOpts = e),
          this[i(1180)]()
      }
      [t(2578)](e) {
        const i = t
        return Cesium[i(1637)][i(733)]
      }
      [t(1078)]() {
        const e = t
        this[e(1070)] = Cesium[e(1154)][e(2008)](this[e(1035)][e(1070)])
      }
      [t(1474)](e, i) {
        const s = t
        return (
          ((i = Cesium[s(1960)](i, {}))[s(202)] = Cesium.Property[s(2133)](this[s(1070)], e)),
          (i[s(228)] = this.speed),
          i
        )
      }
      [t(1525)](e) {
        const i = t
        return (
          this === e ||
          (e instanceof Yt &&
            Cesium[i(1169)][i(1525)](this[i(1070)], e[i(1070)]) &&
            this[i(2535)] == e[i(2535)])
        )
      }
    }
    Object[t(610)](Yt[t(1727)], { color: Cesium.createPropertyDescriptor(t(1070)) })
    class Xt extends Ot {
      constructor(e = {}) {
        const i = t
        super(e), (this._name = i(1076)), this[i(589)](e)
      }
      [t(589)](e) {
        const i = t
        ;(e.color = Cesium[i(1960)](e[i(1070)], i(456))),
          (e[i(2535)] = Cesium[i(1960)](e[i(2535)], 1)),
          (e[i(159)] = Cesium.defaultValue(e[i(159)], 0.01)),
          (e[i(1166)] = Cesium.defaultValue(e[i(1166)], 0.1)),
          (this[i(1035)] = e),
          this.mergeOpts()
      }
      getType(e) {
        const i = t
        return Cesium[i(1637)][i(1066)]
      }
      [t(1078)]() {
        const e = t
        this[e(1070)] = Cesium[e(1154)][e(2008)](this[e(1035)][e(1070)])
      }
      [t(1474)](e, i) {
        const s = t
        return (
          ((i = Cesium[s(1960)](i, {})).color = Cesium[s(1169)][s(2133)](this[s(1070)], e)),
          (i[s(2535)] = this[s(2535)]),
          (i[s(159)] = this[s(159)]),
          (i[s(1166)] = this[s(1166)]),
          i
        )
      }
      [t(1525)](e) {
        const i = t
        return (
          this === e ||
          (e instanceof Xt &&
            Cesium.Property.equals(this[i(1070)], e[i(1070)]) &&
            this[i(2535)] == e[i(2535)] &&
            this.gradient == e[i(159)] &&
            this[i(1166)] == e[i(1166)])
        )
      }
    }
    Object[t(610)](Xt[t(1727)], { color: Cesium[t(1328)]('color') })
    class Qt extends Ot {
      constructor(e = {}) {
        const i = t
        super(e), (this[i(1165)] = '流动线'), this[i(589)](e)
      }
      [t(589)](e) {
        const i = t
        ;(e.color = Cesium[i(1960)](e.color, 'rgba(255, 0, 0, 0.5)')),
          (e.speed = Cesium[i(1960)](e.speed, 1)),
          (e.image = Cesium[i(1960)](e.image, '')),
          (e[i(638)] = Cesium[i(1960)](e[i(638)], 3)),
          (e[i(2190)] = Cesium.defaultValue(e[i(2190)], 0)),
          (this[i(1035)] = e),
          this[i(1180)]()
      }
      [t(2578)](e) {
        const i = t
        return Cesium[i(1637)][i(799)]
      }
      _mergeOpts() {
        const e = t
        this.color = Cesium[e(1154)][e(2008)](this[e(1035)].color)
      }
      [t(1474)](e, i) {
        const s = t
        return (
          ((i = Cesium[s(1960)](i, {}))[s(1070)] = Cesium[s(1169)][s(2133)](this[s(1070)], e)),
          (i[s(2535)] = this.speed),
          (i[s(638)] = this.repeat),
          i[s(259)] != this[s(259)] && (i[s(259)] = this[s(259)]),
          (i.sampling = this[s(2190)]),
          i
        )
      }
      [t(1525)](e) {
        const i = t
        return (
          this === e ||
          (e instanceof Qt &&
            Cesium[i(1169)][i(1525)](this[i(1070)], e.color) &&
            this[i(2535)] == e[i(2535)] &&
            this[i(259)] == e[i(259)] &&
            this[i(638)] == e[i(638)] &&
            this[i(2190)] == e[i(2190)])
        )
      }
    }
    Object[t(610)](Qt[t(1727)], { color: Cesium[t(1328)](t(1070)) })
    class Zt extends Ot {
      constructor(e = {}) {
        const i = t
        super(e), (this[i(1165)] = '渐变面'), this[i(589)](e)
      }
      [t(589)](e) {
        const i = t,
          s = { x: 0.5, y: 0.5 }
        ;(e[i(1070)] = Cesium[i(1960)](e[i(1070)], i(456))),
          (e[i(365)] = Cesium.defaultValue(e[i(365)], s)),
          (this._materialOpts = e),
          this.mergeOpts()
      }
      [t(2578)](e) {
        const i = t
        return Cesium[i(1637)][i(382)]
      }
      [t(1078)]() {
        const e = t
        this[e(1070)] = Cesium[e(1154)][e(2008)](this._materialOpts[e(1070)])
      }
      [t(1474)](e, i) {
        const s = t
        return (
          ((i = Cesium[s(1960)](i, {}))[s(1070)] = Cesium[s(1169)][s(2133)](this[s(1070)], e)),
          (i[s(365)] = this[s(365)]),
          i
        )
      }
      [t(1525)](e) {
        const i = t
        return (
          this === e ||
          (e instanceof Zt &&
            Cesium.Property[i(1525)](this[i(1070)], e[i(1070)]) &&
            this[i(365)].x == e.center.x &&
            this.center.y == e[i(365)].y)
        )
      }
    }
    Object[t(610)](Zt[t(1727)], { color: Cesium[t(1328)](t(1070)) })
    class Kt extends Ot {
      constructor(e = {}) {
        const i = t
        super(e), (this[i(1165)] = '水面'), this[i(589)](e)
      }
      [t(589)](e) {
        const i = t
        ;(e.baseWaterColor = Cesium[i(1960)](e.color, i(836))),
          (e.blendColor = Cesium[i(1960)](e[i(1794)], i(1475))),
          (e.frequency = Cesium[i(1960)](e[i(553)], 1e3)),
          (e[i(1769)] = Cesium[i(1960)](e[i(1769)], 0.01)),
          (e[i(2550)] = Cesium.defaultValue(e[i(2550)], 10)),
          (e[i(849)] = Cesium[i(1960)](e[i(849)], 0.5)),
          (e[i(378)] = Cesium[i(1960)](e.fadeFactor, 1)),
          (e[i(2345)] = Cesium[i(1960)](e[i(2345)], 1)),
          (this[i(1035)] = e),
          this.mergeOpts()
      }
      [t(2578)](e) {
        return t(2562)
      }
      [t(1078)]() {
        const e = t
        ;(this[e(2095)] = Cesium[e(1154)][e(2008)](this[e(1035)][e(1070)])),
          (this[e(1794)] = Cesium.Color.fromCssColorString(this._materialOpts[e(1794)]))
      }
      [t(1474)](e, i) {
        const s = t
        return (
          (i.normalMap = this.normalMap),
          (i.bgMap = this[s(2204)]),
          (i.frequency = this[s(553)]),
          (i.animationSpeed = this.animationSpeed),
          (i[s(2550)] = this[s(2550)]),
          (i[s(849)] = this[s(849)]),
          (i[s(2095)] = this[s(2095)]),
          (i[s(1794)] = this[s(1794)]),
          (i.fadeFactor = this[s(378)]),
          (i[s(2345)] = this[s(2345)]),
          i
        )
      }
      [t(1525)](e) {
        const i = t
        return (
          this === e ||
          (e instanceof Kt &&
            Cesium.Property[i(1525)](this[i(1070)], e[i(1070)]) &&
            Cesium[i(1169)][i(1525)](this[i(1794)], e[i(1794)]) &&
            this[i(553)] == e[i(553)] &&
            this[i(1769)] == e.animationSpeed &&
            this[i(849)] == e[i(849)] &&
            this[i(2550)] == e[i(2550)] &&
            this[i(2345)] == e[i(2345)])
        )
      }
    }
    Object[t(610)](Kt[t(1727)], {
      color: Cesium[t(1328)]('color'),
      color: Cesium.createPropertyDescriptor(t(1794))
    })
    class Jt extends Ot {
      constructor(e = {}) {
        const i = t
        super(e), (this[i(1165)] = i(396)), this[i(589)](e)
      }
      [t(589)](e) {
        const i = t
        ;(e[i(1070)] = Cesium[i(1960)](e[i(1070)], i(2235))),
          (e[i(2254)] = Cesium.defaultValue(e.spacing, 200)),
          (e.width = Cesium[i(1960)](e[i(575)], 1)),
          (this[i(1035)] = e),
          this[i(1180)]()
      }
      [t(2578)](e) {
        const i = t
        return Cesium[i(1637)][i(734)]
      }
      [t(1078)]() {
        const e = t
        this[e(1070)] = Cesium[e(1154)][e(2008)](this[e(1035)].color)
      }
      getValue(e, i) {
        const s = t
        return (
          ((i = Cesium[s(1960)](i, {}))[s(1070)] = Cesium[s(1169)][s(2133)](this[s(1070)], e)),
          (i[s(2254)] = this[s(2254)]),
          (i[s(575)] = this[s(575)]),
          i
        )
      }
      equals(e) {
        const i = t
        return (
          this === e ||
          (e instanceof Jt &&
            Cesium[i(1169)][i(1525)](this.color, e[i(1070)]) &&
            this[i(2254)] == e[i(2254)] &&
            this[i(575)] == e[i(575)])
        )
      }
    }
    Object[t(610)](Jt[t(1727)], { color: Cesium[t(1328)](t(1070)) })
    class $t extends Ot {
      constructor(e = {}) {
        const i = t
        super(e), (this._name = i(2181)), this[i(589)](e)
      }
      setOpts(e) {
        const i = t
        ;(e.color = Cesium[i(1960)](e.color, i(456))),
          (e.speed = Cesium[i(1960)](e[i(2535)], 3)),
          (this._materialOpts = e),
          this.mergeOpts()
      }
      [t(1078)]() {
        const e = t
        this.color = Cesium[e(1154)][e(2008)](this[e(1035)].color)
      }
      [t(2578)](e) {
        const i = t
        return Cesium[i(1637)][i(1966)]
      }
      [t(1474)](e, i) {
        const s = t
        return (
          ((i = Cesium[s(1960)](i, {})).color = Cesium[s(1169)].getValueOrUndefined(
            this[s(1070)],
            e
          )),
          (i[s(2535)] = this[s(2535)]),
          i
        )
      }
      [t(1525)](e) {
        const i = t
        return (
          this === e ||
          (e instanceof $t &&
            Cesium.Property.equals(this.color, e[i(1070)]) &&
            Cesium.Property.equals(this[i(2535)], e[i(2535)]))
        )
      }
    }
    Object[t(610)]($t[t(1727)], {
      color: Cesium.createPropertyDescriptor(t(1070)),
      speed: Cesium[t(1328)](t(2535))
    })
    class te extends Ot {
      constructor(e = {}) {
        const i = t
        super(e), (this._name = '球体扫描材质'), this[i(589)](e)
      }
      setOpts(e) {
        const i = t
        ;(e[i(1070)] = Cesium[i(1960)](e[i(1070)], i(456))),
          (e[i(2535)] = Cesium[i(1960)](e[i(2535)], 3)),
          (this[i(1035)] = e),
          this.mergeOpts()
      }
      [t(1078)]() {
        const e = t
        this.color = Cesium.Color[e(2008)](this[e(1035)][e(1070)])
      }
      [t(2578)](e) {
        const i = t
        return Cesium[i(1637)][i(2099)]
      }
      [t(1474)](e, i) {
        const s = t
        return (
          ((i = Cesium.defaultValue(i, {}))[s(1070)] = Cesium.Property[s(2133)](this[s(1070)], e)),
          (i.speed = this[s(2535)]),
          i
        )
      }
      [t(1525)](e) {
        const i = t
        return (
          this === e ||
          (e instanceof te &&
            Cesium[i(1169)][i(1525)](this.color, e[i(1070)]) &&
            Cesium.Property[i(1525)](this[i(2535)], e[i(2535)]))
        )
      }
    }
    Object[t(610)](te[t(1727)], {
      color: Cesium[t(1328)]('color'),
      speed: Cesium[t(1328)](t(2535))
    })
    class ee extends Ot {
      constructor(e = {}) {
        const i = t
        super(e), (this[i(1165)] = i(2289)), this[i(589)](e)
      }
      setOpts(e) {
        const i = t
        ;(e[i(1070)] = Cesium.defaultValue(e.color, i(456))),
          (e[i(2535)] = Cesium[i(1960)](e[i(2535)], 3)),
          (this[i(1035)] = e),
          this.mergeOpts()
      }
      [t(1078)]() {
        const e = t
        this[e(1070)] = Cesium[e(1154)][e(2008)](this[e(1035)][e(1070)])
      }
      [t(2578)](e) {
        const i = t
        return Cesium[i(1637)][i(1535)]
      }
      [t(1474)](e, i) {
        const s = t
        return (
          ((i = Cesium[s(1960)](i, {})).color = Cesium[s(1169)][s(2133)](this[s(1070)], e)),
          (i[s(2535)] = this[s(2535)]),
          i
        )
      }
      [t(1525)](e) {
        const i = t
        return (
          this === e ||
          (e instanceof ee &&
            Cesium[i(1169)][i(1525)](this[i(1070)], e[i(1070)]) &&
            Cesium.Property.equals(this[i(2535)], e[i(2535)]))
        )
      }
    }
    Object.defineProperties(ee[t(1727)], {
      color: Cesium.createPropertyDescriptor(t(1070)),
      speed: Cesium[t(1328)](t(2535))
    })
    class ie extends Ot {
      constructor(e = {}) {
        const i = t
        super(e), (this._name = i(1612)), this.setOpts(e)
      }
      [t(589)](e) {
        const i = t
        ;(e[i(1070)] = Cesium[i(1960)](e[i(1070)], i(456))),
          (e[i(2535)] = Cesium.defaultValue(e.speed, 3)),
          (e[i(638)] = Cesium[i(1960)](e[i(638)], 2)),
          (e[i(1657)] = Cesium[i(1960)](e[i(1657)], 1)),
          (this[i(1035)] = e),
          this[i(1180)]()
      }
      [t(1078)]() {
        const e = t
        this[e(1070)] = Cesium[e(1154)][e(2008)](this._materialOpts.color)
      }
      getType(e) {
        const i = t
        return Cesium.Material[i(1890)]
      }
      getValue(e, i) {
        const s = t
        return (
          ((i = Cesium[s(1960)](i, {})).color = Cesium[s(1169)][s(2133)](this[s(1070)], e)),
          (i[s(2535)] = this[s(2535)]),
          (i[s(1657)] = this[s(1657)]),
          (i[s(638)] = this[s(638)]),
          (i[s(259)] = this[s(259)]),
          i
        )
      }
      [t(1525)](e) {
        const i = t
        return (
          this === e ||
          (e instanceof ie &&
            Cesium[i(1169)][i(1525)](this[i(1070)], e[i(1070)]) &&
            Cesium[i(1169)][i(1525)](this[i(2535)], e[i(2535)]) &&
            this[i(259)] == e[i(259)] &&
            this[i(638)] == e.repeat &&
            this[i(1657)] == e[i(1657)])
        )
      }
    }
    Object[t(610)](ie.prototype, {
      color: Cesium[t(1328)]('color'),
      speed: Cesium[t(1328)](t(2535))
    })
    class se extends Ot {
      constructor(e = {}) {
        const i = t
        super(e), (this[i(1165)] = i(1547)), this[i(589)](e)
      }
      [t(589)](e) {
        const i = t
        ;(e[i(1070)] = Cesium[i(1960)](e[i(1070)], i(456))),
          (e.breathe = Cesium[i(1960)](e[i(963)], !1)),
          (this[i(1035)] = e),
          this[i(1180)]()
      }
      [t(1078)]() {
        const e = t
        this[e(1070)] = Cesium[e(1154)][e(2008)](this._materialOpts[e(1070)])
      }
      [t(2578)](e) {
        const i = t
        return Cesium.Material[i(1220)]
      }
      getValue(e, i) {
        const s = t
        return (
          ((i = Cesium[s(1960)](i, {}))[s(1070)] = Cesium[s(1169)].getValueOrUndefined(
            this[s(1070)],
            e
          )),
          (i[s(963)] = this[s(963)]),
          i
        )
      }
      [t(1525)](e) {
        const i = t
        return (
          this === e ||
          (e instanceof se &&
            Cesium[i(1169)][i(1525)](this[i(1070)], e[i(1070)]) &&
            this.breathe == e.breathe)
        )
      }
    }
    Object[t(610)](se[t(1727)], { color: Cesium[t(1328)](t(1070)) })
    class ne extends Ot {
      constructor(e = {}) {
        const i = t
        super(e), (this._name = '墙体跑马灯材质'), this[i(589)](e)
      }
      [t(589)](e) {
        const i = t
        ;(e[i(1070)] = Cesium.defaultValue(e[i(1070)], i(456))),
          (e[i(2535)] = Cesium.defaultValue(e.speed, 3)),
          (e[i(638)] = Cesium[i(1960)](e[i(638)], 3)),
          (this._materialOpts = e),
          this[i(1180)]()
      }
      _mergeOpts() {
        const e = t
        this[e(1070)] = Cesium[e(1154)][e(2008)](this[e(1035)][e(1070)])
      }
      getType(e) {
        const i = t
        return Cesium[i(1637)][i(480)]
      }
      [t(1474)](e, i) {
        const s = t
        return (
          ((i = Cesium.defaultValue(i, {}))[s(1070)] = Cesium.Property[s(2133)](this.color, e)),
          (i[s(638)] = this[s(638)]),
          (i[s(2535)] = this[s(2535)]),
          i
        )
      }
      [t(1525)](e) {
        const i = t
        return (
          this === e ||
          (e instanceof ne &&
            Cesium[i(1169)].equals(this[i(1070)], e[i(1070)]) &&
            Cesium.Property[i(1525)](this.speed, e[i(2535)]) &&
            Cesium[i(1169)][i(1525)](this[i(638)], e[i(638)]))
        )
      }
    }
    Object[t(610)](ne[t(1727)], {
      color: Cesium.createPropertyDescriptor('color'),
      color: Cesium[t(1328)](t(638)),
      speed: Cesium[t(1328)](t(2535))
    })
    class oe extends Ot {
      constructor(e = {}) {
        const i = t
        super(e), (this._name = '墙体科幻材质'), this[i(589)](e)
      }
      [t(589)](e) {
        const i = t
        ;(e[i(1070)] = Cesium[i(1960)](e.color, i(456))),
          (e[i(2535)] = Cesium[i(1960)](e[i(2535)], 1)),
          (e[i(638)] = Cesium[i(1960)](e[i(638)], 5)),
          (this[i(1035)] = e),
          this[i(1180)]()
      }
      [t(1078)]() {
        const e = t
        this[e(1070)] = Cesium.Color.fromCssColorString(this[e(1035)][e(1070)])
      }
      [t(2578)](e) {
        const i = t
        return Cesium[i(1637)][i(2326)]
      }
      [t(1474)](e, i) {
        const s = t
        return (
          ((i = Cesium[s(1960)](i, {}))[s(1070)] = Cesium[s(1169)][s(2133)](this[s(1070)], e)),
          (i[s(2535)] = this[s(2535)]),
          (i[s(638)] = this[s(638)]),
          (i[s(259)] = this[s(1660)] + s(1228)),
          (i[s(2160)] = this.texture + s(2268)),
          (i[s(269)] = this[s(1660)] + s(1933)),
          i
        )
      }
      [t(1525)](e) {
        const i = t
        return (
          this === e ||
          (e instanceof oe &&
            Cesium[i(1169)][i(1525)](this[i(1070)], e[i(1070)]) &&
            Cesium[i(1169)][i(1525)](this[i(2535)], e[i(2535)]) &&
            this[i(638)] == e[i(638)])
        )
      }
    }
    Object.defineProperties(oe.prototype, {
      color: Cesium[t(1328)](t(1070)),
      speed: Cesium[t(1328)](t(2535))
    })
    const re = {
      ColorMaterialProperty: class extends Ot {
        constructor(e = {}) {
          const i = t
          super(e),
            (this[i(1296)] = new Cesium[i(808)]()),
            (this[i(1265)] = this._material),
            (this[i(1165)] = '纯色'),
            this[i(589)](e)
        }
        get proxy() {
          return this._proxy
        }
        [t(589)](e) {
          const i = t
          ;(e[i(1070)] = Cesium[i(1960)](e[i(1070)], 'rgba(255, 0, 0, 0.5)')),
            (this[i(1035)] = e),
            this[i(1180)]()
        }
        _mergeOpts() {
          const e = t
          this[e(1296)][e(1070)] = Cesium.Color[e(2008)](this[e(1035)][e(1070)])
        }
      },
      ImageMaterialProperty: class extends Ot {
        constructor(e = {}) {
          const i = t
          super(e),
            (this._material = new Cesium.ImageMaterialProperty()),
            (this[i(1265)] = this._material),
            (this._name = '图片'),
            this[i(589)](e)
        }
        get proxy() {
          return this[t(1265)]
        }
        [t(589)](e) {
          const i = t,
            s = { x: 1, y: 1 }
          ;(e.color = Cesium[i(1960)](e[i(1070)], i(2217))),
            (e[i(1983)] = Cesium.defaultValue(e.transparent, !1)),
            (e.repeat = Cesium[i(1960)](e[i(638)], s)),
            (this[i(1035)] = e),
            this[i(1180)]()
        }
        [t(1078)]() {
          const e = t
          ;(this[e(1296)][e(1070)] = Cesium[e(1154)][e(2008)](this[e(1035)][e(1070)])),
            (this[e(1296)][e(1983)] = this[e(1035)][e(1983)]),
            (this[e(1296)][e(638)] = this[e(1035)][e(638)]),
            (this._material[e(259)] = this[e(1035)][e(259)])
        }
      },
      CircleColorfulMaterialProperty: Bt,
      CircleWaveMaterialProperty: Vt,
      CircleLightRingMaterialProperty: Nt,
      CircleSpiralMaterialProperty: Ht,
      CircleImageDiffuseMaterialProperty: Gt,
      CircleImageRotateMaterialProperty: Wt,
      CircleScanMaterialProperty_1: Ut,
      CircleScanMaterialProperty_2: jt,
      CircleScanMaterialProperty_3: qt,
      CircleScanMaterialProperty_4: Yt,
      PolylineTrailMaterialProperty: Xt,
      PolylineFlowMaterialProperty: Qt,
      PolylineGlowMaterialProperty: class extends Ot {
        constructor(e = {}) {
          const i = t
          super(e),
            (this[i(1296)] = new Cesium.PolylineGlowMaterialProperty()),
            (this[i(1265)] = this[i(1296)]),
            (this[i(1165)] = i(2252)),
            this.setOpts(e)
        }
        get [t(892)]() {
          return this[t(1265)]
        }
        [t(589)](e) {
          const i = t
          ;(e.color = Cesium.defaultValue(e.color, 'rgba(255, 0, 0, 0.5)')),
            (e.glowPower = Cesium[i(1960)](e[i(1523)], 0.25)),
            (e.taperPower = Cesium.defaultValue(e.taperPower, 1)),
            (this._materialOpts = e),
            this.mergeOpts()
        }
        [t(1078)]() {
          const e = t
          this[e(1296)][e(1070)] = Cesium.Color[e(2008)](this[e(1035)][e(1070)])
        }
      },
      PolyGradientMaterialProperty: Zt,
      PolyWaterMaterialProperty: Kt,
      PolyElevationContourMaterialProperty: Jt,
      SphereElectricMaterialProperty: $t,
      SphereScanMaterialProperty: te,
      WallFlowMaterialProperty: ee,
      WallGradientMaterialProperty: class extends Ot {
        constructor(e = {}) {
          const i = t
          super(e), (this._name = i(936)), this[i(589)](e)
        }
        [t(589)](e) {
          const i = t,
            s = {}
          s[0] = i(1880)
          const n = { 1: 'rgba(32,153,230,0.5)' }
          let o = Cesium[i(1960)](e[i(1751)], [s, n])
          ;(e[i(1751)] = o),
            (this[i(259)] = (function (t) {
              const e = i
              let s = document[e(1945)](e(1493))
              ;(s[e(2306)] = 100), (s.width = 200)
              let n = s[e(1019)]('2d'),
                o = n[e(1624)](200, 50, 0, 50)
              return (
                t[e(1602)]((t) => {
                  const i = e
                  for (let e in t) o[i(803)](e, t[e])
                }),
                (n.fillStyle = o),
                n[e(412)](0, 0, 200, 100),
                n[e(555)](),
                s[e(1184)](e(1942)).replace(e(1942), e(1662))
              )
            })(o)),
            (this[i(1035)] = e),
            this[i(1180)]()
        }
        [t(1078)]() {}
        [t(2578)](e) {
          const i = t
          return Cesium[i(1637)][i(2558)]
        }
        [t(1474)](e, i) {
          const s = t
          return ((i = Cesium[s(1960)](i, {}))[s(259)] = this.image), i
        }
        [t(1525)](t) {
          return !1
        }
      },
      WallImageFlowMaterialProperty: ie,
      WallLightMaterialProperty: se,
      WallScrollMaterialProperty: ne,
      WallCoolMaterialProperty: oe
    }
    class ae {
      static create(e = 'color', i = {}) {
        const s = t
        switch (e) {
          case C.color:
            return new re.ColorMaterialProperty(i)
          case C.image:
            return new re[s(1100)](i)
          case C[s(899)]:
            return new re[s(700)](i)
          case C[s(2037)]:
            return new re.PolyWaterMaterialProperty(i)
          case C[s(1911)]:
            return new re[s(210)](i)
          case C[s(1103)]:
            return new re[s(798)](i)
          case C.trailLine:
            return new re[s(1006)](i)
          case C[s(1915)]:
            return new re[s(1847)](i)
          case C[s(326)]:
            return new re[s(1132)](i)
          case C[s(1294)]:
            return new re[s(1732)](i)
          case C[s(1823)]:
            return new re[s(1279)](i)
          case C[s(584)]:
            return new re.CircleLightRingMaterialProperty(i)
          case C[s(1064)]:
            return new re[s(569)](i)
          case C[s(1170)]:
            return new re[s(2012)](i)
          case C[s(479)]:
            return new re[s(2589)](i)
          case C.circleScan_4:
            return new re[s(1261)](i)
          case C[s(2481)]:
            return new re[s(1395)](i)
          case C.circleSpiral:
            return new re[s(629)](i)
          case C[s(1422)]:
            return new re[s(1389)](i)
          case C[s(1845)]:
            return new re[s(991)](i)
          case C[s(1677)]:
            return new re[s(952)](i)
          case C.wallGradient:
            return new re[s(1042)](i)
          case C.wallLight:
            return new re[s(2354)](i)
          case C[s(533)]:
            return new re[s(578)](i)
          case C[s(1376)]:
            return new re[s(1157)](i)
          case C[s(497)]:
            return new re[s(305)](i)
        }
      }
    }
    class he extends b {
      constructor(e) {
        const i = t
        super(e),
          (this[i(343)][i(506)] = Cesium.defaultValue(
            this[i(343)].selectedColor,
            'rgba(255,255,0,0.8)'
          )),
          (this[i(343)][i(224)] = Cesium[i(1960)](this._style[i(224)], i(1070))),
          (this[i(343)][i(1713)] = Cesium[i(1960)](this[i(343)][i(1713)], {})),
          (this._style.materialOpts[i(1070)] = Cesium.defaultValue(
            this[i(343)][i(1713)][i(1070)],
            i(2431)
          )),
          (this._materialType = this[i(343)].materialType),
          (this[i(1296)] = null),
          (this[i(1020)] = null)
      }
      [t(2239)]() {
        const e = t
        this._materialAttachTarget && (this[e(1296)] || this[e(1414)](), this[e(609)]())
      }
      [t(1414)]() {
        const e = t,
          i = this._style
        ;(this._material = ae[e(1951)](i[e(224)], i.materialOpts)),
          this[e(1296)][e(892)]
            ? (this[e(1020)][e(322)] = this[e(1296)].proxy)
            : (this[e(1020)][e(322)] = this._material)
      }
      [t(609)]() {
        const e = t,
          i = this[e(343)]
        if (i.material)
          return (
            (this._material = i.material),
            (this._materialAttachTarget[e(322)] = this[e(1296)]),
            this[e(1496)]
              ? ((this[e(2390)] = this[e(1296)][e(1070)] && this[e(1296)][e(1070)][e(2143)]),
                this._color && (this[e(2390)] = this[e(2390)][e(1886)]()),
                void (this[e(1296)][e(1070)] = Cesium.Color[e(2008)](this[e(343)].selectedColor)))
              : void (
                  this._color &&
                  (this[e(1296)][e(1070)] = Cesium[e(1154)].fromCssColorString(this[e(2390)]))
                )
          )
        if (this[e(1476)] != i[e(224)])
          return this[e(1414)](i), void (this._materialType = i.materialType)
        let s = { ...this[e(343)][e(1713)] }
        this[e(1496)] && (s[e(1070)] = this._style.selectedColor), this[e(1296)][e(589)](s)
      }
    }
    class le extends he {
      constructor(e) {
        const i = t
        super(e),
          (this[i(1115)] = M),
          (this[i(2336)] = Q[i(1351)]),
          (this[i(528)] = '线'),
          (this[i(314)] = i(2310)),
          (this[i(324)] = 2),
          (this[i(343)].show = Cesium.defaultValue(this[i(343)].show, !0)),
          (this[i(343)][i(575)] = Cesium[i(1960)](this[i(343)].width, 2)),
          (this._style[i(481)] = Cesium[i(1960)](this[i(343)][i(481)], Cesium[i(1597)][i(1722)])),
          (this[i(343)][i(2080)] = Cesium[i(1960)](this[i(343)][i(2080)], !1)),
          (this[i(2469)] = Cesium[i(1960)](e[i(2333)], [])),
          this[i(2515)](this[i(2469)]),
          (this[i(987)] = this[i(2547)]()),
          this[i(2231)]()
      }
      [t(415)](e) {
        const i = t
        this._entity[i(1482)] = e
      }
      [t(2547)]() {
        const e = t,
          i = {}
        return (
          (i.graphicId = this._id),
          (i[e(811)] = {}),
          (i[e(811)][e(1482)] = this._style[e(1482)]),
          (i[e(811)].positions = this[e(1378)]),
          (i[e(811)][e(575)] = 20),
          new Cesium[e(1380)](i)
        )
      }
      [t(2515)](e = []) {
        const i = t
        ;(this[i(2469)] = e),
          (this[i(1912)] = e[0]),
          (this[i(794)] = this._positions),
          this._positions.length < 2 ||
            ((this[i(1378)] = this[i(1264)](e)),
            this._cartesian3Array &&
              (this._style[i(160)] && this[i(1378)][i(2553)](this[i(1378)][0]),
              !this[i(2520)] && this[i(987)] && (this[i(987)].polyline[i(2333)] = this[i(1378)])))
      }
      [t(200)]() {
        const e = t
        let i = this[e(1378)]
        this._editMode &&
          (i = new Cesium[e(2569)](
            (t) => (this[e(1378)] && this[e(1378)][e(277)] > 1 ? this._cartesian3Array : []),
            !1
          )),
          (this[e(987)][e(811)][e(2333)] = i)
      }
      _generatePositions(e) {
        const i = t
        if (!(e.length < 2)) return this[i(1690)](e)
      }
      _convertPositions(e) {
        const i = t
        let s = Cesium[i(310)][i(774)]([][i(1500)][i(560)]([], e))
        return (this[i(1024)] = Cesium[i(1242)][i(951)](s)), s
      }
      [t(548)]() {
        const e = t
        for (const t in this[e(343)])
          if (Object[e(782)].call(this[e(343)], t)) {
            const i = this[e(343)][t]
            null != i &&
              null != i &&
              this[e(987)][e(811)][e(782)](t) &&
              (this[e(987)].polyline[t] = i)
          }
      }
      [t(2231)]() {
        const e = t
        this[e(548)]()
        const i = this._style
        ;(this[e(987)][e(811)][e(575)] = i[e(575)]),
          (this[e(987)][e(811)].clampToGround = i[e(2080)]),
          (this[e(987)][e(811)][e(481)] = i[e(481)]),
          i[e(978)] && (this[e(987)][e(811)][e(978)] = i[e(978)]),
          (this._entity[e(811)][e(1482)] = i.show),
          (this[e(1378)] = this[e(1264)](this[e(2469)])),
          this._editMode && this[e(1920)] && this[e(1920)](),
          this[e(2336)] != Q[e(975)] ||
            this[e(2520)] ||
            (this._entity[e(811)][e(2333)] = this._cartesian3Array),
          (this[e(1020)] = this[e(987)][e(811)]),
          this[e(2239)](),
          this[e(1762)]()
      }
      _setStyleHook() {}
      [t(1545)](e) {
        const i = t
        ;(this[i(987)][i(2122)] = e.id), (this._layer = e), e[i(245)][i(118)][i(1861)](this[i(987)])
      }
      _removeHook(e) {
        const i = t
        e[i(245)][i(118)][i(1896)](this[i(987)])
      }
    }
    let ce = {
      readFeature(e) {
        const i = t
        let s = e.geometry[i(2365)],
          n = e.properties
        return (n.positions = s), this.create(n)
      },
      create(e) {
        const i = t
        switch (e[i(1158)]) {
          case Q[i(1351)]:
            return new le(e)
          case Q[i(975)]:
            return new (class extends le {
              constructor(e) {
                const i = t
                super(e),
                  (this._graphicType = Q[i(975)]),
                  (this._typeName = i(2024)),
                  (this[i(350)] = 2),
                  (this[i(324)] = null),
                  (this[i(343)][i(381)] = Cesium[i(1960)](this[i(343)][i(381)], 1e4)),
                  (this[i(343)][i(1250)] = Cesium[i(1960)](this[i(343)][i(1250)], 50)),
                  (this[i(343)][i(2080)] = !1),
                  this._setStyle()
              }
              [t(1264)](e) {
                const i = t
                if (!(e[i(277)] < 2))
                  return (
                    (e = (function (t, e, s, n) {
                      const o = i
                      ;(t = Cesium[o(310)].fromDegrees(t[0], t[1], t[2])),
                        (e = Cesium[o(310)][o(667)](e[0], e[1], e[2]))
                      let r = []
                      if (0 === Cesium[o(310)][o(1849)](t, e)) return r
                      const a =
                          (function (t, e) {
                            const i = o,
                              s = Cesium[i(2285)].fromCartesian(t),
                              n = Cesium.Cartographic[i(2579)](e),
                              r = (180 * s[i(2106)]) / Math.PI,
                              a = (180 * s[i(199)]) / Math.PI,
                              h = (180 * n[i(2106)]) / Math.PI,
                              l = (180 * n[i(199)]) / Math.PI
                            return Math.sqrt((r - h) * (r - h) + (a - l) * (a - l))
                          })(t, e) * s,
                        h = Cesium.Cartesian3[o(1902)](t),
                        l = Cesium.Cartesian3.clone(e),
                        c = Cesium.Cartesian3.distance(h, Cesium.Cartesian3.ZERO),
                        u = Cesium[o(310)][o(1849)](l, Cesium[o(310)][o(1738)])
                      Cesium.Cartesian3[o(379)](h, h), Cesium.Cartesian3[o(379)](l, l)
                      const m = Cesium[o(310)][o(354)](h, l)
                      r[o(2553)](t)
                      for (let t = 1; t < n - 1; t++) {
                        const e = (1 * t) / (n - 1),
                          i = 1 - e,
                          s = Math[o(884)](i * m) / Math[o(884)](m),
                          p = Math.sin(e * m) / Math[o(884)](m),
                          d = Cesium[o(310)].multiplyByScalar(h, s, new Cesium.Cartesian3()),
                          f = Cesium[o(310)][o(606)](l, p, new Cesium.Cartesian3())
                        let C = Cesium[o(310)].add(d, f, new Cesium[o(310)]())
                        const v = e * Math.PI,
                          _ = c * i + u * e + Math[o(884)](v) * a
                        ;(C = Cesium.Cartesian3[o(606)](C, _, C)), r[o(2553)](C)
                      }
                      r[o(2553)](e)
                      for (let t = 0; t < r[o(277)]; t++) {
                        const e = r[t],
                          i = Cesium[o(2285)][o(2579)](e)
                        r[t] = [
                          Cesium[o(475)][o(363)](i[o(2106)]),
                          Cesium[o(475)][o(363)](i.latitude),
                          i[o(2306)]
                        ]
                      }
                      return r
                    })(e[0], e[e[i(277)] - 1], this[i(343)][i(381)], this[i(343)].count)),
                    this[i(1690)](e)
                  )
              }
            })(e)
          case Q.uprightLine:
            return new (class extends le {
              constructor(e = {}) {
                const i = t
                ;(e[i(1679)] = Cesium[i(1960)](e[i(1679)], {})),
                  (e[i(1679)][i(2080)] = !1),
                  (e[i(1679)][i(2306)] = Cesium[i(1960)](e[i(1679)][i(2306)], 100)),
                  e[i(2333)] && (e.style[i(2306)] = e[i(2333)][1][2] - e[i(2333)][0][2]),
                  super(e),
                  (this[i(2336)] = Q[i(1596)]),
                  (this[i(528)] = '竖线'),
                  (this[i(350)] = 1),
                  e[i(2251)] && this[i(1253)](e[i(2251)])
              }
              [t(1253)](e) {
                const i = t
                this._position = e
                let s = [e, [e[0], e[1], e[2] + this[i(343)][i(2306)]]]
                this[i(2515)](s)
              }
              [t(1762)]() {
                const e = t
                this[e(1912)] && this[e(1253)](this._position)
              }
            })(e)
        }
      }
    }
    let ue = {
      readFeature(e) {
        const i = t
        let s = e[i(138)][i(2365)],
          n = e[i(1004)]
        return (n.positions = s), this.create(n)
      },
      create(e) {
        const i = t
        if (e[i(1158)] === Q[i(1121)])
          return new (class extends he {
            constructor(e) {
              const i = t
              super(e),
                (this._graphicClassType = A),
                (this[i(2336)] = Q[i(1121)]),
                (this._typeName = '管道'),
                (this[i(314)] = i(2310)),
                (this[i(350)] = 2),
                (this[i(343)][i(1482)] = Cesium.defaultValue(this[i(343)][i(1482)], !0)),
                (this._style[i(1981)] = Cesium.defaultValue(this[i(343)].radius, 2)),
                (this[i(343)][i(1640)] = Cesium[i(1960)](this[i(343)][i(1640)], 0)),
                (this[i(343)].repeatDis = Cesium[i(1960)](this[i(343)][i(1324)], 0)),
                (this._positions = Cesium.defaultValue(e.positions, [])),
                this[i(2515)](this[i(2469)]),
                (this[i(987)] = this._createEntiy()),
                this._setStyle()
            }
            _setVisible(e) {
              const i = t
              this[i(987)][i(1482)] = e
            }
            [t(2547)]() {
              const e = t
              return new Cesium.Entity({
                graphicId: this[e(1570)],
                polylineVolume: {
                  show: this._style[e(1482)],
                  positions: this[e(1378)],
                  cornerType: Cesium[e(657)][e(824)]
                }
              })
            }
            [t(2515)](e = []) {
              const i = t
              ;(this[i(2469)] = e),
                (this[i(794)] = this._positions),
                this[i(2469)][i(277)] < 2 ||
                  ((this[i(1378)] = this[i(1264)](e)),
                  this[i(1378)] &&
                    (this._editMode &&
                      this[i(987)] &&
                      (this._entity.polylineVolume[i(2333)] = this[i(1378)]),
                    this[i(343)][i(1324)] > 0 &&
                      (this[i(2422)](),
                      this[i(1296)] &&
                        (this[i(1296)][i(638)] = this[i(343)].materialOpts[i(638)]))))
            }
            _generatePositions(e) {
              const i = t
              if (!(e[i(277)] < 2)) return this[i(1690)](e)
            }
            [t(1690)](e) {
              const i = t
              let s = Cesium[i(310)].fromDegreesArrayHeights([][i(1500)][i(560)]([], e))
              return (this[i(1024)] = Cesium.BoundingSphere[i(951)](s)), s
            }
            _computeShape() {
              const e = t
              let i = this[e(2491)]()
              return 1 == this._style[e(1640)] && (i = this[e(2297)]()), i
            }
            _computeSquareShape() {
              const e = t,
                i = this[e(343)].radius,
                s = []
              return (
                s.push(new Cesium[e(194)](-i, -i)),
                s[e(2553)](new Cesium[e(194)](i, -i)),
                s[e(2553)](new Cesium[e(194)](i, i)),
                s[e(2553)](new Cesium[e(194)](-i, i)),
                s
              )
            }
            [t(2491)]() {
              const e = t,
                i = this[e(343)][e(1981)],
                s = []
              for (let t = 0; t < 360; t += 1) {
                const n = Cesium.Math.toRadians(t)
                s[e(2553)](new Cesium.Cartesian2(i * Math.cos(n), i * Math[e(884)](n)))
              }
              return s
            }
            [t(548)]() {
              const e = t
              for (const t in this[e(343)])
                if (Object[e(782)].call(this[e(343)], t)) {
                  const i = this[e(343)][t]
                  null != i &&
                    null != i &&
                    this[e(987)][e(2482)][e(782)](t) &&
                    (this[e(987)][e(2482)][t] = i)
                }
            }
            [t(2231)]() {
              const e = t
              this[e(548)]()
              const i = this._style
              i[e(978)] && (this[e(987)][e(2482)][e(978)] = i.distanceDisplayCondition)
              let s = this[e(1413)]()
              ;(this[e(987)][e(2482)][e(1101)] = s),
                (this._entity[e(2482)][e(1482)] = i[e(1482)]),
                (this[e(1020)] = this._entity[e(2482)]),
                this[e(2422)](),
                this[e(2239)](),
                this[e(1762)]()
            }
            [t(2422)]() {
              const e = t
              if (
                this._style.repeatDis > 0 &&
                this[e(1378)] &&
                2 == this._cartesian3Array[e(277)]
              ) {
                let t =
                  Cesium[e(310)].distance(this[e(1378)][0], this._cartesian3Array[1]) /
                  this[e(343)][e(1324)]
                this[e(343)].materialOpts[e(638)] = parseInt(t)
              }
            }
            _setStyleHook() {}
            _addHook(e) {
              const i = t
              ;(this[i(987)][i(2122)] = e.id),
                (this[i(2462)] = e),
                e[i(245)][i(118)][i(1861)](this[i(987)])
            }
            _removeHook(e) {
              const i = t
              e._viewer[i(118)][i(1896)](this[i(987)])
            }
          })(e)
      }
    }
    class me extends he {
      constructor(e) {
        const i = t
        super(e),
          (this[i(1115)] = T),
          (this[i(2336)] = Q[i(2087)]),
          (this[i(314)] = 'Polygon'),
          (this[i(528)] = '面'),
          (this[i(324)] = 3)
        const s = Cesium.defaultValue(this._style[i(529)], {})
        ;(s[i(1482)] = Cesium[i(1960)](s[i(1482)], !0)),
          (s[i(575)] = Cesium[i(1960)](s.width, 1)),
          (s[i(1070)] = Cesium.defaultValue(s[i(1070)], i(741))),
          (s[i(506)] = Cesium[i(1960)](s[i(506)], s.color)),
          (this._style[i(529)] = s),
          (this[i(343)][i(2080)] = Cesium[i(1960)](this[i(343)][i(2080)], !0)),
          (this[i(1497)] = this[i(343)][i(2080)]),
          (this[i(343)].extrudedHeight = Cesium[i(1960)](this._style[i(398)], 0)),
          (this[i(343)][i(481)] = Cesium[i(1960)](this[i(343)][i(481)], Cesium[i(1597)][i(1722)])),
          this[i(1200)](),
          (this[i(2469)] = Cesium.defaultValue(e[i(2333)], [])),
          (this[i(261)] = []),
          (this[i(1516)] = null),
          this[i(2515)](this[i(2469)]),
          (this[i(1255)] = new Cesium[i(2229)](this[i(1378)])),
          (this[i(987)] = this[i(2547)]()),
          this[i(2231)]()
      }
      [t(1200)]() {}
      [t(415)](e) {
        const i = t
        this._entity[i(1482)] = e
      }
      [t(2547)]() {
        const e = t
        return new Cesium[e(1380)]({
          graphicId: this[e(1570)],
          polygon: {
            hierarchy: this[e(1255)],
            classificationType: Cesium.ClassificationType.BOTH,
            zIndex: 1
          },
          polyline: { show: this[e(343)][e(529)].show, positions: this[e(1378)] }
        })
      }
      [t(1733)]() {
        const e = t
        ;(this[e(987)][e(2087)][e(603)] = this._hierarchy),
          this[e(2462)]._viewer[e(118)].remove(this[e(987)]),
          (this._entity = this[e(2547)]()),
          (this[e(987)] = this._layer[e(245)][e(118)][e(1861)](this[e(987)])),
          (this[e(987)][e(2122)] = this[e(2462)].id),
          (this[e(987)][e(2087)][e(322)] = this[e(1296)][e(892)]
            ? this._material[e(892)]
            : this[e(1296)]),
          this.setPositions(this[e(2469)]),
          this[e(200)]()
      }
      setPositions(e = []) {
        const i = t
        if (((this[i(2469)] = e), (this[i(794)] = this[i(2469)]), this[i(2469)].length < 2)) return
        let s = []
        e[i(1602)]((t) => {
          s[i(2553)]([t[0], t[1]])
        }),
          (this[i(261)] = s),
          this[i(261)][i(2553)](this._points[0])
        const n = this._style
        !n.clampToGround && n.extrudedHeight > 0 && this[i(2309)](e),
          (this[i(1378)] = this[i(1264)](e)),
          this[i(1378)] &&
            (this[i(1378)][i(2553)](this[i(1378)][0]),
            (this[i(1255)] = new Cesium[i(2229)](this[i(1378)])))
      }
      [t(2309)](e) {
        const i = t
        ;(this._height = e[0][2]),
          this._height < 0 && (this[i(1516)] = 0),
          e[i(1602)]((t) => {
            const e = i
            t[2] = this[e(1516)]
          }),
          this[i(987)] && (this._entity[i(2087)][i(2306)] = this[i(1516)]),
          this[i(987)] && (this._entity[i(2087)][i(398)] = this[i(343)][i(398)] + this[i(1516)])
      }
      [t(200)]() {
        const e = t
        let i = this[e(1255)],
          s = this[e(1378)]
        this[e(2520)] &&
          ((i = new Cesium.CallbackProperty((t) => this[e(1255)], !1)),
          (s = new Cesium[e(2569)](
            (t) => (this[e(1378)] && this[e(1378)][e(277)] > 1 ? this[e(1378)] : []),
            !1
          ))),
          (this[e(987)].polygon.hierarchy = i),
          this[e(987)][e(811)] && (this[e(987)][e(811)].positions = s)
      }
      _generatePositions(e) {
        const i = t
        if (!(e[i(277)] < 2)) return this[i(1690)](e)
      }
      [t(1690)](e) {
        const i = t
        let s = Cesium[i(310)].fromDegreesArrayHeights([][i(1500)][i(560)]([], e))
        return (this[i(1024)] = Cesium.BoundingSphere[i(951)](s)), s
      }
      [t(2231)]() {
        const e = t,
          i = this._style
        this._clampToGround != this[e(343)].clampToGround &&
          (this[e(1733)](), (this[e(1497)] = this._style[e(2080)]))
        const s = i[e(529)],
          n = this._isSelected ? s[e(506)] : s[e(1070)]
        if (
          ((this[e(987)][e(811)][e(322)] = Cesium.Color.fromCssColorString(n)),
          (this[e(987)][e(811)][e(481)] = this[e(343)].classificationType),
          (this[e(987)].polyline[e(575)] = s[e(575)]),
          (this[e(987)].polyline[e(2080)] = i.clampToGround),
          (this[e(987)][e(2087)][e(1826)] = !1),
          (this[e(987)][e(811)][e(1482)] = s[e(1482)]),
          (this[e(987)].polygon[e(481)] = this[e(343)].classificationType),
          this[e(343)].clampToGround
            ? (this[e(987)][e(2087)].perPositionHeight = !1)
            : ((this[e(987)][e(2087)][e(1819)] = !0),
              (this[e(987)][e(2087)][e(1826)] = s.show && s[e(575)] > 0),
              (this[e(987)][e(2087)][e(725)] = Cesium.Color[e(2008)](n))),
          !this[e(343)][e(2080)] && this._style[e(398)] > 0)
        ) {
          let t = this._positions[e(277)] > 0 ? this[e(2469)][0][2] : 0
          t < 0 && (t = 0),
            (this[e(987)][e(2087)][e(2306)] = t),
            (this._entity[e(811)][e(1482)] = !1),
            (this[e(987)][e(2087)].extrudedHeight = t + this[e(343)][e(398)])
        }
        this[e(2520)] && this[e(1920)] && this._reAddHook(),
          (this[e(1020)] = this[e(987)][e(2087)]),
          this[e(2239)]()
      }
      [t(1545)](e) {
        const i = t
        ;(this[i(987)][i(2122)] = e.id),
          (this[i(2462)] = e),
          e._viewer[i(118)].add(this[i(987)]),
          (this[i(987)].polygon[i(1836)] = e[i(245)][i(118)][i(736)][i(277)])
      }
      [t(2389)](e) {
        const i = t
        e[i(245)].entities[i(1896)](this[i(987)])
      }
    }
    let pe = {
      readFeature(e) {
        const i = t
        let s = e[i(138)][i(2365)],
          n = e[i(1004)]
        return (n[i(2333)] = s), this[i(1951)](n)
      },
      create(e) {
        if (e[t(1158)] === Q.polygon) return new me(e)
      }
    }
    class de extends he {
      constructor(e) {
        const i = t
        super(e),
          (this._graphicClassType = E),
          (this._graphicType = Q[i(2584)]),
          (this[i(528)] = '圆'),
          (this[i(350)] = 2),
          (this[i(314)] = i(1865)),
          (this[i(343)][i(1981)] = Cesium[i(1960)](this._style.radius, 10)),
          (this[i(343)][i(2080)] = Cesium[i(1960)](this._style[i(2080)], !0)),
          (this[i(1497)] = this._style.clampToGround),
          (this[i(343)].classificationType = Cesium[i(1960)](
            this[i(343)][i(481)],
            Cesium[i(1597)][i(1722)]
          )),
          (this[i(343)][i(398)] = Cesium[i(1960)](this[i(343)][i(398)], 0)),
          (this[i(343)].outline = Cesium.defaultValue(this[i(343)][i(1826)], !1)),
          (this[i(343)][i(376)] = Cesium[i(1960)](this[i(343)][i(376)], 1)),
          (this._style[i(725)] = Cesium.defaultValue(this[i(343)].outlineColor, i(429))),
          (this[i(987)] = this._createEntiy())
        let s = Cesium[i(1960)](e[i(2251)], [111, 28, 0])
        this[i(1253)](s),
          this[i(2231)](),
          (this[i(2119)] = new Dt(e[i(2240)])),
          this[i(2119)][i(1567)](this[i(987)][i(2240)])
      }
      get [t(2240)]() {
        return this._label
      }
      set label(e) {
        const i = t
        ;(this[i(2119)] = e), this[i(2119)][i(1567)](this._entity[i(2240)])
      }
      _setVisible(e) {
        const i = t
        this[i(987)][i(1482)] = e
      }
      [t(2547)]() {
        const e = t,
          i = {}
        i[e(1836)] = 1
        const s = {}
        return (
          (s.graphicId = this[e(1570)]),
          (s[e(2251)] = this[e(2238)]),
          (s[e(1111)] = i),
          (s[e(2240)] = {}),
          new Cesium[e(1380)](s)
        )
      }
      [t(2515)](e) {
        const i = t
        ;(this[i(2469)] = e), (this[i(2469)][1][2] = this._positions[0][2])
        let s = Cesium[i(310)][i(774)]([].concat[i(560)]([], e))
        this.setPosition(this[i(2469)][0])
        let n = Cesium[i(310)][i(1849)](s[0], s[1])
        this[i(343)].radius = Number(n[i(1268)](2))
      }
      setPosition(e) {
        const i = t
        ;(this[i(1912)] = e),
          (this[i(2238)] = Cesium[i(310)][i(667)](e[0], e[1], e[2])),
          (this[i(794)] = this[i(1912)]),
          !this._editMode && (this[i(987)].position = this._cartesian3),
          this[i(343)][i(2080)] ||
            ((this._entity[i(1111)][i(2306)] = this[i(1912)][2]),
            (this[i(987)][i(1111)][i(398)] = this._style.extrudedHeight + this[i(1912)][2])),
          (this._boundingSphere = new Cesium[i(1242)](this[i(2238)], this[i(343)].radius)),
          this[i(343)][i(2080)] ||
            (this._style.extrudedHeight > this[i(343)][i(1981)] &&
              (this[i(1024)] = new Cesium[i(1242)](this[i(2238)], this[i(343)].extrudedHeight)))
      }
      [t(200)]() {
        const e = t
        let i = this[e(2238)],
          s = this[e(343)][e(1981)]
        this[e(2520)] &&
          ((i = new Cesium[e(2569)]((t) => this[e(2238)], !1)),
          (s = new Cesium[e(2569)]((t) => this[e(343)][e(1981)], !1))),
          (this[e(987)][e(1111)][e(2404)] = this[e(987)][e(1111)][e(2471)] = s),
          (this[e(987)][e(2251)] = i)
      }
      [t(548)]() {
        const e = t
        let i = { ...this[e(343)] }
        delete i[e(2306)], delete i[e(398)]
        const s = this[e(987)][e(1111)]
        for (const t in i)
          if (Object[e(782)].call(i, t)) {
            const e = i[t]
            null != e && null != e && (s[t] = e)
          }
        i[e(2080)] || (s[e(398)] = this[e(343)][e(398)] + this[e(1912)][2]),
          (s[e(725)] = Cesium.Color.fromCssColorString(i[e(725)]))
      }
      [t(2231)]() {
        const e = t,
          i = this._style
        this[e(1497)] != this._style[e(2080)]
          ? (this[e(1733)](), (this[e(1497)] = this[e(343)][e(2080)]))
          : this[e(1497)] || (this[e(987)][e(1111)].height = this[e(1912)][2]),
          this[e(548)](),
          this._editMode && this[e(1920)] && this._reAddHook(),
          this[e(2520)] ||
            (this[e(987)][e(1111)][e(2404)] = this[e(987)][e(1111)][e(2471)] = i[e(1981)]),
          (this._materialAttachTarget = this[e(987)][e(1111)]),
          this.updateMaterial()
      }
      _reAdd() {
        const e = t
        ;(this[e(987)][e(1111)][e(2251)] = this[e(2238)]),
          (this[e(987)][e(1111)].semiMajorAxis = this[e(987)][e(1111)][e(2471)] = 10),
          this[e(2462)][e(245)][e(118)][e(1896)](this[e(987)]),
          (this._entity = this._createEntiy()),
          (this[e(987)] = this._layer[e(245)][e(118)][e(1861)](this._entity)),
          (this[e(987)][e(2122)] = this._layer.id),
          (this._entity[e(1111)][e(322)] = this[e(1296)][e(892)]
            ? this._material[e(892)]
            : this[e(1296)]),
          this[e(1253)](this._position),
          this[e(200)](),
          this[e(2119)][e(1567)](this._entity[e(2240)]),
          this._style.clampToGround || (this[e(987)][e(1111)][e(2306)] = this[e(1912)][2])
      }
      [t(1545)](e) {
        const i = t
        ;(this[i(987)].layerId = e.id),
          (this[i(2462)] = e),
          (this[i(987)] = e[i(245)][i(118)].add(this[i(987)]))
      }
      [t(2389)](e) {
        const i = t
        e[i(245)].entities[i(1896)](this._entity)
      }
    }
    let fe = {
      readFeature(e) {
        const i = t
        let s = e.geometry.coordinates,
          n = e[i(1004)]
        return (n[i(2251)] = s), this.create(n)
      },
      create(e) {
        const i = t
        if (e[i(1158)] === Q[i(2584)]) return new de(e)
      }
    }
    class Ce extends he {
      constructor(e) {
        const i = t
        super(e),
          (this[i(1115)] = L),
          (this[i(2336)] = Q[i(1074)]),
          (this[i(314)] = i(1865)),
          (this._typeName = '球体'),
          (this[i(350)] = 2),
          (this._style[i(1981)] = Cesium[i(1960)](this[i(343)].radius, 10)),
          (this[i(343)].minimumClock = Cesium[i(1960)](this[i(343)][i(1319)], 0)),
          (this[i(343)].maximumClock = Cesium[i(1960)](this[i(343)].maximumClock, 360)),
          (this[i(343)][i(1574)] = Cesium.defaultValue(this[i(343)][i(1574)], 0)),
          (this[i(343)].maximumCone = Cesium[i(1960)](this._style[i(1663)], 180)),
          (this[i(1912)] = Cesium.defaultValue(e[i(2251)], [111, 28, 0])),
          (this[i(987)] = this._createEntiy()),
          this[i(1912)] && this[i(1253)](this[i(1912)]),
          this[i(2231)](),
          (this[i(2119)] = new Dt(e[i(2240)])),
          this[i(2119)][i(1567)](this[i(987)][i(2240)])
      }
      get [t(2240)]() {
        return this[t(2119)]
      }
      set [t(2240)](e) {
        const i = t
        ;(this[i(2119)] = e), this._label.merge(this[i(987)].label)
      }
      [t(415)](e) {
        const i = t
        this[i(987)][i(1482)] = e
      }
      [t(2547)]() {
        const e = t,
          i = {}
        return (
          (i[e(2366)] = this[e(1570)]), (i[e(1565)] = {}), (i[e(2240)] = {}), new Cesium[e(1380)](i)
        )
      }
      [t(2515)](e) {
        const i = t
        ;(this[i(2469)] = e), (this[i(2469)][1][2] = this._positions[0][2])
        let s = Cesium.Cartesian3[i(774)]([].concat[i(560)]([], e))
        this[i(1253)](this._positions[0])
        let n = Cesium[i(310)][i(1849)](s[0], s[1])
        this._style[i(1981)] = Number(n[i(1268)](2))
      }
      setPosition(e) {
        const i = t
        ;(this._position = e),
          (this[i(2238)] = Cesium.Cartesian3[i(667)](e[0], e[1], e[2])),
          this[i(987)] && !this._editMode && (this[i(987)][i(2251)] = this._cartesian3),
          (this[i(794)] = this._position),
          (this[i(1024)] = new Cesium[i(1242)](this[i(2238)], this[i(343)][i(1981)]))
      }
      _setEditMode() {
        const e = t,
          i = this[e(343)]
        let s = this[e(2238)],
          n = new Cesium.Cartesian3(i[e(1981)], i[e(1981)], i.radius)
        this[e(2520)] &&
          ((s = new Cesium[e(2569)]((t) => this[e(2238)], !1)),
          (n = new Cesium.CallbackProperty(
            (t) => new Cesium[e(310)](i[e(1981)], i[e(1981)], i[e(1981)]),
            !1
          ))),
          (this[e(987)][e(1565)].radii = n),
          (this[e(987)][e(2251)] = s)
      }
      [t(2231)]() {
        const e = t,
          i = this._style
        this[e(1804)] ||
          (this[e(987)][e(1565)].radii = new Cesium[e(310)](i.radius, i[e(1981)], i[e(1981)])),
          (this[e(987)][e(1565)][e(1319)] = Cesium[e(475)][e(1149)](i.minimumClock)),
          (this._entity[e(1565)].maximumClock = Cesium[e(475)].toRadians(i.maximumClock)),
          (this[e(987)][e(1565)][e(1574)] = Cesium[e(475)].toRadians(i[e(1574)])),
          (this[e(987)].ellipsoid[e(1663)] = Cesium[e(475)][e(1149)](i[e(1663)])),
          (this[e(1024)] = new Cesium[e(1242)](this._cartesian3, this[e(343)][e(1981)])),
          (this[e(1020)] = this[e(987)][e(1565)]),
          this[e(2239)]()
      }
      [t(1545)](e) {
        const i = t
        ;(this[i(987)][i(2122)] = e.id),
          (this[i(2462)] = e),
          e[i(245)][i(118)][i(1861)](this[i(987)])
      }
      [t(2389)](e) {
        const i = t
        e._viewer.entities[i(1896)](this._entity)
      }
    }
    let ve = {
      readFeature(e) {
        const i = t
        let s = e[i(138)][i(2365)],
          n = e[i(1004)]
        return (n[i(2251)] = s), this.create(n)
      },
      create(t) {
        if (t.graphicType === Q.sphere) return new Ce(t)
      }
    }
    class _e extends he {
      constructor(e) {
        const i = t
        super(e),
          (this[i(1115)] = z),
          (this._typeName = '墙体'),
          (this._style[i(1482)] = Cesium.defaultValue(this[i(343)][i(1482)], !0)),
          (this[i(343)][i(2306)] = Cesium[i(1960)](this[i(343)][i(2306)], 10)),
          (this[i(2469)] = e[i(2333)]),
          (this[i(2427)] = e[i(2332)]),
          (this[i(1326)] = e[i(1089)]),
          (this[i(987)] = this[i(2547)]()),
          (this._geometryType = i(2310))
      }
      _setVisible(t) {
        this._entity.show = t
      }
      [t(2547)]() {
        const e = t,
          i = {}
        return (i[e(2366)] = this[e(1570)]), (i[e(2467)] = {}), new Cesium.Entity(i)
      }
      [t(2515)](e = []) {
        const i = t
        ;(this[i(2469)] = e),
          this[i(2469)][i(277)] < 2 ||
            ((this[i(1378)] = this[i(1264)](e)),
            this[i(1378)] &&
              ((this[i(794)] = this[i(2469)]),
              this[i(1804)] ||
                ((this[i(987)].wall[i(2333)] = this[i(1378)]),
                (this[i(987)].wall[i(2332)] = this[i(487)]()),
                (this[i(987)][i(2467)].maximumHeights = this[i(2447)]()))))
      }
      [t(1253)](t) {}
      [t(487)]() {
        const e = t
        if (this[e(2427)]) return this[e(2427)]
        let i = []
        for (let t = 0; t < this._positions[e(277)]; t++) i[e(2553)](this._positions[t][2])
        return i
      }
      _covertMaximumHeights() {
        const e = t
        if (this[e(1326)]) return this._maximumHeights
        let i = []
        for (let t = 0; t < this[e(2469)][e(277)]; t++)
          i.push(this._positions[t][2] + this[e(343)].height)
        return i
      }
      [t(200)]() {
        const e = t
        let i = this._cartesian3Array,
          s = this[e(2447)](),
          n = this[e(487)]()
        this[e(2520)] &&
          ((i = new Cesium.CallbackProperty(
            (t) => (this[e(1378)] && this[e(1378)][e(277)] > 1 ? this[e(1378)] : []),
            !1
          )),
          (s = new Cesium[e(2569)]((t) => this[e(2447)](), !1)),
          (n = new Cesium[e(2569)]((t) => this[e(487)](), !1))),
          (this[e(987)].wall.positions = i),
          (this[e(987)][e(2467)][e(2332)] = n),
          (this._entity[e(2467)][e(1089)] = s)
      }
      _generatePositions(e) {
        if (!(e[t(277)] < 2)) return this._convertPositions(e)
      }
      [t(1690)](e) {
        const i = t
        let s = Cesium[i(310)][i(774)]([].concat[i(560)]([], e))
        return (this[i(1024)] = Cesium[i(1242)][i(951)](s)), s
      }
      [t(2231)]() {
        const e = t
        this._style,
          (this[e(1020)] = this[e(987)][e(2467)]),
          this[e(2239)](),
          this._graphicType == Q[e(738)]
            ? this[e(2515)](this._positions)
            : this[e(1253)](this[e(1912)]),
          this[e(2520)] && this._reAddHook && this[e(1920)]()
      }
      [t(1545)](e) {
        const i = t
        ;(this[i(987)][i(2122)] = e.id),
          (this[i(2462)] = e),
          e._viewer[i(118)][i(1861)](this._entity)
      }
      [t(2389)](e) {
        const i = t
        e[i(245)][i(118)][i(1896)](this._entity)
      }
      _getProperties() {
        const e = t,
          i = {}
        return (i[e(2332)] = this._minimumHeights), (i.maximumHeights = this[e(1326)]), i
      }
    }
    let ge = {
      readFeature(e) {
        const i = t
        let s = e[i(138)].coordinates,
          n = e[i(1004)]
        return (n[i(2333)] = s), (n.position = s), this.create(n)
      },
      create(e) {
        const i = t
        switch (e[i(1158)]) {
          case Q[i(738)]:
            return new (class extends _e {
              constructor(e) {
                const i = t
                super(e),
                  (this._minPointCount = 2),
                  (this[i(2336)] = Q[i(738)]),
                  (this._typeName = '简单墙体'),
                  this[i(2231)]()
              }
            })(e)
          case Q[i(1921)]:
            return new (class extends _e {
              constructor(e) {
                const i = t
                super(e),
                  (this._typeName = i(658)),
                  (this._fixPointCount = 2),
                  (this[i(2336)] = Q[i(1921)]),
                  (this[i(343)][i(212)] = Cesium.defaultValue(this[i(343)][i(212)], 120)),
                  this[i(343)][i(212)] < 3 && (this[i(343)][i(212)] = 3),
                  (this._style[i(1981)] = Cesium[i(1960)](this[i(343)].radius, 10)),
                  (this[i(1912)] = e[i(2251)]),
                  (this[i(314)] = i(1865)),
                  this._setStyle()
              }
              setPositions(e = []) {
                const i = t
                if (((this[i(2469)] = e), this[i(1253)](this._positions[0]), e[i(277)] < 2)) return
                let s = Cesium[i(310)].fromDegreesArrayHeights([][i(1500)].apply([], e)),
                  n = Cesium[i(310)][i(1849)](s[0], s[1])
                this[i(343)][i(1981)] = Number(n[i(1268)](2))
              }
              [t(1253)](e) {
                const i = t
                ;(this[i(1912)] = e || [0, 0, 0]),
                  (this[i(1378)] = this[i(1264)](this._position)),
                  (this[i(1024)] = Cesium[i(1242)].fromPoints(this._cartesian3Array)),
                  (this._coordinates = this._position),
                  this.editMode ||
                    ((this[i(987)].wall[i(2333)] = this[i(1378)]),
                    (this[i(987)][i(2467)][i(2332)] = this._covertMinimumHeights()),
                    (this[i(987)][i(2467)].maximumHeights = this[i(2447)]()))
              }
              [t(1264)](e) {
                const i = t
                e = Cesium[i(310)][i(667)](e[0], e[1], e[2])
                let s = Cesium[i(2058)][i(1648)](e),
                  n = []
                const o = this[i(343)][i(212)],
                  r = this[i(343)].radius
                for (let t = 0; t < o; t++) {
                  let e = (t / o) * Cesium[i(475)].TWO_PI,
                    a = Math[i(1272)](e),
                    h = Math.sin(e),
                    l = new Cesium[i(310)](a * r, h * r, 0)
                  n[i(2553)](Cesium[i(2066)][i(2483)](s, l, new Cesium[i(310)]()))
                }
                return n[i(2553)](n[0]), n
              }
              _covertMinimumHeights() {
                const e = t
                return new Array(this[e(343)][e(212)] + 1)[e(2327)](this._position[2])
              }
              [t(2447)]() {
                const e = t
                return new Array(this[e(343)][e(212)] + 1)[e(2327)](
                  this._position[2] + this._style[e(2306)]
                )
              }
            })(e)
          case Q[i(549)]:
            return new (class extends _e {
              constructor(e) {
                const i = t
                super(e),
                  (this[i(528)] = i(658)),
                  (this[i(350)] = 2),
                  (this[i(2336)] = Q[i(549)]),
                  (this[i(343)].slices = Cesium[i(1960)](this._style[i(212)], 120)),
                  this[i(343)][i(212)] < 3 && (this[i(343)][i(212)] = 3),
                  (this[i(343)][i(1981)] = Cesium.defaultValue(this[i(343)][i(1981)], 10)),
                  (this[i(343)][i(2535)] = Cesium[i(1960)](this[i(343)].speed, 0.5)),
                  (this[i(343)][i(2237)] = Cesium[i(1960)](this._style.dynamicHeight, !0)),
                  (this[i(2438)] = Cesium[i(1960)](this[i(343)][i(771)], 0.1)),
                  (this[i(2110)] = Cesium[i(1960)](this._style[i(397)], 0.01)),
                  (this[i(1397)] = this[i(343)].height),
                  (this[i(894)] = this[i(2110)]),
                  (this[i(314)] = i(1865)),
                  (this[i(1912)] = e.position),
                  (this[i(722)] = {}),
                  this[i(1253)](this._position),
                  this[i(2231)]()
              }
              [t(2515)](e = []) {
                const i = t
                if (((this[i(2469)] = e), this[i(1253)](this[i(2469)][0]), e[i(277)] < 2)) return
                let s = Cesium[i(310)].fromDegreesArrayHeights([].concat[i(560)]([], e)),
                  n = Cesium[i(310)].distance(s[0], s[1])
                this[i(343)].radius = Number(n[i(1268)](2))
              }
              setPosition(e) {
                const i = t
                ;(this[i(1912)] = e || [0, 0, 0]),
                  (e = Cesium[i(310)][i(667)](
                    this[i(1912)][0],
                    this[i(1912)][1],
                    this[i(1912)][2]
                  )),
                  (this[i(2238)] = e),
                  (this[i(1636)] = Cesium[i(2058)][i(1648)](e)),
                  (this._cartesian3Array = this[i(1264)](this[i(1912)])),
                  (this[i(1024)] = new Cesium[i(1242)](this[i(2238)], this[i(343)][i(1981)])),
                  (this[i(794)] = this[i(1912)])
                let s = new Cesium.CallbackProperty(
                    (t) =>
                      this[i(1378)] && this._cartesian3Array[i(277)] > 1
                        ? this[i(1264)](this._position)
                        : [],
                    !1
                  ),
                  n = new Cesium.CallbackProperty((t) => this[i(2447)](), !1),
                  o = new Cesium[i(2569)]((t) => this._covertMinimumHeights(), !1)
                ;(this[i(987)][i(2467)][i(2333)] = s),
                  (this[i(987)][i(2467)][i(2332)] = o),
                  (this[i(987)][i(2467)][i(1089)] = n)
              }
              [t(200)]() {}
              [t(1264)]() {
                const e = t
                let i = []
                ;(this[e(1397)] -= (this[e(343)][e(2306)] * this._style[e(2535)]) / 100),
                  (this[e(894)] += (this[e(343)].radius * this[e(343)].speed) / 100),
                  (this[e(894)] > this[e(343)].radius || this[e(1397)] < this._minHeight) &&
                    ((this._currentRadius = this._minRadius), (this[e(1397)] = this[e(343)].height))
                const s = this[e(343)][e(212)],
                  n = this[e(894)]
                if (this[e(722)]['_' + n]) return this[e(722)]['_' + n]
                for (let t = 0; t < s; t++) {
                  let o = (t / s) * Cesium[e(475)][e(271)],
                    r = Math.cos(o),
                    a = Math[e(884)](o),
                    h = new Cesium.Cartesian3(r * n, a * n, 0)
                  i.push(Cesium[e(2066)].multiplyByPoint(this[e(1636)], h, new Cesium[e(310)]()))
                }
                return i[e(2553)](i[0]), (this[e(722)]['_' + n] = i), i
              }
              [t(487)]() {
                const e = t
                return new Array(this[e(343)][e(212)] + 1)[e(2327)](this._position[2])
              }
              [t(2447)]() {
                const e = t
                let i = this[e(343)][e(2237)] ? this[e(1397)] : this[e(343)][e(2306)]
                return (i += this[e(1912)][2]), new Array(this[e(343)].slices + 1).fill(i)
              }
            })(e)
        }
      }
    }
    class ye extends b {
      constructor(e = {}) {
        const i = t
        super(e),
          (this[i(528)] = i(2504)),
          (this[i(1115)] = R),
          (this[i(314)] = i(1865)),
          (this[i(2299)] = !1),
          (this[i(1844)] = e[i(2171)])
        const s = this[i(1844)]
        if (!s || s.length < 4) throw new Cesium.DeveloperError(i(428))
        let n = s[i(1041)]('.')
        if (
          ((this._videoExt = s[i(347)](n + 1)), [i(451), i(1950), i(554)][i(877)](this[i(699)]) < 0)
        )
          throw new Cesium[i(2360)](i(2220))
        ;(this[i(1959)] = this[i(2423)]()), this[i(2073)]()
      }
      get [t(2171)]() {
        return this[t(1844)]
      }
      get [t(1109)]() {
        return this[t(1959)]
      }
      get playing() {
        return this[t(2299)]
      }
      [t(2073)]() {
        const e = t
        this._videoPlayer[e(2073)](), (this[e(2299)] = !0)
      }
      [t(2528)]() {
        const e = t
        this[e(1409)][e(2528)](), (this[e(2299)] = !1)
      }
      [t(2423)]() {
        const e = t,
          i = this[e(1844)]
        return e(451) == this[e(699)]
          ? this[e(1329)](i)
          : e(1950) == this._videoExt
            ? this[e(864)](i)
            : e(554) == this[e(699)]
              ? this._createFlvVideoEle(i)
              : void 0
      }
      [t(1329)](e) {
        const i = t
        let s = document[i(1945)]('video')
        s[i(793)]('control', !1),
          s[i(793)]('autoplay', i(2445)),
          s[i(793)]('preload', i(2363)),
          s.setAttribute(i(370), 'true'),
          s[i(793)](i(1947), !0),
          document.body[i(1621)](s)
        let n = i(2337) + new Date().getTime()
        return (
          s[i(793)]('id', n), s[i(793)](i(1744), e), (this[i(1409)] = s), document.getElementById(n)
        )
      }
      [t(864)](e) {
        const i = t
        let s = document[i(1945)](i(316))
        s[i(2453)][i(1861)](i(367)),
          s[i(2453)][i(1861)](i(2156)),
          s[i(793)]('controls', !0),
          s.setAttribute(i(2445), i(2445)),
          s[i(793)](i(2375), i(2363)),
          s[i(793)](i(1947), !0),
          document[i(879)][i(1621)](s)
        let n = i(2337) + new Date().getTime()
        s[i(793)]('id', n)
        let o = videojs(n)
        return (
          o.ready((t) => {
            const s = i,
              n = {}
            ;(n[s(1744)] = e), (n[s(1640)] = 'application/x-mpegURL')
            var r = [n]
            o[s(1744)](r), o.load()
          }),
          (this[i(1409)] = o),
          document[i(2081)](n + i(2187))
        )
      }
      [t(247)]() {
        const e = t
        let i = document[e(1945)](e(316))
        ;(i[e(1947)] = !0),
          i[e(793)](e(1947), e(1947)),
          i[e(793)](e(2445), !0),
          i[e(793)](e(370), 'loop'),
          i[e(793)](e(2408), ''),
          i[e(793)]('controls', !0),
          i[e(793)](e(746), -1),
          (i[e(1679)][e(2582)] = e(1772)),
          document[e(879)][e(1621)](i)
        const s = {}
        ;(s[e(1640)] = 'flv'), (s[e(2171)] = this[e(1844)])
        const n = flvjs[e(2587)](s)
        return n[e(1489)](i), n[e(2338)](), (this[e(1409)] = n), i
      }
      _getProperties() {
        const e = t,
          i = {}
        return (i[e(2171)] = this[e(1844)]), i
      }
    }
    class we {
      constructor(e) {
        const i = t
        ;(this[i(245)] = e[i(395)]),
          (this[i(343)] = e[i(1679)]),
          (this[i(2238)] = e[i(1131)]),
          (this[i(1365)] = this[i(1357)]()),
          this[i(2214)]()
      }
      [t(1357)]() {
        const e = t
        let i = new Cesium[e(2162)](this[e(245)].scene)
        i.setView({
          destination: this[e(2238)],
          orientation: {
            heading: Cesium[e(475)][e(1149)](this[e(343)][e(1198)]),
            pitch: Cesium.Math[e(1149)](this[e(343)][e(711)]),
            roll: 0
          }
        })
        const s = this._style.hAngle / this[e(343)].vAngle
        let n = this[e(343)][e(2519)]
        return (
          s < 1 && (n = this[e(343)][e(2493)]),
          (i[e(1850)] = new Cesium.PerspectiveFrustum({
            fov: Cesium[e(475)][e(1149)](n),
            aspectRatio: s,
            near: this._style.near,
            far: this[e(343)][e(348)]
          })),
          i
        )
      }
      _createGeometry() {
        const e = t
        let i = this[e(1365)],
          s = i[e(1258)],
          n = i.upWC,
          o = i[e(586)],
          r = new Cesium.Cartesian3(),
          a = new Cesium[e(819)](),
          h = new Cesium[e(1472)]()
        ;(o = Cesium.Cartesian3[e(875)](o, r)),
          Cesium.Matrix3[e(573)](a, 0, o, a),
          Cesium[e(819)][e(573)](a, 1, n, a),
          Cesium[e(819)][e(573)](a, 2, s, a)
        let l = Cesium.Quaternion.fromRotationMatrix(a, h),
          c = this[e(1365)][e(1850)],
          u = new Cesium[e(885)]({
            origin: this._cartesian3,
            orientation: l,
            frustum: c,
            _drawNearPlane: !0
          })
        this._frustumOutlineGeometry = u
        let m = new Cesium[e(1615)]({
          origin: this._cartesian3,
          orientation: l,
          frustum: c,
          vertexFormat: Cesium[e(1788)][e(2050)],
          _drawNearPlane: !0
        })
        this[e(375)] = m
      }
      [t(342)]() {
        const e = t
        let i = this[e(545)],
          s = new Float64Array(24)
        Cesium[e(1615)]._computeNearFarPlanes(i[e(1957)], i[e(291)], i._frustumType, i[e(2004)], s)
        let n = []
        for (let t = 12; t < 24; t += 3) n[e(2553)](new Cesium[e(310)](s[t], s[t + 1], s[t + 2]))
        return n
      }
    }
    class xe extends ye {
      constructor(e = {}) {
        const i = t
        super(e),
          (this._typeName = i(2250)),
          (this._fixPointCount = 2),
          (this[i(343)][i(1292)] = Cesium.defaultValue(this._style[i(1292)], 0.01)),
          (this[i(343)].far = Cesium[i(1960)](this._style[i(348)], 50)),
          (this[i(343)][i(2519)] = Cesium[i(1960)](this[i(343)][i(2519)], 10)),
          (this._style[i(2493)] = Cesium.defaultValue(this[i(343)][i(2493)], 10)),
          (this._style.alpha = Cesium.defaultValue(this[i(343)][i(2315)], 1)),
          (this[i(343)].debugFrustum = Cesium.defaultValue(this[i(343)][i(1146)], !0)),
          (this[i(343)][i(1968)] = Cesium[i(1960)](this[i(343)][i(1968)], !0)),
          (this[i(343)].debugFrustumColor = Cesium.defaultValue(this[i(343)][i(1818)], i(741))),
          (this[i(343)][i(506)] = Cesium[i(1960)](this._style.selectedColor, '#00FFFF')),
          (this[i(343)].heading = Cesium[i(1960)](this._style.heading, 0)),
          (this[i(343)][i(711)] = Cesium[i(1960)](this._style[i(711)], -90)),
          (this[i(1912)] = Cesium[i(1960)](e[i(2251)], [120, 40, 0])),
          (this[i(794)] = this[i(1912)]),
          (this._cartesian3 = Cesium[i(310)].fromDegrees(
            this[i(1912)][0],
            this[i(1912)][1],
            this._position[2]
          )),
          (this[i(1024)] = new Cesium.BoundingSphere(this[i(2238)], this[i(343)][i(348)]))
      }
      [t(2515)](e) {
        const i = t
        let s = []
        e.forEach((t) => {
          const e = a0_0x3b79
          let i = Cesium[e(310)][e(667)](t[0], t[1], t[2])
          s.push(i)
        }),
          (this._style[i(348)] = Cesium[i(310)].distance(s[0], s[1])),
          (this[i(343)][i(1198)] = $(s[0], s[1])),
          (this[i(343)][i(711)] = (function (t, e) {
            const s = i
            let n = new Cesium[s(310)](),
              o = Cesium[s(2058)][s(1648)](t)
            return (
              Cesium.Matrix4[s(814)](o, o),
              Cesium.Matrix4[s(2483)](o, e, n),
              Cesium.Cartesian3[s(379)](n, n),
              Cesium[s(475)][s(363)](Math[s(2463)](n.z))
            )
          })(s[0], s[1])),
          this.setPosition(e[0])
      }
      [t(1253)](e) {
        const i = t
        ;(this._position = Cesium[i(1960)](e, [0, 0, 0])),
          (this[i(2238)] = Cesium[i(310)][i(667)](
            this[i(1912)][0],
            this[i(1912)][1],
            this._position[2]
          )),
          (this[i(987)].position = this._cartesian3),
          (this[i(794)] = this[i(1912)]),
          this[i(2231)]()
      }
      [t(2231)]() {
        const e = t
        this[e(2462)][e(245)][e(696)][e(1346)][e(1896)](this._debugFrustum),
          (this[e(827)] = this._createVideoCamera()),
          this[e(1891)](),
          (this._debugFrustum = this[e(2128)]()),
          (this[e(1699)][e(2122)] = this._layer.id),
          (this[e(1699)][e(2366)] = this[e(1570)]),
          this[e(2462)][e(245)][e(696)].primitives[e(1861)](this._debugFrustum),
          this._entity && (this[e(987)][e(1482)] = this[e(343)][e(1968)]),
          (this._boundingSphere = new Cesium[e(1242)](this[e(2238)], this[e(343)][e(348)]))
      }
      [t(415)](e) {
        const i = t
        this[i(1699)] && (this[i(1699)][i(1482)] = e && this._style[i(1146)]), this[i(1418)](e)
      }
      [t(2128)]() {
        const e = t
        let i = this[e(1496)] ? this[e(343)][e(506)] : this[e(343)][e(1818)]
        return (
          (i = Cesium.Color[e(2008)](i)[e(329)](0.3)),
          new Cesium[e(2541)]({
            camera: this._videoCamera,
            color: i,
            show: this._style.debugFrustum
          })
        )
      }
      [t(2049)]() {
        const e = t,
          i = {}
        ;(i.viewer = this[e(2462)][e(245)]),
          (i[e(1131)] = this[e(2238)]),
          (i[e(1679)] = this[e(343)])
        let s = new we(i)
        return (this[e(1781)] = s), s._camera
      }
      [t(1251)]() {
        const e = t
        this[e(2462)][e(245)][e(1391)][e(2138)]({
          destination: this[e(2238)],
          orientation: { direction: this[e(827)].direction, up: this[e(827)].up },
          duration: 2
        })
      }
      _addOriginPoint() {
        const e = t,
          i = { pixelSize: 10 }
        ;(i.color = Cesium[e(1154)][e(588)]),
          (this._entity = this[e(2462)][e(245)][e(118)][e(1861)]({
            position: this[e(2238)],
            graphicId: this._id,
            show: this[e(343)].showOrigin,
            point: i
          })),
          (this[e(987)].layerId = this._layer.id)
      }
      [t(1545)](e) {
        const i = t
        ;(this[i(2462)] = e),
          (this[i(827)] = this._createVideoCamera()),
          this[i(1891)](),
          (this[i(1699)] = this[i(2128)]()),
          (this[i(1699)][i(2122)] = e.id),
          (this[i(1699)][i(2366)] = this[i(1570)]),
          this._layer[i(245)][i(696)][i(1346)][i(1861)](this[i(1699)]),
          this[i(1964)]()
      }
      [t(2389)](e) {
        const i = t
        e._viewer[i(118)][i(1896)](this._entity),
          e._viewer[i(696)][i(1346)].remove(this[i(1699)]),
          this[i(1959)][i(1896)](),
          this[i(1642)](e),
          this[i(2039)]()
      }
      [t(2039)]() {
        const e = t
        for (const t in this) Object[e(782)][e(1669)](this, t) && delete this[t]
      }
    }
    let be = {
      readFeature(e) {
        const i = t
        let s = e.geometry[i(2365)],
          n = e[i(1004)]
        return (n[i(2251)] = s), this[i(1951)](n)
      },
      create(e) {
        const i = t
        switch (e.graphicType) {
          case Q.videoFusion:
            return new (class extends xe {
              constructor(e = {}) {
                const i = t
                super(e), (this[i(528)] = i(1034)), (this[i(2336)] = Q.videoFusion)
              }
              [t(1891)]() {
                const e = t
                ;(this._viewShadowMap = this[e(234)]()),
                  this[e(1151)] ||
                    ((this[e(1151)] = this[e(419)]()),
                    this[e(2462)]._viewer.scene[e(2165)][e(1861)](this[e(1151)]),
                    this[e(2462)][e(245)][e(696)].primitives[e(1861)](this),
                    this._createTexture())
              }
              _setComponentVisible(e) {
                const i = t
                this[i(1151)] && (this[i(1151)][i(299)] = e)
              }
              [t(1504)]() {
                const e = t
                ;(this.activeVideoListener = () => {
                  const t = a0_0x3b79
                  this._videoTexture && this[t(1191)][t(2039)](),
                    (this[t(1191)] = new Cesium[t(1641)]({
                      context: this._layer[t(245)][t(696)][t(2325)],
                      source: this[t(1959)],
                      width: 1,
                      height: 1,
                      pixelFormat: Cesium[t(2277)].RGBA,
                      pixelDatatype: Cesium[t(547)][t(1382)]
                    })),
                    this._maskTexture && this[t(2291)][t(2039)](),
                    (this[t(2291)] = new Cesium[t(1641)]({
                      context: this[t(2462)][t(245)].scene[t(2325)],
                      source: i
                    }))
                }),
                  this[e(2462)]._viewer.clock.onTick.addEventListener(this[e(330)])
                let i = new Image()
                i[e(1744)] = this[e(343)][e(1617)]
              }
              [t(234)]() {
                const e = t,
                  i = {}
                return (
                  (i[e(2418)] = this[e(827)]),
                  (i[e(1036)] = !1),
                  (i[e(1390)] = 1),
                  (i[e(640)] = !1),
                  (i.isSpotLight = !0),
                  (i[e(364)] = !1),
                  (i.context = this[e(2462)][e(245)].scene[e(2325)]),
                  (i.pointLightRadius = this._style[e(348)]),
                  new Cesium[e(433)](i)
                )
              }
              [t(419)]() {
                const e = t
                let i = this,
                  s = i._viewShadowMap[e(1216)] ? i._viewShadowMap._pointBias : i[e(2373)][e(270)]
                return new Cesium[e(2506)]({
                  fragmentShader:
                    '\n uniform float opacity;\n uniform sampler2D colorTexture;\n uniform sampler2D stcshadow;\n uniform sampler2D videoTexture;\n uniform sampler2D maskTexture;\n uniform sampler2D depthTexture;\n uniform mat4 _shadowMap_matrix;\n uniform vec4 shadowMap_lightPositionEC;\n uniform vec4 shadowMap_normalOffsetScaleDistanceMaxDistanceAndDarkness;\n uniform vec4 shadowMap_texelSizeDepthBiasAndNormalShadingSmooth; \n in vec2 v_textureCoordinates;\n vec4 toEye( in vec2 uv, in float depth) {\n     vec2 xy = vec2((uv.x * 2. - 1.), (uv.y * 2. - 1.));\n     vec4 posInCamera = czm_inverseProjection * vec4(xy, depth, 1.);\n     posInCamera = posInCamera / posInCamera.w;\n     return posInCamera;\n \n }\n float getDepth( in vec4 depth) {\n     float z_window = czm_unpackDepth(depth);\n     z_window = czm_reverseLogDepth(z_window);\n     float n_range = czm_depthRange.near;\n     float f_range = czm_depthRange.far;\n     return (2. * z_window - n_range - f_range) / (f_range - n_range);\n }\n float _czm_sampleShadowMap(sampler2D shadowMap, vec2 uv) {\n     return texture(shadowMap, uv)\n         .r;\n }\n float _czm_shadowDepthCompare(sampler2D shadowMap, vec2 uv, float depth) {\n     return step(depth, _czm_sampleShadowMap(shadowMap, uv));\n }\n float _czm_shadowVisibility(sampler2D shadowMap, czm_shadowParameters shadowParameters) {\n     float depthBias = shadowParameters.depthBias;\n     float depth = shadowParameters.depth;\n     float nDotL = shadowParameters.nDotL;\n     float normalShadingSmooth = shadowParameters.normalShadingSmooth;\n     float darkness = shadowParameters.darkness;\n     vec2 uv = shadowParameters.texCoords;\n     depth -= depthBias;\n     vec2 texelStepSize = shadowParameters.texelStepSize;\n     float radius = 1.;\n     float dx0 = -texelStepSize.x * radius;\n     float dy0 = -texelStepSize.y * radius;\n     float dx1 = texelStepSize.x * radius;\n     float dy1 = texelStepSize.y * radius;\n     float visibility = (_czm_shadowDepthCompare(shadowMap, uv, depth) + _czm_shadowDepthCompare(shadowMap, uv + vec2(dx0, dy0), depth) + _czm_shadowDepthCompare(shadowMap, uv + vec2(0., dy0), depth) + _czm_shadowDepthCompare(shadowMap, uv + vec2(dx1, dy0), depth) + _czm_shadowDepthCompare(shadowMap, uv + vec2(dx0, 0.), depth) + _czm_shadowDepthCompare(shadowMap, uv + vec2(dx1, 0.), depth) + _czm_shadowDepthCompare(shadowMap, uv + vec2(dx0, dy1), depth) + _czm_shadowDepthCompare(shadowMap, uv + vec2(0., dy1), depth) + _czm_shadowDepthCompare(shadowMap, uv + vec2(dx1, dy1), depth)) * (1. / 9.);\n     return visibility;\n }\n vec3 pointProjectOnPlane( in vec3 planeNormal, in vec3 planeOrigin, in vec3 point) {\n     vec3 v01 = point - planeOrigin;\n     float d = dot(planeNormal, v01);\n     return (point - planeNormal * d);\n }\n float ptm(vec3 pt) {\n     return sqrt(pt.x * pt.x + pt.y * pt.y + pt.z * pt.z);\n }\n void main() {\n     const float PI = 3.141592653589793;\n     vec4 color = texture(colorTexture, v_textureCoordinates);\n     vec4 currD = texture(depthTexture, v_textureCoordinates);\n     if (currD.r >= 1.) {\n         out_FragColor = color;\n         return;\n     }\n     float depth = getDepth(currD);\n     vec4 positionEC = toEye(v_textureCoordinates, depth);\n     vec3 normalEC = vec3(1.);\n     czm_shadowParameters shadowParameters;\n     shadowParameters.texelStepSize = shadowMap_texelSizeDepthBiasAndNormalShadingSmooth.xy;\n     shadowParameters.depthBias = shadowMap_texelSizeDepthBiasAndNormalShadingSmooth.z;\n     shadowParameters.normalShadingSmooth = shadowMap_texelSizeDepthBiasAndNormalShadingSmooth.w;\n     shadowParameters.darkness = shadowMap_normalOffsetScaleDistanceMaxDistanceAndDarkness.w;\n     shadowParameters.depthBias *= max(depth * .01, 1.);\n     vec3 directionEC = normalize(positionEC.xyz - shadowMap_lightPositionEC.xyz);\n     float nDotL = clamp(dot(normalEC, -directionEC), 0., 1.);\n     vec4 shadowPosition = _shadowMap_matrix * positionEC;\n     shadowPosition /= shadowPosition.w;\n     if (any(lessThan(shadowPosition.xyz, vec3(0.))) || any(greaterThan(shadowPosition.xyz, vec3(1.)))) {\n         out_FragColor = color;\n         return;\n     }\n     shadowParameters.texCoords = shadowPosition.xy;\n     shadowParameters.depth = shadowPosition.z;\n     shadowParameters.nDotL = nDotL;\n     float visibility = _czm_shadowVisibility(stcshadow, shadowParameters);\n     vec4 videoColor = texture(videoTexture, shadowPosition.xy); \n     vec4 maskColor = texture(maskTexture, shadowPosition.xy);\n     videoColor *= maskColor;  \n     if (visibility == 1.) {\n         out_FragColor = mix(color, vec4(videoColor.xyz, 1.), opacity * videoColor.a);\n     } else {\n         if (abs(shadowPosition.z - 0.) < .01) {\n             return;\n         }\n         out_FragColor = color;\n     }\n }',
                  uniforms: {
                    hiddenAreaColor: function () {
                      return new Cesium[e(2141)](1, 0, 0, 1)
                    },
                    opacity: function () {
                      const t = e
                      return i[t(343)][t(2315)]
                    },
                    stcshadow: function () {
                      const t = e
                      return i[t(2373)][t(1247)]
                    },
                    videoTexture: function () {
                      const t = e
                      return i[t(1191)] || i[t(2373)][t(1247)]
                    },
                    maskTexture: function () {
                      const t = e
                      return i[t(2291)] || i._videoTexture || i[t(2373)][t(1247)]
                    },
                    _shadowMap_matrix: function () {
                      const t = e
                      return i[t(2373)][t(829)]
                    },
                    shadowMap_lightPositionEC: function () {
                      const t = e
                      return i[t(2373)][t(600)]
                    },
                    shadowMap_texelSizeDepthBiasAndNormalShadingSmooth: function () {
                      const t = e
                      let n = new Cesium[t(194)]()
                      return (
                        (n.x = 1 / i._viewShadowMap[t(1188)].x),
                        (n.y = 1 / i[t(2373)]._textureSize.y),
                        Cesium[t(2141)][t(2434)](
                          n.x,
                          n.y,
                          s[t(1658)],
                          s[t(1379)],
                          this.combinedUniforms1
                        )
                      )
                    },
                    shadowMap_normalOffsetScaleDistanceMaxDistanceAndDarkness: function () {
                      const t = e
                      return Cesium.Cartesian4[t(2434)](
                        s[t(184)],
                        i[t(2373)][t(2301)],
                        i._viewShadowMap[t(1288)],
                        i[t(2373)][t(585)],
                        this.combinedUniforms2
                      )
                    }
                  }
                })
              }
              update() {
                const e = t
                this[e(2373)] &&
                  this[e(1144)] &&
                  this[e(2462)][e(245)][e(696)].frameState[e(1747)].push(this[e(2373)])
              }
              _removeComponent(e) {
                const i = t
                e[i(245)].scene[i(2165)][i(1896)](this[i(1151)]),
                  e[i(245)].clock.onTick.removeEventListener(this[i(330)]),
                  this._videoTexture && this[i(1191)][i(2039)](),
                  e[i(245)][i(696)][i(1346)][i(1896)](this)
              }
            })(e)
          case Q[i(608)]:
            return new (class extends xe {
              constructor(e = {}) {
                const i = t
                super(e), (this._typeName = i(135)), (this[i(2336)] = Q[i(608)])
              }
              [t(1891)]() {
                const e = t
                let i = this[e(1969)]()
                this[e(496)](i),
                  this[e(2357)] ||
                    ((this[e(2357)] = this[e(1594)]()),
                    this[e(2462)][e(245)][e(118)][e(1861)](this[e(2357)]))
              }
              [t(496)](e) {
                const i = t
                let s = new Cesium.Cartesian3()
                e[i(1602)]((t) => {
                  const e = i
                  Cesium[e(310)][e(1861)](t, s, s)
                }),
                  (s = Cesium.Cartesian3[i(2317)](s, 4, new Cesium.Cartesian3()))
                let n = Cesium[i(2058)][i(1648)](this._cartesian3),
                  o = Cesium[i(2066)][i(814)](n, new Cesium[i(2066)]()),
                  r = Cesium[i(2066)].multiplyByPoint(o, s, new Cesium[i(310)]())
                ;(r = Cesium[i(310)][i(379)](r, r)),
                  (r = Cesium[i(310)][i(875)](r, r)),
                  (this[i(1580)] = new Cesium[i(2314)](
                    r,
                    Cesium[i(310)][i(1849)](this._cartesian3, s)
                  )),
                  (this._dimensions = new Cesium[i(194)](
                    Cesium[i(310)][i(1849)](e[0], e[1]),
                    Cesium.Cartesian3[i(1849)](e[0], e[3])
                  ))
              }
              _computeHierarchy() {
                const e = t
                return this[e(1781)][e(342)]()
              }
              _createVideoPlane() {
                const e = t
                return new Cesium[e(1380)]({
                  layerId: this._layer.id,
                  graphicId: this.id,
                  position: new Cesium[e(2569)]((t) => this._cartesian3, !1),
                  plane: {
                    plane: new Cesium.CallbackProperty((t) => this[e(1580)], !1),
                    dimensions: new Cesium[e(2569)]((t) => this[e(1613)], !1),
                    material: this[e(1959)]
                  }
                })
              }
              [t(1418)](e) {
                const i = t
                this[i(2357)] && (this._videoPlane[i(1482)] = e)
              }
              [t(1642)](e) {
                const i = t
                e[i(245)][i(118)].remove(this[i(2357)])
              }
            })(e)
        }
      }
    }
    let Se = {
      readFeature(e) {
        const i = t
        let s = e.geometry[i(2365)],
          n = e.properties
        return (n[i(2251)] = s), this[i(1951)](n)
      },
      create(e) {
        if (e.graphicType === Q.svg)
          return new (class extends kt {
            constructor(e = {}) {
              const i = t
              e[i(909)] && (e[i(909)][i(1323)] = !1),
                super(e),
                (this[i(1115)] = k),
                (this[i(528)] = i(2440)),
                (this._graphicType = Q[i(1396)])
              const s = this[i(343)]
              ;(s[i(1837)] = Cesium[i(1960)](s.rotate, 0)),
                (s[i(2313)] = Cesium[i(1960)](s[i(2313)], 1)),
                (s[i(2022)] = Cesium[i(1960)](s.scaleWidth, 1)),
                (s[i(1040)] = Cesium[i(1960)](s[i(1040)], 4)),
                (s[i(2281)] = Cesium.defaultValue(s[i(2281)], 0.2)),
                (s[i(2531)] = Cesium.defaultValue(s[i(2531)], 4)),
                (s.minScaleWidth = Cesium.defaultValue(s.minScaleWidth, 0.2)),
                (s[i(2313)] = Math[i(1113)](s.minScaleHeight, s.scaleHeight)),
                (s[i(2313)] = Math[i(1429)](s[i(1040)], s[i(2313)])),
                (s.scaleWidth = Math[i(1113)](s[i(258)], s[i(2022)])),
                (s[i(2022)] = Math[i(1429)](s.maxScaleWidth, s[i(2022)])),
                (this._style = s),
                (s[i(2327)] = Cesium[i(1960)](s[i(2327)], i(1626))),
                (s[i(189)] = Cesium.defaultValue(s[i(2327)], i(2159))),
                (s[i(1521)] = Cesium[i(1960)](s[i(1521)], 1)),
                (s.selectedColor = Cesium.defaultValue(s.selectedColor, i(1655))),
                (this._image = this._billboard[i(259)]),
                this[i(1444)]()
            }
            [t(1444)]() {
              const e = t,
                i = this[e(556)].image
              if (!i) return
              const s = i.substring(i[e(1041)]('.') + 1)
              if (i && e(1396) != s) throw new Cesium[e(2360)]('不支持的文件类型！')
              this[e(2192)] && this[e(232)] == i
                ? this[e(770)]()
                : (this[e(1740)](i), (this[e(232)] = i))
            }
            _fetchSvg(e) {
              const i = t
              fetch(e)
                [i(687)]((t) => t[i(1982)]())
                [i(687)]((t) => {
                  const e = i
                  let s = new DOMParser()[e(1079)](t, e(574))[e(630)](e(1396))[0]
                  if (!(s instanceof SVGElement)) throw new Cesium.DeveloperError(e(1119))
                  ;(this[e(1462)] = parseFloat(s[e(1548)][e(2306)].value) || 64),
                    (this._originalWidth = parseFloat(s[e(1548)][e(575)][e(922)]) || 64),
                    (this[e(2192)] = s),
                    this[e(770)]()
                })
                .catch((t) => {
                  const e = i
                  throw new Cesium[e(2360)](e(453), t)
                })
            }
            [t(2231)]() {
              this._setSvgDom()
            }
            [t(770)]() {
              const e = t
              let i = this[e(343)][e(2313)],
                s = this[e(343)][e(2022)],
                n = Math.max(i, s),
                o = this[e(1462)] * n,
                r = this[e(1262)] * n
              this[e(2192)][e(793)](e(2306), o),
                this._svgDom[e(793)](e(575), r),
                this[e(230)](this[e(2192)]),
                this._updateBillboard(this[e(2192)])
            }
            [t(230)](e) {
              const i = t
              let s = this._isSelected ? this._style[i(506)] : this[i(343)][i(2327)],
                n = this._isSelected ? this._style[i(506)] : this[i(343)][i(189)]
              for (let t = 0; t < e[i(683)][i(277)]; t++)
                e[i(683)][t][i(1679)] &&
                  (this[i(343)][i(2327)] && (e.childNodes[t].style[i(2327)] = s),
                  this[i(343)].stroke && (e[i(683)][t][i(1679)][i(189)] = n),
                  this[i(343)][i(1521)] &&
                    (e.childNodes[t].style['stroke-width'] = this._style[i(1521)]))
            }
            _updateBillboard(e) {
              const i = t
              let s = Cesium.getTimestamp()
              if (s < this._lastUpdate + 100) return
              this[i(2057)] = s
              let n = new XMLSerializer()[i(1591)](e),
                o = new Image()
              o.src = i(2218) + btoa(unescape(encodeURIComponent(n)))
              let r = this.getSize()
              ;(this[i(987)][i(909)].image = o),
                (this[i(987)].billboard[i(2306)] = r.height),
                (this[i(987)][i(909)].width = r[i(575)]),
                (this[i(987)].billboard[i(2307)] = Cesium.Math[i(1149)](this._style.rotate))
            }
            getSize() {
              const e = t,
                i = {}
              return (
                (i.height = this[e(1462)] * this[e(343)][e(2313)]),
                (i[e(575)] = this[e(1262)] * this[e(343)][e(2022)]),
                i
              )
            }
            [t(1400)](e) {
              const i = t
              e &&
                ((this[i(343)][i(1837)] = e),
                (this[i(987)][i(909)].rotation = Cesium.Math[i(1149)](this[i(343)][i(1837)])))
            }
            setScale(e, i) {
              const s = t
              e &&
                ((this[s(1679)][s(2313)] = Math.max(this[s(1679)][s(2281)], e)),
                (this[s(1679)][s(2313)] = Math[s(1429)](this[s(1679)][s(1040)], e))),
                i &&
                  ((this[s(1679)][s(2022)] = Math[s(1113)](this[s(1679)][s(258)], i)),
                  (this[s(1679)][s(2022)] = Math.min(this[s(1679)][s(2531)], i))),
                this[s(1444)]()
            }
          })(e)
      }
    }
    class Pe extends b {
      constructor(e = {}) {
        const i = t
        super(e), (this._graphicClassType = D), (this[i(314)] = i(1865))
        let s = Cesium[i(1960)](e.position, [111, 28, 0])
        const n = { near: 0.1 }
        ;(n[i(348)] = 1e4),
          (this[i(248)] = new Cesium.Cartesian2()),
          this[i(1253)](s),
          (this._style[i(2274)] = Cesium[i(1960)](this._style[i(2274)], [0, 0])),
          (this._style[i(1491)] = Cesium.defaultValue(
            this._style.horizontalOrigin,
            Cesium[i(2324)][i(2064)]
          )),
          (this[i(343)][i(978)] = Cesium.defaultValue(this._style.distanceDisplayCondition, n)),
          (this[i(343)][i(645)] = Cesium[i(1960)](this[i(343)].pointerEvents, !0)),
          (this[i(343)][i(622)] = Cesium.defaultValue(this._style.userSelect, !1)),
          (this[i(1736)] = this[i(2271)]()),
          (this[i(350)] = 1),
          (this[i(779)] = 0),
          (this._lastTime = 0)
      }
      get domEle() {
        return this._dom
      }
      [t(1253)](e) {
        const i = t
        ;(this[i(1912)] = e),
          (this[i(794)] = this[i(1912)]),
          (this[i(2238)] = Cesium[i(310)][i(667)](
            this[i(2251)][0],
            this[i(2251)][1],
            this.position[2]
          )),
          (this[i(1024)] = new Cesium.BoundingSphere(this[i(2238)], 5))
      }
      [t(2515)](e) {
        const i = t
        Array[i(1108)](e) && this[i(1253)](e[0])
      }
      [t(2271)]() {
        const e = t,
          i = document[e(1945)](e(1120))
        return i.classList.add(e(2533)), i
      }
      _addHook(e) {
        const i = t
        ;(this[i(2462)] = e),
          e[i(245)][i(1392)][i(1621)](this._dom),
          (this[i(1736)][i(1803)] = (t) => {
            const e = i
            this.onDomClick && this[e(974)](this)
          })
      }
      [t(2389)](e) {
        const i = t
        e[i(245)][i(1392)][i(2269)](this._dom)
      }
      [t(415)](e) {
        const i = t
        if ((this[i(1736)][i(2453)].remove(i(1532)), this[i(1736)][i(2453)][i(1896)](i(558)), e))
          return (
            this[i(1736)][i(2453)][i(1861)](i(1532)),
            (this._dom.style[i(645)] = this[i(343)][i(645)] ? i(2363) : i(1772)),
            void (this[i(1736)][i(1679)][i(622)] = this[i(343)][i(622)] ? i(1982) : i(1772))
          )
        this[i(1736)][i(2453)][i(1861)](i(558)), (this[i(1736)][i(1679)][i(645)] = i(1772))
      }
      [t(1550)](e) {
        const i = t
        let s = new Date()[i(1566)]()
        if (s - this[i(1830)] < 10) return
        this[i(1830)] = s
        let n = this[i(248)]
        const o = this[i(2462)]._viewer[i(696)][i(1493)][i(2306)]
        Cesium[i(2176)][i(708)](this[i(2462)][i(245)][i(696)], this[i(1131)], n),
          (this[i(1461)][i(1679)][i(2251)] = i(2076))
        let r = this[i(343)][i(2274)][0],
          a = this[i(343)][i(2274)][1]
        this[i(343)].horizontalOrigin == Cesium.HorizontalOrigin[i(2064)] &&
          (r -= this[i(1461)][i(168)] / 2),
          (this[i(1461)][i(1679)][i(1709)] = n.x + r + 'px'),
          (this[i(1461)].style[i(151)] = o - n.y + a + 'px')
        const h = this[i(779)]
        this[i(1461)][i(1679)].zIndex = 1e3 - Math.trunc((h / e) * 1e3)
        let l = this[i(343)][i(978)]
        l && l[i(348)] && (h < l[i(348)] && h > l.near ? this[i(415)](!0) : this[i(415)](!1))
      }
      _setStyle() {
        const e = t
        ;(this[e(1736)][e(1679)].pointerEvents = this[e(343)].pointerEvents ? e(2363) : e(1772)),
          (this[e(1736)][e(1679)][e(622)] = this._style[e(622)] ? e(1982) : e(1772)),
          this[e(2294)]()
      }
      _setSelectedStyle() {
        const e = t
        this[e(2296)](), this.isSelected && this[e(865)]()
      }
      [t(2296)]() {
        const e = t
        this._seletedDom && (this._seletedDom[e(1896)](), (this[e(1136)] = null))
      }
      _addSelectedDom() {
        const e = t,
          i = document[e(1945)](e(1120))
        ;(i[e(1679)].cssText =
          e(2144) + (this[e(1736)][e(281)] + 8) + e(1199) + (this[e(1736)][e(1481)] + 8) + e(1257)),
          this._dom[e(1621)](i),
          (this[e(1136)] = i)
      }
    }
    class Me extends Pe {
      constructor(e) {
        const i = t
        super(e),
          (this._typeName = '自定义模板'),
          (this[i(343)][i(2274)] = [0, 10]),
          (this[i(2336)] = Q[i(807)]),
          (this[i(992)] = Cesium.defaultValue(e[i(1922)], '')),
          (this[i(343)][i(2274)] = Cesium[i(1960)](this[i(343)].offset, [0, 10])),
          (this[i(1922)] = this[i(992)])
      }
      get [t(1922)]() {
        return this[t(992)]
      }
      set [t(1922)](e) {
        const i = t
        ;(this[i(992)] = e), (this[i(1736)].innerHTML = this[i(992)]), this[i(2231)]()
      }
      [t(1114)]() {
        const e = t,
          i = {}
        return (i.content = this[e(1922)]), i
      }
    }
    function Ae(e, i, s) {
      const n = t
      e && (e[n(493)] ? e[n(493)]('on' + i, s) : e[n(1973)] ? e[n(1973)](i, s) : (e['on' + i] = s))
    }
    function Te(e, i, s) {
      const n = t
      e &&
        (e[n(152)]
          ? e[n(152)]('on' + i, s)
          : e[n(704)]
            ? e.removeEventListener(i, s)
            : (e['on' + i] = null))
    }
    let Ee = {
      readFeature(e) {
        const i = t
        let s = e[i(138)][i(2365)],
          n = e.properties
        return (n[i(2251)] = s), this.create(n)
      },
      create(e) {
        const i = t
        switch (e[i(1158)]) {
          case Q.filletLabel:
            return new (class extends Pe {
              constructor(e) {
                const i = t
                if (
                  (super(e),
                  (this[i(528)] = '圆角标签'),
                  (this[i(2336)] = Q.filletLabel),
                  (this[i(343)][i(2274)] = [-2, 70]),
                  (this[i(343)].fontBold = Cesium[i(1960)](this[i(343)][i(822)], !1)),
                  (this[i(343)][i(1070)] = Cesium[i(1960)](this[i(343)][i(1070)], i(287))),
                  (this[i(343)].fontSize = Cesium[i(1960)](this[i(343)][i(735)], 14)),
                  (this._style[i(821)] = Cesium[i(1960)](this[i(343)][i(821)], [
                    '#070acbbf',
                    i(1450)
                  ])),
                  !Array[i(1108)](this._style[i(821)]) || 2 != this[i(343)][i(821)][i(277)])
                )
                  throw new Cesium[i(2360)](i(2074))
                this[i(431)](), (this[i(992)] = e[i(1922)]), (this.content = this[i(992)])
              }
              get content() {
                return this._content
              }
              set [t(1922)](e) {
                const i = t
                ;(this[i(992)] = e), (this[i(1430)].innerHTML = e), this[i(2231)]()
              }
              [t(431)]() {
                const e = t,
                  i = document[e(1945)]('div')
                ;(this[e(1430)] = i), i[e(2453)][e(1861)](e(1797)), this._dom[e(1621)](i)
                const s = document[e(1945)](e(1120))
                s[e(2453)][e(1861)](e(1507)), (this[e(844)] = s), this._dom.appendChild(s)
              }
              [t(2231)]() {
                const e = t,
                  i = this[e(343)][e(1070)],
                  s = this[e(343)][e(735)],
                  n = this[e(1679)][e(821)]
                let o = this[e(343)][e(822)] ? e(1339) : e(226)
                ;(this._contentDom[e(1679)].cssText =
                  'color:' +
                  i +
                  ';font-size:' +
                  s +
                  'px;font-weight:' +
                  o +
                  e(631) +
                  n[0] +
                  ',' +
                  n[1] +
                  ');'),
                  (this[e(844)][e(1679)][e(1811)] = e(1851) + n[1]),
                  this[e(2294)]()
              }
              [t(1114)]() {
                const t = {}
                return (t.content = this.content), t
              }
            })(e)
          case Q[i(2052)]:
            return new (class extends Pe {
              constructor(e = {}) {
                const i = t
                super(e),
                  (this[i(528)] = i(2379)),
                  (this._graphicType = Q[i(2052)]),
                  (this._style[i(2274)] = [0, 30]),
                  (this[i(343)][i(1070)] = Cesium[i(1960)](this[i(343)].color, i(1678))),
                  (this[i(343)][i(1981)] = Cesium[i(1960)](this[i(343)].radius, 10)),
                  (this._style.horizontalOrigin = Cesium[i(1960)](
                    this[i(343)].horizontalOrigin,
                    Cesium[i(2324)][i(2064)]
                  )),
                  this[i(431)](),
                  this[i(2231)]()
              }
              [t(431)]() {
                const e = t,
                  i = document[e(1945)](e(1120))
                i.classList.add('xt3d-divgraphic-animation')
                const s = document[e(1945)]('p')
                i[e(1621)](s), (this._contentDom = i), this[e(1736)][e(1621)](i)
              }
              [t(2231)]() {
                const e = t,
                  i = this[e(1679)][e(1070)]
                ;(this._contentDom.style[e(1070)] = i),
                  (this[e(1430)][e(1679)][e(575)] = this[e(1430)][e(1679)][e(2306)] =
                    this[e(343)][e(1981)] + 'px'),
                  this[e(2294)]()
              }
              [t(865)]() {
                const e = t,
                  i = document.createElement(e(1120))
                ;(i[e(1679)][e(1811)] =
                  'height: 30px; width: 29px; background: #ff000042; border: 1px dashed red; position: absolute;  left: -15px; top: -15px;pointer-events:none;z-index:-1'),
                  this._dom[e(1621)](i),
                  (this[e(1136)] = i)
              }
            })(e)
          case Q.dynamicBorder:
            return new (class extends Pe {
              constructor(e) {
                const i = t
                super(e),
                  (this[i(528)] = i(2494)),
                  (this[i(2336)] = Q[i(1810)]),
                  (this[i(343)][i(2274)] = [0, 26]),
                  (this[i(343)][i(1070)] = Cesium[i(1960)](this[i(343)][i(1070)], i(287))),
                  (this[i(343)].borderColor = Cesium[i(1960)](this._style[i(2198)], i(2042))),
                  (this._style[i(735)] = Cesium.defaultValue(this[i(343)][i(735)], 14)),
                  (this[i(992)] = Cesium[i(1960)](e[i(1922)], '')),
                  this._createContentDom(),
                  this._setStyle()
              }
              get [t(1922)]() {
                return this[t(992)]
              }
              set [t(1922)](e) {
                const i = t
                ;(this[i(992)] = e), (this[i(1430)][i(176)] = e), this._setStyle()
              }
              _createContentDom() {
                const e = t
                let i = document[e(1945)]('div')
                i[e(2453)].add(e(2429))
                let s = document.createElement(e(1120))
                s[e(2453)].add(e(891))
                let n = document[e(1945)](e(1081))
                n[e(2453)][e(1861)](e(988)),
                  s.appendChild(n),
                  (n[e(176)] = this[e(992)]),
                  i.appendChild(s),
                  this[e(1736)][e(1621)](i),
                  (this._contentDom = n)
              }
              [t(865)]() {
                const e = t,
                  i = document[e(1945)](e(1120))
                let s = this[e(1736)][e(1666)](e(2542))
                ;(i[e(1679)].cssText =
                  'height: ' + (s[e(281)] + 14) + e(1199) + (s[e(1481)] + 14) + e(1807)),
                  this[e(1736)][e(1621)](i),
                  (this._seletedDom = i)
              }
              [t(2231)]() {
                const e = t,
                  i = this._style.color
                ;(this[e(1430)][e(1679)].color = i),
                  (this[e(1430)].style[e(735)] = this._style[e(735)] + 'px')
                let s = this._dom[e(1666)](e(2542)),
                  n =
                    e(871) +
                    2 * this[e(343)].fontSize +
                    e(292) +
                    (this[e(992)][e(277)] * this._style[e(735)] + 40) +
                    'px'
                ;(this[e(1736)][e(1679)][e(1811)] = n),
                  (s[e(1679)].cssText =
                    e(775) +
                    this[e(343)][e(2198)] +
                    '; box-shadow: inset 0 0 0 1px ' +
                    this[e(343)][e(2198)] +
                    ';' +
                    n),
                  this[e(2294)]()
              }
              [t(1114)]() {
                const e = t,
                  i = {}
                return (i[e(1922)] = this[e(1922)]), i
              }
            })(e)
          case Q.brightenDiv:
            return new (class extends Pe {
              constructor(e) {
                const i = t
                super(e),
                  (this[i(528)] = i(1432)),
                  (this[i(2336)] = Q[i(163)]),
                  (this[i(2045)] = Cesium[i(1960)](e[i(1038)], '')),
                  (this[i(992)] = Cesium.defaultValue(e[i(1922)], '')),
                  (this[i(343)][i(2274)] = [0, 10]),
                  (this[i(343)].horizontalOrigin = Cesium[i(2324)][i(2451)]),
                  (this[i(343)][i(2246)] = Cesium[i(1960)](this._style[i(2246)], i(2193))),
                  (this[i(1922)] = this._content)
              }
              [t(2271)]() {
                const e = t
                let i = document[e(1945)](e(1120))
                return i[e(2453)][e(1861)](e(1687)), (this._contentDom = i), i
              }
              get [t(1038)]() {
                return this[t(2045)]
              }
              set [t(1038)](e) {
                ;(this[t(2045)] = e), this._setContent()
              }
              get [t(1922)]() {
                return this._content
              }
              set content(e) {
                const i = t
                ;(this[i(992)] = e), this[i(1332)]()
              }
              _setContent() {
                const e = t
                ;(this._contentDom[e(176)] =
                  e(1419) + this[e(2045)] + e(650) + this._content + e(2032)),
                  this[e(2231)]()
              }
              [t(2231)]() {
                const e = t
                this._dom[e(1666)](e(142))[e(1679)][e(1811)] = e(2361) + this[e(343)].theme + e(180)
                let i = 'background:' + this._style[e(2246)] + e(841) + this[e(343)][e(2246)] + ';'
                ;(this._dom[e(1666)](e(424)).style[e(1811)] = i),
                  (this[e(1736)][e(1666)](e(1388))[e(1679)][e(1811)] = i),
                  (this._dom[e(1666)](e(1388))[e(1679)][e(1811)] = i),
                  (this[e(1736)].querySelector(e(1877)).style[e(1811)] = i),
                  (this._dom[e(1666)]('.b-b-r').style.cssText = i),
                  (this[e(1736)][e(1666)](e(1623))[e(1679)].cssText = i),
                  (this[e(1736)][e(1666)]('.b-r')[e(1679)][e(1811)] = i),
                  (this[e(1736)][e(1666)](e(478))[e(1679)][e(1811)] = i),
                  (this._dom[e(1666)](e(2019)).style[e(1811)] = i),
                  (this._dom[e(1666)](e(1388))[e(1679)][e(1811)] = i)
                const s = Cesium.Color[e(2008)](this[e(343)][e(2246)])[e(329)](0.4)[e(1886)]()
                ;(this[e(1736)][e(1666)](e(1705))[e(1679)].cssText =
                  e(414) +
                  s +
                  ' 30px, ' +
                  s +
                  ' 50%, transparent 50%), linear-gradient( -45deg, transparent 30px, ' +
                  s +
                  e(1190) +
                  s +
                  e(2565)),
                  this._setSelectedStyle()
              }
              [t(1114)]() {
                const e = t,
                  i = {}
                return (i.title = this[e(1038)]), (i.content = this[e(1922)]), i
              }
            })(e)
          case Q.imageLabel:
            return new (class extends Pe {
              constructor(e) {
                const i = t
                super(e),
                  (this._typeName = i(1013)),
                  (this[i(2336)] = Q.imageLabel),
                  (this[i(343)][i(2274)] = [0, 10]),
                  (this[i(1231)] = e[i(185)]),
                  this[i(431)](),
                  this[i(2231)]()
              }
              get [t(185)]() {
                return this._imageUrl
              }
              set imageUrl(e) {
                const i = t
                ;(this[i(1231)] = e), (this[i(976)].src = this._imageUrl)
              }
              _createContentDom() {
                const e = t
                let i = document[e(1945)](e(2575))
                ;(i[e(1744)] = this[e(1231)]), this._dom[e(1621)](i), (this[e(976)] = i)
              }
              _setStyle() {
                const e = t
                this[e(343)][e(2306)] &&
                  (this._conetntDom[e(1679)][e(2306)] = this[e(343)].height + 'px'),
                  this[e(343)].width &&
                    (this[e(976)][e(1679)][e(575)] = this[e(343)][e(575)] + 'px'),
                  this._setSelectedStyle()
              }
              [t(1114)]() {
                const e = t,
                  i = {}
                return (i[e(185)] = this[e(1231)]), i
              }
            })(e)
          case Q[i(835)]:
            return new (class extends Pe {
              constructor(e = {}) {
                const i = t
                super(e),
                  (this._typeName = i(532)),
                  (this[i(2336)] = Q.liquidfill),
                  (this[i(992)] = Cesium[i(1960)](e[i(1922)], 0.5)),
                  (this[i(343)].offset = [0, 10]),
                  (this[i(343)][i(1070)] = Cesium[i(1960)](this[i(343)][i(1070)], '#ff9501')),
                  (this[i(343)].radius = Cesium[i(1960)](this[i(343)][i(1981)], 80)),
                  (this[i(343)][i(735)] = Cesium[i(1960)](this[i(343)][i(735)], 15)),
                  this[i(2231)]()
              }
              [t(2271)]() {
                const e = t,
                  i = document[e(1945)](e(1120))
                return (this[e(1692)] = i), i
              }
              [t(2442)]() {
                const e = t
                this[e(2126)] || (this._chart = echarts[e(1719)](this[e(1692)])), this[e(1583)]()
              }
              get [t(1922)]() {
                return this[t(992)]
              }
              set content(e) {
                const i = t
                ;(this._content = e), this[i(1583)]()
              }
              _setChartOption() {
                const e = t,
                  i = { show: !1 },
                  s = {}
                ;(s[e(1640)] = 'liquidFill'),
                  (s.data = [this[e(992)]]),
                  (s[e(1981)] = e(2101)),
                  (s[e(1826)] = i),
                  (s[e(2240)] = {}),
                  (s[e(1842)] = {}),
                  (s[e(2240)][e(2251)] = ['50%', '65%']),
                  (s[e(2240)][e(1909)] = {}),
                  (s[e(2240)][e(1909)].fontSize = this._style[e(735)]),
                  (s[e(2240)][e(1909)][e(2527)] = e(1293)),
                  (s[e(1842)][e(1070)] = this._style[e(1070)])
                const n = {}
                n[e(1700)] = [s]
                let o = n
                this[e(2126)][e(449)](o), this[e(2126)][e(809)]()
              }
              [t(2389)](e) {
                const i = t
                this[i(2126)][i(335)]()
              }
              [t(2231)]() {
                const e = t
                ;(this[e(1736)][e(1679)][e(1811)] =
                  'height:' + this[e(343)][e(1981)] + e(292) + this[e(343)][e(1981)] + e(2405)),
                  this._initChart(),
                  this[e(2294)]()
              }
              _getProperties() {
                const e = t,
                  i = {}
                return (i[e(1922)] = this[e(1922)]), i
              }
            })(e)
          case Q[i(1944)]:
            return new (class extends Pe {
              constructor(e) {
                const i = t
                super(e),
                  (this[i(528)] = i(2199)),
                  (this[i(2336)] = Q[i(1944)]),
                  (this[i(992)] = Cesium[i(1960)](e[i(1922)], [])),
                  (this[i(343)].type = Cesium[i(1960)](this[i(1679)][i(1640)], 1)),
                  (this[i(343)][i(1070)] = Cesium[i(1960)](this[i(343)][i(1070)], i(287))),
                  (this[i(343)].fontSize = Cesium[i(1960)](this[i(343)][i(735)], 14)),
                  (this[i(2473)] = []),
                  this[i(431)](),
                  this[i(854)](),
                  this[i(2231)]()
              }
              get [t(1922)]() {
                return this[t(992)]
              }
              set content(e) {
                const i = t
                ;(this[i(992)] = e), this[i(2419)](), this[i(854)](), this[i(2231)]()
              }
              _createContentDom() {
                const e = t
                let i = document[e(1945)](e(1120))
                i.classList[e(1861)]('xt3d-divgraphic-multiline'), (this._contentDom = i)
                let s = document[e(1945)](e(1120))
                ;(this[e(2426)] = s), this._dom[e(1621)](i), this[e(1736)][e(1621)](s)
              }
              [t(854)]() {
                const e = t
                for (let t = 0; t < this[e(992)].length; t++) {
                  const i = this[e(992)][t]
                  let s = document[e(1945)]('div')
                  s.classList[e(1861)](e(1179)),
                    (s[e(176)] = i),
                    this[e(2473)].push(s),
                    this[e(1430)][e(1621)](s)
                }
              }
              [t(2419)]() {
                const e = t
                if (this[e(2473)] && !(this._contentDoms[e(277)] < 1)) {
                  for (let t = 0; t < this[e(2473)].length; t++) {
                    const i = this[e(2473)][t]
                    this[e(1430)][e(2269)](i)
                  }
                  this[e(2473)] = []
                }
              }
              [t(2231)]() {
                const e = t
                1 == this._style[e(1640)]
                  ? ((this[e(343)][e(1491)] = Cesium[e(2324)][e(2064)]),
                    (this[e(343)].offset = [0, 15]),
                    this[e(2426)].classList[e(1861)](e(1156)))
                  : ((this[e(343)].horizontalOrigin = Cesium[e(2324)][e(2451)]),
                    (this[e(343)].offset = [16, 48]),
                    this[e(2426)].classList[e(1861)]('xt3d-divgraphic-multiline-img_0'),
                    (this[e(2426)][e(1679)].backgroundImage =
                      e(1748) + this[e(343)][e(311)] + ')')),
                  this[e(2473)].forEach((t) => {
                    const i = e
                    ;(t[i(1679)][i(1070)] = this[i(343)][i(1070)]),
                      (t.style[i(735)] = this[i(343)][i(735)] + 'px')
                  }),
                  this[e(2294)]()
              }
              _getProperties() {
                const t = {}
                return (t.content = this.content), t
              }
            })(e)
          case Q.rectangleLabel:
            return new (class extends Pe {
              constructor(e) {
                const i = t
                if (
                  (super(e),
                  (this[i(528)] = i(611)),
                  (this[i(2336)] = Q[i(665)]),
                  (this._content = Cesium.defaultValue(e[i(1922)], '')),
                  (this[i(343)].color = Cesium[i(1960)](this[i(343)][i(1070)], i(287))),
                  (this[i(343)][i(2274)] = [0, 15]),
                  (this[i(343)].fontSize = Cesium.defaultValue(this[i(343)][i(735)], 14)),
                  (this._style[i(821)] = Cesium[i(1960)](this[i(343)].backgroundColor, [
                    i(1302),
                    i(1450)
                  ])),
                  !Array.isArray(this[i(343)][i(821)]) || 2 != this[i(343)][i(821)][i(277)])
                )
                  throw new Cesium.DeveloperError(i(2074))
                this._createContentDom(), this._setStyle()
              }
              get [t(1922)]() {
                return this[t(992)]
              }
              set [t(1922)](e) {
                const i = t
                ;(this[i(992)] = e), (this[i(1430)][i(176)] = e)
              }
              [t(431)]() {
                const e = t
                let i = document[e(1945)](e(1120))
                i[e(2453)][e(1861)](e(1321)),
                  (i.innerHTML = this[e(992)]),
                  (this._contentDom = i),
                  this[e(1736)].appendChild(i)
              }
              [t(2231)]() {
                const e = t,
                  i = this._style[e(1070)],
                  s = this[e(343)][e(735)],
                  n = this[e(1679)][e(821)]
                ;(this[e(1430)][e(1679)].cssText =
                  'color:' + i + e(571) + s + e(1917) + n[0] + ',' + n[1] + ');'),
                  this[e(2294)]()
              }
              [t(1114)]() {
                const e = {}
                return (e[t(1922)] = this.content), e
              }
            })(e)
          case Q[i(1093)]:
            return new (class extends Pe {
              constructor(e) {
                const i = t
                super(e),
                  (this[i(528)] = i(463)),
                  (this[i(343)][i(2274)] = [0, 10]),
                  (this._graphicType = Q.uprightLabel),
                  (this[i(992)] = Cesium[i(1960)](e[i(1922)], '')),
                  (this[i(343)][i(1070)] = Cesium[i(1960)](this[i(343)][i(1070)], i(741))),
                  (this[i(343)][i(735)] = Cesium[i(1960)](this[i(343)][i(735)], 16)),
                  (this[i(343)].lineHeight = Cesium[i(1960)](this[i(343)][i(349)], 100)),
                  this[i(431)](),
                  this[i(2231)]()
              }
              get [t(1922)]() {
                return this._content
              }
              set content(e) {
                const i = t
                ;(this._content = e), (this[i(1430)][i(176)] = e)
              }
              [t(431)]() {
                const e = t
                let i = document[e(1945)](e(1120))
                i.classList[e(1861)](e(802))
                let s = document[e(1945)](e(1120))
                ;(s[e(176)] = this[e(992)]),
                  s[e(2453)].add(e(906)),
                  i.appendChild(s),
                  (this[e(1430)] = s)
                let n = document[e(1945)]('div')
                n[e(2453)].add(e(2475)), i[e(1621)](n), (this[e(1340)] = n)
                let o = document[e(1945)](e(1120))
                o[e(2453)].add(e(815)), (this._circle = o), i[e(1621)](o), this[e(1736)][e(1621)](i)
              }
              [t(2231)]() {
                const e = t,
                  i = this[e(343)].color
                ;(this[e(1430)][e(1679)].color = i),
                  (this[e(1430)][e(1679)].fontSize = this[e(343)][e(735)] + 'px'),
                  (this[e(1340)].style[e(821)] = i),
                  (this[e(1340)][e(1679)][e(2306)] = this[e(343)][e(349)] + 'px'),
                  (this[e(1162)].style[e(821)] = i),
                  this[e(2294)]()
              }
              _getProperties() {
                const e = t,
                  i = {}
                return (i[e(1922)] = this[e(1922)]), i
              }
            })(e)
          case Q[i(807)]:
            return new Me(e)
          case Q[i(1193)]:
            return new (class extends Me {
              constructor(e) {
                const i = t
                super(e),
                  (this._typeName = '文字模板'),
                  (this[i(343)][i(2274)] = [0, 10]),
                  (this[i(2336)] = Q[i(1193)]),
                  this._dom.classList[i(1861)](this[i(343)][i(1304)])
              }
            })(e)
          case Q[i(1384)]:
            return new (class extends Pe {
              constructor(e = {}) {
                const i = t
                super(e),
                  (this[i(528)] = i(394)),
                  (this._graphicType = Q[i(1384)]),
                  (this._content = Cesium[i(1960)](e[i(1922)], '')),
                  (this[i(343)][i(2274)] = [-8, 20]),
                  (this[i(343)][i(1306)] = Cesium[i(1960)](this[i(343)][i(1306)], 50)),
                  (this[i(343)][i(1466)] = Cesium[i(1960)](this[i(343)].moveDomTop, -100)),
                  (this[i(343)].autoPoistion = Cesium[i(1960)](this[i(343)][i(442)], !0)),
                  (this[i(343)].horizontalPoistion = Cesium[i(1960)](
                    this[i(343)][i(1317)],
                    i(1709)
                  )),
                  (this[i(343)].verticalPoistion = Cesium[i(1960)](this[i(343)][i(1317)], i(151))),
                  (this[i(343)][i(1018)] = Cesium[i(1960)](this._style[i(1018)], i(741))),
                  this._createContentDom()
              }
              set [t(1922)](e) {
                const i = t
                ;(this[i(992)] = e), (this[i(1822)][i(176)] = this[i(992)])
              }
              get [t(1922)]() {
                return this[t(992)]
              }
              [t(431)]() {
                const e = t
                ;(this[e(1822)] = document[e(1945)](e(1120))),
                  this[e(1822)][e(2453)][e(1861)](e(319)),
                  (this[e(1822)][e(176)] = this[e(992)]),
                  (this[e(462)] = document[e(1945)]('div')),
                  this[e(462)][e(2453)][e(1861)]('divIndicator-fixed'),
                  (this[e(1209)] = document[e(1945)](e(1120))),
                  this[e(1209)][e(2453)][e(1861)](e(806)),
                  this[e(1736)][e(1621)](this._container_drag),
                  this[e(1736)][e(1621)](this[e(1209)]),
                  this._dom[e(1621)](this._container_fixed),
                  this[e(1822)][e(1973)](e(2464), this[e(512)][e(1505)](this)),
                  setTimeout(() => {
                    this._updateLineStyle()
                  }, 500)
              }
              _setStyle() {
                const e = t
                this[e(790)](), this[e(2294)]()
              }
              [t(512)](e) {
                const i = t
                var s, n
                e[i(1511)](), e[i(1031)]()
                const o = e[i(1260)] - (null == (s = this[i(1822)]) ? void 0 : s[i(805)]),
                  r = e.clientY - (null == (n = this[i(1822)]) ? void 0 : n[i(1129)])
                Ae(document.documentElement, i(2263), h),
                  Ae(document[i(486)], i(1298), function t(e) {
                    const s = i
                    e[s(1511)](),
                      e[s(1031)](),
                      Te(document.documentElement, s(2263), h),
                      Te(document[s(486)], s(1298), t),
                      Te(a[s(1392)], s(2263), h)
                  }),
                  Ae(this[i(1392)], 'mousemove', h)
                const a = this
                function h(t) {
                  const e = i
                  t.preventDefault(),
                    t[e(1031)](),
                    (a[e(1679)][e(1306)] = t[e(1260)] - o),
                    (a[e(1679)].moveDomTop = t[e(1958)] - r),
                    a[e(790)]()
                }
              }
              [t(790)]() {
                const e = t
                ;(this._container_drag[e(1679)][e(1709)] = this[e(343)][e(1306)] + 'px'),
                  (this._container_drag[e(1679)][e(2223)] = this[e(343)][e(1466)] + 'px')
                const i = this[e(1822)][e(2472)](),
                  s = this._container_fixed.getBoundingClientRect()
                let n, o
                ;(this[e(1209)][e(1679)][e(1811)] = e(1858) + this[e(343)][e(1018)]),
                  this[e(343)][e(442)]
                    ? ((o = i[e(1709)] < s[e(1709)] ? e(1883) : e(1709)),
                      (n = i[e(2223)] < s[e(2223)] ? e(151) : e(2223)))
                    : ((n = this[e(343)][e(781)]), (o = this._style[e(1317)]))
                const r = s.y + s.height / 2,
                  a = s.x + s.width / 2,
                  h = i[n],
                  l = i[o],
                  c = (function (t, i) {
                    const s = e,
                      n = Math[s(636)](t.x - i.x),
                      o = Math.abs(t.y - i.y)
                    return Number(
                      Math[s(1345)](Math[s(2394)](n, 2) + Math[s(2394)](o, 2)).toFixed(2)
                    )
                  })({ x: r, y: a }, { x: h, y: l }),
                  u = (l - a) / 2 - 1,
                  m = (h - r) / 2 - c / 2,
                  p = -Math[e(2003)](a - l, r - h) * (180 / Math.PI),
                  d = {}
                ;(d[e(2306)] = c + 'px'),
                  (d.transform =
                    'translateX(' +
                    u +
                    'px) translateY(' +
                    m +
                    'px) scale(1) rotate(' +
                    p +
                    e(748)),
                  Object[e(250)](this._container_line[e(1679)], d)
              }
            })(e)
        }
      }
    }
    class ze extends me {
      constructor(t) {
        super(t), (this._graphicClassType = I)
      }
      [t(1200)]() {}
      [t(2515)](e = []) {
        const i = t
        ;(this[i(2469)] = e), (this[i(794)] = this[i(2469)])
        let s = []
        e[i(1602)]((t) => {
          s.push([t[0], t[1]])
        })
        const n = this[i(343)]
        !n[i(2080)] && n[i(398)] > 0 && e[i(277)] > 0 && this[i(2309)](e),
          (this._cartesian3Array = this[i(1264)](s)),
          this._cartesian3Array &&
            (this[i(1378)].push(this[i(1378)][0]),
            (this[i(1255)] = new Cesium[i(2229)](this[i(1378)])))
      }
      _generatePositions(t) {}
      _convertPositions(e) {
        const i = t
        e[i(1602)]((t) => {
          const e = i
          t[2] = this[e(1516)]
        })
        let s = Cesium[i(310)][i(774)]([][i(1500)][i(560)]([], e))
        return (this[i(1024)] = Cesium.BoundingSphere.fromPoints(s)), s
      }
    }
    const De = (e, i) =>
        Math[t(1345)](Math[t(2394)](e[0] - i[0], 2) + Math[t(2394)](e[1] - i[1], 2)),
      Ie = (e) => {
        const i = t
        let s = 0
        return (
          e &&
            Array.isArray(e) &&
            e[i(277)] > 0 &&
            e[i(1602)]((t, i) => {
              i < e.length - 1 && (s += De(t, e[i + 1]))
            }),
          s
        )
      },
      ke = (t) => Math.pow(Ie(t), 0.99),
      Fe = (t, e) => [(t[0] + e[0]) / 2, (t[1] + e[1]) / 2],
      Re = (e, i) => {
        const s = t
        let n,
          o = Math[s(2463)](Math[s(636)](i[1] - e[1]) / De(e, i))
        return (
          i[1] >= e[1] && i[0] >= e[0]
            ? (n = o + Math.PI)
            : i[1] >= e[1] && i[0] < e[0]
              ? (n = 2 * Math.PI - o)
              : i[1] < e[1] && i[0] < e[0]
                ? (n = o)
                : i[1] < e[1] && i[0] >= e[0] && (n = Math.PI - o),
          n
        )
      },
      Le = (t, e, i) => {
        let s = Re(e, t) - Re(e, i)
        return s < 0 ? s + 2 * Math.PI : s
      },
      Oe = (t, e, i) => (i[1] - t[1]) * (e[0] - t[0]) > (e[1] - t[1]) * (i[0] - t[0]),
      Be = (e, i, s, n, o) => {
        const r = t
        e = Math[r(1113)](Math[r(1429)](e, 1), 0)
        let [a, h] = [1 - e, e * e],
          l = h * e,
          c = a * a,
          u = c * a
        return [
          u * i[0] + 3 * c * e * s[0] + 3 * a * h * n[0] + l * o[0],
          u * i[1] + 3 * c * e * s[1] + 3 * a * h * n[1] + l * o[1]
        ]
      },
      Ve = (e, i, s, n, o) => {
        const r = t
        let a = Re(e, i),
          h = o ? a + s : a - s,
          l = n * Math[r(1272)](h),
          c = n * Math.sin(h)
        return [i[0] + l, i[1] + c]
      },
      Ne = (e, i, s, o) => {
        const r = t
        let a = He(i, s, o),
          [h, l, c, u, m] = [null, null, null, null, null],
          p = Math[r(1345)](a[0] * a[0] + a[1] * a[1]),
          d = a[0] / p,
          f = a[1] / p,
          C = De(i, s),
          v = De(s, o)
        return (
          p > n
            ? Oe(i, s, o)
              ? ((c = e * C),
                (h = [(u = s[0] - c * f), (m = s[1] + c * d)]),
                (c = e * v),
                (l = [(u = s[0] + c * f), (m = s[1] - c * d)]))
              : ((c = e * C),
                (h = [(u = s[0] + c * f), (m = s[1] - c * d)]),
                (c = e * v),
                (l = [(u = s[0] - c * f), (m = s[1] + c * d)]))
            : ((h = [(u = s[0] + e * (i[0] - s[0])), (m = s[1] + e * (i[1] - s[1]))]),
              (l = [(u = s[0] + e * (o[0] - s[0])), (m = s[1] + e * (o[1] - s[1]))])),
          [h, l]
        )
      },
      He = (e, i, s) => {
        const n = t
        let o = e[0] - i[0],
          r = e[1] - i[1],
          a = Math.sqrt(o * o + r * r)
        ;(o /= a), (r /= a)
        let h = s[0] - i[0],
          l = s[1] - i[1],
          c = Math[n(1345)](h * h + l * l)
        return [o + (h /= c), r + (l /= c)]
      },
      Ge = function (e) {
        const i = t
        if (e[i(277)] <= 2) return e
        {
          let t = [],
            s = e[i(277)] - 1
          for (let n = 0; n <= 1; n += 0.01) {
            let [o, r] = [0, 0]
            for (let t = 0; t <= s; t++) {
              let a = Ue(s, t),
                h = Math[i(2394)](n, t),
                l = Math[i(2394)](1 - n, s - t)
              ;(o += a * h * l * e[t][0]), (r += a * h * l * e[t][1])
            }
            t[i(2553)]([o, r])
          }
          return t[i(2553)](e[s]), t
        }
      },
      We = (t) => {
        let e = 1
        switch (t) {
          case t <= 1:
            e = 1
            break
          case 2 === t:
            e = 2
            break
          case 3 === t:
            e = 6
            break
          case 24 === t:
            e = 24
            break
          case 5 === t:
            e = 120
            break
          default:
            for (let i = 1; i <= t; i++) e *= i
        }
        return e
      },
      Ue = (t, e) => We(t) / (We(e) * We(t - e)),
      je = (e) => {
        const i = t
        if (e.length <= 2) return e
        {
          let [t, s] = [2, []],
            n = e[i(277)] - t - 1
          s[i(2553)](e[0])
          for (let o = 0; o <= n; o++)
            for (let n = 0; n <= 1; n += 0.05) {
              let [r, a] = [0, 0]
              for (let i = 0; i <= t; i++) {
                let t = qe(i, n)
                ;(r += t * e[o + i][0]), (a += t * e[o + i][1])
              }
              s[i(2553)]([r, a])
            }
          return s[i(2553)](e[e.length - 1]), s
        }
      },
      qe = (e, i) => {
        const s = t
        let n = 0
        return (
          0 === e
            ? (n = Math.pow(i - 1, 2) / 2)
            : 1 === e
              ? (n = (-2 * Math[s(2394)](i, 2) + 2 * i + 1) / 2)
              : 2 === e && (n = Math[s(2394)](i, 2) / 2),
          n
        )
      }
    class Ye extends ze {
      constructor(e = {}) {
        const i = t
        super(e), (this[i(2336)] = Q[i(941)]), (this[i(528)] = i(1812))
      }
      _initConsts() {
        const e = t
        ;(this[e(133)] = 0.18),
          (this[e(2358)] = 0.3),
          (this[e(128)] = 0.85),
          (this[e(2343)] = 0.15),
          (this[e(1010)] = 0.8),
          (this[e(324)] = 3)
      }
      [t(1264)](e) {
        const i = t
        let s = e[i(277)]
        if (s < 2) return
        if (2 == s) return this._convertPositions(e)
        let n = e,
          o = n[0],
          r = n[1]
        Oe(n[0], n[1], n[2]) && ((o = n[1]), (r = n[0]))
        let a = [Fe(o, r)][i(1500)](n.slice(2)),
          h = this[i(2321)](a, o, r),
          l = h[0],
          c = h[4],
          u = De(o, r) / ke(a),
          m = this[i(1649)](a, l, c, u)
        s = m[i(277)]
        let p = [o][i(1500)](m[i(123)](0, s / 2))
        p[i(2553)](l)
        let d = [r][i(1500)](m[i(123)](s / 2, s))
        return d[i(2553)](c), (p = je(p)), (d = je(d)), this[i(1690)](p[i(1500)](h, d.reverse()))
      }
      _getArrowHeadPoints(e, i, n) {
        const o = t
        let r = ke(e),
          a = r * this[o(133)],
          h = e[e[o(277)] - 1]
        r = De(h, e[e[o(277)] - 2])
        let l = De(i, n)
        a > l * this[o(1010)] && (a = l * this[o(1010)])
        let c = a * this[o(2358)],
          u = a * this._neckWidthFactor,
          m = (a = a > r ? r : a) * this[o(128)],
          p = Ve(e[e[o(277)] - 2], h, 0, a, !0),
          d = Ve(e[e[o(277)] - 2], h, 0, m, !0),
          f = Ve(h, p, s, c, !1),
          C = Ve(h, p, s, c, !0)
        return [Ve(h, d, s, u, !1), f, h, C, Ve(h, d, s, u, !0)]
      }
      [t(1649)](e, i, s, n) {
        const o = t
        let r = Ie(e),
          a = ke(e) * n,
          h = (a - De(i, s)) / 2,
          l = 0,
          c = [],
          u = []
        for (let t = 1; t < e.length - 1; t++) {
          let i = Le(e[t - 1], e[t], e[t + 1]) / 2,
            s = (a / 2 - ((l += De(e[t - 1], e[t])) / r) * h) / Math[o(884)](i),
            n = Ve(e[t - 1], e[t], Math.PI - i, s, !0),
            m = Ve(e[t - 1], e[t], i, s, !1)
          c[o(2553)](n), u[o(2553)](m)
        }
        return c[o(1500)](u)
      }
    }
    class Xe extends ze {
      constructor(e = {}) {
        const i = t
        super(e), (this[i(2336)] = Q.fineArrow), (this[i(528)] = i(1210))
      }
      [t(1200)]() {
        const e = t
        ;(this._tailWidthFactor = 0.1),
          (this[e(2343)] = 0.2),
          (this[e(2358)] = 0.25),
          (this._headAngle = Math.PI / 8.5),
          (this._neckAngle = Math.PI / 13),
          (this[e(350)] = 2)
      }
      [t(1264)](e) {
        const i = t
        if (e[i(277)] < 2) return
        let n = e,
          o = n[0],
          r = n[1],
          a = ke(n),
          h = a * this[i(1243)],
          l = a * this[i(2343)],
          c = a * this[i(2358)],
          u = Ve(r, o, s, h, !0),
          m = Ve(r, o, s, h, !1),
          p = Ve(o, r, this[i(148)], c, !1),
          d = Ve(o, r, this._headAngle, c, !0),
          f = [u, Ve(o, r, this[i(1229)], l, !1), p, r, d, Ve(o, r, this[i(1229)], l, !0), m]
        return this[i(1690)](f)
      }
    }
    class Qe extends Ye {
      constructor(e = {}) {
        const i = t
        super(e), (this._graphicType = Q[i(1668)]), (this[i(528)] = i(1741))
      }
      [t(1200)]() {
        const e = t
        ;(this[e(133)] = 0.18),
          (this[e(2358)] = 0.3),
          (this[e(128)] = 0.85),
          (this[e(2343)] = 0.15),
          (this._tailWidthFactor = 0.1),
          (this[e(619)] = 1),
          (this._swallowTailPnt = null),
          (this[e(324)] = 2)
      }
      [t(1264)](e) {
        const i = t
        let s = e[i(277)]
        if (!(s < 2)) {
          var n = e,
            o = this[i(752)](n),
            r = this[i(2321)](n, o[0], o[2]),
            a = r[0],
            h = r[4],
            l = this[i(1649)](n, a, h, this[i(1243)])
          s = l[i(277)]
          var c = [o[0]][i(1500)](l[i(123)](0, s / 2))
          c.push(a)
          var u = [o[2]][i(1500)](l[i(123)](s / 2, s))
          return (
            u[i(2553)](h),
            (c = je(c)),
            (u = je(u)),
            this[i(1690)](c[i(1500)](r, u[i(1027)](), [o[1], c[0]]))
          )
        }
      }
      _getTailPoints(e) {
        const i = t
        var n = ke(e) * this[i(1243)],
          o = Ve(e[1], e[0], s, n, !1),
          r = Ve(e[1], e[0], s, n, !0),
          a = n * this[i(619)]
        return [o, Ve(e[1], e[0], 0, a, !0), r]
      }
    }
    let Ze = {
      readFeature(e) {
        const i = t
        let s = e.geometry.coordinates,
          n = e[i(1004)]
        return (n.positions = s), this[i(1951)](n)
      },
      create(e) {
        const n = t
        switch (e[n(1158)]) {
          case Q.doubleArrow:
            return new (class extends ze {
              constructor(e = {}) {
                const i = t
                super(e), (this[i(2336)] = Q[i(682)]), (this[i(528)] = '钳击箭头')
              }
              [t(1200)]() {
                const e = t
                ;(this[e(133)] = 0.25),
                  (this[e(2358)] = 0.3),
                  (this[e(128)] = 0.85),
                  (this[e(2343)] = 0.15),
                  (this[e(2029)] = null),
                  (this._tempPoint4 = null),
                  (this._fixPointCount = 5)
              }
              [t(1264)](e) {
                const i = t
                let s = e[i(277)]
                if (s < 2) return
                if (2 == s) return this[i(1690)](e)
                let n,
                  o,
                  r = e[0],
                  a = e[1],
                  h = e[2]
                ;(this[i(413)] = 3 == s ? this[i(227)](r, a, h) : e[3]),
                  (this[i(2029)] = 3 == s || 4 == s ? Fe(r, a) : e[4]),
                  Oe(r, a, h)
                    ? ((n = this[i(2027)](r, this[i(2029)], this._tempPoint4, !1)),
                      (o = this[i(2027)](this[i(2029)], a, h, !0)))
                    : ((n = this[i(2027)](a, this[i(2029)], h, !1)),
                      (o = this[i(2027)](this[i(2029)], r, this[i(413)], !0)))
                let l = n[i(277)],
                  c = (l - 5) / 2,
                  u = n[i(123)](0, c),
                  m = n[i(123)](c, c + 5),
                  p = n[i(123)](c + 5, l),
                  d = o.slice(0, c),
                  f = o.slice(c, c + 5),
                  C = o[i(123)](c + 5, l)
                d = Ge(d)
                let v = Ge(C[i(1500)](u[i(123)](1)))
                p = Ge(p)
                let _ = d[i(1500)](f, v, m, p)
                return this._convertPositions(_)
              }
              [t(2027)](e, i, n, o) {
                const r = t
                let a = Fe(e, i),
                  h = De(a, n),
                  l = Ve(n, a, 0, 0.3 * h, !0),
                  c = Ve(n, a, 0, 0.5 * h, !0),
                  u = [a, (l = Ve(a, l, s, h / 5, o)), (c = Ve(a, c, s, h / 4, o)), n],
                  m = this[r(2321)](u, this[r(133)], this[r(2358)], this[r(128)], this[r(2343)]),
                  p = m[0],
                  d = m[4],
                  f = De(e, i) / ke(u) / 2,
                  C = this._getArrowBodyPoints(u, p, d, f),
                  v = C[r(277)],
                  _ = C[r(123)](0, v / 2),
                  g = C[r(123)](v / 2, v)
                return (
                  _.push(p),
                  g[r(2553)](d),
                  (_ = _[r(1027)]())[r(2553)](i),
                  (g = g[r(1027)]())[r(2553)](e),
                  _.reverse().concat(m, g)
                )
              }
              _getArrowHeadPoints(e, i, n) {
                const o = t
                let r = ke(e) * this[o(133)],
                  a = e[e.length - 1]
                De(i, n)
                let h = r * this[o(2358)],
                  l = r * this[o(2343)],
                  c = r * this[o(128)],
                  u = Ve(e[e[o(277)] - 2], a, 0, r, !0),
                  m = Ve(e[e[o(277)] - 2], a, 0, c, !0),
                  p = Ve(a, u, s, h, !1),
                  d = Ve(a, u, s, h, !0)
                return [Ve(a, m, s, l, !1), p, a, d, Ve(a, m, s, l, !0)]
              }
              [t(1649)](e, i, s, n) {
                const o = t
                let r = Ie(e),
                  a = ke(e) * n,
                  h = (a - De(i, s)) / 2,
                  l = 0,
                  c = [],
                  u = []
                for (let t = 1; t < e[o(277)] - 1; t++) {
                  let i = Le(e[t - 1], e[t], e[t + 1]) / 2,
                    s = (a / 2 - ((l += De(e[t - 1], e[t])) / r) * h) / Math[o(884)](i),
                    n = Ve(e[t - 1], e[t], Math.PI - i, s, !0),
                    m = Ve(e[t - 1], e[t], i, s, !1)
                  c[o(2553)](n), u.push(m)
                }
                return c[o(1500)](u)
              }
              [t(227)](e, i, n) {
                const o = t
                let r,
                  a,
                  h,
                  l,
                  c = Fe(e, i),
                  u = De(c, n),
                  m = Le(e, c, n)
                return (
                  m < s
                    ? ((a = u * Math[o(884)](m)),
                      (h = u * Math[o(1272)](m)),
                      (l = Ve(e, c, s, a, !1)),
                      (r = Ve(c, l, s, h, !0)))
                    : m >= s && m < Math.PI
                      ? ((a = u * Math[o(884)](Math.PI - m)),
                        (h = u * Math[o(1272)](Math.PI - m)),
                        (l = Ve(e, c, s, a, !1)),
                        (r = Ve(c, l, s, h, !1)))
                      : m >= Math.PI && m < 1.5 * Math.PI
                        ? ((a = u * Math.sin(m - Math.PI)),
                          (h = u * Math.cos(m - Math.PI)),
                          (l = Ve(e, c, s, a, !0)),
                          (r = Ve(c, l, s, h, !0)))
                        : ((a = u * Math[o(884)](2 * Math.PI - m)),
                          (h = u * Math[o(1272)](2 * Math.PI - m)),
                          (l = Ve(e, c, s, a, !0)),
                          (r = Ve(c, l, s, h, !1))),
                  r
                )
              }
            })(e)
          case Q[n(941)]:
            return new Ye(e)
          case Q[n(2286)]:
            return new (class extends ze {
              constructor(e = {}) {
                const i = t
                super(e), (this[i(2336)] = Q[i(2286)]), (this[i(528)] = i(1999))
              }
              _initConsts() {
                ;(this[t(324)] = 3), (this.t = 0.3)
              }
              [t(1264)](e) {
                const s = t
                let n = e.length
                if (n < 2) return
                if (2 == n) return this[s(1690)](e)
                let o = e
                o[s(2553)](o[0], o[1])
                let r = []
                for (let t = 0; t < o[s(277)] - 2; t++) {
                  let e = Ne(this.t, o[t], o[t + 1], o[t + 2])
                  r = r[s(1500)](e)
                }
                r = [r[(n = r[s(277)]) - 1]].concat(r[s(123)](0, n - 1))
                let a = []
                for (let t = 0; t < o[s(277)] - 2; t++) {
                  let e = o[t],
                    n = o[t + 1]
                  a[s(2553)](e)
                  for (let o = 0; o <= i; o++) {
                    let h = Be(o / i, e, r[2 * t], r[2 * t + 1], n)
                    a[s(2553)](h)
                  }
                  a[s(2553)](n)
                }
                return this[s(1690)](a)
              }
            })(e)
          case Q[n(1181)]:
            return new Xe(e)
          case Q[n(2568)]:
            return new (class extends Xe {
              constructor(e = {}) {
                const i = t
                super(e), (this[i(2336)] = Q[i(2568)]), (this._typeName = i(1307))
              }
              [t(1200)]() {
                const e = t
                ;(this._tailWidthFactor = 0.2),
                  (this[e(2343)] = 0.25),
                  (this[e(2358)] = 0.3),
                  (this[e(148)] = Math.PI / 4),
                  (this[e(1229)] = 0.17741 * Math.PI),
                  (this[e(350)] = 2)
              }
            })(e)
          case Q.gatheringPlace:
            return new (class extends ze {
              constructor(e = {}) {
                const i = t
                super(e), (this[i(2336)] = Q.gatheringPlace), (this[i(528)] = i(2150))
              }
              [t(1200)]() {
                const e = t
                ;(this._t = 0.4), (this[e(350)] = 3)
              }
              [t(1264)](e) {
                const n = t
                let o = e.length,
                  r = e
                if (o < 2) return
                if (2 == o) {
                  let t = Fe(r[0], r[1]),
                    e = De(r[0], t) / 0.9,
                    i = Ve(r[0], t, s, e, !0)
                  r = [r[0], i, r[1]]
                }
                let a = Fe(r[0], r[2])
                r[n(2553)](a, r[0], r[1])
                let h,
                  l,
                  c,
                  u = []
                for (let t = 0; t < r[n(277)] - 2; t++) {
                  ;(h = r[t]), (l = r[t + 1]), (c = r[t + 2])
                  let e = Ne(this._t, h, l, c)
                  u = u[n(1500)](e)
                }
                u = [u[(o = u.length) - 1]].concat(u[n(123)](0, o - 1))
                let m = []
                for (let t = 0; t < r[n(277)] - 2; t++) {
                  ;(h = r[t]), (l = r[t + 1]), m.push(h)
                  for (let e = 0; e <= i; e++) {
                    let s = Be(e / i, h, u[2 * t], u[2 * t + 1], l)
                    m[n(2553)](s)
                  }
                  m.push(l)
                }
                return this[n(1690)](m)
              }
            })(e)
          case Q[n(2197)]:
            return new (class extends Ye {
              constructor(e = {}) {
                const i = t
                super(e), (this._graphicType = Q[i(2197)]), (this[i(528)] = i(2207))
              }
              _initConsts() {
                const e = t
                ;(this._headHeightFactor = 0.18),
                  (this[e(2358)] = 0.3),
                  (this._neckHeightFactor = 0.85),
                  (this[e(2343)] = 0.15),
                  (this._headTailFactor = 0.8),
                  (this._tailWidthFactor = 0.1),
                  (this[e(324)] = 2)
              }
              _generatePositions(e) {
                const i = t
                let s = e[i(277)]
                if (s < 2) return
                let n = e,
                  o = this[i(752)](n),
                  r = this[i(2321)](n, o[0], o[1]),
                  a = r[0],
                  h = r[4],
                  l = this[i(1649)](n, a, h, this._tailWidthFactor)
                s = l[i(277)]
                let c = [o[0]][i(1500)](l.slice(0, s / 2))
                c[i(2553)](a)
                let u = [o[1]][i(1500)](l[i(123)](s / 2, s))
                return (
                  u[i(2553)](h),
                  (c = je(c)),
                  (u = je(u)),
                  this._convertPositions(c.concat(r, u[i(1027)]()))
                )
              }
              [t(752)](e) {
                const i = t
                let n = ke(e) * this[i(1243)]
                return [Ve(e[1], e[0], s, n, !1), Ve(e[1], e[0], s, n, !0)]
              }
            })(e)
          case Q[n(1668)]:
            return new Qe(e)
          case Q[n(1394)]:
            return new (class extends Ye {
              constructor(e = {}) {
                const i = t
                super(e), (this[i(2336)] = Q.tailedAttackArrow), (this[i(528)] = i(1368))
              }
              [t(1200)]() {
                const e = t
                ;(this[e(133)] = 0.18),
                  (this._headWidthFactor = 0.3),
                  (this._neckHeightFactor = 0.85),
                  (this[e(2343)] = 0.15),
                  (this[e(1010)] = 0.8),
                  (this[e(1243)] = 0.1),
                  (this[e(619)] = 1),
                  (this[e(1953)] = null),
                  (this[e(324)] = 2)
              }
              _generatePositions(e) {
                const i = t
                let s = e[i(277)]
                if (s < 2) return
                if (2 == s) return this[i(1690)](e)
                let n = e,
                  o = n[0],
                  r = n[1]
                Oe(n[0], n[1], n[2]) && ((o = n[1]), (r = n[0]))
                let a = [Fe(o, r)][i(1500)](n.slice(2)),
                  h = this[i(2321)](a, o, r),
                  l = h[0],
                  c = h[4],
                  u = De(o, r),
                  m = ke(a),
                  p = m * this._tailWidthFactor * this[i(619)]
                this[i(1953)] = Ve(a[1], a[0], 0, p, !0)
                let d = u / m,
                  f = this._getArrowBodyPoints(a, l, c, d)
                s = f.length
                let C = [o][i(1500)](f.slice(0, s / 2))
                C.push(l)
                let v = [r][i(1500)](f.slice(s / 2, s))
                return (
                  v[i(2553)](c),
                  (C = je(C)),
                  (v = je(v)),
                  this[i(1690)](C[i(1500)](h, v[i(1027)](), [this[i(1953)], C[0]]))
                )
              }
            })(e)
          case Q.tailedFineArrow:
            return new (class extends Qe {
              constructor(e = {}) {
                const i = t
                super(e), (this._graphicType = Q.tailedFineArrow), (this[i(528)] = i(2284))
              }
              [t(1200)]() {
                const e = t
                ;(this[e(133)] = 0.18),
                  (this[e(2358)] = 0.3),
                  (this[e(128)] = 0.85),
                  (this[e(2343)] = 0.15),
                  (this[e(1010)] = 0.8),
                  (this[e(1243)] = 0.1),
                  (this[e(619)] = 1),
                  (this[e(1953)] = null),
                  (this[e(350)] = 2)
              }
            })(e)
        }
      }
    }
    class Ke extends b {
      constructor(e) {
        const i = t
        super(e),
          (this._graphicClassType = O),
          (this._typeName = '光源'),
          (this[i(350)] = 1),
          (this[i(314)] = i(1865)),
          (this[i(987)] = this[i(2547)]()),
          (this[i(343)][i(1981)] = Cesium[i(1960)](this[i(343)].radius, 1e3)),
          (this[i(343)][i(1070)] = Cesium.defaultValue(this[i(343)][i(1070)], 'rgba(255,0,0,1)')),
          (this[i(343)][i(301)] = Cesium[i(1154)][i(2008)](this[i(343)][i(1070)])),
          (this._style[i(2480)] = Cesium[i(1960)](this[i(343)][i(2480)], 1e3)),
          (this[i(343)][i(1968)] = Cesium[i(1960)](this[i(343)].showOrigin, !0)),
          (this[i(1912)] = Cesium[i(1960)](e.position, [111, 28, 0])),
          this[i(1253)](this._position),
          this[i(2231)]()
      }
      [t(415)](e) {
        const i = t
        this._entity[i(1482)] = e
      }
      [t(2547)]() {
        const e = t
        return new Cesium[e(1380)]({
          graphicId: this[e(1570)],
          point: { pixelSize: 10, color: Cesium.Color[e(588)] }
        })
      }
      [t(1253)](e) {
        const i = t
        ;(this[i(1912)] = e),
          (this[i(2238)] = Cesium[i(310)].fromDegrees(e[0], e[1], e[2])),
          this[i(987)] && (this[i(987)][i(2251)] = this._cartesian3),
          this[i(1365)] && (this[i(1365)].position = this[i(2238)]),
          (this[i(794)] = this[i(1912)]),
          (this[i(1024)] = new Cesium[i(1242)](this._cartesian3, 3))
      }
      createShadowMap(t) {}
      [t(304)](t) {}
      [t(139)]() {
        const e = t
        return new Cesium[e(2506)]({
          sampleMode: Cesium[e(1051)][e(766)],
          fragmentShader: this[e(1764)](),
          uniforms: this[e(266)]()
        })
      }
      [t(1764)]() {}
      getUniform() {}
      [t(2231)]() {
        const e = t
        ;(this[e(343)].lightColor = Cesium[e(1154)][e(2008)](this._style[e(1070)])),
          this._shadowMap && (this._shadowMap[e(426)] = this[e(343)][e(1981)]),
          (this[e(987)].show = this._show && this._style[e(1968)]),
          this[e(1737)](),
          (this[e(1024)] = new Cesium[e(1242)](this[e(2238)], 3))
      }
      [t(1737)]() {}
      [t(415)](e) {
        const i = t
        ;(this[i(987)][i(1482)] = this[i(1144)] && this[i(343)][i(1968)]),
          this._postProcessStage && (this[i(2364)][i(299)] = e)
      }
      _init(e) {
        const i = t
        ;(this[i(1365)] = this[i(304)](e)),
          (this[i(2444)] = this[i(637)](e)),
          (this[i(2364)] = this[i(139)]())
      }
      [t(1545)](e) {
        const i = t
        ;(this._entity[i(2122)] = e.id),
          (this._layer = e),
          this[i(144)](e._viewer),
          e[i(245)].entities[i(1861)](this._entity),
          e[i(245)][i(696)][i(1346)][i(1861)](this),
          e[i(245)][i(696)][i(2165)][i(1861)](this[i(2364)])
      }
      [t(2389)](e) {
        const i = t
        e[i(245)][i(696)][i(2165)][i(1896)](this[i(2364)]),
          e._viewer[i(696)].primitives[i(1896)](this),
          e[i(245)].entities[i(1896)](this[i(987)])
      }
      [t(720)](e) {
        const i = t
        this[i(2444)] && e.shadowMaps.push(this[i(2444)])
      }
      [t(2039)]() {}
    }
    let Je = {
      readFeature(e) {
        const i = t
        let s = e.geometry[i(2365)],
          n = e[i(1004)]
        return (n[i(2251)] = s), this[i(1951)](n)
      },
      create(e) {
        const i = t
        switch (e[i(1158)]) {
          case Q[i(1421)]:
            return new (class extends Ke {
              constructor(e) {
                const i = t
                super(e), (this[i(528)] = i(956)), (this[i(2336)] = Q[i(1421)])
              }
              [t(637)](e) {
                const i = t,
                  s = {}
                return (
                  (s.lightCamera = this._camera),
                  (s[i(1036)] = !1),
                  (s[i(1390)] = 1),
                  (s[i(640)] = !0),
                  (s[i(364)] = !1),
                  (s.context = e.scene[i(2325)]),
                  (s[i(1386)] = this._style.radius),
                  (s.fromLightSource = !1),
                  new Cesium[i(433)](s)
                )
              }
              [t(304)](e) {
                const i = t
                let s = new Cesium[i(2162)](e.scene)
                const n = this[i(2251)]
                let o = Cesium[i(310)][i(667)](n[0], n[1], n[2])
                return (
                  (s[i(2251)] = o), (s.up = Cesium[i(310)].normalize(o, new Cesium[i(310)]())), s
                )
              }
              getUniform() {
                const e = t
                let i = this[e(2444)][e(2537)],
                  s = Cesium[e(2141)][e(2434)](
                    1 / this._shadowMap[e(1188)].x,
                    1 / this[e(2444)]._textureSize.y,
                    i.depthBias,
                    i[e(1379)]
                  ),
                  n = new Cesium.Cartesian2(this[e(2444)].darkness, 1)
                const o = {}
                return (
                  (o[e(162)] = () => this[e(2444)][e(600)]),
                  (o.intensity = () => this._style[e(2480)]),
                  (o.lightColor = () => this[e(343)][e(301)]),
                  (o[e(2398)] = () => this._shadowMap[e(829)]),
                  (o[e(220)] = () => s),
                  (o.shadowMapDarknessType = () => n),
                  (o.direction = () => Cesium[e(310)].ZERO),
                  (o.outerConeCos = () => 0),
                  (o[e(670)] = () => 0),
                  (o[e(773)] = () => this[e(2444)][e(1247)]),
                  o
                )
              }
              [t(1764)]() {
                return '\n        uniform sampler2D colorTexture;\n        uniform sampler2D depthTexture; \n\n        uniform vec4 lightPositionEC;\n        uniform float intensity;\n        uniform vec3 lightColor;\n        uniform vec3 direction;\n        uniform float outerConeCos;\n        uniform float innerConeCos;\n        uniform mat4 shadowMapMatrix;\n        uniform vec4 shadowMapTexelSizeDepthBiasAndNormalShadingSmooth;\n        uniform vec2 shadowMapDarknessType;\n\n        uniform samplerCube lightShadowMapCube;\n        in vec2 v_textureCoordinates;\n\n        const float M_PI = 3.141592653589793;\n\n        vec3 getEyeCoordinate3FromWindowCoordinate(vec2 fragCoord, float logDepthOrDepth) {\n        vec4 eyeCoordinate = czm_windowToEyeCoordinates(fragCoord, logDepthOrDepth);\n        return eyeCoordinate.xyz / eyeCoordinate.w;\n        }\n\n        vec3 vectorFromOffset(vec4 eyeCoordinate, vec2 positiveOffset) {\n        vec2 glFragCoordXY = v_textureCoordinates.xy * czm_viewport.zw;\n        float upOrRightLogDepth = czm_unpackDepth(texture(depthTexture, (glFragCoordXY + positiveOffset) / czm_viewport.zw));\n        float downOrLeftLogDepth = czm_unpackDepth(texture(depthTexture, (glFragCoordXY - positiveOffset) / czm_viewport.zw));\n\n        bvec2 upOrRightInBounds = lessThan(glFragCoordXY + positiveOffset, czm_viewport.zw);\n        float useUpOrRight = float(upOrRightLogDepth > 0.0 && upOrRightInBounds.x && upOrRightInBounds.y);\n        float useDownOrLeft = float(useUpOrRight == 0.0);\n        vec3 upOrRightEC = getEyeCoordinate3FromWindowCoordinate(glFragCoordXY + positiveOffset, upOrRightLogDepth);\n        vec3 downOrLeftEC = getEyeCoordinate3FromWindowCoordinate(glFragCoordXY - positiveOffset, downOrLeftLogDepth);\n        return (upOrRightEC - (eyeCoordinate.xyz / eyeCoordinate.w)) * useUpOrRight + ((eyeCoordinate.xyz / eyeCoordinate.w) - downOrLeftEC) * useDownOrLeft;\n        }\n\n        float getRangeAttenuation(float range, float d) {\n        if(range <= 0.0) {\n        return 1.0 / pow(d, 2.0);\n        }\n        return max(min(1.0 - pow(d / range, 4.0), 1.0), 0.0) / pow(d, 2.0);\n        } \n\n        vec3 getLightIntensity(vec3 color, float intensity, float type, float range, vec3 pointToLight, vec3 direction, float outerConeCos, float innerConeCos) {\n        float rangeAttenuation = 1.0; \n        rangeAttenuation = getRangeAttenuation(range, length(pointToLight)); \n        return rangeAttenuation   * intensity * color;\n        }\n\n        float czm_private_shadowVisibility(float visibility, float nDotL, float normalShadingSmooth, float darkness) {\n        float strength = clamp(nDotL / normalShadingSmooth, 0.0, 1.0);\n        visibility *= strength;\n        visibility = max(visibility, darkness);\n        return visibility;\n        }\n\n        struct _shadowParameters {\n        vec3 texCoordsCube;\n        vec2 texCoords;\n        float depthBias;\n        float depth;\n        float nDotL;\n        vec2 texelStepSize;\n        float normalShadingSmooth;\n        float darkness;\n        };\n\n        float shadowVisibilityCube(samplerCube shadowMap, _shadowParameters shadowParameters) {\n        float depthBias = shadowParameters.depthBias;\n        float depth = shadowParameters.depth;\n        float nDotL = shadowParameters.nDotL;\n        float normalShadingSmooth = shadowParameters.normalShadingSmooth;\n        float darkness = shadowParameters.darkness;\n        vec3 uvw = shadowParameters.texCoordsCube;\n\n        depth -= depthBias;\n        return czm_shadowDepthCompare(shadowMap, uvw, depth);\n        }\n        float shadowVisibility2D(sampler2D shadowMap, _shadowParameters shadowParameters) {\n        float depthBias = shadowParameters.depthBias;\n        float depth = shadowParameters.depth;\n        float nDotL = shadowParameters.nDotL;\n        float normalShadingSmooth = shadowParameters.normalShadingSmooth;\n        float darkness = shadowParameters.darkness;\n        vec2 uv = shadowParameters.texCoords;\n\n        depth -= depthBias;\n        return czm_shadowDepthCompare(shadowMap, uv, depth);\n        }\n\n        vec3 getPointLightColor(vec3 normalEC, vec3 positionEC,   samplerCube lightShadowMapCube) {\n        vec4 lightPEC = lightPositionEC ;\n        vec2 shadowMapDT = shadowMapDarknessType ;\n        vec3 pointToLightEC = positionEC - lightPEC.xyz;\n        float pointToLightECLength = length(pointToLightEC);\n        vec3 l = normalize(pointToLightEC);\n        float NdotL = clamp(dot(- normalEC, l), 0.0, 1.0);\n\n        float visibility = 0.0;\n        float radius = lightPEC.w;\n        float type = shadowMapDT.y;\n\n        if(pointToLightECLength <= radius) {\n        vec4 shadowMapTSDBANSS = shadowMapTexelSizeDepthBiasAndNormalShadingSmooth;\n\n        _shadowParameters shadowParameters;\n        shadowParameters.texelStepSize = shadowMapTSDBANSS.xy;\n        shadowParameters.depthBias = shadowMapTSDBANSS.z;\n        shadowParameters.normalShadingSmooth = shadowMapTSDBANSS.w;\n        shadowParameters.darkness = shadowMapDT.x;\n        shadowParameters.depth = pointToLightECLength / radius;\n        shadowParameters.nDotL = NdotL;\n        shadowParameters.texCoordsCube = czm_inverseViewRotation * l;\n        visibility = shadowVisibilityCube(lightShadowMapCube, shadowParameters);\n        }\n\n        if(visibility == 1.0) {\n        vec3 colorIntensity = getLightIntensity(lightColor, intensity, type, lightPEC.w, pointToLightEC, direction, outerConeCos, innerConeCos);\n        return NdotL * colorIntensity;\n        }\n        return vec3(0.0);\n        } \n\n        void main() {\n        vec4 color = texture(colorTexture, v_textureCoordinates);\n        float logDepthOrDepth = czm_unpackDepth(texture(depthTexture, v_textureCoordinates));\n        if(logDepthOrDepth >= 1.0) {\n        out_FragColor = color;\n        return;\n        }\n\n        vec4 eyeCoordinate = czm_windowToEyeCoordinates(v_textureCoordinates.xy * czm_viewport.zw, logDepthOrDepth);\n        vec3 downUp = vectorFromOffset(eyeCoordinate, vec2(0.0, 1.0));\n        vec3 leftRight = vectorFromOffset(eyeCoordinate, vec2(1.0, 0.0));\n        vec3 normalEC = normalize(cross(leftRight, downUp));\n        vec3 positionEC = eyeCoordinate.xyz / eyeCoordinate.w; \n        vec3 lightColor= getPointLightColor(normalEC, positionEC,lightShadowMapCube);\n        out_FragColor = vec4(color.xyz + lightColor, 1.0);\n        }\n        '
              }
            })(e)
          case Q.spotLight:
            return new (class extends Ke {
              constructor(e) {
                const i = t
                super(e),
                  (this._typeName = '聚光灯'),
                  (this._graphicType = Q.spotLight),
                  (this[i(343)][i(1198)] = Cesium[i(1960)](this._style[i(1198)], 0)),
                  (this[i(343)][i(358)] = Cesium.defaultValue(this[i(343)][i(358)], 0)),
                  (this._style.roll = Cesium[i(1960)](this[i(343)].roll, 0)),
                  (this[i(343)][i(982)] = Cesium[i(1960)](this[i(343)].outerConeCos, 45)),
                  (this[i(343)][i(670)] = Cesium.defaultValue(this[i(343)][i(670)], 10)),
                  this[i(1737)]()
              }
              [t(637)](e) {
                const i = t,
                  s = {}
                return (
                  (s.lightCamera = this[i(1365)]),
                  (s[i(1036)] = !1),
                  (s[i(1390)] = 1),
                  (s[i(640)] = !1),
                  (s[i(1600)] = !0),
                  (s[i(364)] = !1),
                  (s.context = e[i(696)][i(2325)]),
                  (s[i(1386)] = this[i(343)][i(1981)]),
                  (s.fromLightSource = !1),
                  new Cesium.ShadowMap(s)
                )
              }
              [t(304)](e) {
                const i = t
                let s = new Cesium[i(2162)](e[i(696)])
                const n = this.position
                let o = Cesium[i(310)][i(667)](n[0], n[1], n[2])
                return (s[i(2251)] = o), (this[i(1365)] = s), this[i(1737)](), s
              }
              [t(266)]() {
                const e = t
                let i = this[e(2444)][e(270)],
                  s = Cesium.Cartesian4[e(2434)](
                    1 / this._shadowMap[e(1188)].x,
                    1 / this[e(2444)][e(1188)].y,
                    i.depthBias,
                    i.normalShadingSmooth
                  ),
                  n = new Cesium.Cartesian2(this._shadowMap.darkness, 2)
                return {
                  lightPositionEC: () => this._shadowMap[e(600)],
                  intensity: () => this[e(343)][e(2480)],
                  lightColor: () => this._style[e(301)],
                  shadowMapMatrix: () => this[e(2444)][e(829)],
                  shadowMapTexelSizeDepthBiasAndNormalShadingSmooth: () => s,
                  shadowMapDarknessType: () => n,
                  direction: () => this[e(2444)]._lightDirectionEC,
                  outerConeCos: () => Math.cos(Cesium.Math.toRadians(this._style[e(982)])),
                  innerConeCos: () =>
                    Math[e(1272)](Cesium[e(475)][e(1149)](this[e(343)].innerConeCos)),
                  lightShadowMapCube: () => this[e(2444)][e(1247)]
                }
              }
              [t(1764)]() {
                return t(912)
              }
              _setCamera() {
                const e = t
                this._camera &&
                  (this[e(1365)][e(1206)]({
                    destination: this[e(2238)],
                    orientation: {
                      heading: Cesium[e(475)][e(1149)](this[e(343)][e(1198)]),
                      pitch: Cesium.Math.toRadians(this[e(343)][e(358)]),
                      roll: Cesium[e(475)][e(1149)](this[e(343)][e(1143)])
                    }
                  }),
                  (this[e(1365)][e(1850)].fov = Cesium.Math[e(1149)](2 * this[e(343)][e(982)])),
                  (this[e(1365)].frustum[e(348)] = this[e(343)][e(1981)]),
                  (this[e(1365)][e(1850)][e(1292)] = 0.1),
                  (this[e(1365)][e(1850)][e(124)] = 1))
              }
            })(e)
        }
      }
    }
    class $e extends b {
      constructor(e = {}) {
        const i = t
        super(e),
          (this[i(1115)] = F),
          (this._geometryType = i(1632)),
          (this[i(943)] = Cesium[i(1960)](this._style[i(820)], {})),
          (this[i(252)] = 1036800),
          (this[i(2469)] = Cesium[i(1960)](e[i(2333)], [])),
          (this[i(794)] = this._positions)
      }
      [t(1504)]() {
        const e = t
        this[e(484)]()
        let i = this._convertData(this[e(2469)])
        this._setCartesian3Array()
        let s = Math[e(1427)](i[e(2135)]),
          n = Math[e(1427)](i.canvasWidth)
        ;(this[e(1331)] = this[e(1985)](s, n)), this[e(2023)]()[e(1621)](this._heatContainer)
        let o = i.heatmapData
        const r = {}
        ;(r.data = o),
          (r.max = i[e(1113)]),
          (this[e(2142)] = h337[e(1951)]({ container: this[e(1331)], ...this._heatStyle })),
          this[e(2142)][e(373)](r),
          (this._rect = i[e(1342)]),
          (this[e(1325)] = this[e(2142)]._renderer[e(1493)]),
          this[e(2287)]()
      }
      [t(2231)]() {
        this[t(1504)]()
      }
      [t(2287)]() {}
      _createHeatContainer(e, i) {
        const s = t
        let n = document[s(1945)](s(1120))
        return (n[s(1679)][s(2306)] = e + 'px'), (n[s(1679)][s(575)] = i + 'px'), n
      }
      [t(484)]() {
        const e = t
        this[e(1331)] && (this._heatContainer[e(1896)](), (this[e(1331)] = null))
      }
      [t(2023)]() {
        const e = t
        let i = document.getElementById(e(977))
        return (
          i ||
            ((i = document[e(1945)](e(1120))).setAttribute('id', 'xt3d-heatmap-container'),
            (i.style[e(1811)] = e(2529)),
            document[e(879)][e(1621)](i)),
          i
        )
      }
      [t(391)](e) {
        const i = t
        let s = 180,
          n = 90,
          o = 0,
          r = 0
        e[i(1602)]((t) => {
          const e = i
          ;(s = Math.min(t[0], s)),
            (n = Math[e(1429)](t[1], n)),
            (o = Math[e(1113)](t[0], o)),
            (r = Math[e(1113)](t[1], r))
        })
        let a = o - s,
          h = r - n,
          l = Math[i(1427)](Math[i(1345)](h * this[i(252)])),
          c = Math[i(1427)]((l * a) / h),
          u = [],
          m = -1e4
        return (
          e[i(1602)]((t) => {
            const e = i
            let n = ((t[0] - s) / a) * c,
              o = (-(t[1] - r) / h) * l
            const p = {}
            ;(p.x = n), (p.y = o), (p.value = t[2]), u[e(2553)](p), (m = Math[e(1113)](t[2], m))
          }),
          {
            max: m,
            heatmapData: u,
            canvasWidth: c,
            canvasHeight: l,
            rect: new Cesium.Rectangle(
              Cesium[i(475)][i(1149)](s),
              Cesium[i(475)].toRadians(n),
              Cesium[i(475)][i(1149)](o),
              Cesium[i(475)].toRadians(r)
            )
          }
        )
      }
      [t(947)]() {
        const e = t
        this[e(2469)][e(1602)]((t) => {
          const i = e
          this[i(1378)][i(2553)](Cesium[i(310)][i(667)](t[0], t[1], 0))
        }),
          (this[e(1024)] = Cesium[e(1242)][e(951)](this[e(1378)]))
      }
      [t(2515)](e) {
        const i = t
        ;(this[i(2469)] = e), (this[i(794)] = this[i(2469)]), this[i(1504)]()
      }
    }
    let ti = {
      readFeature(e) {
        const i = t
        let s = e[i(138)][i(2365)],
          n = e[i(1004)]
        return (n.positions = s), this[i(1951)](n)
      },
      create(e) {
        const i = t
        switch (e.graphicType) {
          case Q[i(2102)]:
            return new (class extends $e {
              constructor(e) {
                const i = t
                super(e),
                  (this[i(2336)] = Q[i(2102)]),
                  (this[i(998)] = Cesium[i(515)][i(2407)]({
                    cull: { enabled: !0 },
                    depthTest: { enabled: !0 },
                    stencilTest: {
                      enabled: !0,
                      frontFunction: Cesium.StencilFunction.ALWAYS,
                      frontOperation: {
                        fail: Cesium.StencilOperation[i(2424)],
                        zFail: Cesium[i(2340)][i(2424)],
                        zPass: Cesium[i(2340)][i(2168)]
                      },
                      backFunction: Cesium[i(878)][i(1881)],
                      backOperation: {
                        fail: Cesium[i(2340)][i(2424)],
                        zFail: Cesium[i(2340)][i(2424)],
                        zPass: Cesium[i(2340)][i(2168)]
                      },
                      reference: 2,
                      mask: 2
                    },
                    blending: Cesium[i(1870)][i(1785)]
                  })),
                  (this[i(1679)][i(2306)] = Cesium.defaultValue(this.style.height, 10)),
                  (this[i(1679)][i(666)] = Cesium.defaultValue(this[i(1679)][i(666)], 0)),
                  (this[i(1679)][i(2556)] = Cesium[i(1960)](this.style[i(2556)], 2)),
                  (this[i(1679)][i(2315)] = Cesium[i(1960)](this[i(1679)][i(2315)], 1)),
                  (this[i(1934)] = this._createPrimitive()),
                  this._createTexture()
              }
              [t(2287)]() {
                const e = t
                let i = this[e(2214)]()
                const s = {}
                ;(s[e(138)] = i), (this[e(1934)][e(1122)] = new Cesium[e(2155)](s))
                const n = {}
                ;(n[e(599)] = 'texture(image, fract(materialInput.st)).rgb'), (n[e(2315)] = e(832))
                let o = new Cesium[e(985)]({
                  aboveGround: !0,
                  renderState: this[e(998)],
                  material: new Cesium[e(1637)]({
                    fabric: {
                      uniforms: {
                        image: this[e(1325)],
                        maxHeight: this.style[e(2306)],
                        heightRatio: this[e(1679)].heightRatio,
                        alpha: this[e(1679)].alpha
                      },
                      components: n
                    },
                    translucent: function (t) {
                      return !0
                    }
                  }),
                  vertexShaderSource: e(251)
                })
                this[e(1934)].appearance = o
              }
              [t(438)]() {
                const e = t,
                  i = {}
                return (i[e(1122)] = []), (i[e(644)] = !1), new Cesium.Primitive(i)
              }
              [t(2214)]() {
                const e = t
                return new Cesium[e(492)]({
                  rectangle: this[e(476)],
                  vertexFormat: Cesium[e(1788)].ALL,
                  granularity: this[e(343)].granularity || 5e-5,
                  height: this[e(1679)][e(666)]
                })
              }
              _setCartesian3Array() {
                const e = t
                this[e(2469)][e(1602)]((t) => {
                  const i = e
                  this[i(1378)][i(2553)](
                    Cesium.Cartesian3.fromDegrees(t[0], t[1], this[i(1679)][i(666)])
                  )
                }),
                  (this[e(1024)] = Cesium[e(1242)][e(951)](this._cartesian3Array))
              }
              [t(415)](e) {
                const i = t
                this[i(1934)][i(1482)] = e
              }
              _addHook(e) {
                const i = t
                ;(this[i(1934)].layerId = e.id), e[i(245)].scene[i(1346)][i(1861)](this[i(1934)])
              }
              [t(2389)](e) {
                const i = t
                e[i(245)][i(696)].primitives[i(1896)](this[i(1934)])
              }
            })(e)
          case Q.planeHeat:
            return new (class extends $e {
              constructor(e) {
                const i = t
                super(e),
                  (this[i(2336)] = Q.planeHeat),
                  (this[i(987)] = this[i(1377)]()),
                  this[i(1504)]()
              }
              [t(2287)]() {
                const e = t
                let i = new Cesium[e(1100)]({
                  image: new Cesium[e(2569)]((t) => this[e(1325)], !1),
                  transparent: !0
                })
                this._entity[e(1412)] = { coordinates: this[e(476)], ...this[e(1679)], material: i }
              }
              [t(1377)]() {
                return new Cesium[t(1380)]({})
              }
              [t(415)](t) {
                this._entity.show = t
              }
              _addHook(e) {
                const i = t
                ;(this[i(987)][i(2122)] = e.id), e[i(245)].entities.add(this[i(987)])
              }
              [t(2389)](e) {
                const i = t
                e[i(245)][i(118)].remove(this._entity)
              }
            })(e)
        }
      }
    }
    const ei = {}
    ;(ei[t(2080)] = t(2080)), (ei[t(1053)] = 'upright'), (ei[t(659)] = t(659))
    const ii = ei
    let si = {
      readFeature(e) {
        const i = t
        let s = e.geometry.coordinates,
          n = e[i(1004)]
        return (n.position = s), this.create(n)
      },
      create(e) {
        const i = t
        if (e.graphicType === Q[i(2316)])
          return new (class extends b {
            constructor(e) {
              const i = t
              super(e),
                (this[i(1115)] = B),
                (this[i(2336)] = Q[i(2316)]),
                (this[i(528)] = '站牌'),
                (this._fixPointCount = 1),
                (this[i(314)] = i(1865)),
                (this[i(343)][i(1640)] = Cesium.defaultValue(this._style[i(1640)], ii[i(1053)])),
                ii[this[i(343)][i(1640)]] || (this._style[i(1640)] = ii[i(1053)]),
                (this[i(343)][i(1241)] = Cesium.defaultValue(this[i(343)][i(1241)], 0)),
                (this[i(343)][i(2306)] = Cesium.defaultValue(this[i(343)][i(2306)], 20)),
                (this[i(343)].width = Cesium[i(1960)](this[i(343)][i(575)], 20)),
                (this[i(343)][i(735)] = Cesium[i(1960)](this[i(343)][i(735)], 20)),
                (this[i(343)][i(1318)] = Cesium.defaultValue(this[i(343)].fontColor, i(1655))),
                (this[i(343)][i(1982)] = Cesium[i(1960)](this[i(343)].text, '站牌图元'))
              let s = Cesium[i(1960)](e.position, [111, 28, 0])
              this.setPosition(s), (this[i(987)] = this._createEntiy()), this[i(2231)]()
            }
            [t(415)](e) {
              const i = t
              this._entity[i(1482)] = e
            }
            [t(2547)]() {
              const e = t,
                i = {}
              return (
                (i.graphicId = this[e(1570)]),
                (i[e(2467)] = {}),
                (i[e(1412)] = {}),
                new Cesium[e(1380)](i)
              )
            }
            _loadMaterial() {
              const e = t
              if (this[e(343)][e(259)]) {
                let t = new Image()
                ;(t[e(1744)] = this[e(343)][e(259)]),
                  (t[e(1701)] = (e) => {
                    this._setMaterial(t)
                  }),
                  (t[e(1978)] = (t) => {
                    const i = e
                    throw new Cesium.DeveloperError(i(1843))
                  })
              } else {
                let t = this._style
                ;(t[e(735)] = t[e(735)]), (t[e(1318)] = t[e(1318)])
                const i = t.text,
                  s = document[e(1945)](e(1493)),
                  n = (i + '')[e(277)] * t[e(735)]
                ;(s[e(575)] = n + 10), (s[e(2306)] = t[e(735)])
                const o = s.getContext('2d')
                ;(o[e(537)] = t[e(1318)]),
                  (o[e(1291)] = e(886) + t[e(735)] + 'px 微软雅黑'),
                  (o[e(504)] = e(520)),
                  o[e(1104)](i, 5, 0, s.width, s.height),
                  this[e(609)](s)
              }
            }
            _setMaterial(e) {
              const i = t,
                s = {}
              ;(s[i(259)] = e), (s.transparent = !0)
              let n = new Cesium.ImageMaterialProperty(s)
              ;(this[i(987)].wall.material = n), (this[i(987)][i(1412)].material = n)
            }
            [t(1253)](e) {
              const i = t
              ;(this[i(1912)] = e),
                (this[i(794)] = this[i(1912)]),
                (this[i(2238)] = Cesium[i(310)][i(667)](e[0], e[1], e[2])),
                this._entity && this._updateGeometry(),
                (this._boundingSphere = new Cesium[i(1242)](this[i(2238)], this[i(343)].width))
            }
            [t(2231)]() {
              const e = t
              this[e(1092)](),
                this[e(2270)](),
                (this[e(1024)] = new Cesium.BoundingSphere(this[e(2238)], this[e(343)][e(575)]))
            }
            [t(2270)]() {
              const e = t
              let i = Cesium.Cartesian3[e(667)](
                this.position[0],
                this[e(2251)][1],
                this[e(2251)][2]
              )
              if (this[e(343)][e(1640)] == ii[e(1053)])
                return (
                  this[e(464)](i),
                  (this[e(987)].rectangle.show = !1),
                  void (this[e(987)][e(2467)][e(1482)] = !0)
                )
              this[e(2182)](i),
                (this[e(987)][e(1412)][e(1482)] = !0),
                (this._entity[e(2467)].show = !1)
            }
            [t(464)](e) {
              const i = t
              let s = this._style.width / 2,
                n = this[i(1426)](e, -this._style[i(1241)], s),
                o = this[i(1426)](e, -this[i(343)][i(1241)] - 180, s),
                r = Cesium.Cartographic.fromCartesian(e),
                a = [o, n],
                h = new Array(2).fill(r[i(2306)]),
                l = new Array(2).fill(r[i(2306)] + this._style[i(2306)])
              this[i(2520)] &&
                ((a = new Cesium[i(2569)]((t) => [o, n])),
                (h = new Cesium.CallbackProperty((t) => new Array(2)[i(2327)](r[i(2306)]))),
                (l = new Cesium[i(2569)]((t) =>
                  new Array(2)[i(2327)](r[i(2306)] + this[i(343)][i(2306)])
                ))),
                (this[i(987)][i(2467)][i(2333)] = a),
                (this._entity.wall[i(2332)] = h),
                (this[i(987)][i(2467)][i(1089)] = l)
            }
            _udpateRectangleGeometry(e) {
              const i = t
              let s = this[i(343)][i(575)] / 2,
                n = this._getPositoinByAngle(e, 0, s),
                o = this._getPositoinByAngle(e, -180, s)
              const r = this[i(343)].height / 2
              let a = Cesium[i(310)][i(1676)](o, n, new Cesium.Cartesian3()),
                h = [],
                l = Cesium[i(310)].cross(a, n, new Cesium.Cartesian3())
              l = Cesium[i(310)][i(379)](l, l)
              let c = new Cesium.Ray(o, l)
              e = Cesium[i(2007)].getPoint(c, r, new Cesium[i(310)]())
              let u = Cesium[i(2285)].fromCartesian(e)
              h[i(2553)](u[i(2106)], u.latitude),
                Cesium.Cartesian3.negate(l, l),
                (c = new Cesium.Ray(n, l)),
                (e = Cesium[i(2007)][i(1867)](c, r, new Cesium[i(310)]())),
                (u = Cesium.Cartographic.fromCartesian(e)),
                h.push(u[i(2106)], u[i(199)])
              let m = new Cesium[i(613)](h[0], h[1], h[2], h[3])
              this._style[i(1640)] == ii[i(2080)]
                ? delete this[i(987)][i(1412)].height
                : (this._entity[i(1412)].height = this[i(2251)][2])
              let p = m
              this[i(2520)] && (p = new Cesium[i(2569)]((t) => m)),
                (this[i(987)].rectangle[i(2365)] = p),
                (this[i(987)][i(1412)][i(2307)] = Cesium.Math[i(1149)](this[i(343)].angle)),
                (this[i(987)][i(1412)].stRotation = Cesium[i(475)][i(1149)](
                  this[i(343)][i(1241)] + 90
                ))
            }
            [t(1426)](e, i, s) {
              const n = t
              let o = Cesium.Cartographic[n(2579)](e)
              ;(o[n(2106)] = Cesium[n(475)][n(363)](o[n(2106)])),
                (o[n(199)] = Cesium[n(475)].toDegrees(o[n(199)]))
              let r = s * Math[n(884)]((i * Math.PI) / 180),
                a = s * Math[n(1272)]((i * Math.PI) / 180),
                h = 6356725 + (21412 * (90 - o[n(199)])) / 90,
                l = h * Math[n(1272)]((o[n(199)] * Math.PI) / 180)
              return (
                (o[n(2106)] = (180 * (r / l + (o[n(2106)] * Math.PI) / 180)) / Math.PI),
                (o[n(199)] = (180 * (a / h + (o[n(199)] * Math.PI) / 180)) / Math.PI),
                Cesium[n(310)].fromDegrees(o[n(2106)], o[n(199)], o[n(2306)])
              )
            }
            [t(1545)](e) {
              const i = t
              ;(this._entity[i(2122)] = e.id),
                (this[i(2462)] = e),
                e._viewer[i(118)][i(1861)](this._entity)
            }
            [t(2389)](e) {
              const i = t
              e[i(245)].entities[i(1896)](this[i(987)])
            }
          })(e)
      }
    }
    class ni {
      constructor(e = {}) {
        const i = t
        ;(this[i(1482)] = Cesium[i(1960)](e[i(1482)], !0)),
          (this[i(634)] = Cesium[i(1960)](e[i(634)], '')),
          (this.scale = Cesium[i(1960)](e[i(2341)], 1)),
          (this[i(328)] = Cesium[i(1960)](e[i(328)], 0)),
          (this[i(206)] = Cesium.defaultValue(e.maximumScale, null)),
          (this[i(1512)] = Cesium[i(1960)](e[i(1512)], !0)),
          (this[i(283)] = Cesium[i(1960)](e[i(283)], !0)),
          (this[i(208)] = Cesium[i(1960)](e[i(208)], Cesium[i(862)].ENABLED)),
          (this[i(697)] = Cesium[i(1960)](e[i(697)], Cesium[i(893)][i(2184)])),
          (this[i(1234)] = Cesium[i(1960)](e[i(1234)], i(2396))),
          (this[i(351)] = Cesium[i(1960)](e.silhouetteSize, 0)),
          (this[i(1070)] = Cesium[i(1960)](e[i(1070)], i(287))),
          (this.colorBlendMode = Cesium.defaultValue(e[i(217)], Cesium.ColorBlendMode[i(1161)])),
          (this[i(1353)] = Cesium[i(1960)](e[i(1353)], 0.5)),
          (this.distanceDisplayCondition = Cesium[i(1960)](e[i(978)], null))
      }
      [t(1567)](e) {
        const i = t
        for (const t in this) {
          const i = this[t]
          null != i && null != i && (e[t] = i)
        }
        ;(e[i(1070)] = Cesium[i(1154)][i(2008)](this[i(1070)])),
          (e[i(1234)] = Cesium.Color[i(2008)](this[i(1234)])),
          (e.show = this[i(1482)] && this[i(634)])
      }
    }
    let oi = {
      readFeature(e) {
        const i = t
        let s = e[i(138)].coordinates,
          n = e[i(1004)]
        return (n.positions = s), this[i(1951)](n)
      },
      create(e) {
        const i = t
        if (e[i(1158)] === Q[i(265)])
          return new (class extends b {
            constructor(e) {
              const i = t
              super(e),
                (this[i(1115)] = N),
                (this[i(2336)] = Q[i(265)]),
                (this[i(528)] = i(620)),
                (this[i(324)] = 2),
                (this[i(968)] = Cesium[i(1960)](e[i(2535)], 0.5)),
                (this[i(1102)] = Cesium[i(1960)](e[i(2328)], !1)),
                (this[i(1055)] = Cesium[i(1960)](e.trackedView, !1)),
                (this[i(2572)] = Cesium[i(1960)](e.headingPitchRange, {
                  heading: 90,
                  pitch: -45,
                  range: 100
                })),
                (this[i(343)].show = Cesium[i(1960)](this[i(343)][i(1482)], !0)),
                (this[i(343)][i(575)] = Cesium.defaultValue(this[i(343)].width, 1)),
                (this[i(343)].color = Cesium[i(1960)](this[i(343)].color, i(615))),
                (this[i(343)][i(2256)] = Cesium[i(1960)](this._style[i(2256)], !0)),
                (this._style[i(990)] = Cesium[i(1960)](this._style[i(990)], 4)),
                (this[i(343)][i(1029)] = Cesium[i(1960)](this._style[i(1029)], i(1599))),
                (this._positionProperty = null),
                (this[i(291)] = null),
                (this[i(1118)] = null),
                (this[i(1828)] = 0),
                (this._drivenDistance = 0),
                (this[i(737)] = 0),
                (this[i(263)] = 0),
                (this[i(2469)] = Cesium[i(1960)](e.positions, [])),
                (this[i(987)] = this._createEntiy()),
                this[i(2515)](this[i(2469)]),
                this[i(2231)](),
                (this[i(303)] = []),
                (this[i(280)] = []),
                (this._runing = !1),
                (this._geometryType = i(2310)),
                (this[i(2119)] = new Dt(e[i(2240)])),
                (this._model = new ni(e.model)),
                (this[i(556)] = new It(e[i(909)]))
            }
            get totalTime() {
              return Number(this._totalTime.toFixed(1))
            }
            get [t(2353)]() {
              return this[t(263)]
            }
            get [t(1367)]() {
              const e = t
              return Number(this[e(1828)][e(1268)](2))
            }
            get [t(804)]() {
              return this[t(1782)]
            }
            get [t(2154)]() {
              return this[t(2323)]
            }
            get [t(925)]() {
              return this[t(291)]
            }
            get [t(904)]() {
              return this._availability
            }
            get [t(1106)]() {
              return this[t(303)]
            }
            set [t(2328)](e) {
              const i = t
              ;(this[i(1102)] = e),
                !e && this._layer[i(245)].camera[i(1610)](Cesium[i(2066)][i(1393)])
            }
            get [t(2328)]() {
              return this._lockView
            }
            set trackedView(e) {
              const i = t
              ;(this[i(1055)] = e),
                e
                  ? (this[i(2462)][i(245)][i(1889)] = this[i(1467)])
                  : (this[i(2462)][i(245)].trackedEntity = null)
            }
            get [t(2123)]() {
              return this._trackedView
            }
            set [t(380)](t) {
              this._headingPitchRange = t
            }
            get [t(380)]() {
              return this[t(2572)]
            }
            get [t(2535)]() {
              return this[t(968)]
            }
            set [t(2535)](e) {
              this[t(968)] = e
            }
            [t(498)]() {
              const e = t
              if (this[e(2350)]) return
              let i = this[e(2462)][e(245)]
              this[e(690)](),
                i[e(2391)][e(654)][e(1973)](this._tickEventHandler, this),
                (i[e(2391)][e(205)] = this[e(692)][e(205)].clone()),
                (i.clock[e(2072)] = this[e(692)][e(2072)][e(1902)]()),
                (i.clock[e(2503)] = this[e(692)][e(205)][e(1902)]()),
                (i[e(2391)][e(604)] = Cesium[e(1846)][e(2435)]),
                (i[e(2391)][e(1433)] = !0),
                (this[e(2350)] = !0),
                this[e(450)](p.routeStart),
                this[e(1729)]()
            }
            [t(1133)]() {
              const e = t,
                i = this._layer[e(245)],
                s = i[e(2391)][e(2503)][e(1902)](),
                n = this._stopsTimes[this[e(795)]]
              let o = this._positionProperty[e(1474)](s)
              if (
                (o &&
                  ((this[e(303)] = this[e(280)].concat(o)),
                  (this[e(1782)] = this[e(1014)](this[e(303)]))),
                (this._drivenTime = Cesium.JulianDate[e(2479)](s, i[e(2391)].startTime)),
                this[e(450)](p[e(1017)], o),
                Cesium[e(1665)][e(944)](n, s) &&
                  (this[e(450)](p[e(1205)], this[e(795)]),
                  this[e(280)].push(this[e(1378)][this[e(795)]]),
                  this._nextStopsIndex++,
                  this[e(795)] == this._stopsTimes[e(277)]))
              )
                this[e(519)]()
              else if (this[e(1102)] && !this[e(1055)]) {
                let t = this[e(291)].getValue(s),
                  n = Cesium[e(2058)][e(1648)](o)
                n = Cesium.Matrix4[e(942)](Cesium.Matrix3.fromQuaternion(t), o)
                const r = Cesium[e(475)].toRadians(this[e(380)].heading),
                  a = Cesium.Math.toRadians(this[e(380)][e(711)])
                let h = this[e(380)][e(812)]
                h < 0.1 && (h = 0.1),
                  i[e(1391)].lookAtTransform(n, new Cesium.HeadingPitchRange(r, a, h))
              }
            }
            stop() {
              const e = t
              this._runing = !1
              let i = this[e(2462)]._viewer
              this[e(450)](p[e(2414)]),
                i[e(2391)][e(654)][e(704)](this[e(1133)], this),
                i[e(1391)][e(1610)](Cesium.Matrix4.IDENTITY),
                this[e(410)]()
            }
            _createAnimateEntity() {
              const e = t
              let i = this[e(2462)]._viewer[e(118)][e(1861)]({
                position: this[e(2154)],
                orientation: this[e(925)],
                label: {},
                billboard: {},
                model: {},
                polyline: {
                  show: this[e(343)][e(2256)],
                  positions: new Cesium[e(2569)]((t) => this.drivenPositions, !1),
                  material: Cesium[e(1154)].fromCssColorString(this[e(343)][e(1029)]),
                  width: this[e(343)].drivenLineWidth
                }
              })
              ;(i.viewFrom = new Cesium[e(310)](-180, -180, 200)),
                this[e(1055)] && (this[e(2462)][e(245)][e(1889)] = i),
                this[e(2119)][e(1567)](i.label),
                this[e(556)][e(1567)](i.billboard),
                this[e(2226)][e(1567)](i.model),
                (this[e(1467)] = i)
            }
            [t(410)]() {
              const e = t
              this[e(1467)] &&
                (this[e(2462)][e(245)][e(118)][e(1896)](this[e(1467)]), (this._tempEntity = null))
            }
            [t(200)]() {
              const e = t
              let i = this._cartesian3Array
              this._editMode && (i = new Cesium.CallbackProperty((t) => this[e(1378)], !1)),
                (this._entity[e(811)].positions = i)
            }
            [t(415)](e) {
              const i = t
              this[i(987)][i(1482)] = e
            }
            _createEntiy() {
              const e = t,
                i = {}
              return (
                (i[e(2366)] = this[e(1570)]),
                (i[e(811)] = {}),
                (i[e(811)].show = this._style[e(1482)]),
                (i[e(811)][e(939)] = 3),
                new Cesium[e(1380)](i)
              )
            }
            setPositions(e = []) {
              const i = t
              ;(this[i(2469)] = e),
                this[i(2469)].length < 2 ||
                  ((this[i(1378)] = this._generatePositions(e)),
                  this[i(1378)] &&
                    (this[i(987)] &&
                      !this[i(2520)] &&
                      (this[i(987)][i(811)][i(2333)] = this[i(1378)]),
                    this[i(690)](),
                    (this[i(794)] = this[i(2469)]),
                    (this[i(1024)] = Cesium.BoundingSphere.fromPoints(
                      this[i(1378)],
                      new Cesium.BoundingSphere()
                    ))))
            }
            [t(690)]() {
              const e = t
              let i = this._computeLineInfo(this[e(1378)], this[e(968)])
              ;(this._lineInfo = i),
                (this._positionProperty = this[e(1226)](
                  this[e(1378)],
                  i[e(205)],
                  i[e(1688)][e(1926)]
                )),
                (this._availability = new Cesium[e(1314)]([
                  new Cesium[e(760)]({ start: i[e(205)], stop: i[e(2072)] })
                ])),
                (this._orientation = new Cesium[e(1559)](this[e(2323)])),
                (this[e(1730)] = []),
                (this[e(795)] = 0),
                (this[e(303)] = []),
                (this._drivenPositionsTemp = []),
                (this._totalDistance = i[e(1688)][e(1367)]),
                (this[e(737)] = i[e(1688)].timeSum),
                (this._drivenDistance = 0),
                (this[e(263)] = 0),
                i[e(1688)].siteTimes[e(1602)]((t) => {
                  const s = e,
                    n = Cesium.JulianDate[s(470)](i[s(205)], t, new Cesium[s(1665)]())
                  this[s(1730)].push(n)
                })
            }
            [t(1226)](e, i, s) {
              const n = t
              let o = new Cesium[n(765)]()
              for (let t = 0; t < e.length; t++) {
                const r = Cesium[n(1665)].addSeconds(i, s[t], new Cesium[n(1665)]())
                o[n(460)](r, e[t])
              }
              return o
            }
            [t(2505)](e, i) {
              const s = t
              let n = {}
              return (
                (n.timeInfo = this[s(614)](e, i)),
                (n[s(205)] = Cesium[s(1665)][s(1447)](new Date())),
                (n[s(2072)] = Cesium[s(1665)].addSeconds(
                  n[s(205)],
                  n[s(1688)].timeSum,
                  new Cesium[s(1665)]()
                )),
                n
              )
            }
            [t(614)](e, i) {
              const s = t
              let n,
                o = [0],
                r = 0,
                a = 0
              for (let t = 1; t < e.length; t++)
                (a += n = this[s(1014)]([e[t - 1], e[t]])), (r += n / i), o[s(2553)](r)
              const h = {}
              return (h[s(126)] = r), (h[s(1926)] = o), (h.totalDistance = a), h
            }
            [t(1014)](e) {
              const i = t
              if (!e || e[i(277)] < 2) return 0
              let s = 0
              for (let t = 1; t < e[i(277)]; t++) s += Cesium.Cartesian3[i(1849)](e[t - 1], e[t])
              return s
            }
            _generatePositions(e) {
              const i = t
              if (!(e[i(277)] < 2)) return this[i(1690)](e)
            }
            [t(1690)](e) {
              const i = t
              return Cesium[i(310)][i(774)]([][i(1500)][i(560)]([], e))
            }
            [t(2231)]() {
              const e = t,
                i = this._style
              ;(this[e(987)][e(811)].show = i.show),
                (this[e(987)].polyline[e(322)] = Cesium.Color[e(2008)](this._style.color))
            }
            [t(1545)](e) {
              const i = t
              ;(this[i(987)].layerId = e.id),
                (this._layer = e),
                e[i(245)][i(118)][i(1861)](this[i(987)])
            }
            [t(2389)](e) {
              const i = t
              this[i(519)](), e[i(245)][i(118)][i(1896)](this._entity)
            }
            [t(1114)]() {
              const e = t,
                i = {}
              return (
                (i[e(2535)] = this[e(2535)]),
                (i.lockView = this[e(2328)]),
                (i[e(380)] = this[e(380)]),
                i
              )
            }
          })(e)
      }
    }
    class ri extends b {
      constructor(e) {
        const i = t
        super(e),
          (this[i(1115)] = G),
          (this._typeName = '粒子'),
          (this[i(350)] = 1),
          (this._geometryType = i(1865)),
          (this[i(987)] = this[i(2547)]()),
          (this._particle = this[i(1084)]()),
          (this[i(1912)] = Cesium[i(1960)](e[i(2251)], [0, 0, 0])),
          this[i(1253)](this._position),
          (this[i(343)][i(1198)] = Cesium[i(1960)](this._style.heading, 0)),
          (this[i(343)][i(711)] = Cesium[i(1960)](this[i(343)][i(711)], 0)),
          (this[i(742)] = new Cesium[i(2066)]()),
          (this[i(1359)] = new Cesium.Cartesian3()),
          (this[i(219)] = new Cesium[i(1472)]()),
          (this._hpr = new Cesium.HeadingPitchRoll()),
          (this[i(579)] = new Cesium[i(236)]()),
          this._addWindowEvent()
      }
      [t(2205)]() {
        const e = t
        let i =
            'hidden' in document
              ? e(186)
              : e(2420) in document
                ? e(2420)
                : 'mozHidden' in document
                  ? e(1758)
                  : null,
          s = i.replace(/hidden/i, e(677))
        ;(this[e(713)] = i), (this[e(1533)] = s)
        let n = this,
          o = function () {
            const t = e
            n[t(2462)] &&
              (n[t(2462)]._viewer[t(696)][t(1346)].remove(n._particle),
              document[i] ||
                ((n._particle = n[t(1084)]()),
                (n[t(2496)][t(1482)] = n._show),
                n[t(2462)][t(245)][t(696)].primitives[t(1861)](n[t(2496)])))
          }
        ;(this[e(1460)] = o), document.addEventListener(s, o)
      }
      [t(2547)]() {
        const e = t
        return new Cesium.Entity({
          graphicId: this._id,
          point: {
            pixelSize: 15,
            color: Cesium[e(1154)][e(2474)],
            disableDepthTestDistance: 100,
            distanceDisplayCondition: new Cesium.DistanceDisplayCondition(0, 1e3),
            show: !1
          }
        })
      }
      [t(1084)]() {
        return new Cesium[t(2145)]({})
      }
      setPosition(e) {
        const i = t
        ;(this._position = e),
          (this._cartesian3 = Cesium.Cartesian3[i(667)](e[0], e[1], e[2])),
          (this[i(987)][i(2251)] = this[i(2238)]),
          (this._coordinates = this[i(1912)]),
          (this[i(1024)] = new Cesium.BoundingSphere(this[i(2238)], 3))
      }
      _setStyle() {
        const e = t
        let i = this._style,
          s = this[e(2496)]
        delete i[e(1573)], delete i.endColor
        for (const t in i)
          if (Object.hasOwnProperty[e(1669)](i, t)) {
            const e = i[t]
            null != e && null != e && (s[t] = e)
          }
        this[e(987)][e(1542)][e(1482)] = this._isSelected
      }
      [t(415)](e) {
        const i = t
        this[i(2462)] &&
          (this._layer[i(245)].scene.primitives.remove(this[i(2496)]),
          e &&
            ((this[i(2496)] = this[i(1084)]()),
            (this._particle[i(1482)] = e),
            this[i(2462)]._viewer.scene[i(1346)][i(1861)](this[i(2496)])))
      }
      [t(1545)](e) {
        const i = t
        ;(this._entity[i(2122)] = e.id),
          (this[i(2462)] = e),
          e[i(245)][i(118)].add(this[i(987)]),
          e[i(245)][i(696)].primitives[i(1861)](this[i(2496)]),
          e[i(245)][i(696)][i(647)][i(1973)](this[i(1501)], this)
      }
      [t(1501)](e, i) {
        const s = t,
          n = this[s(2496)]
        let o = this._entity[s(2577)](i, new Cesium.Matrix4()),
          r = Cesium[s(2066)][s(1420)](new Cesium[s(310)](0, 0, 3.5), new Cesium[s(2066)]())
        ;(n.modelMatrix = Cesium[s(2066)][s(999)](o, r, new Cesium[s(2066)]())),
          Cesium[s(183)].fromDegrees(this[s(343)][s(1198)], this[s(343)][s(711)], 0, this[s(2412)]),
          (this[s(579)][s(570)] = Cesium[s(310)].fromElements(0, 0, 0, this._translation)),
          (this[s(579)][s(2307)] = Cesium[s(1472)][s(1245)](this[s(2412)], this[s(219)])),
          (n[s(1423)] = Cesium[s(2066)][s(918)](this[s(579)], this[s(742)]))
      }
      _removeHook(e) {
        const i = t
        e[i(245)][i(118)][i(1896)](this[i(987)]),
          e[i(245)][i(696)][i(647)][i(704)](this[i(1501)], this),
          e[i(245)][i(696)].primitives[i(1896)](this[i(2496)]),
          document[i(704)](this[i(1533)], this[i(1460)])
      }
    }
    let ai = {
        readFeature(e) {
          const i = t
          let s = e[i(138)][i(2365)],
            n = e[i(1004)]
          return (n.position = s), this[i(1951)](n)
        },
        create(e) {
          const i = t
          switch (e[i(1158)]) {
            case Q[i(1946)]:
              return new (class extends ri {
                constructor(e) {
                  const i = t
                  super(e),
                    (this._typeName = i(2586)),
                    (this[i(2336)] = Q[i(1946)]),
                    (this[i(343)][i(1573)] = i(2217)),
                    (this[i(343)][i(1072)] = i(1655)),
                    (this._style[i(195)] = 2),
                    (this[i(343)][i(507)] = 3),
                    (this._style[i(119)] = 1.5),
                    (this._style[i(2104)] = 1.5),
                    (this[i(343)][i(1211)] = 1.8),
                    (this[i(343)].minimumSpeed = 7),
                    (this[i(343)][i(1854)] = 9),
                    (this[i(343)][i(1080)] = new Cesium.Cartesian2(2, 2)),
                    (this[i(343)].emissionRate = 200),
                    (this[i(343)][i(1116)] = 16)
                }
                _createParticle() {
                  const e = t
                  return new Cesium[e(2145)]({
                    image: this[e(343)][e(259)],
                    startColor: new Cesium[e(1154)](1, 1, 1, 1),
                    endColor: new Cesium[e(1154)](0.5, 0, 0, 0),
                    particleSize: 2,
                    startScale: 3,
                    endScale: 1.5,
                    minimumParticleLife: 1.5,
                    maximumParticleLife: 1.8,
                    minimumSpeed: 7,
                    maximumSpeed: 9,
                    imageSize: new Cesium[e(194)](2, 2),
                    emissionRate: 200,
                    lifetime: 16,
                    loop: !0,
                    emitter: new Cesium[e(1575)](Cesium[e(475)][e(1149)](45)),
                    sizeInMeters: !0
                  })
                }
              })(e)
            case Q[i(726)]:
              return new (class extends ri {
                constructor(e) {
                  const i = t
                  super(e),
                    (this._typeName = i(417)),
                    (this[i(2336)] = Q[i(726)]),
                    (this[i(1498)] = new Cesium[i(310)]()),
                    (this._style[i(1286)] = Cesium.defaultValue(this._style[i(1286)], -3.5)),
                    (this[i(343)][i(507)] = Cesium[i(1960)](this[i(343)][i(507)], 3)),
                    (this[i(343)][i(1080)] = Cesium[i(1960)](
                      this._style[i(1080)],
                      new Cesium.Cartesian2(2, 2)
                    )),
                    (this[i(343)][i(1116)] = Cesium.defaultValue(this[i(343)].lifetime, 16)),
                    (this[i(343)][i(370)] = Cesium[i(1960)](this._style[i(370)], !0)),
                    (this[i(343)][i(2177)] = Cesium.defaultValue(this[i(343)][i(2177)], !0)),
                    (this[i(343)][i(1509)] = Cesium[i(1960)](this[i(343)][i(1509)], 40)),
                    (this[i(343)][i(2104)] = Cesium[i(1960)](this[i(343)][i(2104)], 6)),
                    (this[i(343)][i(1211)] = Cesium[i(1960)](this[i(343)][i(1211)], 7)),
                    (this[i(343)][i(1278)] = Cesium[i(1960)](this[i(343)][i(1278)], 9)),
                    (this._style.maximumSpeed = Cesium[i(1960)](this._style[i(1854)], 9.5)),
                    (this[i(343)][i(507)] = Cesium[i(1960)](this._style[i(507)], 1)),
                    (this._style[i(119)] = Cesium.defaultValue(this[i(343)][i(119)], 7)),
                    (this[i(343)][i(195)] = Cesium.defaultValue(this[i(343)][i(195)], 1)),
                    (this[i(343)][i(1116)] = Cesium[i(1960)](this._style.lifetime, 16)),
                    this[i(2231)]()
                }
                [t(1084)]() {
                  const e = t
                  return new Cesium[e(2145)]({
                    image: this._style[e(259)],
                    startColor: new Cesium[e(1154)](1, 1, 1, 0.6),
                    endColor: new Cesium[e(1154)](0.8, 0.86, 1, 0.4),
                    imageSize: new Cesium.Cartesian2(1, 1),
                    lifetime: 16,
                    loop: !0,
                    emitter: new Cesium.CircleEmitter(0.2),
                    updateCallback: (t, i) => this[e(2178)](t, i),
                    sizeInMeters: !0,
                    performance: !1,
                    emissionRate: 60,
                    gravity: -5,
                    minimumParticleLife: 1,
                    maximumParticleLife: 2,
                    minimumSpeed: 9,
                    maximumSpeed: 9.5,
                    startScale: 2,
                    endScale: 4,
                    particleSize: 1,
                    heading: 110,
                    pitch: -45
                  })
                }
                [t(2178)](e, i) {
                  const s = t
                  Cesium[s(310)][s(379)](e[s(2251)], this[s(1498)]),
                    Cesium[s(310)][s(606)](
                      this._gravityScratch,
                      this[s(343)][s(1286)] * i,
                      this._gravityScratch
                    ),
                    (e.velocity = Cesium.Cartesian3[s(1861)](e[s(1009)], this[s(1498)], e[s(1009)]))
                }
              })(e)
          }
        }
      },
      hi = {
        readFeature(e) {
          const i = t
          let s = e[i(138)][i(2365)],
            n = e[i(1004)]
          return (n[i(2251)] = s), this[i(1951)](n)
        },
        create(e) {
          const i = t
          if (e[i(1158)] === Q[i(868)])
            return new (class extends he {
              constructor(e) {
                const i = t
                super(e),
                  (this[i(1115)] = j),
                  (this._graphicType = Q[i(868)]),
                  (this[i(528)] = '平面'),
                  (this._fixPointCount = 1),
                  (this[i(314)] = i(1865)),
                  (this[i(343)][i(1370)] = Cesium[i(1960)](this[i(343)].dimensions, {
                    x: 10,
                    y: 10
                  })),
                  (this[i(343)].outline = Cesium[i(1960)](this[i(343)].outline, !1)),
                  (this._style[i(376)] = Cesium[i(1960)](this[i(343)].outlineWidth, 1)),
                  (this[i(343)][i(725)] = Cesium[i(1960)](this[i(343)][i(725)], i(429))),
                  (this[i(987)] = this[i(2547)]())
                let s = Cesium.defaultValue(e[i(2251)], [111, 28, 0])
                this[i(1253)](s),
                  this[i(2231)](),
                  (this[i(2119)] = new Dt(e[i(2240)])),
                  this._label[i(1567)](this[i(987)][i(2240)])
              }
              get label() {
                return this._label
              }
              set [t(2240)](e) {
                const i = t
                ;(this._label = e), this[i(2119)][i(1567)](this._entity[i(2240)])
              }
              [t(415)](e) {
                this[t(987)].show = e
              }
              [t(2547)]() {
                const e = t,
                  i = new Cesium[e(2314)](Cesium[e(310)][e(1246)], 0),
                  s = {}
                s[e(868)] = i
                const n = {}
                return (
                  (n.graphicId = this[e(1570)]),
                  (n[e(2251)] = this._cartesian3),
                  (n[e(868)] = s),
                  (n.label = {}),
                  new Cesium.Entity(n)
                )
              }
              setPosition(e) {
                const i = t
                ;(this._position = e),
                  (this[i(2238)] = Cesium[i(310)][i(667)](e[0], e[1], e[2])),
                  (this[i(794)] = this[i(1912)]),
                  !this._editMode && (this._entity.position = this._cartesian3)
                let s =
                  this._style[i(1370)].x > this[i(343)][i(1370)].y
                    ? this[i(343)].dimensions.x
                    : this._style[i(1370)].y
                ;(this[i(1024)] = new Cesium[i(1242)](this[i(2238)], s)), this[i(2455)]()
              }
              [t(2455)]() {
                this[t(2462)]
              }
              [t(200)]() {
                const e = t
                let i = this[e(2238)]
                this[e(2520)] && (i = new Cesium[e(2569)]((t) => this[e(2238)], !1)),
                  (this[e(987)][e(2251)] = i)
              }
              _merge() {
                const e = t
                let i = { ...this[e(343)] }
                const s = this[e(987)].plane
                for (const t in i)
                  if (Object[e(782)].call(i, t)) {
                    const e = i[t]
                    null != e && null != e && (s[t] = e)
                  }
                s[e(725)] = Cesium[e(1154)][e(2008)](i[e(725)])
              }
              [t(2231)]() {
                const e = t
                this[e(548)](), (this[e(1020)] = this[e(987)][e(868)])
                let i =
                  this._style[e(1370)].x > this[e(343)][e(1370)].y
                    ? this[e(343)][e(1370)].x
                    : this[e(343)][e(1370)].y
                ;(this[e(1024)] = new Cesium[e(1242)](this._cartesian3, i)), this[e(2239)]()
              }
              [t(1545)](e) {
                const i = t
                ;(this._entity[i(2122)] = e.id),
                  (this[i(2462)] = e),
                  (this._entity = e[i(245)][i(118)][i(1861)](this._entity)),
                  this[i(2455)]()
              }
              [t(2389)](e) {
                const i = t
                e[i(245)][i(118)][i(1896)](this[i(987)])
              }
            })(e)
        }
      },
      li = {
        readFeature(e) {
          const i = t
          let s = e[i(138)][i(2365)],
            n = e[i(1004)]
          return (n[i(2251)] = s), this.create(n)
        },
        create(e) {
          const i = t
          switch (e[i(1158)]) {
            case Q.rectPyramid:
              return new tt(e)
            case Q[i(1625)]:
              return new gt(e)
            case Q[i(1270)]:
              return new pt(e)
            case Q[i(427)]:
              return new ot(e)
            case Q.probeRadar:
              return new rt(e)
            case Q[i(1650)]:
              return new nt(e)
          }
        }
      }
    class ci {
      constructor(e) {
        const i = t
        ;(this[i(2273)] = new Cesium[i(525)]()),
          (this[i(2390)] = void 0),
          (this[i(1056)] = void 0),
          (this.color = e[i(1070)]),
          (this[i(259)] = e[i(259)])
      }
    }
    const ui = {
        get: function () {
          return !1
        }
      },
      mi = {}
    ;(mi[t(902)] = function () {
      return this[t(2273)]
    }),
      Object.defineProperties(ci[t(1727)], {
        isConstant: ui,
        definitionChanged: mi,
        color: Cesium[t(1328)](t(1070))
      }),
      (ci.prototype[t(2578)] = function (e) {
        return t(1805)
      }),
      (ci.prototype[t(1474)] = function (e, i) {
        const s = t
        return (
          Cesium[s(2330)](i) || (i = {}),
          (i[s(1070)] = Cesium.Property[s(1659)](
            this[s(2390)],
            e,
            Cesium[s(1154)][s(588)],
            i[s(1070)]
          )),
          (i[s(259)] = this[s(259)]),
          i
        )
      }),
      (ci[t(1727)][t(1525)] = function (e) {
        const i = t
        return this === e || (e instanceof ci && Cesium[i(1169)][i(1525)](this[i(2390)], e._color))
      }),
      (Cesium[t(1637)][t(823)] = t(1805)),
      (Cesium[t(1637)][t(2146)] = t(1587)),
      Cesium[t(1637)][t(2257)][t(1371)](Cesium.Material[t(823)], {
        fabric: {
          type: Cesium[t(1637)][t(823)],
          uniforms: { color: new Cesium.Color(1, 0, 0, 0.5), image: '', time: 0 },
          source: Cesium[t(1637)].CylinderGlowGradientWallSource
        },
        translucent: function (t) {
          return !0
        }
      })
    class pi {
      constructor(e) {
        const i = t
        ;(this[i(2273)] = new Cesium.Event()),
          (this[i(2390)] = void 0),
          (this[i(1056)] = void 0),
          (this[i(1070)] = e[i(1070)]),
          (this[i(642)] = e[i(642)]),
          (this[i(1233)] = new Date()[i(1566)]()),
          (this[i(259)] = e[i(259)])
      }
    }
    const di = {}
    di[t(902)] = function () {
      return !1
    }
    const fi = {}
    ;(fi[t(902)] = function () {
      return this[t(2273)]
    }),
      Object[t(610)](pi[t(1727)], {
        isConstant: di,
        definitionChanged: fi,
        color: Cesium[t(1328)](t(1070))
      }),
      (pi[t(1727)][t(2578)] = function (e) {
        return t(1763)
      }),
      (pi.prototype[t(1474)] = function (e, i) {
        const s = t
        return (
          Cesium[s(2330)](i) || (i = {}),
          (i[s(1070)] = Cesium[s(1169)][s(1659)](
            this._color,
            e,
            Cesium[s(1154)].WHITE,
            i[s(1070)]
          )),
          (i[s(416)] = ((new Date()[s(1566)]() - this[s(1233)]) % this[s(642)]) / this[s(642)]),
          (i[s(259)] = this[s(259)]),
          i
        )
      }),
      (pi[t(1727)][t(1525)] = function (e) {
        const i = t
        return (
          this === e ||
          (e instanceof pi &&
            Cesium[i(1169)][i(1525)](this._color, e[i(2390)]) &&
            this[i(642)] == e[i(642)])
        )
      }),
      (Cesium.Material.CylinderGlowFlowWallType = t(1763)),
      (Cesium[t(1637)][t(1021)] = t(1776)),
      Cesium[t(1637)][t(2257)][t(1371)](Cesium[t(1637)][t(2242)], {
        fabric: {
          type: Cesium.Material[t(2242)],
          uniforms: { color: new Cesium[t(1154)](1, 0, 0, 0.5), image: '', time: 0 },
          source: Cesium[t(1637)][t(1021)]
        },
        translucent: function (t) {
          return !0
        }
      })
    class Ci {
      constructor(e) {
        const i = t
        ;(this._definitionChanged = new Cesium.Event()),
          (this[i(1070)] = void 0),
          (this[i(1056)] = void 0),
          (this[i(1070)] = e[i(1070)]),
          (this[i(259)] = e.image),
          (this.speed = e[i(2535)])
      }
    }
    const vi = {}
    vi[t(902)] = function () {
      return !1
    }
    const _i = {}
    ;(_i[t(902)] = function () {
      return this[t(2273)]
    }),
      Object.defineProperties(Ci[t(1727)], {
        isConstant: vi,
        definitionChanged: _i,
        color: Cesium[t(1328)](t(1070))
      }),
      (Ci[t(1727)].getType = function (e) {
        return t(173)
      }),
      (Ci[t(1727)][t(1474)] = function (e, i) {
        const s = t
        return (
          Cesium[s(2330)](i) || (i = {}),
          (i[s(1070)] = Cesium[s(1169)][s(1659)](
            this[s(1070)],
            e,
            Cesium[s(1154)][s(588)],
            i[s(1070)]
          )),
          (i[s(259)] = this[s(259)]),
          (i[s(2535)] = this[s(2535)]),
          i
        )
      }),
      (Ci.prototype.equals = function (e) {
        const i = t
        return (
          this === e ||
          (e instanceof Ci &&
            Cesium[i(1169)][i(1525)](this[i(1070)], e[i(1070)]) &&
            this[i(2535)] == e[i(2535)])
        )
      }),
      (Cesium[t(1637)][t(1840)] = t(173)),
      (Cesium[t(1637)][t(2206)] = t(643)),
      Cesium[t(1637)][t(2257)][t(1371)](Cesium.Material[t(1840)], {
        fabric: {
          type: Cesium[t(1637)][t(1840)],
          uniforms: { color: new Cesium[t(1154)](1, 0, 0, 0.5), image: '', speed: 0 },
          source: Cesium[t(1637)][t(2206)]
        },
        translucent: function (t) {
          return !0
        }
      })
    let gi = {
      readFeature(e) {
        const i = t
        let s = e[i(138)][i(2365)],
          n = e[i(1004)]
        return (n[i(2251)] = s), this.create(n)
      },
      create(e) {
        const i = t
        switch (e[i(1158)]) {
          case Q.glowCylinder:
            return new (class extends b {
              constructor(e = {}) {
                const i = t
                super(e),
                  (this[i(1912)] = e[i(2251)]),
                  (this[i(528)] = i(1702)),
                  (this[i(314)] = i(1865)),
                  (this[i(1115)] = U),
                  (this._graphicType = Q[i(203)]),
                  (this[i(343)][i(1070)] = Cesium.defaultValue(
                    e[i(1679)][i(1070)],
                    'rgba(255,0,0,1)'
                  )),
                  (this[i(343)][i(506)] = Cesium.defaultValue(e.style[i(506)], i(615))),
                  (this[i(343)][i(1981)] = Cesium[i(1960)](e[i(1679)][i(1981)], 100)),
                  (this._style[i(2306)] = Cesium[i(1960)](e[i(1679)][i(2306)], 1500)),
                  (this._style.texture = Cesium[i(1960)](e.style[i(1660)], '/')),
                  (this[i(343)][i(1513)] = Cesium.defaultValue(e[i(1679)][i(1513)], !0)),
                  (this._style[i(295)] = Cesium[i(1960)](e[i(1679)][i(295)], !0)),
                  (this[i(1516)] = this[i(343)][i(2306)]),
                  (this[i(920)] = this[i(343)].radius)
                let s = this[i(1496)] ? this[i(343)].selectedColor : this[i(343)][i(1070)]
                this[i(2390)] = Cesium[i(1154)].fromCssColorString(s)
                let n = Cesium.defaultValue(e.position, [111, 28, 0])
                ;(this[i(2238)] = Cesium[i(310)][i(667)](n[0], n[1], n[2])),
                  (this[i(794)] = this._position = n),
                  (this[i(1024)] = new Cesium[i(1242)](this[i(2238)], this[i(343)][i(2306)])),
                  (this[i(1244)] = []),
                  (this._isAdd = !1)
              }
              [t(415)](e) {
                const i = t
                this[i(1244)].forEach((t) => {
                  t[i(1482)] = e
                })
              }
              [t(2231)]() {
                const e = t
                let i = this[e(1496)] ? this[e(343)][e(506)] : this[e(343)][e(1070)]
                ;(i = Cesium.Color[e(2008)](i)),
                  (this[e(920)] = this[e(343)][e(1981)]),
                  (this[e(116)].ellipse[e(2471)] = this[e(343)][e(1981)]),
                  (this._bottomRing[e(1111)][e(2404)] = this._style.radius),
                  (this[e(1417)].ellipse.semiMinorAxis = this._style.radius),
                  (this[e(1417)][e(1111)][e(2404)] = this[e(343)].radius)
                let s = this._getCirclePosition(this[e(343)][e(1981)])
                ;(this[e(455)][e(2467)][e(2333)] = s),
                  (this[e(1970)][e(2467)].positions = s),
                  (this._height = this[e(343)][e(2306)])
                let n = new Array(s[e(277)])[e(2327)](this[e(1912)][2]),
                  o = new Array(s[e(277)])[e(2327)](this[e(1912)][2] + this[e(343)].height)
                ;(this[e(455)][e(2467)].minimumHeights = n),
                  (this[e(455)][e(2467)][e(1089)] = o),
                  (this[e(1970)].wall[e(2332)] = n),
                  (this._outerCylinder[e(2467)][e(1089)] = o),
                  (this[e(1970)][e(2467)].material.color = i),
                  (this[e(455)].wall[e(322)][e(1070)] = i),
                  (this[e(116)].ellipse[e(322)][e(1070)] = i),
                  (this[e(1417)].ellipse.material[e(1070)] = i),
                  (this[e(1417)].show = this[e(343)].showBottomCircle),
                  (this[e(1970)][e(1482)] = this[e(343)][e(295)]),
                  (this._boundingSphere = new Cesium[e(1242)](this[e(2238)], this._style.height))
              }
              setPosition(e) {
                const i = t
                ;(this[i(1912)] = e),
                  (this[i(2238)] = Cesium[i(310)][i(667)](e[0], e[1], e[2])),
                  (this[i(794)] = this._position)
                let s = this[i(122)](this._style[i(1981)])
                s[i(2553)](s[0]), (this[i(455)][i(2467)][i(2333)] = s)
                let n = new Array(s[i(277)])[i(2327)](this._position[2]),
                  o = new Array(s[i(277)])[i(2327)](this._position[2] + this[i(343)][i(2306)])
                ;(this[i(455)][i(2467)][i(2332)] = n),
                  (this[i(455)][i(2467)].maximumHeights = o),
                  (this._outerCylinder.wall[i(2333)] = s),
                  (this[i(1970)].wall[i(2332)] = n),
                  (this[i(1970)].wall.maximumHeights = o),
                  (this[i(1417)][i(2251)] = this[i(2238)]),
                  (this._bottomCircle[i(1111)][i(2306)] = e[2]),
                  (this[i(116)][i(1111)][i(2306)] = e[2]),
                  (this._bottomRing[i(2251)] = this._cartesian3),
                  (this[i(1024)] = new Cesium[i(1242)](this[i(2238)], this[i(343)][i(2306)]))
              }
              [t(1510)]() {
                const e = t
                let i = [],
                  s = this[e(461)]()
                return (
                  i.push(s),
                  (s = this._createBottomRing()),
                  i.push(s),
                  (s = this[e(1724)]()),
                  i.push(s),
                  (s = this.addOuterCylinder()),
                  i[e(2553)](s),
                  i
                )
              }
              [t(1724)]() {
                const e = t
                let i = this[e(122)](this[e(343)][e(1981)])
                return (
                  i[e(2553)](i[0]),
                  (this._innerCylinder = new Cesium[e(1380)]({
                    position: this._cartesian3,
                    wall: {
                      positions: i,
                      minimumHeights: new Array(i[e(277)]).fill(this[e(1912)][2]),
                      maximumHeights: new Array(i[e(277)])[e(2327)](
                        this[e(1912)][2] + this[e(343)][e(2306)]
                      ),
                      material: new ci({
                        color: this[e(2390)],
                        image: this[e(343)][e(1660)] + '/042DDE230DC541238F3ECA5FE7DFC79E.png'
                      })
                    }
                  })),
                  this[e(455)]
                )
              }
              [t(1164)]() {
                const e = t
                let i = this[e(122)](this[e(343)][e(1981)])
                return (
                  i[e(2553)](i[0]),
                  (this[e(1970)] = new Cesium.Entity({
                    position: this._cartesian3,
                    show: this[e(343)][e(295)],
                    wall: {
                      positions: i,
                      minimumHeights: new Array(i[e(277)])[e(2327)](this[e(1912)][2]),
                      maximumHeights: new Array(i[e(277)])[e(2327)](
                        this[e(1912)][2] + this[e(343)].height
                      ),
                      material: new pi({
                        color: this[e(2390)],
                        duration: 1e3,
                        image: this[e(343)][e(1660)] + '/B3D3B3374F584346AB74A9AFA9127259.png'
                      })
                    }
                  })),
                  this[e(1970)]
                )
              }
              [t(122)](e) {
                const i = t
                e /= 1.6
                let s = [],
                  n = Cesium[i(2058)][i(1648)](this[i(2238)]),
                  o = (2 * Math.PI) / 120,
                  r = (2 * Math.PI * 270) / 360
                for (let t = 0; t < 120; t++) {
                  let a = r - o * t,
                    h = new Cesium[i(310)](Math.sin(a) * e, Math[i(1272)](a) * e, 0)
                  s.push(Cesium[i(2066)][i(2483)](n, h, new Cesium.Cartesian3()))
                }
                return s[i(2553)](s[0]), s
              }
              [t(461)]() {
                const e = t
                return (
                  (this[e(1417)] = new Cesium[e(1380)]({
                    position: this._cartesian3,
                    show: this[e(343)][e(1513)],
                    ellipse: {
                      semiMinorAxis: this[e(343)].radius,
                      semiMajorAxis: this[e(343)].radius,
                      height: this[e(1912)][2],
                      material: new Ci({
                        image: this[e(343)][e(1660)] + e(1217),
                        color: this[e(2390)],
                        speed: 10
                      })
                    }
                  })),
                  this[e(1417)]
                )
              }
              [t(392)]() {
                const e = t
                return (
                  (this[e(116)] = new Cesium[e(1380)]({
                    position: this[e(2238)],
                    ellipse: {
                      semiMinorAxis: this[e(343)].radius,
                      semiMajorAxis: this[e(343)][e(1981)],
                      height: this._position[2],
                      material: new Ci({
                        image: this[e(343)][e(1660)] + e(422),
                        color: this[e(2390)],
                        speed: 0
                      })
                    }
                  })),
                  this._bottomRing
                )
              }
              [t(1545)](e) {
                const i = t
                ;(this[i(1244)] = this[i(1510)]()),
                  this[i(1244)][i(1602)]((t) => {
                    const s = i
                    e[s(245)][s(118)].add(t), (t.layerId = e.id), (t[s(2366)] = this[s(1570)])
                  }),
                  (this[i(2462)] = e),
                  (this[i(2549)] = !0)
              }
              _removeHook(e) {
                const i = t
                this._entities[i(1602)]((t) => {
                  const s = i
                  e[s(245)][s(118)].remove(t)
                }),
                  (this[i(2549)] = !1),
                  (this[i(1244)] = [])
              }
            })(e)
          case Q[i(1458)]:
            return new (class extends b {
              constructor(e = {}) {
                const i = t
                super(e),
                  (this[i(1912)] = e[i(2251)]),
                  (this[i(528)] = i(851)),
                  (this[i(314)] = i(1865)),
                  (this._graphicClassType = U),
                  (this[i(2336)] = Q[i(1458)]),
                  (this._style[i(1070)] = Cesium[i(1960)](e[i(1679)].color, 'rgba(255,0,0,1)')),
                  (this[i(343)].selectedColor = Cesium.defaultValue(e[i(1679)][i(506)], i(615))),
                  (this._style.radius = Cesium.defaultValue(e[i(1679)].radius, 100)),
                  (this[i(343)][i(2306)] = Cesium[i(1960)](e[i(1679)][i(2306)], 1500))
                let s = this[i(1496)] ? this._style[i(506)] : this[i(343)][i(1070)]
                this[i(2390)] = Cesium[i(1154)][i(2008)](s)
                let n = Cesium[i(1960)](e[i(2251)], [111, 28, 0])
                this[i(1253)](n), (this[i(1955)] = []), (this[i(2549)] = !1)
              }
              _setVisible(e) {
                const i = t
                this[i(1955)][i(1602)]((t) => {
                  t.show = e
                })
              }
              _setStyle() {
                const e = t
                let i = this[e(1496)] ? this[e(343)][e(506)] : this._style[e(1070)]
                ;(this[e(2390)] = Cesium.Color[e(2008)](i)),
                  (this[e(1024)] = new Cesium[e(1242)](this[e(2238)], this[e(343)][e(2306)])),
                  this._isAdd && (this[e(2389)](this[e(2462)]), this[e(1545)](this[e(2462)]))
              }
              setPosition(e) {
                const i = t
                ;(this._position = e),
                  (this[i(2238)] = Cesium.Cartesian3.fromDegrees(e[0], e[1], e[2])),
                  (this._coordinates = this[i(1912)]),
                  this[i(2231)]()
              }
              _createPrimitives() {
                const e = t
                let i = [],
                  s = this[e(392)]()
                return (
                  i[e(2553)](s),
                  (s = this._createBottomCircle()),
                  i.push(s),
                  (s = this.addInnereCylinder()),
                  i[e(2553)](s),
                  (s = this[e(1164)]()),
                  i[e(2553)](s),
                  i
                )
              }
              [t(1724)]() {
                const e = t
                let i = this[e(122)](0.5 * this[e(343)].radius * 0.7),
                  s = this._getCirclePosition(1),
                  n = this._color,
                  o = this[e(157)](i, s, n),
                  r = new Cesium[e(2305)]({
                    material: new Cesium[e(1637)]({
                      fabric: { uniforms: { u_color: n }, source: e(313) }
                    })
                  })
                const a = {}
                return (
                  (a.geometryInstances = o),
                  (a[e(1994)] = r),
                  (a[e(644)] = !1),
                  new Cesium[e(2186)](a)
                )
              }
              [t(1164)]() {
                const e = t
                let i = this[e(122)](0.5 * this._style[e(1981)]),
                  s = this[e(122)](1),
                  n = this[e(2390)],
                  o = this[e(157)](i, s, n),
                  r = new Cesium.Primitive({ geometryInstances: o, asynchronous: !1 })
                if (this._outImage) this[e(2561)](r, this._outImage)
                else {
                  let t = new Image()
                  ;(t[e(1744)] =
                    'data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAACAAAAEACAYAAADSoXR2AAAAGXRFWHRTb2Z0d2FyZQBBZG9iZSBJbWFnZVJlYWR5ccllPAAAAyZpVFh0WE1MOmNvbS5hZG9iZS54bXAAAAAAADw/eHBhY2tldCBiZWdpbj0i77u/IiBpZD0iVzVNME1wQ2VoaUh6cmVTek5UY3prYzlkIj8+IDx4OnhtcG1ldGEgeG1sbnM6eD0iYWRvYmU6bnM6bWV0YS8iIHg6eG1wdGs9IkFkb2JlIFhNUCBDb3JlIDUuNi1jMTQ1IDc5LjE2MzQ5OSwgMjAxOC8wOC8xMy0xNjo0MDoyMiAgICAgICAgIj4gPHJkZjpSREYgeG1sbnM6cmRmPSJodHRwOi8vd3d3LnczLm9yZy8xOTk5LzAyLzIyLXJkZi1zeW50YXgtbnMjIj4gPHJkZjpEZXNjcmlwdGlvbiByZGY6YWJvdXQ9IiIgeG1sbnM6eG1wPSJodHRwOi8vbnMuYWRvYmUuY29tL3hhcC8xLjAvIiB4bWxuczp4bXBNTT0iaHR0cDovL25zLmFkb2JlLmNvbS94YXAvMS4wL21tLyIgeG1sbnM6c3RSZWY9Imh0dHA6Ly9ucy5hZG9iZS5jb20veGFwLzEuMC9zVHlwZS9SZXNvdXJjZVJlZiMiIHhtcDpDcmVhdG9yVG9vbD0iQWRvYmUgUGhvdG9zaG9wIENDIDIwMTkgKFdpbmRvd3MpIiB4bXBNTTpJbnN0YW5jZUlEPSJ4bXAuaWlkOjExQTg0NDEyMDEzQjExRUFBNDhBRjhGMUMzOUUyNTU0IiB4bXBNTTpEb2N1bWVudElEPSJ4bXAuZGlkOjExQTg0NDEzMDEzQjExRUFBNDhBRjhGMUMzOUUyNTU0Ij4gPHhtcE1NOkRlcml2ZWRGcm9tIHN0UmVmOmluc3RhbmNlSUQ9InhtcC5paWQ6MTFBODQ0MTAwMTNCMTFFQUE0OEFGOEYxQzM5RTI1NTQiIHN0UmVmOmRvY3VtZW50SUQ9InhtcC5kaWQ6MTFBODQ0MTEwMTNCMTFFQUE0OEFGOEYxQzM5RTI1NTQiLz4gPC9yZGY6RGVzY3JpcHRpb24+IDwvcmRmOlJERj4gPC94OnhtcG1ldGE+IDw/eHBhY2tldCBlbmQ9InIiPz41vRwAAAAE90lEQVR42uydyW4UMRCG3T2dgYSAEGs4sp44cCJBcGUJbwCvALwWPAI8ABwAiUVwgLBdkEikJEiAGMhkZqhfU1aa1sy0g+yaJPyWSupOpPjz0uVyucrJer2eG2fJ3ZgLAQhAgC0PgN8XIlkqgGLE75oih0WmRVZEvop0rHog18rnRe6IzInsthyCXFt+TuSKyGmRXZZDgK5eFrkv8l7kiUhrxN/JSo3pigSvcNmI1bCh3b5LK2+NqHyvyEF9x3z5HgqRRViOAXhW5JrCoNdei/y20gMYxhmRSyIXRY6JTMT4DENLW+SdyD19x/NajDmwmYIW79Hnn+MA4GJEAAIQgABJbMJB+n5Sl9zWZvR9DABUfkJkVt8fi3zUldAEoKl24Y2S1fPZEmBdZFHkkb4vxTLRQ5djbyUf0ncYrD/UADUB8MZno2Q19yx7gIqIAAQgAAEIQIAtYZT+85LvRnjQCoPKR3rQUgN4Y3ZeK30g8qps0qeeAzDn4TWD9+ySG+BBS90DaGnVg9a2tgnR4il9/jkOACoiAhCAAASIZpAAFl6yulO0JACo/IDrn6CedP1zxOciv6wAYNnAP3RVBeWtJQCWzG8iL/X9g4vkqNzMctzUXtintt2KM/aS+bmQaY90rb8CF6tSKiICEIAABCAAAQiwbfcFIWVCbUcc9bZDLKiYAKj4uOv7Bb+ILLj+8W7PCgC7Jhzv44Qdh9yrru8V61jPgd645gD2io+15RiC5ZA5ENNPmGmDmtrta9YAVEQE2FlrQUM/J6efU8cSAD87ojodBaccS6kgiiErGiq/qe93XT+3oGU5BJnbyKrIUs6BQZoQUEcrQ7CYagiGqeKiMgnXLXsg5uJUm/RQJKy8fGS7PMw6SgVQPrJFGZr0kCfUL/7IFjIzrLGpemBQ0kPbehL6I9vMjUh6oEVEAAIQgAAEIAABdpyj0qcK71YrqGUJANNrv+sHOcAiRpDDM1cTY1BEHk7kqyP/2Ac5vLEEQEG8aDnIoTb5ObZRWg5ywG5o1dX4ClNYxT7IISgJPsXGpLut9MCWU0S5bqkmdAa3LQEwcZDAfkZ3tthQfkwNkVeeD6kiuaUabcp6DpS9YyalrAdQ8bT1EFQVUabarND9fNuyB2gREYAABCAAAQhAAAIQoM770UgBXARCIr/koDogVlzgrWuxemBKd0m4JnBed0+5ZQ9gn3BK5LK+P4y5ewoB8DFiKO91COLFfwVuTCZV1hSoYw1APZCsFIkbN1Haaa87wwiKQe6eT27AEX6esPVVd8/kOOdAz/oz9FE0p+uGIKUeqE5C8xgSKiICEIAABCAAAQhAAAIQgAAEIAABCECA6AA+sbkREyDUUYlKD6ggZwzRkr8tewBhusisv+0in5qE/hFchoZTE0TL4p8sTbtIpyahQ4Ag5fKpSVBmfdDECvQTZjoM0U9N6KgkwH9xY7PXpPiCmq5yuaLVldHIO7jgNi5XfOEhCqNh9udHV/RnC5YAUDTftdV4/ivvwEoRVS9XXPWa1FIT5ird8jpSB+BDN3rO8AaGMnGy0I0QRYTvdk6NkOsucuhGCMAgI8Q0isaHbmAevNMhiD4P6iZhstANGiQE2PEGCRo2NcgAsQDwKnxWFdlT17/duWU1BJmuHz6A5bwbEsCScgh+qAHit3Jr1oooOPExdRBLbeJjSpuwuy30AAH+CDAAPH5ltESNYl4AAAAASUVORK5CYII='),
                    (t.onload = () => {
                      const i = e
                      let s = document[i(1945)](i(1493))
                      ;(s[i(575)] = 64), (s.height = 256)
                      let n = s[i(1019)]('2d')
                      n[i(2264)](0, 0, 64, 256),
                        n[i(723)](t, 0, 0),
                        n[i(723)](t, 33, 0),
                        (this[i(489)] = s.toDataURL()),
                        this._setOuterCylinderAppearence(r, this._outImage)
                    })
                }
                return r
              }
              [t(2561)](e, i) {
                const s = t,
                  n = {}
                ;(n[s(202)] = this[s(2390)]), (n[s(259)] = i)
                const o = {}
                ;(o[s(1285)] = n),
                  (o[s(660)] =
                    'uniform vec4 u_color;\nczm_material czm_getMaterial(czm_materialInput materialInput){\n    czm_material material = czm_getDefaultMaterial(materialInput);\n    vec2 st = materialInput.st;\n    float time = fract(czm_frameNumber / 90.) ;\n    vec2 new_st = fract(st-vec2(time,time));\n    vec4 color = texture(image,new_st);\n\n    vec3 diffuse = color.rgb;\n    float alpha = color.a;\n    diffuse *= u_color.rgb;\n    alpha *= u_color.a;\n    alpha *= u_color.a;\n    material.diffuse = diffuse;\n    material.alpha = alpha * pow(1. - st.t,u_color.a);\n    return material;\n}')
                const r = {}
                ;(r[s(1061)] = o),
                  (e.appearance = new Cesium[s(2305)]({ material: new Cesium.Material(r) }))
              }
              _createCylinderInstance(e, i, s) {
                const n = t
                let o = this[n(343)].height,
                  r = e[n(123)](),
                  a = e[n(277)],
                  h = 2 * a,
                  l = [],
                  c = 1 / (a - 1),
                  u = [],
                  m = []
                for (let t = 0; t < a; t++) {
                  let e = Cesium[n(2285)][n(2579)](i[t]),
                    s = Cesium.Cartesian3.fromRadians(e[n(2106)], e[n(199)], e[n(2306)] + o)
                  m[n(2553)](s), l.push(t * c), l[n(2553)](0)
                  let r = (t + 1) % a,
                    p = h - (t + 1)
                  u.push[n(560)](u, [p - 1, p, t]), u[n(2553)].apply(u, [t, r, p - 1])
                }
                for (let t = 0; t < m.length; t++)
                  r[n(2553)](m[a - t - 1]), l[n(2553)](1 - t * c), l.push(1)
                let p = new Cesium[n(1765)]({
                  polygonHierarchy: new Cesium.PolygonHierarchy(r),
                  perPositionHeight: !0
                })
                return (
                  ((p = Cesium[n(1765)].createGeometry(p))[n(495)] = u),
                  (p[n(1548)].st[n(736)] = l),
                  new Cesium.GeometryInstance({
                    geometry: p,
                    attributes: { color: Cesium[n(1127)].fromColor(s) }
                  })
                )
              }
              _createBottomRing() {
                const e = t
                let i = this[e(792)](this[e(343)][e(1981)]),
                  s = this._createRingCanvas()
                const n = {}
                ;(n[e(202)] = this[e(2390)]), (n[e(259)] = s)
                const o = {}
                ;(o[e(1285)] = n), (o[e(660)] = e(1907))
                const r = {}
                return (
                  (r[e(1061)] = o),
                  new Cesium[e(2186)]({
                    geometryInstances: i,
                    appearance: new Cesium.MaterialAppearance({ material: new Cesium.Material(r) }),
                    asynchronous: !1
                  })
                )
              }
              [t(499)]() {
                const e = t
                if (this[e(1207)]) return this[e(1207)]
                let i = document.createElement('canvas')
                ;(i[e(575)] = 512), (i[e(2306)] = 512)
                let s = i[e(1019)]('2d'),
                  n = s[e(980)](256, 256, 0, 256, 256, 256)
                return (
                  n.addColorStop(0.1, 'rgba(255, 255, 255, 1.0)'),
                  n.addColorStop(0.2, e(744)),
                  n.addColorStop(0.3, e(1448)),
                  n[e(803)](0.5, e(744)),
                  n.addColorStop(0.9, e(2522)),
                  n.addColorStop(1, e(882)),
                  s.clearRect(0, 0, 512, 512),
                  s.beginPath(),
                  s[e(1855)](256, 256, 256, 0, 2 * Math.PI, !0),
                  (s[e(537)] = n),
                  s[e(2327)](),
                  s[e(555)](),
                  (this._circleImage = i),
                  i
                )
              }
              _createCircleInstance() {
                const e = t,
                  i = new Cesium[e(2320)]({
                    center: Cesium[e(310)][e(667)](
                      this[e(1912)][0],
                      this[e(1912)][1],
                      this[e(1912)][2]
                    ),
                    semiMajorAxis: this[e(343)][e(1981)],
                    semiMinorAxis: this[e(343)].radius,
                    height: this[e(1912)][2]
                  }),
                  s = Cesium[e(2320)][e(933)](i),
                  n = {}
                return (n[e(138)] = s), new Cesium[e(2155)](n)
              }
              [t(461)]() {
                const e = t
                let i = this[e(792)](this[e(343)][e(1981)]),
                  s = this[e(996)]()
                const n = {}
                ;(n[e(202)] = this[e(2390)]), (n[e(259)] = s)
                const o = {}
                ;(o.uniforms = n), (o[e(660)] = e(1631))
                const r = {}
                return (
                  (r[e(1061)] = o),
                  new Cesium[e(2186)]({
                    geometryInstances: i,
                    appearance: new Cesium[e(2305)]({ material: new Cesium.Material(r) })
                  })
                )
              }
              _createRingCanvas() {
                const e = t
                if (this._ringImage) return this[e(1007)]
                let i = document.createElement('canvas')
                ;(i[e(575)] = 512), (i[e(2306)] = 512)
                let s = i[e(1019)]('2d')
                return (
                  (s.fillStyle = e(635)),
                  (s[e(1354)] = e(1049)),
                  s[e(693)]([50, 50]),
                  (s[e(839)] = 30),
                  s[e(1735)](),
                  s[e(1855)](256, 256, 150, 0, 2 * Math.PI, !0),
                  s[e(189)](),
                  s[e(555)](),
                  (this[e(1007)] = i),
                  i[e(1184)]()
                )
              }
              [t(122)](e) {
                const i = t
                let s = [],
                  n = Cesium[i(2058)][i(1648)](this[i(2238)]),
                  o = (2 * Math.PI) / 120,
                  r = (2 * Math.PI * 270) / 360
                for (let t = 0; t < 120; t++) {
                  let a = r - o * t,
                    h = new Cesium[i(310)](Math[i(884)](a) * e, Math[i(1272)](a) * e, 0)
                  s[i(2553)](Cesium[i(2066)][i(2483)](n, h, new Cesium[i(310)]()))
                }
                return s[i(2553)](s[0]), s
              }
              [t(1545)](e) {
                const i = t
                ;(this[i(1955)] = this[i(768)]()),
                  this[i(1955)][i(1602)]((t) => {
                    const s = i
                    e[s(245)][s(696)].primitives[s(1861)](t),
                      (t[s(2122)] = e.id),
                      (t[s(2366)] = this[s(1570)])
                  }),
                  (this._layer = e),
                  (this._isAdd = !0)
              }
              [t(2389)](e) {
                const i = t
                this[i(1955)][i(1602)]((t) => {
                  const s = i
                  e[s(245)][s(696)][s(1346)][s(1896)](t)
                }),
                  (this._isAdd = !1),
                  (this[i(1955)] = [])
              }
            })(e)
        }
      }
    }
    const yi = 1,
      wi = 2,
      xi = 3,
      bi = 1,
      Si = 2
    class Pi {
      constructor(e) {
        const i = t
        ;(this._viewer = e[i(395)]),
          (this[i(406)] = e.shaderType || yi),
          (this.cullFaceType = e[i(2489)] || xi),
          (this[i(592)] = e[i(1929)] || null),
          (this[i(968)] = e[i(2535)] || 1),
          (this[i(1233)] = e[i(416)] || 1),
          (this[i(2390)] = e[i(1070)] || new Cesium.Color(1, 1, 1, 1)),
          (this[i(1227)] = e[i(1529)] || 1),
          (this[i(1630)] = e[i(2359)] || new Cesium.Cartesian4(1, 1, 1, 1)),
          (this[i(2236)] = e[i(1996)] || new Cesium[i(2141)](1, 1, 1, 1)),
          (this[i(1720)] = e[i(1720)] || null),
          (this[i(995)] = e.vertexShader),
          (this[i(2458)] = e.fragmentShader),
          (this.images = e[i(1434)] || []),
          (this[i(432)] = e[i(432)] || []),
          (this[i(1406)] =
            e[i(664)] ||
            new Cesium[i(2203)]({
              wrapS: Cesium[i(1576)][i(1558)],
              wrapT: Cesium[i(1576)][i(1558)],
              minificationFilter: Cesium[i(1576)][i(2523)],
              magnificationFilter: Cesium[i(1576)][i(2523)]
            })),
          (this[i(1473)] = e.depthTest || !0),
          (this.depthMask = e[i(898)] || !0),
          (this._renderStateOptions = {}),
          this.init()
      }
      init() {
        const e = t
        let i = this,
          s = []
        0 === this[e(432)][e(277)] &&
          (this[e(1434)].forEach(function (t) {
            const i = e
            s.push(Cesium[i(1283)][i(1030)](t).fetchImage())
          }),
          Promise.all(s).then((t) => {
            i.textures = t.map(function (t) {
              const e = a0_0x3b79,
                s = {}
              return (
                (s[e(2325)] = i[e(245)][e(696)].context),
                (s[e(660)] = t),
                (s[e(664)] = i[e(1406)]),
                new Cesium[e(1641)](s)
              )
            })
          })),
          (this[e(2458)] = this[e(2267)]()),
          (this[e(1720)] = this[e(2564)]())
      }
      [t(2267)]() {
        const e = t
        let i,
          s = this
        switch (this[e(406)]) {
          case yi:
            return this[e(2458)]
          case wi:
            return (
              (i = i =
                (i = (i = Cesium[e(231)][e(1152)](s[e(2458)], e(2005)))[e(2225)](
                  new RegExp('fragColor', 'gm'),
                  e(318)
                ))[e(2225)](new RegExp('fragCoord', 'gm'), 'gl_FragCoord')),
              (s[e(592)] = s._outFragColor || e(1681)),
              '\n                       uniform float u_speed;\n                       uniform vec3 u_color;\n                       uniform float u_time;\n                       uniform float u_glow;\n                       uniform vec4 u_param;\n                       uniform vec4 iDate; // 日期（年，月，日，时）\n                       uniform float iTimeDelta; // 渲染时间，单位秒\n                       uniform vec4 iMouse; // xy分量取得当鼠标浮动的坐标位置，zw获得鼠标点击时的坐标位置。\n                       uniform sampler2D iChannel0; // 获得第一个四个通道对应的输入。\n                       uniform sampler2D iChannel1; // 获得第二个四个通道对应的输入。\n                       uniform sampler2D iChannel2; // 获得第三个四个通道对应的输入。\n                       uniform float iChannelTime[4]; // 针对于我们四个输入通道，每个通道的回放时间，以秒计时。\n                       uniform vec3 iChannelResolution[4]; // 四个输入通道的分辨率，以像素为单位。\n                       uniform sampler2D colorTexture; \n                       in vec2 v_st; \n                       vec2 iResolution;\n                       float iTime;  \n                       ' +
                i +
                e(1128) +
                s[e(592)] +
                e(278)
            )
        }
      }
      createUniformMap() {
        const e = t
        let i = this
        if (this[e(1720)]) return this[e(1720)]
        switch (this[e(406)]) {
          case yi:
            return {
              u_color: function () {
                return i[e(2390)]
              },
              u_speed: function () {
                return i[e(968)]
              },
              u_time: function () {
                return i[e(1233)]
              },
              u_glow: function () {
                return i[e(1227)]
              },
              u_param: function () {
                return i[e(1630)]
              },
              image: function () {
                const t = e
                return Cesium[t(2330)](i[t(432)][0])
                  ? i.textures[0]
                  : i[t(245)][t(696)][t(2325)][t(1627)]
              }
            }
          case wi:
            return {
              u_color: function () {
                return i[e(2390)]
              },
              u_speed: function () {
                return i._speed
              },
              u_time: function () {
                return i[e(1233)]
              },
              u_glow: function () {
                return i[e(1227)]
              },
              u_param: function () {
                return i[e(1630)]
              },
              iChannel0: function () {
                const t = e
                return Cesium.defined(i[t(432)][0])
                  ? i[t(432)][0]
                  : i[t(245)][t(696)].context.defaultTexture
              },
              iChannel1: function () {
                const t = e
                return Cesium.defined(i[t(432)][1])
                  ? i[t(432)][1]
                  : i[t(245)][t(696)][t(2325)][t(1627)]
              },
              iChannel2: function () {
                const t = e
                return Cesium[t(2330)](i[t(432)][2])
                  ? i[t(432)][2]
                  : i[t(245)][t(696)].context.defaultTexture
              },
              iDate: function () {
                const t = e
                let i = new Date()
                return { x: i[t(1290)]() + 1, y: i[t(2517)](), z: i[t(1503)](), w: i.getSeconds() }
              },
              iMouse: function () {
                return i[e(2236)]
              },
              iChannelResolution: function () {
                const t = { x: 100, y: 100, z: 0, w: 0 },
                  e = { x: 100, y: 100, z: 0, w: 0 }
                return [t, e]
              }
            }
        }
      }
    }
    const Mi = {}
    ;(Mi.box = t(843)), (Mi[t(1565)] = t(1565)), (Mi.plane = t(868))
    const Ai = Mi
    class Ti {
      constructor(e) {
        const i = t
        ;(this[i(568)] = e), (this[i(1265)] = this[i(933)]())
      }
      [t(933)]() {}
      get [t(892)]() {
        return this[t(1265)]
      }
    }
    class Ei extends Ti {
      constructor(t = {}) {
        super(t)
      }
      [t(933)]() {
        const e = t
        return Cesium[e(1607)][e(933)](
          new Cesium[e(1607)]({ vertexFormat: Cesium.VertexFormat[e(2050)] })
        )
      }
    }
    const zi = {}
    ;(zi[t(214)] = t(2276)),
      (zi[t(1579)] = t(729)),
      (zi.RealFlameFS = t(1742)),
      (zi.GlowPyramidFS = t(2351)),
      (zi[t(2056)] = t(1527)),
      (zi[t(1833)] =
        '\nfloat glow = 0.0;\nfloat t;\nconst float EPSILON = 0.0001;\n\nmat3 makeRotationMatrix(vec3 a)\n{\n    return mat3(\n    cos(a.x) * cos(a.z) - sin(a.x) * cos(a.y) * sin(a.z),\n        -cos(a.x) * sin(a.z) - sin(a.x) * cos(a.y) * cos(a.z),\n        sin(a.x) * sin(a.y),\n        sin(a.x) * cos(a.z) + cos(a.x) * cos(a.y) * sin(a.z),\n        -sin(a.x) * sin(a.z) + cos(a.x) * cos(a.y) * cos(a.z),\n        -cos(a.x) * sin(a.y),\n        sin(a.y) * sin(a.z),\n        sin(a.y) * cos(a.z),\n        cos(a.y)\n    );\n}\n\n#define sph(p, r) (length(p) - r)\n#define cube(p, b) length(max(abs(p) - vec3(b), 0.))\n\n// http://iquilezles.org/www/articles/distfunctions/distfunctions.htm\nfloat sdCapsule( vec3 p, vec3 a, vec3 b, float r )\n{\n    vec3 pa = p - a, ba = b - a;\n    float h = clamp( dot(pa,ba)/dot(ba,ba), 0.0, 1.0 );\n    return length( pa - ba*h ) - r;\n}\n\nfloat capsules(vec3 q) {\n    vec3 a = vec3(0.05);\n    vec3 b = vec3(0.2);\n    float r = 0.05 * abs(t);\n    \n    float c1 = sdCapsule(q, a, b, r);\n    float c2 = sdCapsule(q, -a, -b, r);\n    float c3 = sdCapsule(vec3(-q.x, -q.y, q.z), a, b, r);\n    float c4 = sdCapsule(vec3(-q.x, -q.y, q.z), -a, -b, r);\n    \n    float c5 = sdCapsule(vec3(-q.x, q.y, q.z), a, b, r);\n    float c6 = sdCapsule(vec3(-q.x, q.y, q.z), -a, -b, r);\n    float c7 = sdCapsule(vec3(q.x, -q.y, q.z), a, b, r);\n    float c8 = sdCapsule(vec3(q.x, -q.y, q.z), -a, -b, r);\n    \n    return min(min(min(min(min(min(min(c1, c2), c3), c4), c5), c6), c7), c8);\n}\n\nfloat sceneSDF(vec3 p) {\n    mat3 rot = makeRotationMatrix(vec3(iTime, iTime * 0.3, iTime * 0.6));\n    vec3 q = rot * p.xyz;\n    t = sin(iTime * 4.);\nfloat d = max(sph(p, mix(.125, .15, t)), cube(q, .1)); \n    \n    d = max(-capsules(q), d);\n    \n    d = max(-sph(p.xyz, mix(.1, .125, t)), d);\n    \n\n// https://www.shadertoy.com/view/4t2yW1\nglow += 0.0001 / (.01 + d * d);\nreturn d;\n}\n\nvec3 estimateNormal(vec3 p) {\n    return normalize(vec3(\n        sceneSDF(vec3(p.x + EPSILON, p.y, p.z)) - sceneSDF(vec3(p.x - EPSILON, p.y, p.z)),\n        sceneSDF(vec3(p.x, p.y + EPSILON, p.z)) - sceneSDF(vec3(p.x, p.y - EPSILON, p.z)),\n        sceneSDF(vec3(p.x, p.y, p.z  + EPSILON)) - sceneSDF(vec3(p.x, p.y, p.z - EPSILON))\n    ));\n}\n\nvoid main()\n{\n// vec2 uv = ( fragCoord - .5*iResolution.xy ) / iResolution.y;\nvec2 uv = 2.*v_st.xy - vec2(1., 1.);\nvec3 ro = vec3(0., 0., 0.5), p;\nvec3 rd = normalize(vec3(uv, -1));\np = ro;\n\nfloat t = 0.;\n    float hit = 0.;\n//  vec3 normal = vec3(0.0);\nfor (float i = 0.; i < 2.0; i += .01) {\n    p = ro + rd * t;\n    float d = sceneSDF(p);\n        if (d < .001) {\n            hit = 1.0;\n            // normal = estimateNormal(p);\n            break;\n        }\n        \n    t += d * 0.2; // avoid clipping, enhance the glow\n        if (t > 10.) {\n            hit = 0.0;\n            break;\n        }\n}\n    \nvec3 c = vec3(glow, glow, sqrt(glow)); // normal;\nfragColor = vec4(c *u_glow, 1.);\n}\n'),
      (zi.FlameCloudFS = t(2370)),
      (zi[t(694)] = t(2234)),
      (zi[t(1564)] = t(2443)),
      (zi[t(273)] = t(1096)),
      (zi[t(671)] =
        '\nfloat rand(vec2 co){\n    // https://stackoverflow.com/questions/4200224/random-noise-functions-for-glsl\n    return fract(sin(dot(co.xy ,vec2(12.9898,78.233))) * 43758.5453);\n}\n\n/*\nfloat rand(float v){\n    return fract(sin(12.9898*v)*43758.5453);\n}\n*/\n\nfloat rand(float p)\n{\n    // Hash function by Dave Hoskins\n    // https://www.shadertoy.com/view/4djSRW\n    p = fract(p * .1031);\n    p *= p + 33.33;\n    p *= p + p;\n    return fract(p);\n}\n\n\nvec3 lastExplosion(float time)\n{\n    // vec3(time since last explosion,\n    //      index of last explosion,\n    //      time until next explosion)\n    float t = mod(time, 10.);\n    float interval = floor(time/10.);\n    float t0max = 0., imax=-1.;\n    float t0next = 10.;\n    for(float i=0.; i<10.; i++)\n    {\n        float t0 = rand(vec2(interval, i)) * 10.;\n        if(t > t0 && t0 > t0max)\n        {\n            t0max = t0;\n            imax = i;\n        }\n        if(t < t0 && t0 < t0next)\n        {\n            t0next = t0;\n        }\n    }\n    return vec3(t-t0max, 10.*interval+imax, t0next-t);\n}\nvec3 glow(vec2 p, vec2 lpos)\n{\n    vec2 q = p - lpos;\n    float atten = 1./dot(q,q);\n    //atten *= (1. + atten*1e-4); // Make the inside slightly sharper\n    return vec3(1.0) * atten;\n} \n\nvoid main()\n{\n    vec2 p =  2.*v_st.xy - vec2(1., 1.);//(2.*fragCoord-iResolution.xy)/iResolution.y;\n\n    vec3 col = vec3(0);\n    \n    vec3 lastExpl = lastExplosion(iTime);\n    float t = lastExpl.x, explNum = lastExpl.y, tFadeout = lastExpl.z;\n    \n    // Fireworks base color\n    vec3 baseCol =  0.4*sin(vec3(1.)*explNum+vec3(0.,2.1,-2.1));\n    \n    // Number of particles\n    float N_LIGHTS = 100.;\n    for(float i=0.; i<N_LIGHTS; i++)\n    {\n        \n        // Generate points uniformly on hemisphere\n        // (see Total Compendium eq. (34))\n        float f = i/N_LIGHTS;\n        float r = sqrt(1. - f*f);\n        float th = 2.*0.618033*3.14159*i; // Use Golden Ratio for a quasirandom sequence\n        float hash = sin(explNum+i*85412.243);\n        float weight = (1.-0.2*hash);\n        th += hash *3.* 6.28/N_LIGHTS;\n        // Only take x and y coordinates\n        vec2 lpos = vec2(cos(th), sin(th)) * r;\n        // Add some physics\n        lpos.xy *= (1.-exp(-3.*t/weight)) * weight; // explosion, easing out\n        lpos.y += t*0.3*weight - t*(1.-exp(-t*weight)) * 0.6 * weight; // vertical free-fall motion\n        float intensity = 2e-4;\n        intensity *= exp(-2.*t); // Fade out with time\n        intensity *= (1.-0.5*hash); // Randomize per particle\n        intensity *= (1.+10.*exp(-20.*t)); // Intensity burst at explosion\n        intensity *= clamp(3.*tFadeout, 0., 1.); // Fade out before next explosion\n        col += glow(p, lpos) * intensity * baseCol;\n    }\n    \n    \n    col = max(col, 0.);\n    //col = 1.-exp(-col); // Tone mapping\n    col = (col*(2.51*col+0.03))/(col*(2.43*col+0.59)+0.14); // Tone mapping\n    //col = col/(1.+col);\n    col = sqrt(col); // gamma correction\n    \n    fragColor = vec4(col,1.0);\n}\n'),
      (zi[t(663)] = t(1778)),
      (zi[t(1976)] =
        '\n        #define TWO_PI 6.28318530718\n        #define PI 3.14159265359\n\n        float noise(vec2 p) {\n            return texture(iChannel0, p * 0.05 ).x;\n        }\n\n        float fbm(vec2 p) {\n            float a =1.;\n            float f = 1.;\n            return a*noise(p) \n                + a*0.5 * noise(p*f*2. ) \n                + a*0.25 * noise(p*f*4. )\n                + a*0.1 * noise(p*f*8. );\n        }\n\n        float circle(vec2 p) {\n            float r = length(p);\n            float radius = 0.4;\n            float height = 1.; \n            float width = 150.;\n            \n          return height - pow(r - radius, 2.) *width ;\n        }\n\n        void main()\n        {\n            \n            vec2 uv = 2.0*v_st.xy - vec2(1., 1.);\n            uv /=2.05;\n            vec2 st  = vec2(\n                    atan(uv.y, uv.x) ,\n                    length(uv) * 1. + iTime * 0.1\n                );\n\n            \n            st.x += st.y * 1.1;// - iTime * 0.3;\n            st.x = mod(st.x , TWO_PI);\n          \n            \n            float n = fbm(st ) * 1.5 -1. ;\n            n = max(n, 0.1);\n            float circle =  max(1.- circle(uv ), 0.)  ;\n            \n            float color = n/circle;\n            float mask = smoothstep(0.48, 0.4, length(uv));\n            \n            color *= mask;\n            vec3 rez = vec3(1., 0.5, 0.25) * color;\n            // Output to screen\n            fragColor = vec4(rez * u_glow,1.0);\n        }\n    '),
      (zi[t(2428)] =
        '\n        float stepping(float t){\n            if(t<0.)return -1.+pow(1.+t,2.);\n            else return 1.-pow(1.-t,2.);\n        }\n        void main()\n        {\n            vec2 uv = 2.*v_st.xy - vec2(1., 1.);\n            uv *= 2.0;\n            // vec2 uv = (fragCoord*2.-iResolution.xy)/iResolution.y;\n            fragColor = vec4(0);\n            uv = normalize(uv) * length(uv);\n            for(int i=0;i<12;i++){\n                float t = iTime + float(i)*3.141592/12.*(5.+1.*stepping(sin(iTime*3.)));\n                vec2 p = vec2(cos(t),sin(t));\n                p *= cos(iTime + float(i)*3.141592*cos(iTime/8.));\n                vec3 col = cos(vec3(0,1,-1)*3.141592*2./3.+3.141925*(iTime/2.+float(i)/5.)) * 0.5 + 0.5;\n                fragColor += vec4(0.05/length(uv-p*0.9)*col,1.0);\n            }\n            fragColor.xyz = pow(fragColor.xyz ,vec3(3.));\n            fragColor.w = 1.0;\n        }\n    '),
      (zi.ShieldFS = t(2546)),
      (zi.GlowBoxFS = t(1590)),
      (zi[t(1223)] =
        '\n            precision highp float;\n            precision highp int;\n            uniform float u_speed;\n            uniform vec3 u_color;\n            uniform float u_time;\n            uniform float u_glow;\n            in vec3 v_positionEC;\n            in vec3 v_normalEC;\n            in vec2 v_st;\n            mat4 mat  = mat4 ( vec4 ( 1.0 , 0.0 , 0.0 , 0.0 ),\n                vec4 ( 0.0 , 1.0 , 0.0 , 0.0 ),\n                vec4 ( 0.0 , 0.0 , 1.0 , 0.0 ),\n                vec4 ( 0.0 , 0.0 , 0.0 , 1.0 ) );\n\n            vec2 pos;\n\n            vec4 col = vec4 ( 0., 0., 0., 1000.0 );\n            // vec4 col = vec4 ( .3, .1, 0.1, 1000.0 );\n            \n            void Line2 ( vec2 a, vec2 b );\n            void Line2 ( vec2 a, vec2 b ) {\n            float d = distance ( pos , a ) + distance ( pos , b ) - distance ( a , b ) + 1e-5;\n            col += max ( 1. - pow ( d * 14. , 0.1 ) , -0.01 );\n            }\n\n            void Line4 ( vec4 a, vec4 b );\n            void Line4 ( vec4 a, vec4 b ) {\n            a = mat * a;\n            a.xyz /= 1.5 + a.w * 2.;\n            b = mat * b;\n            b.xyz /= 1.5 + b.w * 2.;\n            Line2 ( a.xy , b.xy );\n            }\n\n            void Point ( vec4 p );\n            void Point ( vec4 p ) {\n            p = mat * p;\n            p.xyz /= 1.5 + p.w * 2.;\n            \n            float d = distance ( pos , p.xy );\n            \n            if ( d < .3 )\n            if ( p.z < col.a ) {\n                col.b += max ( 1.0 - pow ( d * 5.0 , .1 ) , 0.0 );\n            }\n            }\n\n            void Rotate ( float angle, float d1, float d2, float d3, float d4);\n            void Rotate ( float angle, float d1, float d2, float d3, float d4) {\n            float c = cos (angle), s = sin (angle);\n            mat *= mat4 ( vec4 (  c*d1+(1.-d1),  s * d2 * d1 , -s * d3 * d1 ,  s * d4 * d1 ),\n                    vec4 ( -s * d1 * d2 ,  c*d2+(1.-d2),  s * d3 * d2 , -s * d4 * d2 ),\n                    vec4 (  s * d1 * d3 , -s * d2 * d3 ,  c*d3+(1.-d3),  s * d4 * d3 ),\n                    vec4 ( -s * d1 * d4 ,  s * d2 * d4 , -s * d3 * d4 ,  c*d4+(1.-d4)) );\n            }\n\n            void main( void ) {\n\n            float time = czm_frameNumber / 60.0;\n            time = time * u_speed ;\n            pos = v_st - 0.5;\n            float pi = 3.141592;\n\n            Rotate ( pi / 3.,      0.0, 1.0, 1.0, 0.0 );\n            Rotate ( time,      1.0, 1.0, 0.0, 0.0 );\n            \n            vec4 point1 = vec4 ( 0., 0., .2, -.4 );\n            vec4 point2 = vec4 ( 0., 0., -.6, .4 );\n\n            vec4 point3 = vec4 ( -.2, 0., 0., -0.2 );\n            vec4 point4 = vec4 ( .2, 0., 0., -0.2 );\n            vec4 point5 = vec4 ( 0., .2, 0., -0.2 );\n            vec4 point6 = vec4 ( 0., -.2, 0., -0.2 );\n\n            Line4 ( point1, point3 );\n            Line4 ( point1, point4 );\n            Line4 ( point1, point5 );\n            Line4 ( point1, point6 );\n\n            Line4 ( point2, point3 );\n            Line4 ( point2, point4 );\n            Line4 ( point2, point5 );\n            Line4 ( point2, point6 );\n\n            Line4 ( point3, point5 );\n            Line4 ( point3, point6 );\n            Line4 ( point4, point5 );\n            Line4 ( point4, point6 );\n\n            Point ( point1 );\n            Point ( point2 );\n            Point ( point3 );\n            Point ( point4 );\n            Point ( point5 );\n            Point ( point6 );\n\n            out_FragColor = vec4( col.xyz * u_glow, 1.0 );\n            out_FragColor = vec4(out_FragColor.rgb + u_color, out_FragColor.b* out_FragColor.b);\n            }\n        '),
      (zi[t(1025)] = t(633)),
      (zi[t(129)] = t(404)),
      (zi[t(2502)] = t(1443)),
      (zi[t(801)] = t(1159)),
      (zi.SwimFS = t(1670)),
      (zi.DoubleLoopFS = t(593)),
      (zi[t(953)] =
        '\nprecision mediump float;\n\nuniform float time;\nuniform vec2 mouse;\nuniform vec2 resolution;\n\n#define MAX_ITER \t100\n#define MAX_DIST \t10.0\n#define EPS\t\t0.001\n#define PI\t\t3.14159\n#define RADIUS\t\t1.5\n\nconst vec3 ZERO\t\t= vec3(0.0);\nconst vec3 X \t\t= vec3(1.0,0.0,0.0);\nconst vec3 Y \t\t= vec3(0.0,1.0,0.0);\nconst vec3 Z \t\t= vec3(0.0,0.0,1.0);\n\nconst vec3 AMBIENT\t= vec3(0.2);\nconst vec3 DIFFUSE\t= vec3(1.0);\nconst vec3 SPECULAR\t= vec3(1.0);\nconst float SHINESS\t= 32.0;\n\nvec3 sphereColor\t= vec3(1.0);\nvec3 lightPosition\t= vec3(3.0,2.0,4.0);\nvec3 camPosition\t= vec3(0.0,0.0,3.0);\nvec3 camUp \t\t= vec3(0.0,1.0,0.0);\nvec3 camLookat\t\t= vec3(0.0);\n\nvec3 bgColor  \t\t= vec3(0.);\n\nmat2 mm2(in float deg){\n\tfloat c = cos(deg); \n\tfloat s = sin(deg);\n\treturn mat2(c,-s,s,c);\n}\n\nfloat rand(float x,float y){\n \treturn fract(sin(x*12.9898+y*78.233) * 43758.5453);\n}\n\nfloat interpolate(float a, float b, float x){\n  \tfloat ft = x * PI;\n  \tfloat f = (1.0 - cos(ft)) * 0.5;\n  \treturn a*(1.0-f) + b*f;\n}\n\nfloat fbm1(float r,float index){\n\tfloat x = floor(r);\n\tfloat y = fract(r);\n\tfloat a = rand(x,2.1*index);\n\tfloat b = rand(x+1.0,2.1*index);\n\treturn interpolate(a,b,y)-0.5;\n\t//return 0.0;\n}\n\nvec2 line(vec3 pos, vec3 from, vec3 to,float index){\n\tvec3 p = pos-from;\n\tfloat lp = length(p);\n\tfloat len = 1.0*length(to-from);\n\tvec3 dir = normalize(to-from);\n\tfloat k = clamp(dot(p,dir),0.0,len);\n\tfloat ins = smoothstep(0.12,.56,lp);\n    \tfloat outs = 0.22+smoothstep(.0,.022,abs(lp-len));\n    \tfloat id = ins*outs;\n\tfloat rr = 0.3*smoothstep(0.1,0.4,k)*(1.0-smoothstep(len-0.3,len,k));\n\tvec3 offset = rr * (Y*fbm1(1.*k,index)+Z*fbm1(1.*k,index+1.0));\n\treturn vec2(0.4*distance(p,from+k*dir+offset)-0.009/id,id);\n}\n\nvec2 map(vec3 p, float i){\n    \tvec3 en = X*RADIUS;\n\tfloat k = mod(i,6.0);\n\tp.yz *= mm2(PI/3.0*k);\n\tfloat j = step(1.0,i)+step(7.0,i)+step(13.0,i);\n\tp.xy *= mm2(PI/3.0*j);\n    \treturn line(p,ZERO,en,i);\n}\n\n\nvec2 raymarch(vec3 ro, vec3 rd, float index){\n\tvec2 dist = vec2(EPS);\n\tfloat k = 0.0;\n\tfor(int i = 0; i < MAX_ITER; i++){\n\t\tdist = map(ro+k*rd,index);\n\t\tif(abs(dist.x)<EPS || dist.x>MAX_DIST){\n\t\t\tbreak;\n\t\t} \n\t\tk += dist.x; \t\n\t}\n\treturn vec2(k,dist.y);\n}\n\nvec2 sphere(vec3 pos,vec3 center, float radius){\n\treturn vec2(distance(pos,center) - radius,1.0);\n}\n\nvec2 sphereDist(vec3 pos){\n\treturn sphere(pos, ZERO, RADIUS);\n}\n\nvec2 raymarchSphere(vec3 ro, vec3 rd){\n\tvec2 dist = vec2(EPS);\n\tfloat k = 0.0;\n\tfor(int i = 0; i < MAX_ITER; i++){\n\t\tdist = sphereDist(ro+k*rd);\n\t\tif(abs(dist.x)<EPS || dist.x>MAX_DIST){\n\t\t\tbreak;\n\t\t} \n\t\tk += dist.x; \t\n\t}\n\treturn vec2(k,dist.y);\n}\n\nvec3 setupCamera(vec3 pos, vec3 up, vec3 lookat, vec2 uv){\n\tvec3 camDir = normalize(lookat-pos);\n\tvec3 camUp = normalize(up);\n\tvec3 camRight = normalize(cross(camDir,camUp));\n\treturn normalize(camRight*uv.x + camUp*uv.y + camDir);\n}\n\nvec3 calcNormal(vec3 p, float d){\n\treturn normalize(vec3(\n\t\t(sphereDist(p+EPS*X) - sphereDist(p-EPS*X)).x,\n\t\t(sphereDist(p+EPS*Y) - sphereDist(p-EPS*Y)).x,\n\t\t(sphereDist(p+EPS*Z) - sphereDist(p-EPS*Z)).x\n\t    )\n\t);\n}\n\nvec3 calcLight(vec3 p, vec3 lightPos){\n\tvec3 dir = normalize(lightPos-p);\n\tvec3 nor = calcNormal(p,0.0);\n\tvec3 dif = DIFFUSE * max(dot(dir,nor),0.0);\n\tvec3 spe = 0.5*SPECULAR * pow(max(dot(reflect(-dir,nor),dir),0.0),SHINESS);\n\treturn AMBIENT + dif + spe;\n}\n\nvec3 calcDpDx(vec3 ro,vec3 rd,vec3 rdx,vec3 rdy, vec3 nor, float t){\n\treturn t*(rdx*dot(rd,nor)/dot(rdx,nor) - rd);\n}\n\nvec3 calcDpDy(vec3 ro,vec3 rd,vec3 rdx,vec3 rdy, vec3 nor, float t){\n\treturn t*(rdy*dot(rd,nor)/dot(rdy,nor) - rd);\n}\n\nvec3 textureMapping(vec3 pos,vec3 dposdx,vec3 dposdy,float oid ){\n\treturn vec3(1.0);\n}\n\nin vec2 v_st;\nvoid main( void ) {\n    float time = czm_frameNumber /60.0;\n    vec2 resolution=czm_viewport.zw;\n    vec2  st=   v_st * 2.0-1.0;\n\tvec2 uv =st; //(2.0 * gl_FragCoord.xy - resolution.xy ) / resolution.y;\n\tvec3 ro = camPosition;\n\tmat2 mm = mm2(time);\n       \tro.xy *= mm;\n        ro.xz *= mm;\n\t//ro.yz *= mm;\n\tvec3 rd = setupCamera(ro, camUp, camLookat, uv);\n\tvec3 color = bgColor;\n\tfor(float i = 0.0 ; i < 14.0; i++){\n\t\tvec2 rz = raymarch(ro,rd,i);\n\t\tif(rz.x > MAX_DIST || rz.x < EPS) continue;\n    \t\tvec3 pos = ro+rz.x*rd;\n\t\tfloat thin = 2.5*map(pos,i).y;\n\t\tcolor = max(color, vec3(1.5,1.2,1.2) * vec3(1.0, thin, 1.0));\n\t}\n\tvec2 sphere = raymarchSphere(ro,rd);\n\tvec3 pos = ro + sphere.x*rd;\n\tif(sphere.x <= MAX_DIST && sphere.x > EPS){\n\t\tcolor = mix(color,sphereColor*calcLight(pos,lightPosition),0.125);\n\t}\t\n\n   \n    // float alpha=(color.x+ color.y + color.z); \n\n    // if(alpha < 0.005){\n    //     // alpha=0.0;\n    //     discard; \n    // }\n\tout_FragColor = vec4(color,color.x* color.y); \n}\n'),
      (zi[t(2374)] = t(2036)),
      (zi[t(505)] =
        '\n#extension GL_OES_standard_derivatives : enable\n\nprecision mediump float; \nin vec2 v_st;\n\nvoid main( void ) {\n  vec2 mouse=vec2(1.0,1.0);\n  float time = czm_frameNumber /100.0;\n  vec2 resolution = czm_viewport.zw;\n  vec2  st =v_st  * 2.5 - 1.25;\n\n    vec3   col =vec3(0.,0.,0.);\n\tvec2 p =st;// ( gl_FragCoord.xy*2.- resolution ) / min(resolution.x,resolution.y );\n\tvec2 m =  vec2(mouse.x*2.0-1.0,-mouse.y*2.0+1.0);\n\n\tfloat u =0.6*sin((atan(p.y,p.x)-time*0.2)*16.0);\n\tfloat t =0.0405/abs(0.15+u-length(p));\n\t  col +=vec3(0., t,t);  \n\t\n\tfloat u1 =0.6*sin((atan(p.y,p.x)+time*0.2)*16.0);\n\tfloat t1 =0.06/abs(0.25+u1-length(p));\t \n\t\n\t  col +=vec3(t1,t1/3.,0.); \n\t\n\t float t2 =0.06/length(p);\n   //col +=vec3(t2,0.,0.);\n\tfloat alpha=(col.x+ col.y + col.z) / 5.0; \n  if(alpha < 0.1){ \n      discard; \n  } \n\tout_FragColor = vec4( col,alpha );\n\n}\n'),
      (zi[t(1571)] =
        '\n#ifdef GL_ES\nprecision mediump float;\n#endif\n\nuniform float time;\nuniform vec2 mouse;\nuniform vec2 resolution;\n\nfloat line( vec2 a, vec2 b, vec2 p )\n{\n\tvec2 aTob = b - a;\n\tvec2 aTop = p - a;\n\t\n\tfloat t = dot( aTop, aTob ) / dot( aTob, aTob);\n\t\n\tt = clamp( t, 0.0, 1.0);\n\t\n\tfloat d = length( p - (a + aTob * t) );\n\td = 1.0 / d;\n\t\n\treturn clamp( d, 0.0, 1.0 );\n}\n\nin vec2 v_st;\n\nvoid main( )\n{\n    vec2 mouse=vec2(1.0,1.0);\n   float time = czm_frameNumber /100.0;\n   vec2 resolution = czm_viewport.zw;\n\tvec2 uv = v_st ;// * 2.5 - 1.25;//( gl_FragCoord.xy / resolution.xy );\n\t\n\tvec2 signedUV = uv * 2.0 - 1.0; //centers\n\tfloat aspectRatio = resolution.x / resolution.y;\n\tsignedUV.x *= aspectRatio;\n\tsignedUV.y *= -1.0;\n\t\n\tfloat scale = 40.0;\n\tconst float v = 90.0;\n\tvec3 finalColor = vec3( 0.0 );\n\t\n\tfloat timePulse = clamp(sin(time*2.0), 0.5, 2.4);\n\tfloat t = line( vec2(0.0, -4.0 + v * 0.4), vec2(0.0, -4.0 -v * 0.3), signedUV * scale );\n\tfinalColor += vec3( 8.0 * t, 4.0 * t, 2.0 * t) * 0.5 * timePulse;\n\t\n\tt = line( vec2(-v * 0.2, -v*0.1), vec2(v * 0.2, -v*0.1), signedUV * scale );\n\tfinalColor += vec3( 8.0 * t, 4.0 * t, 2.0 * t) * 0.4 * timePulse;\t\n\n    float alpha=(finalColor.x+ finalColor.y + finalColor.z) / 10.0; \n    if(alpha < 0.05){ \n        discard; \n    } \n\tout_FragColor = vec4( finalColor, alpha);\n\n}'),
      (zi[t(587)] = t(430)),
      (zi[t(211)] = t(1651)),
      (zi[t(2311)] = t(1001)),
      (zi[t(1618)] =
        '\n  precision mediump float; \n  float time; \n  vec2 resolution;\n  in vec2 v_st;\n  void main( void ) {\n      time = czm_frameNumber /100.0;\n       resolution=czm_viewport.zw;\n      vec2 p =v_st/2.0 - .25 ;\n\t//vec2 p = ( gl_FragCoord.xy - 0.5 * resolution.xy ) / resolution.y;\n\tfloat dd = length(p)*1.25;\n\tp.x += sin(p.y*+time)*0.01;\n\t\n\tfloat pi = 3.141592;\n\t\n\tfloat x = cos(sin(p.y*24.0+time))*4.0;\n\tfloat y = sin(time*1.3+dd+time);\n\t\n\tfloat d0 = length(p + vec2(0, 0.0));\n\t\n\tfloat color0 = sin(2.25 * pi * d0*d0/0.3 - time);\n\tfloat color = 0.65+ (sin(x+time*3.4+color0)*cos(y+time+color0))*0.25;\n\n\n\tvec4 col = vec4( vec3( sin(2.1/color*0.5),color*0.8,color*(0.9+sin(p.y*2.0+time*4.0)*0.3) ), 1.0 )*0.15/dd;\n    col *=0.8-dd;\n    if((col.r+col.g+col.b)<1.)discard;\n    col.a/=10.;\n    out_FragColor=col;\n}\n'),
      (zi[t(524)] = t(1815)),
      (zi[t(616)] =
        '\n\n\n\n\n// 水体\n\n// /*\n\n// A quick experiment with rain drop ripples.\n\n// This effect was written for and used in the launch scene of the\n// 64kB intro "H - Immersion", by Ctrl-Alt-Test.\n\n//  > http://www.ctrl-alt-test.fr/productions/h-immersion/\n//  > https://www.youtube.com/watch?v=27PN1SsXbjM\n\n// -- \n// Zavie / Ctrl-Alt-Test\n\n// */\n\n// // Maximum number of cells a ripple can cross.\n// #define MAX_RADIUS 2\n\n// // Set to 1 to hash twice. Slower, but less patterns.\n// #define DOUBLE_HASH 0\n\n// // Hash functions shamefully stolen from:\n// // https://www.shadertoy.com/view/4djSRW\n// #define HASHSCALE1 .1031\n// #define HASHSCALE3 vec3(.1031, .1030, .0973)\n\n// float hash12(vec2 p)\n// {\n// \tvec3 p3  = fract(vec3(p.xyx) * HASHSCALE1);\n//     p3 += dot(p3, p3.yzx + 19.19);\n//     return fract((p3.x + p3.y) * p3.z);\n// }\n\n// vec2 hash22(vec2 p)\n// {\n// \tvec3 p3 = fract(vec3(p.xyx) * HASHSCALE3);\n//     p3 += dot(p3, p3.yzx+19.19);\n//     return fract((p3.xx+p3.yz)*p3.zy);\n\n// }\n\n// void main( )\n// {\n//    // float resolution = //10. * exp2(-3.*iMouse.x/iResolution.x);\n// \tvec2 uv = v_st;\n//     vec2 p0 = floor(uv);\n\n//     vec2 circles = vec2(0.);\n//     for (int j = -MAX_RADIUS; j <= MAX_RADIUS; ++j)\n//     {\n//         for (int i = -MAX_RADIUS; i <= MAX_RADIUS; ++i)\n//         {\n// \t\t\tvec2 pi = p0 + vec2(i, j);\n//             #if DOUBLE_HASH\n//             vec2 hsh = hash22(pi);\n//             #else\n//             vec2 hsh = pi;\n//             #endif\n//             vec2 p = pi + hash22(hsh);\n\n//             float t = fract(0.3*iTime + hash12(hsh));\n//             vec2 v = p - uv;\n//             float d = length(v) - (float(MAX_RADIUS) + 1.)*t;\n\n//             float h = 1e-3;\n//             float d1 = d - h;\n//             float d2 = d + h;\n//             float p1 = sin(31.*d1) * smoothstep(-0.6, -0.3, d1) * smoothstep(0., -0.3, d1);\n//             float p2 = sin(31.*d2) * smoothstep(-0.6, -0.3, d2) * smoothstep(0., -0.3, d2);\n//             circles += 0.5 * normalize(v) * ((p2 - p1) / (2. * h) * (1. - t) * (1. - t));\n//         }\n//     }\n//     circles /= float((MAX_RADIUS*2+1)*(MAX_RADIUS*2+1));\n\n//     float intensity = mix(0.01, 0.15, smoothstep(0.1, 0.6, abs(fract(0.05*iTime + 0.5)*2.-1.)));\n//     vec3 n = vec3(circles, sqrt(1. - dot(circles, circles)));\n//     vec3 color = texture(iChannel0, uv - intensity*n.xy).rgb + 5.*pow(clamp(dot(n, normalize(vec3(1., 0.7, 0.5))), 0., 1.), 6.);\n// \tfragColor = vec4(color, 1.0);\n// }\n\n// 分解\n// float random(float x) {\n \n//     return fract(sin(x) * 10000.);\n          \n// }\n\n// float noise(vec2 p) {\n\n//     return random(p.x + p.y * 10000.);\n            \n// }\n\n// vec2 sw(vec2 p) { return vec2(floor(p.x), floor(p.y)); }\n// vec2 se(vec2 p) { return vec2(ceil(p.x), floor(p.y)); }\n// vec2 nw(vec2 p) { return vec2(floor(p.x), ceil(p.y)); }\n// vec2 ne(vec2 p) { return vec2(ceil(p.x), ceil(p.y)); }\n\n// float smoothNoise(vec2 p) {\n\n//     vec2 interp = smoothstep(0., 1., fract(p));\n//     float s = mix(noise(sw(p)), noise(se(p)), interp.x);\n//     float n = mix(noise(nw(p)), noise(ne(p)), interp.x);\n//     return mix(s, n, interp.y);\n        \n// }\n\n// float fractalNoise(vec2 p) {\n\n//     float x = 0.;\n//     x += smoothNoise(p      );\n//     x += smoothNoise(p * 2. ) / 2.;\n//     x += smoothNoise(p * 4. ) / 4.;\n//     x += smoothNoise(p * 8. ) / 8.;\n//     x += smoothNoise(p * 16.) / 16.;\n//     x /= 1. + 1./2. + 1./4. + 1./8. + 1./16.;\n//     return x;\n            \n// }\n\n// float movingNoise(vec2 p) {\n \n//     float x = fractalNoise(p + iTime);\n//     float y = fractalNoise(p - iTime);\n//     return fractalNoise(p + vec2(x, y));   \n    \n// }\n\n// // call this for water noise function\n// float nestedNoise(vec2 p) {\n    \n//     float x = movingNoise(p);\n//     float y = movingNoise(p + 100.);\n//     return movingNoise(p + vec2(x, y));\n    \n// }\n// void main(   )\n// {\n// \tvec2 uv = v_st;\n//     float n = nestedNoise(uv * 6.);\n    \n// \tfragColor = vec4(mix(vec3(.4, .6, 1.), vec3(.1, .2, 1.), n), 1.);\n// }\n\n\n// 水\n// A documented, altered, recolored version of "Seascape".\n// The famous original at:\n// https://www.shadertoy.com/view/Ms2SD1\n\n// "Seascape" by Alexander Alekseev aka TDM - 2014\n// Commenting added by bteitler\n//  HSV/color adjustments and additional commenting by CaliCoastReplay - 2016\n\n// License Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.\n\n// PI is a mathematical constant relating the ratio of a circle\'s circumference (distance around\n// the edge) to its diameter (distance between two points opposite on the edge).  \n// Change pi at your own peril, with your own apologies to God.\nconst float PI\t \t= 3.14159265358;\n\n// Can you explain these epsilons to a wide graphics audience?  YOUR comment could go here.\nconst float EPSILON\t= 1e-3;\n#define  EPSILON_NRM\t(0.5 / iResolution.x)\n\n// Constant indicaing the number of steps taken while marching the light ray.  \nconst int NUM_STEPS = 6;\n\n//Constants relating to the iteration of the heightmap for the wave, another part of the rendering\n//process.\nconst int ITER_GEOMETRY = 2;\nconst int ITER_FRAGMENT =5;\n\n// Constants that represent physical characteristics of the sea, can and should be changed and \n//  played with\nconst float SEA_HEIGHT = 0.5;\nconst float SEA_CHOPPY = 3.0;\nconst float SEA_SPEED = 1.9;\nconst float SEA_FREQ = 0.24;\nconst vec3 SEA_BASE = vec3(0.11,0.19,0.22);\nconst vec3 SEA_WATER_COLOR = vec3(0.55,0.9,0.7);\n#define SEA_TIME (iTime * SEA_SPEED)\n\n//Matrix to permute the water surface into a complex, realistic form\nmat2 octave_m = mat2(1.7,1.2,-1.2,1.4);\n\n//Space bar key constant\nconst float KEY_SP    = 32.5/256.0;\n\n//CaliCoastReplay :  These HSV/RGB translation functions are\n//from http://gamedev.stackexchange.com/questions/59797/glsl-shader-change-hue-saturation-brightness\n//This one converts red-green-blue color to hue-saturation-value color\nvec3 rgb2hsv(vec3 c)\n{\n    vec4 K = vec4(0.0, -1.0 / 3.0, 2.0 / 3.0, -1.0);\n    vec4 p = mix(vec4(c.bg, K.wz), vec4(c.gb, K.xy), step(c.b, c.g));\n    vec4 q = mix(vec4(p.xyw, c.r), vec4(c.r, p.yzx), step(p.x, c.r));\n\n    float d = q.x - min(q.w, q.y);\n    float e = 1.0e-10;\n    return vec3(abs(q.z + (q.w - q.y) / (6.0 * d + e)), d / (q.x + e), q.x);\n}\n\n//CaliCoastReplay :  These HSV/RGB translation functions are\n//from http://gamedev.stackexchange.com/questions/59797/glsl-shader-change-hue-saturation-brightness\n//This one converts hue-saturation-value color to red-green-blue color\nvec3 hsv2rgb(vec3 c)\n{\n    vec4 K = vec4(1.0, 2.0 / 3.0, 1.0 / 3.0, 3.0);\n    vec3 p = abs(fract(c.xxx + K.xyz) * 6.0 - K.www);\n    return c.z * mix(K.xxx, clamp(p - K.xxx, 0.0, 1.0), c.y);\n}\n\n// math\n// bteitler: Turn a vector of Euler angles into a rotation matrix\nmat3 fromEuler(vec3 ang) {\n\tvec2 a1 = vec2(sin(ang.x),cos(ang.x));\n    vec2 a2 = vec2(sin(ang.y),cos(ang.y));\n    vec2 a3 = vec2(sin(ang.z),cos(ang.z));\n    mat3 m;\n    m[0] = vec3(a1.y*a3.y+a1.x*a2.x*a3.x,a1.y*a2.x*a3.x+a3.y*a1.x,-a2.y*a3.x);\n\tm[1] = vec3(-a2.y*a1.x,a1.y*a2.y,a2.x);\n\tm[2] = vec3(a3.y*a1.x*a2.x+a1.y*a3.x,a1.x*a3.x-a1.y*a3.y*a2.x,a2.y*a3.y);\n\treturn m;\n}\n\n// bteitler: A 2D hash function for use in noise generation that returns range [0 .. 1].  You could\n// use any hash function of choice, just needs to deterministic and return\n// between 0 and 1, and also behave randomly.  Googling "GLSL hash function" returns almost exactly \n// this function: http://stackoverflow.com/questions/4200224/random-noise-functions-for-glsl\n// Performance is a real consideration of hash functions since ray-marching is already so heavy.\nfloat hash( vec2 p ) {\n    float h = dot(p,vec2(127.1,311.7));\t\n    return fract(sin(h)*83758.5453123);\n}\n\n// bteitler: A 2D psuedo-random wave / terrain function.  This is actually a poor name in my opinion,\n// since its the "hash" function that is really the noise, and this function is smoothly interpolating\n// between noisy points to create a continuous surface.\nfloat noise( in vec2 p ) {\n    vec2 i = floor( p );\n    vec2 f = fract( p );\t\n\n    // bteitler: This is equivalent to the "smoothstep" interpolation function.\n    // This is a smooth wave function with input between 0 and 1\n    // (since it is taking the fractional part of <p>) and gives an output\n    // between 0 and 1 that behaves and looks like a wave.  This is far from obvious, but we can graph it to see\n    // Wolfram link: http://www.wolframalpha.com/input/?i=plot+x*x*%283.0-2.0*x%29+from+x%3D0+to+1\n    // This is used to interpolate between random points.  Any smooth wave function that ramps up from 0 and\n    // and hit 1.0 over the domain 0 to 1 would work.  For instance, sin(f * PI / 2.0) gives similar visuals.\n    // This function is nice however because it does not require an expensive sine calculation.\n    vec2 u = f*f*(3.0-2.0*f);\n\n    // bteitler: This very confusing looking mish-mash is simply pulling deterministic random values (between 0 and 1)\n    // for 4 corners of the grid square that <p> is inside, and doing 2D interpolation using the <u> function\n    // (remember it looks like a nice wave!) \n    // The grid square has points defined at integer boundaries.  For example, if <p> is (4.3, 2.1), we will \n    // evaluate at points (4, 2), (5, 2), (4, 3), (5, 3), and then interpolate x using u(.3) and y using u(.1).\n    return -1.0+2.0*mix( \n                mix( hash( i + vec2(0.0,0.0) ), \n                     hash( i + vec2(1.0,0.0) ), \n                        u.x),\n                mix( hash( i + vec2(0.0,1.0) ), \n                     hash( i + vec2(1.0,1.0) ), \n                        u.x), \n                u.y);\n}\n\n// bteitler: diffuse lighting calculation - could be tweaked to taste\n// lighting\nfloat diffuse(vec3 n,vec3 l,float p) {\n    return pow(dot(n,l) * 0.4 + 0.6,p);\n}\n\n// bteitler: specular lighting calculation - could be tweaked taste\nfloat specular(vec3 n,vec3 l,vec3 e,float s) {    \n    float nrm = (s + 8.0) / (3.1415 * 8.0);\n    return pow(max(dot(reflect(e,n),l),0.0),s) * nrm;\n}\n\n// bteitler: Generate a smooth sky gradient color based on ray direction\'s Y value\n// sky\nvec3 getSkyColor(vec3 e) {\n    e.y = max(e.y,0.0);\n    vec3 ret;\n    ret.x = pow(1.0-e.y,2.0);\n    ret.y = 1.0-e.y;\n    ret.z = 0.6+(1.0-e.y)*0.4;\n    return ret;\n}\n\n// sea\n// bteitler: TLDR is that this passes a low frequency random terrain through a 2D symmetric wave function that looks like this:\n// http://www.wolframalpha.com/input/?i=%7B1-%7B%7B%7BAbs%5BCos%5B0.16x%5D%5D+%2B+Abs%5BCos%5B0.16x%5D%5D+%28%281.+-+Abs%5BSin%5B0.16x%5D%5D%29+-+Abs%5BCos%5B0.16x%5D%5D%29%7D+*+%7BAbs%5BCos%5B0.16y%5D%5D+%2B+Abs%5BCos%5B0.16y%5D%5D+%28%281.+-+Abs%5BSin%5B0.16y%5D%5D%29+-+Abs%5BCos%5B0.16y%5D%5D%29%7D%7D%5E0.65%7D%7D%5E4+from+-20+to+20\n// The <choppy> parameter affects the wave shape.\nfloat sea_octave(vec2 uv, float choppy) {\n    // bteitler: Add the smoothed 2D terrain / wave function to the input coordinates\n    // which are going to be our X and Z world coordinates.  It may be unclear why we are doing this.\n    // This value is about to be passed through a wave function.  So we have a smoothed psuedo random height\n    // field being added to our (X, Z) coordinates, and then fed through yet another wav function below.\n    uv += noise(uv);\n    // Note that you could simply return noise(uv) here and it would take on the characteristics of our \n    // noise interpolation function u and would be a reasonable heightmap for terrain.  \n    // However, that isn\'t the shape we want in the end for an ocean with waves, so it will be fed through\n    // a more wave like function.  Note that although both x and y channels of <uv> have the same value added, there is a \n    // symmetry break because <uv>.x and <uv>.y will typically be different values.\n\n    // bteitler: This is a wave function with pointy peaks and curved troughs:\n    // http://www.wolframalpha.com/input/?i=1-abs%28cos%28x%29%29%3B\n    vec2 wv = 1.0-abs(sin(uv)); \n\n    // bteitler: This is a wave function with curved peaks and pointy troughs:\n    // http://www.wolframalpha.com/input/?i=abs%28cos%28x%29%29%3B\n    vec2 swv = abs(cos(uv));  \n  \n    // bteitler: Blending both wave functions gets us a new, cooler wave function (output between 0 and 1):\n    // http://www.wolframalpha.com/input/?i=abs%28cos%28x%29%29+%2B+abs%28cos%28x%29%29+*+%28%281.0-abs%28sin%28x%29%29%29+-+abs%28cos%28x%29%29%29\n    wv = mix(wv,swv,wv);\n\n    // bteitler: Finally, compose both of the wave functions for X and Y channels into a final \n    // 1D height value, shaping it a bit along the way.  First, there is the composition (multiplication) of\n    // the wave functions: wv.x * wv.y.  Wolfram will give us a cute 2D height graph for this!:\n    // http://www.wolframalpha.com/input/?i=%7BAbs%5BCos%5Bx%5D%5D+%2B+Abs%5BCos%5Bx%5D%5D+%28%281.+-+Abs%5BSin%5Bx%5D%5D%29+-+Abs%5BCos%5Bx%5D%5D%29%7D+*+%7BAbs%5BCos%5By%5D%5D+%2B+Abs%5BCos%5By%5D%5D+%28%281.+-+Abs%5BSin%5By%5D%5D%29+-+Abs%5BCos%5By%5D%5D%29%7D\n    // Next, we reshape the 2D wave function by exponentiation: (wv.x * wv.y)^0.65.  This slightly rounds the base of the wave:\n    // http://www.wolframalpha.com/input/?i=%7B%7BAbs%5BCos%5Bx%5D%5D+%2B+Abs%5BCos%5Bx%5D%5D+%28%281.+-+Abs%5BSin%5Bx%5D%5D%29+-+Abs%5BCos%5Bx%5D%5D%29%7D+*+%7BAbs%5BCos%5By%5D%5D+%2B+Abs%5BCos%5By%5D%5D+%28%281.+-+Abs%5BSin%5By%5D%5D%29+-+Abs%5BCos%5By%5D%5D%29%7D%7D%5E0.65\n    // one last final transform (with choppy = 4) results in this which resembles a recognizable ocean wave shape in 2D:\n    // http://www.wolframalpha.com/input/?i=%7B1-%7B%7B%7BAbs%5BCos%5Bx%5D%5D+%2B+Abs%5BCos%5Bx%5D%5D+%28%281.+-+Abs%5BSin%5Bx%5D%5D%29+-+Abs%5BCos%5Bx%5D%5D%29%7D+*+%7BAbs%5BCos%5By%5D%5D+%2B+Abs%5BCos%5By%5D%5D+%28%281.+-+Abs%5BSin%5By%5D%5D%29+-+Abs%5BCos%5By%5D%5D%29%7D%7D%5E0.65%7D%7D%5E4\n    // Note that this function is called with a specific frequency multiplier which will stretch out the wave.  Here is the graph\n    // with the base frequency used by map and map_detailed (0.16):\n    // http://www.wolframalpha.com/input/?i=%7B1-%7B%7B%7BAbs%5BCos%5B0.16x%5D%5D+%2B+Abs%5BCos%5B0.16x%5D%5D+%28%281.+-+Abs%5BSin%5B0.16x%5D%5D%29+-+Abs%5BCos%5B0.16x%5D%5D%29%7D+*+%7BAbs%5BCos%5B0.16y%5D%5D+%2B+Abs%5BCos%5B0.16y%5D%5D+%28%281.+-+Abs%5BSin%5B0.16y%5D%5D%29+-+Abs%5BCos%5B0.16y%5D%5D%29%7D%7D%5E0.65%7D%7D%5E4+from+-20+to+20\n    return pow(1.0-pow(wv.x * wv.y,0.65),choppy);\n}\n\n// bteitler: Compute the distance along Y axis of a point to the surface of the ocean\n// using a low(er) resolution ocean height composition function (less iterations).\nfloat map(vec3 p) {\n    float freq = SEA_FREQ;\n    float amp = SEA_HEIGHT;\n    float choppy = SEA_CHOPPY;\n    vec2 uv = p.xz; uv.x *= 0.75;\n    \n    // bteitler: Compose our wave noise generation ("sea_octave") with different frequencies\n    // and offsets to achieve a final height map that looks like an ocean.  Likely lots\n    // of black magic / trial and error here to get it to look right.  Each sea_octave has this shape:\n    // http://www.wolframalpha.com/input/?i=%7B1-%7B%7B%7BAbs%5BCos%5B0.16x%5D%5D+%2B+Abs%5BCos%5B0.16x%5D%5D+%28%281.+-+Abs%5BSin%5B0.16x%5D%5D%29+-+Abs%5BCos%5B0.16x%5D%5D%29%7D+*+%7BAbs%5BCos%5B0.16y%5D%5D+%2B+Abs%5BCos%5B0.16y%5D%5D+%28%281.+-+Abs%5BSin%5B0.16y%5D%5D%29+-+Abs%5BCos%5B0.16y%5D%5D%29%7D%7D%5E0.65%7D%7D%5E4+from+-20+to+20\n    // which should give you an idea of what is going.  You don\'t need to graph this function because it\n    // appears to your left :)\n    float d, h = 0.0;    \n    for(int i = 0; i < ITER_GEOMETRY; i++) {\n        // bteitler: start out with our 2D symmetric wave at the current frequency\n    \td = sea_octave((uv+SEA_TIME)*freq,choppy);\n        // bteitler: stack wave ontop of itself at an offset that varies over time for more height and wave pattern variance\n    \t//d += sea_octave((uv-SEA_TIME)*freq,choppy);\n\n        h += d * amp; // bteitler: Bump our height by the current wave function\n        \n        // bteitler: "Twist" our domain input into a different space based on a permutation matrix\n        // The scales of the matrix values affect the frequency of the wave at this iteration, but more importantly\n        // it is responsible for the realistic assymetry since the domain is shiftly differently.\n        // This is likely the most important parameter for wave topology.\n    \tuv *=  octave_m;\n        \n        freq *= 1.9; // bteitler: Exponentially increase frequency every iteration (on top of our permutation)\n        amp *= 0.22; // bteitler: Lower the amplitude every frequency, since we are adding finer and finer detail\n        // bteitler: finally, adjust the choppy parameter which will effect our base 2D sea_octave shape a bit.  This makes\n        // the "waves within waves" have different looking shapes, not just frequency and offset\n        choppy = mix(choppy,1.0,0.2);\n    }\n    return p.y - h;\n}\n\n// bteitler: Compute the distance along Y axis of a point to the surface of the ocean\n// using a high(er) resolution ocean height composition function (more iterations).\nfloat map_detailed(vec3 p) {\n    float freq = SEA_FREQ;\n    float amp = SEA_HEIGHT;\n    float choppy = SEA_CHOPPY;\n    vec2 uv = p.xz; uv.x *= 0.75;\n    \n    // bteitler: Compose our wave noise generation ("sea_octave") with different frequencies\n    // and offsets to achieve a final height map that looks like an ocean.  Likely lots\n    // of black magic / trial and error here to get it to look right.  Each sea_octave has this shape:\n    // http://www.wolframalpha.com/input/?i=%7B1-%7B%7B%7BAbs%5BCos%5B0.16x%5D%5D+%2B+Abs%5BCos%5B0.16x%5D%5D+%28%281.+-+Abs%5BSin%5B0.16x%5D%5D%29+-+Abs%5BCos%5B0.16x%5D%5D%29%7D+*+%7BAbs%5BCos%5B0.16y%5D%5D+%2B+Abs%5BCos%5B0.16y%5D%5D+%28%281.+-+Abs%5BSin%5B0.16y%5D%5D%29+-+Abs%5BCos%5B0.16y%5D%5D%29%7D%7D%5E0.65%7D%7D%5E4+from+-20+to+20\n    // which should give you an idea of what is going.  You don\'t need to graph this function because it\n    // appears to your left :)\n    float d, h = 0.0;    \n    for(int i = 0; i < ITER_FRAGMENT; i++) {\n        // bteitler: start out with our 2D symmetric wave at the current frequency\n    \td = sea_octave((uv+SEA_TIME)*freq,choppy);\n        // bteitler: stack wave ontop of itself at an offset that varies over time for more height and wave pattern variance\n    \td += sea_octave((uv-SEA_TIME)*freq,choppy);\n        \n        h += d * amp; // bteitler: Bump our height by the current wave function\n        \n        // bteitler: "Twist" our domain input into a different space based on a permutation matrix\n        // The scales of the matrix values affect the frequency of the wave at this iteration, but more importantly\n        // it is responsible for the realistic assymetry since the domain is shiftly differently.\n        // This is likely the most important parameter for wave topology.\n    \tuv *= octave_m/1.2;\n        \n        freq *= 1.9; // bteitler: Exponentially increase frequency every iteration (on top of our permutation)\n        amp *= 0.22; // bteitler: Lower the amplitude every frequency, since we are adding finer and finer detail\n        // bteitler: finally, adjust the choppy parameter which will effect our base 2D sea_octave shape a bit.  This makes\n        // the "waves within waves" have different looking shapes, not just frequency and offset\n        choppy = mix(choppy,1.0,0.2);\n    }\n    return p.y - h;\n}\n\n// bteitler:\n// p: point on ocean surface to get color for\n// n: normal on ocean surface at <p>\n// l: light (sun) direction\n// eye: ray direction from camera position for this pixel\n// dist: distance from camera to point <p> on ocean surface\nvec3 getSeaColor(vec3 p, vec3 n, vec3 l, vec3 eye, vec3 dist) {  \n    // bteitler: Fresnel is an exponential that gets bigger when the angle between ocean\n    // surface normal and eye ray is smaller\n    float fresnel = 1.0 - max(dot(n,-eye),0.0);\n    fresnel = pow(fresnel,3.0) * 0.45;\n        \n    // bteitler: Bounce eye ray off ocean towards sky, and get the color of the sky\n    vec3 reflected = getSkyColor(reflect(eye,n))*0.99;    \n    \n    // bteitler: refraction effect based on angle between light surface normal\n    vec3 refracted = SEA_BASE + diffuse(n,l,80.0) * SEA_WATER_COLOR * 0.27; \n    \n    // bteitler: blend the refracted color with the reflected color based on our fresnel term\n    vec3 color = mix(refracted,reflected,fresnel);\n    \n    // bteitler: Apply a distance based attenuation factor which is stronger\n    // at peaks\n    float atten = max(1.0 - dot(dist,dist) * 0.001, 0.0);\n    color += SEA_WATER_COLOR * (p.y - SEA_HEIGHT) * 0.15 * atten;\n    \n    // bteitler: Apply specular highlight\n    color += vec3(specular(n,l,eye,90.0))*0.5;\n    \n    return color;\n}\n\n// bteitler: Estimate the normal at a point <p> on the ocean surface using a slight more detailed\n// ocean mapping function (using more noise octaves).\n// Takes an argument <eps> (stands for epsilon) which is the resolution to use\n// for the gradient.  See here for more info on gradients: https://en.wikipedia.org/wiki/Gradient\n// tracing\nvec3 getNormal(vec3 p, float eps) {\n    // bteitler: Approximate gradient.  An exact gradient would need the "map" / "map_detailed" functions\n    // to return x, y, and z, but it only computes height relative to surface along Y axis.  I\'m assuming\n    // for simplicity and / or optimization reasons we approximate the gradient by the change in ocean\n    // height for all axis.\n    vec3 n;\n    n.y = map_detailed(p); // bteitler: Detailed height relative to surface, temporarily here to save a variable?\n    n.x = map_detailed(vec3(p.x+eps,p.y,p.z)) - n.y; // bteitler approximate X gradient as change in height along X axis delta\n    n.z = map_detailed(vec3(p.x,p.y,p.z+eps)) - n.y; // bteitler approximate Z gradient as change in height along Z axis delta\n    // bteitler: Taking advantage of the fact that we know we won\'t have really steep waves, we expect\n    // the Y normal component to be fairly large always.  Sacrifices yet more accurately to avoid some calculation.\n    n.y = eps; \n    return normalize(n);\n\n    // bteitler: A more naive and easy to understand version could look like this and\n    // produces almost the same visuals and is a little more expensive.\n    // vec3 n;\n    // float h = map_detailed(p);\n    // n.y = map_detailed(vec3(p.x,p.y+eps,p.z)) - h;\n    // n.x = map_detailed(vec3(p.x+eps,p.y,p.z)) - h;\n    // n.z = map_detailed(vec3(p.x,p.y,p.z+eps)) - h;\n    // return normalize(n);\n}\n\n\n//CaliCoastReplay :  Keyboard checking function from the iChannel representing keyboard input\nfloat isKeyPressed(float key)\n{\n\treturn texture( iChannel1, vec2(key, 1.0) ).x;\n}\n\n// bteitler: Find out where a ray intersects the current ocean\nfloat heightMapTracing(vec3 ori, vec3 dir, out vec3 p) {  \n    float tm = 0.0;\n    float tx = 500.0; // bteitler: a really far distance, this could likely be tweaked a bit as desired\n\n    // bteitler: At a really far away distance along the ray, what is it\'s height relative\n    // to the ocean in ONLY the Y direction?\n    float hx = map(ori + dir * tx);\n    \n    // bteitler: A positive height relative to the ocean surface (in Y direction) at a really far distance means\n    // this pixel is pure sky.  Quit early and return the far distance constant.\n    if(hx > 0.0) return tx;   \n\n    // bteitler: hm starts out as the height of the camera position relative to ocean.\n    float hm = map(ori + dir * tm); \n   \n    // bteitler: This is the main ray marching logic.  This is probably the single most confusing part of the shader\n    // since height mapping is not an exact distance field (tells you distance to surface if you drop a line down to ocean\n    // surface in the Y direction, but there could have been a peak at a very close point along the x and z \n    // directions that is closer).  Therefore, it would be possible/easy to overshoot the surface using the raw height field\n    // as the march distance.  The author uses a trick to compensate for this.\n    float tmid = 0.0;\n    for(int i = 0; i < NUM_STEPS; i++) { // bteitler: Constant number of ray marches per ray that hits the water\n        // bteitler: Move forward along ray in such a way that has the following properties:\n        // 1. If our current height relative to ocean is higher, move forward more\n        // 2. If the height relative to ocean floor very far along the ray is much lower\n        //    below the ocean surface, move forward less\n        // Idea behind 1. is that if we are far above the ocean floor we can risk jumping\n        // forward more without shooting under ocean, because the ocean is mostly level.\n        // The idea behind 2. is that if extruding the ray goes farther under the ocean, then \n        // you are looking more orthgonal to ocean surface (as opposed to looking towards horizon), and therefore\n        // movement along the ray gets closer to ocean faster, so we need to move forward less to reduce risk\n        // of overshooting.\n        tmid = mix(tm,tx, hm/(hm-hx));\n        p = ori + dir * tmid; \n                  \n    \tfloat hmid = map(p); // bteitler: Re-evaluate height relative to ocean surface in Y axis\n\n        if(hmid < 0.0) { // bteitler: We went through the ocean surface if we are negative relative to surface now\n            // bteitler: So instead of actually marching forward to cross the surface, we instead\n            // assign our really far distance and height to be where we just evaluated that crossed the surface.\n            // Next iteration will attempt to go forward more and is less likely to cross the boundary.\n            // A naive implementation might have returned <tmid> immediately here, which\n            // results in a much poorer / somewhat indeterministic quality rendering.\n            tx = tmid;\n            hx = hmid;\n        } else {\n            // Haven\'t hit surface yet, easy case, just march forward\n            tm = tmid;\n            hm = hmid;\n        }\n    }\n\n    // bteitler: Return the distance, which should be really close to the height map without going under the ocean\n    return tmid;\n}\n\n// main\nvoid main() {\n    // bteitler: 2D Pixel location passed in as raw pixel, let\'s divide by resolution\n    // to convert to coordinates between 0 and 1\n    iResolution=czm_viewport.zw;\n    vec2 uv = v_st;//fragCoord.xy / iResolution.xy;\n  \n    uv = uv * 2.0  - 1.0; //  bteitler: Shift pixel coordinates from 0 to 1 to between -1 and 1\n    uv.x *= iResolution.x / iResolution.y; // bteitler: Aspect ratio correction - if you don\'t do this your rays will be distorted\n    float time = iTime * 2.7; // bteitler: Animation is based on time, but allows you to scrub the animation based on mouse movement\n        \n    // ray\n\n    // bteitler: Calculated a vector that smoothly changes over time in a sinusoidal (wave) pattern.  \n    // This will be used to drive where the user is looking in world space.\n   // vec3 ang = vec3(sin(time*3.0)*0.1,sin(time)*0.2+0.3,time);\n    float roll = PI + sin(iTime)/14.0 + cos(iTime/2.0)/14.0 ;\n    float pitch = PI*1.021 + (sin(iTime/2.0)+ cos(iTime))/40.0 \n        + (iMouse.y/iResolution.y - .8)*PI/3.0  ;\n    float yaw = iMouse.x/iResolution.x * PI * 4.0;\n    vec3 ang = vec3(roll,pitch,yaw);\n   // vec3 ang = vec3(roll,pitch,0);\n    \n    // bteitler: Calculate the "origin" of the camera in world space based on time.  Camera is located\n    // at height 3.5 atx 0 (zero), and flies over the ocean in the z axis over time.\n    vec3 ori = vec3(0.0,3.5,time*3.0);\n   \n    // bteitler: This is the ray direction we are shooting from the camera location ("ori") that we need to light\n    // for this pixel.  The -2.0 indicates we are using a focal length of 2.0 - this is just an artistic choice and\n    // results in about a 90 degree field of view.\n    //  CaliCoastReplay :  Adjusted slightly to a lower focal length.  Seems to dramatize the scene.\n    vec3 dir = normalize(vec3(uv.xy,-1.6)); \n\n    // bteitler: Distort the ray a bit for a fish eye effect (if you remove this line, it will remove\n    // the fish eye effect and look like a realistic perspective).\n   //  dir.z += length(uv) * 0.15;\n\n    // bteitler: Renormalize the ray direction, and then rotate it based on the previously calculated\n    // animation angle "ang".  "fromEuler" just calculates a rotation matrix from a vector of angles.\n    // if you remove the " * fromEuler(ang)" part, you will disable the camera rotation animation.\n    dir = normalize(dir) * fromEuler(ang);\n    \n    // tracing\n\n    // bteitler: ray-march to the ocean surface (which can be thought of as a randomly generated height map)\n    // and store in p\n    vec3 p;\n    heightMapTracing(ori,dir,p);\n\n    vec3 dist = p - ori; // bteitler: distance vector to ocean surface for this pixel\'s ray\n\n    // bteitler: Calculate the normal on the ocean surface where we intersected (p), using\n    // different "resolution" (in a sense) based on how far away the ray traveled.  Normals close to\n    // the camera should be calculated with high resolution, and normals far from the camera should be calculated with low resolution\n    // The reason to do this is that specular effects (or non linear normal based lighting effects) become fairly random at\n    // far distances and low resolutions and can cause unpleasant shimmering during motion.\n    vec3 n = getNormal(p, \n             dot(dist,dist)   // bteitler: Think of this as inverse resolution, so far distances get bigger at an expnential rate\n                * EPSILON_NRM // bteitler: Just a resolution constant.. could easily be tweaked to artistic content\n           );\n\n    // bteitler: direction of the infinitely far away directional light.  Changing this will change\n    // the sunlight direction.\n    vec3 light = normalize(vec3(0.0,1.0,0.8)); \n             \n    // CaliCoastReplay:  Get the sky and sea colors\n\tvec3 skyColor = getSkyColor(dir);\n    vec3 seaColor = getSeaColor(p,n,light,dir,dist);\n    \n    //Sea/sky preprocessing\n    \n    //CaliCoastReplay:  A distance falloff for the sea color.   Drastically darkens the sea, \n    //this will be reversed later based on day/night.\n    seaColor /= sqrt(sqrt(length(dist))) ;\n    \n    \n    //CaliCoastReplay:  Day/night mode\n    bool night; \t \n    if( isKeyPressed(KEY_SP) > 0.0 )    //night mode!\n    {\n        //Brighten the sea up again, but not too bright at night\n    \tseaColor *= seaColor * 8.5;\n        \n        //Turn down the sky \n    \tskyColor /= 1.69;\n        \n        //Store that it\'s night mode for later HSV calcc\n        night = true;\n    }\n    else  //day mode!\n    {\n        //Brighten the sea up again - bright and beautiful blue at day\n    \tseaColor *= sqrt(sqrt(seaColor)) * 4.0;\n        skyColor *= 1.05;\n        skyColor -= 0.03;\n        night = false;\n    }\n\n    \n    //CaliCoastReplay:  A slight "constrasting" for the sky to match the more contrasted ocean\n    skyColor *= skyColor;\n    \n    \n    //CaliCoastReplay:  A rather hacky manipulation of the high-value regions in the image that seems\n    //to add a subtle charm and "sheen" and foamy effect to high value regions through subtle darkening,\n    //but it is hacky, and not physically modeled at all.  \n    vec3 seaHsv = rgb2hsv(seaColor);\n    if (seaHsv.z > .75 && length(dist) < 50.0)\n       seaHsv.z -= (0.9 - seaHsv.z) * 1.3;\n    seaColor = hsv2rgb(seaHsv);\n    \n    // bteitler: Mix (linear interpolate) a color calculated for the sky (based solely on ray direction) and a sea color \n    // which contains a realistic lighting model.  This is basically doing a fog calculation: weighing more the sky color\n    // in the distance in an exponential manner.\n    \n    vec3 color = mix(\n        skyColor,\n        seaColor,\n    \tpow(smoothstep(0.0,-0.05,dir.y), 0.3) // bteitler: Can be thought of as "fog" that gets thicker in the distance\n    );\n        \n    // Postprocessing\n    \n    // bteitler: Apply an overall image brightness factor as the final color for this pixel.  Can be\n    // tweaked artistically.\n    fragColor = vec4(pow(color,vec3(0.75)), 1.0);\n    \n    // CaliCoastReplay:  Adjust hue, saturation, and value adjustment for an even more processed look\n    // hsv.x is hue, hsv.y is saturation, and hsv.z is value\n    vec3 hsv = rgb2hsv(fragColor.xyz);    \n    //CaliCoastReplay: Increase saturation slightly\n    hsv.y += 0.131;\n    //CaliCoastReplay:\n    //A pseudo-multiplicative adjustment of value, increasing intensity near 1 and decreasing it near\n    //0 to achieve a more contrasted, real-world look\n    hsv.z *= sqrt(hsv.z) * 1.1; \n    \n    if (night)    \n    {\n    ///CaliCoastReplay:\n    //Slight value adjustment at night to turn down global intensity\n        hsv.z -= 0.045;\n        hsv*=0.8;\n        hsv.x += 0.12 + hsv.z/100.0;\n        //Highly increased saturation at night op, oddly.  Nights appear to be very colorful\n        //within their ranges.\n        hsv.y *= 2.87;\n    }\n    else\n    {\n      //CaliCoastReplay:\n        //Add green tinge to the high range\n      //Turn down intensity in day in a different way     \n        \n        hsv.z *= 0.9;\n        \n        //CaliCoastReplay:  Hue alteration \n        hsv.x -= hsv.z/10.0;\n        hsv.x += 0.02 + hsv.z/50.0;\n        //Final brightening\n        hsv.z *= 1.01;\n        //This really "cinemafies" it for the day -\n        //puts the saturation on a squared, highly magnified footing.\n        //Worth looking into more as to exactly why.\n       // hsv.y *= 5.10 * hsv.y * sqrt(hsv.y);\n        hsv.y += 0.07;\n    }\n    \n    //CaliCoastReplay:    \n    //Replace the final color with the adjusted, translated HSV values\n    fragColor.xyz = hsv2rgb(hsv);\n}\n\n')
    const Di = {
        create(e, i) {
          const s = t
          switch (e) {
            case Ai[s(1565)]:
              return new (class extends Ti {
                constructor(t = {}) {
                  super(t)
                }
                [t(933)]() {
                  const e = t
                  return (
                    (this._radii = Cesium.defaultValue(
                      this[e(568)][e(1141)],
                      new Cesium[e(310)](10, 10, 10)
                    )),
                    Cesium[e(2487)][e(933)](
                      new Cesium.EllipsoidGeometry({
                        vertexFormat: Cesium[e(1788)][e(1315)],
                        radii: this._radii
                      })
                    )
                  )
                }
              })(i)
            case Ai[s(868)]:
              return new Ei(i)
            case Ai[s(1111)]:
              return new (class extends Ti {
                constructor(t = {}) {
                  super(t)
                }
                [t(933)]() {
                  const e = t,
                    i = new Cesium[e(2320)]({
                      vertexFormat: Cesium[e(1788)][e(1315)],
                      center: new Cesium[e(310)](0, 0, 0.1),
                      height: 1e3,
                      semiMajorAxis: 1e4,
                      semiMinorAxis: 1e4
                    })
                  return Cesium[e(2320)][e(933)](i)
                }
              })(i)
            case Ai[s(843)]:
              return new (class extends Ti {
                constructor(t = {}) {
                  super(t)
                }
                [t(933)]() {
                  const e = t
                  this[e(1613)] = Cesium.defaultValue(
                    this[e(568)].dimensions,
                    new Cesium.Cartesian3(100, 100, 100)
                  )
                  const i = Cesium[e(1435)][e(1254)]({
                    vertexFormat: Cesium[e(1788)][e(2050)],
                    dimensions: this._dimensions
                  })
                  return Cesium.BoxGeometry.createGeometry(i)
                }
              })(i)
            default:
              return new Ei(i)
          }
        }
      },
      Ii = zi
    class ki extends b {
      constructor(e) {
        const i = t
        super(e),
          (this._geometryType = i(1865)),
          (this[i(1115)] = q),
          (this[i(2336)] = Q[i(2513)]),
          (this[i(528)] = i(817)),
          (this._fixPointCount = 1),
          (this[i(1526)] = Cesium.defaultValue(e[i(2017)], Ai.plane)),
          (this[i(763)] = Cesium[i(1960)](e[i(1173)], {})),
          (this._outFragColorBody = Cesium[i(1960)](e[i(1929)], '')),
          (this._style = Cesium[i(1960)](e[i(1679)], {})),
          (this[i(343)][i(1070)] = Cesium[i(1960)](this._style[i(1070)], i(2217))),
          (this._style[i(506)] = Cesium[i(1960)](this[i(343)][i(506)], i(1655))),
          (this._style[i(2535)] = Cesium[i(1960)](this[i(343)][i(2535)], 6)),
          (this._style[i(1529)] = Cesium[i(1960)](this[i(343)][i(1529)], 1)),
          (this._style[i(1434)] = Cesium.defaultValue(this[i(343)][i(1434)], [])),
          (this[i(343)][i(2341)] = Cesium[i(1960)](
            this._style.scale,
            new Cesium.Cartesian3(1, 1, 1)
          )),
          (this._style[i(2307)] = Cesium[i(1960)](
            this[i(343)][i(2307)],
            new Cesium[i(310)](0, 0, 0)
          )),
          (this[i(425)] = Cesium[i(1960)](e[i(2489)], xi)),
          (this[i(1364)] = Cesium[i(1960)](e[i(406)], wi)),
          (this[i(2530)] = Cesium[i(1960)](e[i(995)], Ii[i(214)])),
          (this._fragmentShader = Cesium[i(1960)](e[i(2458)], Ii[i(1579)])),
          (this[i(995)] = this[i(2530)]),
          (this.fragmentShader = this[i(932)]),
          (this._modelMatrixNeedsUpdate = !0),
          (this[i(577)] = Cesium.defaultValue(e.attributeLocations, {
            position: 0,
            st: 1,
            normal: 2,
            color: 3
          })),
          (this[i(1912)] = Cesium[i(1960)](e[i(2251)], [111, 28, 0])),
          this[i(1253)](this[i(1912)]),
          (this._shaderGeometryOpts[i(2251)] = this[i(2238)]),
          (this[i(1636)] = new Cesium[i(2066)]()),
          Cesium[i(2066)][i(1902)](Cesium[i(2066)][i(1393)], this[i(1636)]),
          (this[i(2179)] = null),
          (this[i(2249)] = null),
          (this[i(1939)] = e[i(1333)] || !1),
          (this[i(357)] = !1 !== e[i(2163)]),
          (this[i(979)] = Cesium[i(1960)](e[i(2088)], !0)),
          (this[i(2259)] = new Cesium.Quaternion())
      }
      [t(1253)](e) {
        const i = t
        e &&
          ((this._position = e),
          (this[i(794)] = this[i(2251)]),
          (this._cartesian3 = Cesium[i(310)][i(667)](
            this[i(1912)][0],
            this[i(1912)][1],
            this[i(1912)][2]
          )),
          (this[i(1024)] = new Cesium[i(1242)](this[i(2238)], 2)),
          (this[i(826)] = !0),
          this[i(245)] && this[i(245)][i(1391)][i(1499)](1e-6))
      }
      _createShaderProgram(e) {
        const i = t,
          s = {}
        s[i(1177)] = [this._material.vertexShader]
        const n = {}
        n.sources = [this[i(1296)][i(2458)]]
        let o = new Cesium[i(231)](s),
          r = new Cesium.ShaderSource(n)
        const a = {}
        return (
          (a.context = e),
          (a[i(2387)] = o),
          (a[i(1363)] = r),
          (a[i(2079)] = this._attributeLocations),
          Cesium[i(178)][i(2407)](a)
        )
      }
      [t(1952)](e) {
        const i = t
        let s = e[i(2325)],
          n = (function (t) {
            const e = i,
              s = {}
            ;(s[e(1292)] = 0), (s[e(348)] = 1)
            const n = {}
            ;(n[e(1626)] = !0), (n[e(1361)] = !0), (n.blue = !0), (n[e(2315)] = !0)
            const o = {}
            ;(o[e(1878)] = t[e(1878)]
              ? Cesium.BlendingState[e(1785)]
              : Cesium.BlendingState[e(1192)]),
              (o[e(1473)] = {}),
              (o[e(1961)] = {}),
              (o.depthRange = s),
              (o.colorMask = n),
              (o[e(898)] = t[e(898)]),
              (o[e(1473)][e(299)] = t[e(1473)]),
              (o[e(1473)][e(2304)] = Cesium[e(1919)][e(689)]),
              (o[e(1961)][e(299)] = !1),
              (o[e(1961)][e(2563)] = Cesium[e(1749)].FRONT)
            let r = o
            switch (((r[e(1961)][e(299)] = !0), t[e(2489)])) {
              case xi:
                r.cull[e(2563)] = Cesium[e(1749)].BACK
                break
              case bi:
                r.cull[e(2563)] = Cesium[e(1749)][e(1402)]
                break
              default:
                r[e(1961)][e(299)] = !1
            }
            return (t[e(1721)] = r), Cesium[e(515)][e(2407)](r)
          })(this[i(1296)]),
          o = Cesium[i(1671)][i(2545)](this[i(1398)])
        o = Cesium[i(1057)][i(1225)]({
          context: s,
          geometry: o,
          attributeLocations: this._attributeLocations,
          bufferUsage: Cesium.BufferUsage[i(1047)]
        })
        let r,
          a = this[i(1943)]
        const h = {}
        return (
          (h[i(2416)] = this),
          (h.id = this[i(1570)]),
          a || ((a = e[i(2325)][i(896)](h)), (this[i(1943)] = a)),
          (r = a[i(1070)]),
          (this._material[i(1720)][i(1208)] = function () {
            return r
          }),
          (this._drawCommand = new Cesium[i(905)]({
            vertexArray: o,
            primitiveType: Cesium[i(1099)][i(2367)],
            renderState: this[i(998)] || n,
            shaderProgram: this._createShaderProgram(s),
            uniformMap: this[i(1296)][i(1720)],
            owner: this,
            pass: Cesium[i(1864)][i(2403)],
            modelMatrix: this[i(357)] ? this[i(1636)] : Cesium.Matrix4[i(1393)][i(1902)](),
            boundingVolume: this[i(357)]
              ? Cesium[i(1242)][i(1863)](
                  this._geometry[i(1611)],
                  this[i(1636)],
                  new Cesium[i(1242)]()
                )
              : this[i(1398)][i(1611)]
          })),
          (this[i(1024)] = this[i(2179)].boundingVolume),
          (s = Cesium[i(231)].replaceMain(this[i(1296)][i(2458)], i(1908))),
          (o = Cesium[i(178)].replaceCache({
            context: e.context,
            shaderProgram: void 0,
            vertexShaderSource: this[i(1296)][i(995)],
            fragmentShaderSource: i(2232) + s + i(1167),
            attributeLocations: this._attributeLocations
          })),
          (this[i(2249)] = new Cesium[i(905)]({
            owner: this,
            pickOnly: !0,
            instanceCount: void 0,
            modelMatrix: this[i(2179)][i(223)],
            primitiveType: this[i(2179)][i(2283)],
            cull: !1,
            pass: Cesium[i(1864)][i(2403)]
          })),
          (this[i(2249)][i(1239)] = this._drawCommand[i(1239)]),
          (this[i(2249)][i(1689)] = this[i(998)] || n),
          (this[i(2249)][i(187)] = o),
          (this[i(2249)][i(1720)] = this[i(2179)][i(1720)]),
          (this[i(2249)][i(1174)] = !1),
          this._drawCommand
        )
      }
      [t(2231)]() {
        const e = t
        let i = this._isSelected ? this[e(343)][e(506)] : this[e(343)].color
        this._material[e(2390)] = Cesium.Color[e(2008)](i)
      }
      [t(720)](e) {
        const i = t
        if (this[i(1144)]) {
          var s, n
          this[i(826)] &&
            ((this[i(1636)] = this._computeModelMatrix()),
            this[i(2179)] &&
              ((this._drawCommand.modelMatrix = this._modelMatrix),
              (this._drawCommand.boundingVolume = this[i(357)]
                ? Cesium.BoundingSphere[i(1863)](
                    this[i(1398)][i(1611)],
                    this._modelMatrix,
                    new Cesium[i(1242)]()
                  )
                : this[i(1398)][i(1611)]),
              (this[i(1024)] = this._drawCommand[i(2172)])),
            (this[i(826)] = !1)),
            Cesium.defined(this[i(2179)])
              ? this[i(2179)] &&
                this._pickCommand &&
                ((this[i(2249)][i(2172)] = this[i(2179)].boundingVolume),
                (this[i(2249)].modelMatrix = this._drawCommand.modelMatrix),
                this[i(300)] &&
                  (this._drawCommand.vertexArray
                    .getAttribute(0)
                    .vertexBuffer.copyFromArrayView(
                      new Float32Array(this[i(1398)][i(1548)][i(2251)][i(736)])
                    ),
                  (this[i(300)] = !1)),
                e.passes[i(1646)]
                  ? (this[i(2249)][i(2172)] ||
                      ((s = o[i(1653)][i(2134)]()),
                      (n = Cesium[i(1671)][i(2545)](s)),
                      (this[i(2249)].vertexArray = Cesium[i(1057)][i(1225)]({
                        context: e.context,
                        geometry: n,
                        attributeLocations: this[i(577)],
                        bufferUsage: Cesium.BufferUsage[i(2200)]
                      })),
                      (this[i(2249)].boundingVolume = Cesium[i(1242)][i(1863)](
                        s[i(1611)],
                        this[i(1636)],
                        s[i(1611)]
                      )),
                      (this[i(2249)][i(223)] = this[i(1636)])),
                    e.commandList[i(2553)](this._pickCommand))
                  : e[i(331)].push(this._drawCommand))
              : (this._drawCommand = this[i(1952)](e)),
            this[i(828)](e)
        }
      }
      [t(828)](e) {
        const i = t
        this[i(245)] &&
          ((this[i(518)] = Cesium.Quaternion[i(814)](
            Cesium[i(1472)][i(2173)](
              Cesium.Matrix4.getMatrix3(
                Cesium[i(2066)][i(739)](this[i(245)].camera),
                new Cesium[i(819)]()
              ),
              new Cesium[i(1472)]()
            ),
            new Cesium[i(1472)]()
          )),
          this._drawCommand &&
            (this[i(1939)] && this[i(518)] && (this[i(1333)] = !0),
            this[i(979)] && this._cameraQ && (this[i(2088)] = !0)))
      }
      [t(1584)]() {
        const e = t,
          i = this[e(343)][e(2307)],
          s = this[e(343)][e(2341)]
        let n = Cesium.Matrix4[e(942)](Cesium[e(819)][e(2377)](Cesium[e(475)][e(1149)](i.x))),
          o = Cesium[e(2066)][e(942)](Cesium[e(819)][e(1693)](Cesium[e(475)].toRadians(i.y))),
          r = Cesium[e(2066)][e(942)](Cesium[e(819)][e(712)](Cesium[e(475)].toRadians(i.z))),
          a = Cesium[e(2058)][e(1648)](this[e(2238)])
        return (
          Cesium[e(2066)][e(999)](a, n, a),
          Cesium[e(2066)][e(999)](a, o, a),
          Cesium[e(2066)][e(999)](a, r, a),
          (n = Cesium[e(2066)].fromScale(new Cesium.Cartesian3(s.x, s.y, s.z))),
          Cesium.Matrix4[e(999)](a, n, new Cesium.Matrix4())
        )
      }
      [t(2039)]() {
        const e = t
        return (
          Cesium[e(2330)](this._drawCommand) &&
            ((this[e(2179)][e(1239)] =
              this._drawCommand.vertexArray && this._drawCommand.vertexArray[e(2039)]()),
            (this[e(2179)][e(187)] = this[e(2179)][e(187)] && this[e(2179)][e(187)][e(2039)]()),
            this[e(2249)] &&
              (this._pickCommand.vertexArray =
                this[e(2249)].vertexArray && this[e(2249)][e(1239)].destroy()),
            this[e(2249)] &&
              (this[e(2249)][e(187)] = this[e(2249)][e(187)] && this[e(2249)][e(187)][e(2039)]())),
          Cesium.destroyObject(this)
        )
      }
      get [t(1333)]() {
        return this[t(1939)]
      }
      set [t(1333)](e) {
        const i = t
        ;(this[i(1939)] = e),
          (this[i(1939)] = e),
          this[i(1939)] &&
            this._height &&
            ((this._angle = this._angle + 0.06),
            0 < Math[i(884)](this[i(366)]) ? (this[i(1516)] = 0.008) : (this._height = -0.008),
            (e = new Cesium[i(310)](0, this[i(1516)], 0)),
            Cesium[i(2066)][i(229)](this[i(2179)][i(1636)], e, this[i(2179)][i(1636)]))
      }
      get [t(2088)]() {
        return this._lookAtState
      }
      set [t(2088)](e) {
        const i = t
        if (!this[i(2179)]) return
        let s,
          n,
          o,
          r,
          a,
          h,
          l,
          c,
          u,
          m,
          p,
          d,
          f,
          C = this._style[i(2341)]
        ;(this[i(979)] = e),
          this._lookAtState &&
            !this[i(2259)][i(1525)](this._cameraQ) &&
            ((e = this[i(2179)]),
            (s = this[i(2238)]),
            (n = this[i(518)]),
            (o = new Cesium[i(310)](C.x, C.y, C.z)),
            (e = e.modelMatrix),
            (l = (r = n.x) * (m = r + r)),
            (c = r * (p = (a = n.y) + a)),
            (r *= d = (h = n.z) + h),
            (u = a * p),
            (a *= d),
            (h *= d),
            (m *= n = n.w),
            (p *= n),
            (n *= d),
            (d = o.x),
            (f = o.y),
            (o = o.z),
            (e[0] = (1 - (u + h)) * d),
            (e[1] = (c + n) * d),
            (e[2] = (r - p) * d),
            (e[3] = 0),
            (e[4] = (c - n) * f),
            (e[5] = (1 - (l + h)) * f),
            (e[6] = (a + m) * f),
            (e[7] = 0),
            (e[8] = (r + p) * o),
            (e[9] = (a - m) * o),
            (e[10] = (1 - (l + u)) * o),
            (e[11] = 0),
            (e[12] = s.x),
            (e[13] = s.y),
            (e[14] = s.z),
            (e[15] = 1),
            (this[i(2259)] = this[i(518)]))
      }
      [t(1414)]() {
        const e = t
        return new Pi({
          color: Cesium[e(1154)].fromCssColorString(this._style[e(1070)]),
          speed: this._style[e(2535)],
          glow: this[e(343)][e(1529)],
          images: this[e(343)][e(1434)],
          viewer: this._viewer,
          shaderType: this[e(1364)],
          cullFaceType: this[e(425)],
          vertexShader: this[e(2530)],
          fragmentShader: this[e(932)],
          outFragColorBody: this[e(140)],
          param: this._param
        })
      }
      [t(1114)]() {
        const e = t,
          i = {}
        return (
          (i[e(995)] = this[e(995)]),
          (i[e(2458)] = this[e(2458)]),
          (i[e(1929)] = this[e(140)]),
          (i[e(406)] = this[e(1364)]),
          (i[e(2088)] = this._lookAtState),
          (i.shaderGeometryType = this[e(1526)]),
          (i.shaderGeometryOpts = this[e(763)]),
          i
        )
      }
      [t(1545)](e) {
        const i = t
        ;(this[i(2462)] = e),
          (this[i(245)] = e[i(245)]),
          (this._geometry = Di[i(1951)](this._shaderGeometryType, this[i(763)])[i(892)]),
          (this._material = this[i(1414)]()),
          (this[i(2122)] = e.id),
          (this[i(2366)] = this._id),
          e[i(245)][i(696)][i(1346)][i(1861)](this)
      }
      [t(2389)](e) {
        const i = t
        e[i(245)][i(696)][i(1346)][i(1896)](this)
      }
    }
    let Fi = {
      readFeature(e) {
        const i = t
        let s = e[i(138)][i(2365)],
          n = e.properties
        return (n.position = s), this[i(1951)](n)
      },
      create(e) {
        const i = t
        switch (e[i(1158)]) {
          case Q[i(2513)]:
            return new ki(e)
          case Q[i(1581)]:
            return new (class extends ki {
              constructor(e) {
                const i = t
                super(e),
                  (this[i(528)] = '着色器特效仿真火焰'),
                  (this._graphicType = Q[i(1581)]),
                  (this[i(1526)] = Ai[i(868)]),
                  (this._shaderType = wi),
                  (this[i(425)] = Si),
                  (this[i(932)] = Ii[i(115)]),
                  (this[i(343)][i(1070)] = Cesium[i(1960)](this[i(343)][i(1070)], 'red')),
                  (this[i(343)][i(2535)] = 1),
                  (this._style[i(1529)] = Cesium.defaultValue(this[i(343)][i(1529)], 1.3)),
                  (this[i(343)][i(2341)] = Cesium[i(1960)](
                    this._style.scale,
                    new Cesium[i(310)](10, 10, 10)
                  ))
              }
            })(e)
          case Q[i(2580)]:
            return new (class extends ki {
              constructor(e) {
                const i = t
                super(e),
                  (this[i(2336)] = Q.shaderEffetColorfulPoint),
                  (this[i(528)] = i(621)),
                  (this[i(1526)] = Ai[i(868)]),
                  (this[i(1364)] = yi),
                  (this[i(425)] = xi),
                  (this[i(932)] = Ii[i(2502)]),
                  (this[i(343)][i(1070)] = Cesium.defaultValue(this[i(343)][i(1070)], i(1626))),
                  (this[i(343)].speed = 1),
                  (this._style[i(1529)] = Cesium[i(1960)](this[i(343)][i(1529)], 1.3)),
                  (this[i(343)][i(2341)] = Cesium[i(1960)](
                    this._style[i(2341)],
                    new Cesium[i(310)](10, 10, 10)
                  ))
              }
            })(e)
          case Q[i(1910)]:
            return new (class extends ki {
              constructor(e) {
                const i = t
                super(e),
                  (this._graphicType = Q[i(1910)]),
                  (this[i(528)] = i(901)),
                  (this[i(1526)] = Ai[i(868)]),
                  (this[i(1364)] = yi),
                  (this[i(425)] = xi),
                  (this[i(932)] = Ii[i(587)]),
                  (this[i(343)][i(1070)] = Cesium[i(1960)](this[i(343)][i(1070)], i(1626))),
                  (this[i(343)][i(2535)] = 1),
                  (this._style.glow = Cesium[i(1960)](this[i(343)][i(1529)], 1.3)),
                  (this[i(343)].scale = Cesium[i(1960)](
                    this[i(343)][i(2341)],
                    new Cesium[i(310)](10, 10, 10)
                  ))
              }
            })(e)
          case Q.shaderEffetConstellationChain:
            return new (class extends ki {
              constructor(e) {
                const i = t
                super(e),
                  (this[i(2336)] = Q.shaderEffetConstellationChain),
                  (this[i(528)] = i(1992)),
                  (this._shaderGeometryType = Ai[i(868)]),
                  (this._shaderType = yi),
                  (this._cullFaceType = xi),
                  (this[i(932)] = Ii[i(211)]),
                  (this[i(343)][i(1070)] = Cesium[i(1960)](this._style[i(1070)], i(1626))),
                  (this._style[i(2535)] = 1),
                  (this._style[i(1529)] = Cesium[i(1960)](this[i(343)].glow, 1.3)),
                  (this[i(343)][i(2341)] = Cesium.defaultValue(
                    this[i(343)][i(2341)],
                    new Cesium.Cartesian3(10, 10, 10)
                  ))
              }
            })(e)
          case Q[i(1629)]:
            return new (class extends ki {
              constructor(e) {
                const i = t
                super(e),
                  (this[i(2336)] = Q[i(1629)]),
                  (this[i(528)] = i(1835)),
                  (this[i(1526)] = Ai.plane),
                  (this[i(1364)] = yi),
                  (this[i(425)] = xi),
                  (this[i(932)] = Ii[i(718)]),
                  (this._style[i(1070)] = Cesium[i(1960)](this._style[i(1070)], i(1626))),
                  (this._style[i(2535)] = 1),
                  (this._style[i(1529)] = Cesium.defaultValue(this._style[i(1529)], 1.3)),
                  (this[i(343)][i(2341)] = Cesium[i(1960)](
                    this[i(343)].scale,
                    new Cesium[i(310)](10, 10, 10)
                  ))
              }
            })(e)
          case Q[i(1820)]:
            return new (class extends ki {
              constructor(e) {
                const i = t
                super(e),
                  (this[i(2336)] = Q[i(1820)]),
                  (this[i(528)] = i(2406)),
                  (this[i(1526)] = Ai[i(868)]),
                  (this[i(1364)] = wi),
                  (this[i(425)] = Si),
                  (this[i(932)] = Ii.CoolBallFS),
                  (this[i(343)].color = Cesium.defaultValue(this[i(343)][i(1070)], i(1626))),
                  (this[i(343)][i(2535)] = 1),
                  (this._style.glow = Cesium.defaultValue(this[i(343)][i(1529)], 1)),
                  (this[i(343)][i(2341)] = Cesium.defaultValue(
                    this[i(343)][i(2341)],
                    new Cesium[i(310)](1, 1, 1)
                  )),
                  (this[i(140)] = i(1123))
              }
            })(e)
          case Q[i(2065)]:
            return new (class extends ki {
              constructor(e) {
                const i = t
                super(e),
                  (this[i(2336)] = Q[i(2065)]),
                  (this[i(528)] = i(873)),
                  (this[i(1526)] = Ai[i(868)]),
                  (this[i(1364)] = wi),
                  (this._cullFaceType = xi),
                  (this[i(932)] = Ii[i(1025)]),
                  (this[i(343)].color = Cesium[i(1960)](this[i(343)].color, i(1626))),
                  (this._style[i(2535)] = 2),
                  (this[i(343)][i(1529)] = Cesium[i(1960)](this._style[i(1529)], 1)),
                  (this[i(343)][i(2341)] = Cesium[i(1960)](
                    this._style[i(2341)],
                    new Cesium.Cartesian3(1, 1, 1)
                  )),
                  (this[i(140)] = i(2484))
              }
            })(e)
          case Q[i(1956)]:
            return new (class extends ki {
              constructor(e) {
                const i = t
                super(e),
                  (this[i(2336)] = Q.shaderEffetFireworks),
                  (this._typeName = i(887)),
                  (this[i(1526)] = Ai.plane),
                  (this._shaderType = wi),
                  (this[i(425)] = Si),
                  (this[i(932)] = Ii[i(671)]),
                  (this[i(343)].color = Cesium.defaultValue(this._style[i(1070)], 'red')),
                  (this[i(343)].speed = 1),
                  (this[i(343)][i(1529)] = Cesium[i(1960)](this[i(343)][i(1529)], 1.3)),
                  (this[i(343)][i(2341)] = Cesium[i(1960)](
                    this[i(343)][i(2341)],
                    new Cesium[i(310)](10, 10, 10)
                  )),
                  (this[i(140)] = i(517))
              }
            })(e)
          case Q.shaderEffetFlameCloud:
            return new (class extends ki {
              constructor(e) {
                const i = t
                super(e),
                  (this[i(2336)] = Q[i(1147)]),
                  (this[i(528)] = i(1050)),
                  (this._shaderGeometryType = Ai[i(868)]),
                  (this._shaderType = wi),
                  (this[i(425)] = xi),
                  (this[i(932)] = Ii.FlameCloudFS),
                  (this[i(1630)] = new Cesium.Cartesian4(0, 1, 0, 0)),
                  (this[i(343)][i(1070)] = Cesium[i(1960)](this[i(343)][i(1070)], i(1153))),
                  (this[i(343)][i(2535)] = 4),
                  (this[i(343)][i(1529)] = Cesium[i(1960)](this[i(343)][i(1529)], 1)),
                  (this._style[i(2341)] = Cesium[i(1960)](
                    this[i(343)].scale,
                    new Cesium[i(310)](10, 10, 10)
                  ))
              }
            })(e)
          case Q[i(2175)]:
            return new (class extends ki {
              constructor(e) {
                const i = t
                super(e),
                  (this[i(2336)] = Q.shaderEffetFlameRing),
                  (this[i(528)] = '着色器特效火圈'),
                  (this[i(1526)] = Ai[i(868)]),
                  (this[i(1364)] = wi),
                  (this[i(425)] = xi),
                  (this[i(932)] = Ii[i(1976)]),
                  (this[i(1630)] = new Cesium[i(2141)](0, 1, 0, 0)),
                  (this._style.color = Cesium.defaultValue(this._style[i(1070)], i(1153))),
                  (this._style.speed = 4),
                  (this._style[i(1529)] = Cesium[i(1960)](this[i(343)].glow, 1)),
                  (this[i(343)].scale = Cesium[i(1960)](
                    this[i(343)].scale,
                    new Cesium[i(310)](10, 10, 10)
                  )),
                  (this[i(979)] = !1)
              }
            })(e)
          case Q[i(2014)]:
            return new (class extends ki {
              constructor(e) {
                const i = t
                super(e),
                  (this[i(2336)] = Q[i(2014)]),
                  (this[i(528)] = i(706)),
                  (this[i(1526)] = Ai[i(868)]),
                  (this[i(1364)] = yi),
                  (this._cullFaceType = xi),
                  (this[i(932)] = Ii[i(724)]),
                  (this[i(1630)] = new Cesium.Cartesian4(0, 1, 0, 0)),
                  (this[i(343)].color = Cesium[i(1960)](this[i(343)][i(1070)], 'rgba(0,0,0,0)')),
                  (this[i(343)][i(2535)] = 1),
                  (this[i(343)][i(1529)] = Cesium[i(1960)](this[i(343)].glow, 1)),
                  (this[i(343)][i(2341)] = Cesium.defaultValue(
                    this[i(343)][i(2341)],
                    new Cesium[i(310)](10, 10, 10)
                  ))
              }
            })(e)
          case Q[i(1553)]:
            return new (class extends ki {
              constructor(e) {
                const i = t
                super(e),
                  (this[i(2336)] = Q[i(1553)]),
                  (this[i(528)] = i(345)),
                  (this._shaderGeometryType = Ai.plane),
                  (this[i(1364)] = yi),
                  (this._cullFaceType = xi),
                  (this._fragmentShader = Ii[i(1571)]),
                  (this[i(343)].color = Cesium[i(1960)](this[i(343)][i(1070)], i(1626))),
                  (this[i(343)].speed = 1),
                  (this[i(343)][i(1529)] = Cesium.defaultValue(this[i(343)].glow, 1.3)),
                  (this[i(343)].scale = Cesium.defaultValue(
                    this[i(343)][i(2341)],
                    new Cesium[i(310)](10, 10, 10)
                  ))
              }
            })(e)
          case Q[i(2498)]:
            return new (class extends ki {
              constructor(e) {
                const i = t
                super(e),
                  (this[i(2336)] = Q.shaderEffetGlowPoint),
                  (this[i(528)] = i(191)),
                  (this[i(1526)] = Ai[i(868)]),
                  (this[i(1364)] = yi),
                  (this[i(425)] = xi),
                  (this[i(932)] = Ii[i(1618)]),
                  (this[i(343)][i(1070)] = Cesium[i(1960)](this[i(343)][i(1070)], i(1626))),
                  (this[i(343)].speed = 1),
                  (this[i(343)][i(1529)] = Cesium[i(1960)](this._style[i(1529)], 1.3)),
                  (this._style[i(2341)] = Cesium[i(1960)](
                    this[i(343)][i(2341)],
                    new Cesium.Cartesian3(10, 10, 10)
                  ))
              }
            })(e)
          case Q[i(1155)]:
            return new (class extends ki {
              constructor(e) {
                const i = t
                super(e),
                  (this[i(528)] = i(1838)),
                  (this[i(2336)] = Q[i(1155)]),
                  (this[i(1526)] = Ai[i(868)]),
                  (this[i(1364)] = yi),
                  (this[i(425)] = xi),
                  (this[i(932)] = Ii[i(524)]),
                  (this._style.color = Cesium.defaultValue(this[i(343)].color, i(1626))),
                  (this._style.speed = 1),
                  (this[i(343)][i(1529)] = Cesium[i(1960)](this[i(343)][i(1529)], 1.3)),
                  (this[i(343)][i(2341)] = Cesium.defaultValue(
                    this[i(343)][i(2341)],
                    new Cesium[i(310)](10, 10, 10)
                  ))
              }
            })(e)
          case Q[i(1301)]:
            return new (class extends ki {
              constructor(e) {
                const i = t
                super(e),
                  (this._graphicType = Q[i(1301)]),
                  (this[i(528)] = '着色器特效渐变环'),
                  (this[i(1526)] = Ai[i(868)]),
                  (this[i(1364)] = yi),
                  (this[i(425)] = xi),
                  (this[i(932)] = Ii[i(129)]),
                  (this[i(343)].color = Cesium[i(1960)](this[i(343)][i(1070)], i(1626))),
                  (this[i(343)].speed = 1),
                  (this[i(343)][i(1529)] = Cesium[i(1960)](this[i(343)].glow, 1.3)),
                  (this._style[i(2341)] = Cesium[i(1960)](
                    this._style[i(2341)],
                    new Cesium.Cartesian3(10, 10, 10)
                  )),
                  (this._lookAtState = !1)
              }
            })(e)
          case Q.shaderEffetMagicBall:
            return new (class extends ki {
              constructor(e) {
                const i = t
                super(e),
                  (this[i(528)] = i(2356)),
                  (this[i(2336)] = Q[i(651)]),
                  (this[i(1526)] = Ai[i(868)]),
                  (this._shaderType = wi),
                  (this[i(425)] = Si),
                  (this[i(932)] = Ii[i(2056)]),
                  (this[i(343)][i(1070)] = Cesium[i(1960)](this._style[i(1070)], i(1626))),
                  (this[i(343)][i(2535)] = 1),
                  (this[i(343)][i(1529)] = Cesium.defaultValue(this[i(343)][i(1529)], 1.3)),
                  (this._style[i(2341)] = Cesium[i(1960)](
                    this[i(343)][i(2341)],
                    new Cesium[i(310)](10, 10, 10)
                  )),
                  (this[i(140)] = i(783))
              }
            })(e)
          case Q[i(1534)]:
            return new (class extends ki {
              constructor(e) {
                const i = t
                super(e),
                  (this[i(528)] = i(1798)),
                  (this._graphicType = Q.shaderEffetMagicBox),
                  (this[i(1526)] = Ai[i(868)]),
                  (this._shaderType = wi),
                  (this._cullFaceType = xi),
                  (this._fragmentShader = Ii.MagicBoxFS),
                  (this[i(343)][i(1070)] = Cesium[i(1960)](this[i(343)][i(1070)], 'red')),
                  (this[i(343)].speed = 1),
                  (this[i(343)][i(1529)] = Cesium.defaultValue(this._style[i(1529)], 1.3)),
                  (this._style[i(2341)] = Cesium[i(1960)](
                    this[i(343)].scale,
                    new Cesium[i(310)](10, 10, 10)
                  )),
                  (this._outFragColorBody = i(1280))
              }
            })(e)
          case Q[i(1171)]:
            return new (class extends ki {
              constructor(e) {
                const i = t
                super(e),
                  (this[i(528)] = '着色器特效魔法环'),
                  (this._graphicType = Q[i(1171)]),
                  (this[i(1526)] = Ai[i(868)]),
                  (this._shaderType = yi),
                  (this[i(425)] = xi),
                  (this._fragmentShader = Ii[i(2374)]),
                  (this[i(343)][i(1070)] = Cesium[i(1960)](this[i(343)][i(1070)], i(1626))),
                  (this[i(343)][i(2535)] = 1),
                  (this._style[i(1529)] = Cesium[i(1960)](this._style[i(1529)], 1.3)),
                  (this._style[i(2341)] = Cesium.defaultValue(
                    this[i(343)].scale,
                    new Cesium[i(310)](10, 10, 10)
                  ))
              }
            })(e)
          case Q[i(797)]:
            return new (class extends ki {
              constructor(e) {
                const i = t
                super(e),
                  (this[i(528)] = i(286)),
                  (this._graphicType = Q[i(797)]),
                  (this[i(1526)] = Ai.plane),
                  (this[i(1364)] = wi),
                  (this[i(425)] = xi),
                  (this[i(932)] = Ii[i(2428)]),
                  (this[i(1630)] = new Cesium.Cartesian4(0, 1, 0, 0)),
                  (this._style.color = Cesium[i(1960)](this[i(343)][i(1070)], i(1153))),
                  (this[i(343)][i(2535)] = 1),
                  (this._style[i(1529)] = Cesium[i(1960)](this[i(343)][i(1529)], 1)),
                  (this[i(343)][i(2341)] = Cesium[i(1960)](
                    this[i(343)][i(2341)],
                    new Cesium[i(310)](10, 10, 10)
                  ))
              }
            })(e)
          case Q[i(628)]:
            return new (class extends ki {
              constructor(e) {
                const i = t
                super(e),
                  (this[i(528)] = i(2348)),
                  (this[i(2336)] = Q[i(628)]),
                  (this[i(1526)] = Ai[i(868)]),
                  (this[i(1364)] = yi),
                  (this[i(425)] = xi),
                  (this._fragmentShader = Ii.NeonPointFS),
                  (this[i(343)][i(1070)] = Cesium.defaultValue(this._style[i(1070)], i(1626))),
                  (this[i(343)][i(2535)] = 1),
                  (this[i(343)][i(1529)] = Cesium.defaultValue(this[i(343)][i(1529)], 1.3)),
                  (this[i(343)].scale = Cesium[i(1960)](
                    this[i(343)][i(2341)],
                    new Cesium[i(310)](10, 10, 10)
                  ))
              }
            })(e)
          case Q[i(243)]:
            return new (class extends ki {
              constructor(e) {
                const i = t
                super(e),
                  (this[i(528)] = '着色器特效花瓣'),
                  (this[i(2336)] = Q.shaderEffetPetal),
                  (this[i(1526)] = Ai[i(868)]),
                  (this[i(1364)] = yi),
                  (this[i(425)] = xi),
                  (this[i(932)] = Ii[i(505)]),
                  (this[i(343)][i(1070)] = Cesium.defaultValue(this._style[i(1070)], i(1626))),
                  (this[i(343)][i(2535)] = 1),
                  (this[i(343)][i(1529)] = Cesium.defaultValue(this[i(343)].glow, 1.3)),
                  (this._style[i(2341)] = Cesium[i(1960)](
                    this[i(343)].scale,
                    new Cesium.Cartesian3(10, 10, 10)
                  ))
              }
            })(e)
          case Q[i(2167)]:
            return new (class extends ki {
              constructor(e) {
                const i = t
                super(e),
                  (this[i(528)] = '着色器特效护盾'),
                  (this[i(2336)] = Q[i(2167)]),
                  (this[i(1526)] = Ai[i(1565)]),
                  (this[i(763)][i(1141)] = Cesium[i(1960)](
                    this._shaderGeometryOpts.radii,
                    new Cesium[i(310)](10, 10, 10)
                  )),
                  (this[i(1364)] = yi),
                  (this[i(425)] = Si),
                  (this[i(932)] = Ii[i(2385)]),
                  (this._param = new Cesium.Cartesian4(0, 1, 0, 0)),
                  (this[i(343)][i(1070)] = Cesium.defaultValue(this[i(343)][i(1070)], i(1153))),
                  (this._style.speed = 1),
                  (this[i(343)][i(1529)] = Cesium[i(1960)](this[i(343)][i(1529)], 1)),
                  (this[i(343)][i(2341)] = Cesium.defaultValue(
                    this._style[i(2341)],
                    new Cesium[i(310)](1, 1, 1)
                  )),
                  (this[i(979)] = !1)
              }
            })(e)
          case Q.shaderEffetRadarScan:
            return new (class extends ki {
              constructor(e) {
                const i = t
                super(e),
                  (this._typeName = i(338)),
                  (this[i(2336)] = Q[i(1140)]),
                  (this[i(1526)] = Ai[i(868)]),
                  (this[i(1364)] = wi),
                  (this[i(425)] = Si),
                  (this[i(932)] = Ii[i(663)]),
                  (this[i(343)][i(1070)] = Cesium[i(1960)](this._style[i(1070)], i(1626))),
                  (this[i(343)][i(2535)] = 1),
                  (this[i(343)][i(1529)] = Cesium[i(1960)](this._style[i(1529)], 1.3)),
                  (this._style[i(2341)] = Cesium[i(1960)](
                    this._style[i(2341)],
                    new Cesium[i(310)](10, 10, 10)
                  )),
                  (this._lookAtState = !1)
              }
            })(e)
          case Q.shaderEffetGlowDiamond:
            return new (class extends ki {
              constructor(e) {
                const i = t
                super(e),
                  (this[i(2336)] = Q[i(1322)]),
                  (this[i(528)] = i(2115)),
                  (this[i(1526)] = Ai[i(868)]),
                  (this[i(1364)] = yi),
                  (this[i(425)] = xi),
                  (this[i(932)] = Ii[i(1223)]),
                  (this[i(1630)] = new Cesium.Cartesian4(0, 1, 0, 0)),
                  (this[i(343)][i(1070)] = Cesium.defaultValue(this[i(343)].color, i(1153))),
                  (this._style[i(2535)] = 1),
                  (this[i(343)].glow = Cesium[i(1960)](this[i(343)][i(1529)], 1)),
                  (this[i(343)][i(2341)] = Cesium[i(1960)](
                    this[i(343)][i(2341)],
                    new Cesium[i(310)](10, 10, 10)
                  ))
              }
            })(e)
          case Q[i(2188)]:
            return new (class extends ki {
              constructor(e) {
                const i = t
                super(e),
                  (this[i(528)] = '着色器特效闪烁环'),
                  (this[i(2336)] = Q[i(2188)]),
                  (this[i(1526)] = Ai[i(868)]),
                  (this._shaderType = yi),
                  (this._cullFaceType = xi),
                  (this[i(932)] = Ii.ShinyRingFS),
                  (this._style[i(1070)] = Cesium[i(1960)](this[i(343)][i(1070)], i(1626))),
                  (this[i(343)][i(2535)] = 1),
                  (this[i(343)][i(1529)] = Cesium.defaultValue(this[i(343)][i(1529)], 1.3)),
                  (this._style.scale = Cesium.defaultValue(
                    this._style[i(2341)],
                    new Cesium[i(310)](10, 10, 10)
                  )),
                  (this._lookAtState = !1)
              }
            })(e)
          case Q[i(928)]:
            return new (class extends ki {
              constructor(e) {
                const i = t
                super(e),
                  (this._typeName = i(1871)),
                  (this[i(2336)] = Q[i(928)]),
                  (this[i(1526)] = Ai[i(868)]),
                  (this[i(1364)] = wi),
                  (this._cullFaceType = Si),
                  (this[i(932)] = Ii[i(1564)]),
                  (this[i(343)].color = Cesium[i(1960)](this._style[i(1070)], i(1626))),
                  (this[i(343)][i(2535)] = 1),
                  (this._style[i(1529)] = Cesium[i(1960)](this._style[i(1529)], 1.3)),
                  (this[i(343)][i(2341)] = Cesium[i(1960)](
                    this[i(343)][i(2341)],
                    new Cesium[i(310)](10, 10, 10)
                  )),
                  (this._outFragColorBody =
                    'out_FragColor = vec4(out_FragColor.rgb + u_color + u_color, out_FragColor.r + out_FragColor.r);')
              }
            })(e)
          case Q[i(1183)]:
            return new (class extends ki {
              constructor(e) {
                const i = t
                super(e),
                  (this[i(528)] = i(2509)),
                  (this._graphicType = Q.shaderEffetSwim),
                  (this[i(1526)] = Ai[i(868)]),
                  (this[i(1364)] = yi),
                  (this[i(425)] = xi),
                  (this[i(932)] = Ii.SwimFS),
                  (this[i(343)].color = Cesium[i(1960)](this[i(343)][i(1070)], 'red')),
                  (this[i(343)][i(2535)] = 1),
                  (this[i(343)][i(1529)] = Cesium[i(1960)](this[i(343)].glow, 1.3)),
                  (this[i(343)][i(2341)] = Cesium[i(1960)](
                    this[i(343)][i(2341)],
                    new Cesium[i(310)](10, 10, 10)
                  ))
              }
            })(e)
          case Q[i(2136)]:
            return new (class extends ki {
              constructor(e) {
                const i = t
                super(e),
                  (this._typeName = i(501)),
                  (this[i(2336)] = Q.shaderEffetVirus),
                  (this[i(1526)] = Ai[i(868)]),
                  (this._shaderType = yi),
                  (this._cullFaceType = xi),
                  (this[i(932)] = Ii.VirusFS),
                  (this[i(343)][i(1070)] = Cesium[i(1960)](this[i(343)].color, i(1626))),
                  (this[i(343)][i(2535)] = 1),
                  (this[i(343)][i(1529)] = Cesium[i(1960)](this[i(343)][i(1529)], 1.3)),
                  (this[i(343)][i(2341)] = Cesium[i(1960)](
                    this[i(343)][i(2341)],
                    new Cesium[i(310)](10, 10, 10)
                  ))
              }
            })(e)
          case Q.shaderEffetWavePetals:
            return new (class extends ki {
              constructor(e) {
                const i = t
                super(e),
                  (this._typeName = i(1037)),
                  (this._graphicType = Q[i(870)]),
                  (this[i(1526)] = Ai[i(868)]),
                  (this._shaderType = wi),
                  (this[i(425)] = Si),
                  (this[i(932)] = Ii.WavePetalsFS),
                  (this._style[i(1070)] = Cesium.defaultValue(this[i(343)].color, i(1626))),
                  (this[i(343)][i(2535)] = 1),
                  (this[i(343)][i(1529)] = Cesium.defaultValue(this._style[i(1529)], 1.3)),
                  (this[i(343)].scale = Cesium.defaultValue(
                    this[i(343)].scale,
                    new Cesium[i(310)](10, 10, 10)
                  )),
                  (this[i(140)] = i(1715))
              }
            })(e)
          case Q.shaderEffetGlowPyramid:
            return new (class extends ki {
              constructor(e) {
                const i = t
                super(e),
                  (this._graphicType = Q.shaderEffetGlowPyramid),
                  (this[i(528)] = i(717)),
                  (this[i(1526)] = Ai[i(868)]),
                  (this[i(1364)] = yi),
                  (this[i(425)] = xi),
                  (this[i(932)] = Ii[i(890)]),
                  (this[i(343)].color = Cesium[i(1960)](this[i(343)][i(1070)], 'red')),
                  (this[i(343)][i(2535)] = 1),
                  (this._style[i(1529)] = Cesium.defaultValue(this[i(343)].glow, 1.3)),
                  (this[i(343)].scale = Cesium[i(1960)](
                    this[i(343)].scale,
                    new Cesium[i(310)](10, 10, 10)
                  ))
              }
            })(e)
        }
      }
    }
    const Ri = 'water_0',
      Li = 'water_1'
    class Oi {
      constructor(e = {}) {
        const i = t
        ;(this[i(2086)] = e), (this._name = i(594)), (this[i(1265)] = null)
      }
      get [t(892)]() {
        return this[t(1265)]
      }
      get [t(981)]() {
        return this[t(2086)]
      }
      get [t(1916)]() {
        return this[t(1165)]
      }
      set [t(1916)](e) {
        this[t(1165)] = e
      }
      [t(1180)]() {
        const e = t
        for (const t in this[e(2086)])
          if (Object[e(782)][e(1669)](this[e(2086)], t)) {
            const e = this._appearanceOpts[t]
            this[t] = e
          }
        ;(this[e(2086)].name = this[e(1165)]), this._mergeOpts()
      }
      [t(1078)]() {}
      [t(589)](t) {}
    }
    const Bi = {
      WaterAppearance_0: class extends Oi {
        constructor(e) {
          const i = t
          super(e), (this[i(1265)] = this._createAppearance()), this[i(589)](e)
        }
        setOpts(e) {
          const i = t
          ;(e[i(2095)] = Cesium.defaultValue(e[i(2095)], i(275))),
            (e[i(1062)] = Cesium.defaultValue(e[i(1062)], Cesium[i(1643)](i(1673)))),
            (e[i(553)] = Cesium[i(1960)](e[i(553)], 1e3)),
            (e[i(1769)] = Cesium[i(1960)](e[i(1769)], 0.03)),
            (e.amplitude = Cesium[i(1960)](e[i(2550)], 10)),
            (e[i(849)] = Cesium[i(1960)](e.specularIntensity, 5)),
            (this[i(2086)] = e),
            this[i(1180)]()
        }
        [t(1078)]() {
          const e = t
          let i = this[e(1265)].material[e(1285)]
          ;(i[e(2095)] = Cesium[e(1154)][e(2008)](this[e(2095)])),
            (i.frequency = this[e(553)]),
            (i.animationSpeed = this[e(1769)]),
            (i[e(2550)] = this[e(2550)]),
            (i[e(849)] = this[e(849)])
        }
        [t(1530)]() {
          const e = t
          return new Cesium[e(2305)]({
            material: new Cesium[e(1637)]({
              fabric: {
                type: e(2562),
                uniforms: {
                  baseWaterColor: Cesium.Color[e(2147)][e(329)](0.3),
                  normalMap: Cesium[e(1643)]('Assets/Textures/waterNormals.jpg'),
                  frequency: 1e3,
                  animationSpeed: 0.03,
                  amplitude: 10,
                  specularIntensity: 5
                }
              }
            }),
            translucent: !0
          })
        }
      },
      WaterAppearance_1: class extends Oi {
        constructor(e) {
          const i = t
          super(e), (this._proxy = this[i(1530)]()), this[i(589)](e)
        }
        [t(589)](e) {
          const i = t
          ;(e[i(849)] = Cesium[i(1960)](e.specularIntensity, 5)),
            (this[i(2086)] = e),
            this[i(1180)]()
        }
        [t(1078)]() {
          const e = t
          this[e(1265)][e(322)][e(1285)][e(849)] = this.specularIntensity
        }
        [t(1530)]() {
          const e = t,
            i = {}
          ;(i[e(1285)] = {}), (i[e(660)] = e(946))
          const s = {}
          return (
            (s[e(1061)] = i),
            new Cesium[e(985)]({
              material: new Cesium[e(1637)](s),
              translucent: !0,
              translucent: !0
            })
          )
        }
      }
    }
    class Vi {
      static create(e = t(1070), i = {}) {
        switch (e) {
          case Ri:
            return new Bi.WaterAppearance_0(i)
          case Li:
            return new Bi.WaterAppearance_1(i)
        }
      }
    }
    class Ni extends b {
      constructor(e) {
        const i = t
        super(e),
          (this[i(314)] = 'Polygon'),
          (this[i(1115)] = Y),
          (this._height = Cesium[i(1960)](e[i(2306)], null)),
          (this[i(343)].height = Cesium[i(1960)](this[i(343)].height, null)),
          (this[i(2306)] = this[i(343)][i(2306)]),
          (this[i(343)][i(1337)] = Cesium[i(1960)](this._style[i(1337)], i(1775))),
          (this[i(343)][i(981)] = Cesium.defaultValue(this[i(343)][i(981)], {})),
          (this[i(343)].flowDir = Cesium[i(1960)](this[i(343)][i(172)], 0))
        let s = Cesium[i(1960)](e[i(2333)], [])
        this[i(2515)](s)
      }
      get [t(2306)]() {
        return this[t(1516)]
      }
      set height(e) {
        const i = t
        if (!e) return
        if (
          ((this[i(1516)] = e),
          this[i(2469)][i(1602)]((t) => {
            t[2] = this._height
          }),
          !this[i(1934)])
        )
          return
        let s = Cesium[i(2285)][i(2579)](this._originalWorldPosition),
          n = Cesium[i(310)][i(1069)](s[i(2106)], s[i(199)], e),
          o = Cesium[i(310)][i(1676)](n, this[i(964)], new Cesium[i(310)]())
        ;(this[i(1934)][i(223)] = Cesium.Matrix4[i(1420)](o, new Cesium[i(2066)]())),
          (this[i(1449)] = n),
          (this[i(2446)] = Cesium[i(1148)][i(276)].geodeticSurfaceNormal(this[i(1449)])),
          (this[i(343)][i(2306)] = e)
      }
      [t(2231)]() {
        const e = t
        this[e(343)].height != this[e(1516)] && this[e(2520)] && this[e(1920)] && this[e(1920)](),
          (this[e(2306)] = this[e(343)][e(2306)])
      }
      _setVisible(e) {
        const i = t
        this[i(1934)] && (this._primitive[i(1482)] = e)
      }
      [t(2515)](e) {
        const i = t
        if (((this[i(2469)] = e), e.length < 3)) return
        ;(this[i(794)] = e),
          (this[i(1378)] = []),
          (this._cartesian3Array = Cesium[i(310)][i(774)]([][i(1500)].apply([], e))),
          (this[i(1024)] = Cesium[i(1242)][i(951)](this[i(1378)]))
        let s = e[i(277)],
          n = new Cesium[i(310)]()
        this._cartesian3Array.forEach((t) => {
          const e = i
          n = Cesium[e(310)][e(1861)](n, t, n)
        }),
          (n = Cesium[i(310)][i(606)](n, 1 / s, n)),
          (this[i(1449)] = n),
          (this[i(964)] = n[i(1902)]()),
          this[i(1733)]()
      }
      _reAdd() {
        const e = t
        this[e(1934)] &&
          this._scene &&
          (this[e(2018)][e(1346)][e(1896)](this[e(1934)]),
          (this[e(1934)] = this._createPrimitive()),
          (this[e(1934)][e(2122)] = this[e(2462)].id),
          (this[e(1934)][e(2366)] = this[e(1570)]),
          this[e(2018)][e(1346)].add(this[e(1934)]))
      }
      _createPrimitive() {
        const e = t
        return new Cesium[e(2186)]({
          geometryInstances: new Cesium.GeometryInstance({
            geometry: new Cesium[e(1765)]({
              polygonHierarchy: new Cesium.PolygonHierarchy(this[e(1378)]),
              perPositionHeight: !0,
              stRotation: Cesium[e(475)][e(1149)](this._style[e(172)]),
              closeTop: !0,
              closeBottom: !1,
              vertexFormat: Cesium[e(1788)][e(1404)]
            })
          }),
          appearance: this[e(1075)][e(892)] ? this[e(1075)][e(892)] : this._appearance,
          asynchronous: !1
        })
      }
      [t(1545)](e) {
        const i = t
        ;(this[i(2462)] = e),
          (this._scene = e._viewer[i(696)]),
          (this[i(1075)] = this[i(1407)]()),
          (this[i(1934)] = this._createPrimitive()),
          this[i(2018)][i(1346)][i(1861)](this[i(1934)]),
          (this[i(1934)][i(2122)] = e.id),
          (this._primitive[i(2366)] = this._id),
          this[i(856)]()
      }
      [t(856)]() {}
      [t(2389)](e) {
        const i = t
        this[i(2018)].primitives[i(1896)](this._primitive), this[i(434)]()
      }
      [t(434)]() {}
      [t(1407)]() {
        const e = t
        return (
          (this[e(595)] = Vi[e(1951)](this[e(343)][e(1337)], this[e(343)][e(981)])),
          this[e(595)][e(892)]
        )
      }
      _beforDestroy() {
        this[t(2018)] = null
      }
    }
    class Hi extends Ni {
      constructor(e) {
        const i = t
        super(e),
          this._overWrite(),
          (this[i(343)][i(1070)] = Cesium.defaultValue(this[i(343)][i(1070)], i(2380))),
          (this[i(343)][i(506)] = Cesium.defaultValue(this._style[i(506)], i(2355)))
        let s = Cesium[i(1960)](e[i(2333)], [])
        this.setPositions(s)
      }
      [t(2515)](e) {
        const i = t
        if (((this[i(2469)] = e), e.length < 3)) return
        ;(this[i(794)] = e), (this[i(1378)] = [])
        let s = e[i(277)],
          n = new Cesium[i(310)]()
        e[i(1602)]((t) => {
          const e = i
          t[2] = this._height
          const s = Cesium[e(310)].fromDegrees(t[0], t[1], this._height)
          this._cartesian3Array[e(2553)](s), (n = Cesium[e(310)].add(n, s, n))
        }),
          (n = Cesium[i(310)].multiplyByScalar(n, 1 / s, n)),
          (this[i(1449)] = n),
          (this[i(964)] = n[i(1902)]()),
          (this[i(2446)] = Cesium[i(1148)][i(276)].geodeticSurfaceNormal(this[i(1449)])),
          (this[i(1024)] = Cesium[i(1242)].fromPoints(this._cartesian3Array)),
          (this[i(2446)] = Cesium[i(1148)][i(276)].geodeticSurfaceNormal(this[i(1449)])),
          (this[i(2280)] = Cesium[i(2314)][i(1163)](this._centerWorldPosition, this[i(2446)])),
          (this[i(997)] = new Cesium[i(2066)](
            -2 * this[i(2280)].normal.x * this[i(2280)].normal.x + 1,
            -2 * this._waterPlane[i(226)].x * this[i(2280)][i(226)].y,
            -2 * this[i(2280)].normal.x * this[i(2280)][i(226)].z,
            -2 * this._waterPlane[i(226)].x * this[i(2280)][i(1849)],
            -2 * this[i(2280)].normal.y * this[i(2280)][i(226)].x,
            -2 * this[i(2280)][i(226)].y * this[i(2280)][i(226)].y + 1,
            -2 * this[i(2280)][i(226)].y * this[i(2280)][i(226)].z,
            -2 * this[i(2280)].normal.y * this[i(2280)][i(1849)],
            -2 * this._waterPlane[i(226)].z * this[i(2280)][i(226)].x,
            -2 * this[i(2280)][i(226)].z * this._waterPlane[i(226)].y,
            -2 * this[i(2280)][i(226)].z * this[i(2280)][i(226)].z + 1,
            -2 * this[i(2280)][i(226)].z * this[i(2280)][i(1849)],
            0,
            0,
            0,
            1
          )),
          (this[i(639)] = Cesium[i(2066)][i(1393)].clone()),
          (this[i(2131)] = Cesium[i(2066)][i(1393)][i(1902)]()),
          this[i(1733)]()
      }
      _createFramebuffer(e, i, s, n) {
        const o = t
        let r = this[o(840)]
        if (Cesium[o(2330)](r) && r[o(575)] === i && r[o(2306)] === s && this._hdr === n) return
        this[o(1515)](), (this[o(867)] = n)
        const a = n
          ? e.halfFloatingPointTexture
            ? Cesium[o(547)].HALF_FLOAT
            : Cesium[o(547)][o(2441)]
          : Cesium[o(547)].UNSIGNED_BYTE
        ;(this[o(840)] = new Cesium[o(1641)]({
          context: e,
          width: i,
          height: s,
          pixelFormat: Cesium.PixelFormat[o(810)],
          pixelDatatype: a,
          sampler: new Cesium[o(2203)]({
            wrapS: Cesium[o(917)][o(1558)],
            wrapT: Cesium.TextureWrap[o(1558)],
            minificationFilter: Cesium.TextureMinificationFilter[o(2523)],
            magnificationFilter: Cesium[o(500)][o(2523)]
          })
        })),
          (this[o(257)] = new Cesium.Texture({
            context: e,
            width: i,
            height: s,
            pixelFormat: Cesium[o(2277)][o(1859)],
            pixelDatatype: Cesium.PixelDatatype.UNSIGNED_INT_24_8
          })),
          (this._colorFramebuffer = new Cesium.Framebuffer({
            context: e,
            colorTextures: [this[o(840)]],
            depthStencilTexture: this[o(257)],
            destroyAttachments: !1
          }))
      }
      [t(1515)]() {
        const e = t
        this[e(840)] && this[e(840)][e(2039)](),
          this._depthStencilTexture && this[e(257)][e(2039)](),
          this[e(2397)] && this._colorFramebuffer[e(2039)](),
          (this._colorTexture = void 0),
          (this._depthStencilTexture = void 0),
          (this[e(2397)] = void 0)
      }
      [t(1568)]() {
        const e = t
        ;(Cesium[e(312)][e(1727)][e(334)] = function (t) {
          const i = e
          Cesium[i(2066)][i(1902)](Cesium[i(1960)](t[i(1297)], t.projectionMatrix), this[i(1311)]),
            (this[i(1726)] = !0),
            (this[i(1514)] = !0),
            (this[i(688)] = !0),
            (this[i(1490)] = !0),
            (this._modelViewProjectionRelativeToEyeDirty = !0),
            Cesium[i(2330)](t[i(1058)]) &&
              (Cesium.Matrix4[i(1902)](t[i(1058)], this[i(2560)]),
              (this._modelViewInfiniteProjectionDirty = !0)),
            (this._currentFrustum.x = t[i(1292)]),
            (this[i(1975)].y = t.far),
            (this[i(1232)] = t[i(348)] - t[i(1292)] + 1),
            (this[i(403)] = Cesium.Math[i(2260)](this[i(1232)])),
            (this[i(559)] = 1 / this[i(403)]),
            Cesium[i(2330)](t[i(935)]) && (t = t._offCenterFrustum),
            (this[i(2209)].x = t.top),
            (this[i(2209)].y = t[i(151)]),
            (this[i(2209)].z = t[i(1709)]),
            (this._frustumPlanes.w = t[i(1883)])
        }),
          (Cesium[e(1484)][e(1727)][e(1902)] = function (t) {
            const i = e
            return (
              !Cesium[i(2330)](t) && (t = new Cesium.PerspectiveFrustum()),
              (t[i(124)] = this.aspectRatio),
              (t[i(1303)] = this[i(1303)]),
              (t[i(1292)] = this[i(1292)]),
              (t[i(348)] = this[i(348)]),
              (t[i(2046)] = void 0),
              (t[i(564)] = void 0),
              (t[i(1235)] = void 0),
              (t[i(2401)] = void 0),
              this._offCenterFrustum[i(1902)](t[i(935)]),
              (t[i(1297)] = this.customProjectionMatrix),
              t
            )
          })
      }
    }
    function Gi(e, i) {
      const s = t
      let n = i.clone(),
        o = e.clone(),
        r = 2 * Cesium[s(310)][s(1872)](e, i)
      return Cesium[s(310)][s(606)](i, r, n), Cesium.Cartesian3[s(1676)](e, n, o)
    }
    const Wi = {}
    Wi[t(134)] = 0
    const Ui = new Cesium[t(2461)](Wi)
    let ji = new Cesium[t(1154)]()
    const qi = {}
    qi[t(134)] = 0
    const Yi = new Cesium[t(2461)](qi)
    let Xi = new Cesium[t(1154)]()
    let Qi = {
      readFeature(e) {
        const i = t
        let s = e[i(138)][i(2365)],
          n = e[i(1004)]
        return (n[i(2333)] = s), this[i(1951)](n)
      },
      create(e) {
        const i = t
        switch (e[i(1158)]) {
          case Q[i(1213)]:
            return new (class extends Ni {
              constructor(e) {
                const i = t
                super(e), (this[i(2336)] = Q[i(1213)])
              }
              [t(2231)]() {
                const e = t
                let i = this[e(343)][e(981)]
                this._appearanceProxy[e(589)](i),
                  this._style[e(2306)] != this[e(1516)] &&
                    this._editMode &&
                    this[e(1920)] &&
                    this[e(1920)](),
                  (this[e(2306)] = this._style.height)
              }
            })(e)
          case Q.reflectWater:
            return new (class extends Hi {
              constructor(e) {
                const i = t
                super(e),
                  (this[i(2336)] = Q[i(2417)]),
                  (this._typeName = i(623)),
                  (this[i(343)][i(757)] = Cesium.defaultValue(this._style.rippleSize, 50)),
                  (this[i(343)].waterAlpha = Cesium[i(1960)](this[i(343)][i(1168)], 0.9)),
                  (this[i(343)][i(386)] = Cesium[i(1960)](this[i(343)][i(386)], 0.1)),
                  (this[i(343)][i(254)] = Cesium[i(1960)](this[i(343)][i(254)], 100)),
                  (this._style[i(1248)] = Cesium[i(1960)](this[i(343)][i(1248)], 3.7)),
                  (this[i(1506)] = Cesium.defaultValue(e[i(1993)], !1))
              }
              [t(2383)](e) {
                const i = t
                let s = new Cesium.Cartesian3(0, 0, -1),
                  n = new Cesium[i(310)]()
                this[i(1557)] = Cesium[i(2162)][i(1902)](e, this[i(1557)])
                const o = e[i(1236)][i(1902)]()
                let r = Cesium[i(310)].subtract(this[i(1449)], o, new Cesium[i(310)]())
                if (Cesium[i(310)][i(1872)](r, this._normal) > 0) return !1
                ;(r = Gi(r, this[i(2446)])),
                  Cesium[i(310)][i(875)](r, r),
                  Cesium.Cartesian3[i(1861)](r, this[i(1449)], r),
                  (this[i(1557)].position = r.clone()),
                  Cesium[i(310)][i(1861)](e[i(1258)], o, s),
                  Cesium[i(310)][i(1676)](this[i(1449)], s, n),
                  (n = Gi(n, this._normal)),
                  Cesium.Cartesian3[i(875)](n, n),
                  Cesium[i(310)].add(n, this[i(1449)], n),
                  (this[i(1557)][i(1657)] = Cesium[i(310)].subtract(
                    n,
                    this[i(1557)][i(2251)],
                    new Cesium.Cartesian3()
                  )),
                  Cesium[i(310)][i(379)](
                    this._virtualCamera.direction,
                    this._virtualCamera[i(1657)]
                  ),
                  Cesium[i(310)][i(1861)](e[i(1334)], o, s),
                  Cesium[i(310)].subtract(this[i(1449)], s, n),
                  (n = Gi(n, this._normal)),
                  Cesium.Cartesian3[i(875)](n, n),
                  Cesium[i(310)][i(1861)](n, this._centerWorldPosition, n),
                  (this[i(1557)].up = Cesium.Cartesian3[i(1676)](
                    n,
                    this._virtualCamera[i(2251)],
                    new Cesium.Cartesian3()
                  )),
                  Cesium[i(310)][i(379)](this[i(1557)].up, this._virtualCamera.up),
                  (this._reflectorProjectionMatrix = this._virtualCamera.frustum.projectionMatrix),
                  (this[i(639)] = this[i(1557)].viewMatrix)
                const a = Cesium[i(2314)][i(1163)](this[i(1449)], this._normal)
                Cesium[i(2314)].transform(a, this._virtualCamera[i(290)], a)
                const h = new Cesium.Cartesian4(a[i(226)].x, a[i(226)].y, a[i(226)].z, a.distance),
                  l = Cesium.Matrix4[i(1902)](this._virtualCamera[i(1850)][i(703)]),
                  c = new Cesium[i(2141)](
                    (Math.sign(h.x) + l[8]) / l[0],
                    (Math[i(1708)](h.y) + l[9]) / l[5],
                    -1,
                    (1 + l[10]) / l[14]
                  )
                return (
                  Cesium[i(2141)][i(606)](h, 2 / Cesium.Cartesian4.dot(h, c), h),
                  (l[2] = h.x),
                  (l[6] = h.y),
                  (l[10] = h.z + 1 - 0),
                  (l[14] = h.w),
                  (this[i(1557)][i(1850)][i(1297)] = Cesium[i(2066)][i(1902)](l)),
                  !0
                )
              }
              [t(1218)](e) {
                const i = t
                if (!this[i(1144)]) return
                let s = e[i(1452)][i(1391)],
                  n = e[i(1694)],
                  o = e[i(2557)].show,
                  r = e[i(2557)][i(255)]
                if (!this[i(2383)](e._defaultView[i(1391)]))
                  return void (this[i(1934)][i(1482)] = !1)
                ;(this._primitive[i(1482)] = !1),
                  (e[i(1452)][i(1391)] = this[i(1557)]),
                  (e[i(1694)] = void 0),
                  (e[i(2557)][i(1482)] = this[i(1506)]),
                  (e[i(2557)][i(255)] = this[i(1506)])
                const a = e[i(2325)],
                  h = a.drawingBufferWidth,
                  l = a[i(237)],
                  c = e.highDynamicRange
                this._createFramebuffer(a, h, l, c),
                  (function (t, e) {
                    const n = i
                    let o = t[n(1138)],
                      r = t.context,
                      a = r[n(656)],
                      h = t[n(1452)]
                    ;(t._view = h),
                      t[n(2021)](),
                      (o[n(1809)][n(1914)] = !0),
                      (o[n(1809)][n(1002)] = t[n(2165)][n(389)]),
                      (o[n(130)] = Ui)
                    let l = Cesium[n(1960)](t.backgroundColor, Cesium[n(1154)].BLACK)
                    t._hdr &&
                      (((l = Cesium[n(1154)][n(1902)](l, ji)).red = Math.pow(l[n(1626)], t.gamma)),
                      (l[n(1361)] = Math[n(2394)](l[n(1361)], t.gamma)),
                      (l[n(2034)] = Math[n(2394)](l[n(2034)], t[n(1488)]))),
                      (o[n(821)] = l),
                      t[n(1927)][n(720)](o),
                      a[n(720)](o)
                    const c = t[n(1694)]
                    Cesium[n(2330)](c) &&
                      c[s(147)] &&
                      (!Cesium[n(2330)](t[n(399)]) || t[n(399)] instanceof Cesium[s(137)]
                        ? Cesium[n(310)].negate(a[n(1160)], t[n(1560)][n(1657)])
                        : Cesium.Cartesian3[n(1902)](t[n(399)][n(1657)], t[n(1560)][n(1657)]),
                      o[n(1747)][n(2553)](c)),
                      (t[n(719)][n(277)] = 0),
                      (t._overlayCommandList.length = 0)
                    const u = h[n(2067)]
                    ;(u.x = 0),
                      (u.y = 0),
                      (u[n(575)] = r.drawingBufferWidth),
                      (u[n(2306)] = r[n(237)])
                    const m = h[n(2212)]
                    ;(m[n(915)] = e),
                      (m[n(2583)] = void 0),
                      (m[n(1703)] = void 0),
                      (m.viewport = Cesium[n(2111)].clone(u, m[n(2067)])),
                      Cesium.defined(t[n(2557)]) && t[n(2557)][n(242)](o),
                      t[n(1273)](),
                      t[n(1645)](m, l),
                      t.resolveFramebuffers(m),
                      Cesium.defined(t[n(2557)]) &&
                        (t.globe[n(137)](o), !t.globe[n(2508)] && (t[n(445)] = !0)),
                      r[n(137)]()
                  })(e, this[i(2397)])
                const u = this._primitive[i(1994)],
                  m = Cesium[i(1641)][i(1814)]({ context: a, framebuffer: this[i(2397)] })
                ;(m[i(1640)] = i(788)),
                  (this[i(1296)][i(1285)].colorTexture = m),
                  (this[i(1296)][i(1285)][i(346)] = Cesium[i(2066)].toArray(
                    this._getFixedFrameToEastNorthUpTransformFromWorldMatrix()
                  )),
                  (u[i(1285)][i(361)] = Cesium[i(2066)][i(1789)](this._reflectMatrix)),
                  (u[i(1285)][i(1044)] = Cesium.Matrix4[i(1789)](this._reflectorProjectionMatrix)),
                  (u[i(1285)][i(448)] = Cesium[i(2066)].toArray(this[i(639)])),
                  (this._primitive.show = !0),
                  (e[i(1452)][i(1391)] = s),
                  (e.shadowMap = n),
                  (e[i(2557)][i(1482)] = o),
                  (e[i(2557)][i(255)] = r)
              }
              [t(1869)]() {
                const e = t,
                  i = Cesium[e(2058)][e(1648)](this._centerWorldPosition)
                return Cesium[e(2066)][e(814)](i, new Cesium.Matrix4())
              }
              [t(856)](e) {
                const i = t
                this[i(2018)][i(1218)][i(1973)](this.preRender, this), (this[i(2018)][i(2581)] = !1)
              }
              [t(1494)](e, i, s, n) {
                const o = t
                let r = this[o(840)]
                if (Cesium[o(2330)](r) && r[o(575)] === i && r[o(2306)] === s && this[o(867)] === n)
                  return
                this[o(1515)](), (this[o(867)] = n)
                const a = n
                  ? e[o(120)]
                    ? Cesium[o(547)][o(1906)]
                    : Cesium.PixelDatatype[o(2441)]
                  : Cesium.PixelDatatype[o(1382)]
                ;(this[o(840)] = new Cesium[o(1641)]({
                  context: e,
                  width: i,
                  height: s,
                  pixelFormat: Cesium.PixelFormat[o(810)],
                  pixelDatatype: a,
                  sampler: new Cesium[o(2203)]({
                    wrapS: Cesium[o(917)][o(1558)],
                    wrapT: Cesium[o(917)][o(1558)],
                    minificationFilter: Cesium[o(485)].LINEAR,
                    magnificationFilter: Cesium[o(500)][o(2523)]
                  })
                })),
                  (this._depthStencilTexture = new Cesium[o(1641)]({
                    context: e,
                    width: i,
                    height: s,
                    pixelFormat: Cesium.PixelFormat.DEPTH_STENCIL,
                    pixelDatatype: Cesium[o(547)].UNSIGNED_INT_24_8
                  })),
                  (this._colorFramebuffer = new Cesium.Framebuffer({
                    context: e,
                    colorTextures: [this[o(840)]],
                    depthStencilTexture: this[o(257)],
                    destroyAttachments: !1
                  }))
              }
              _setStyle() {
                const e = t
                for (const t in this[e(1296)][e(1285)])
                  Object.hasOwnProperty.call(this[e(1296)].uniforms, t) &&
                    this._style[t] &&
                    (this[e(1296)][e(1285)][t] = this[e(343)][t])
                let i = this[e(1496)] ? this._style[e(506)] : this[e(343)][e(1070)]
                ;(this[e(1296)].uniforms[e(2096)] = Cesium[e(1154)][e(2008)](i)),
                  this[e(343)].height != this._height &&
                    ((this[e(2306)] = this[e(343)].height), this[e(1920)] && this[e(1920)]())
              }
              [t(1407)]() {
                const e = t
                this._material = this[e(1414)]()
                const i = {}
                ;(i[e(322)] = this._material), (i[e(2387)] = e(198)), (i[e(2381)] = !0)
                const s = new Cesium[e(985)](i)
                return (
                  (s.uniforms = {}),
                  (s[e(1285)][e(361)] = Cesium[e(2066)][e(1789)](this[e(997)])),
                  (s[e(1285)].reflectorProjectionMatrix = Cesium[e(2066)][e(1789)](this[e(2131)])),
                  (s[e(1285)][e(448)] = Cesium.Matrix4[e(1789)](this[e(639)])),
                  s
                )
              }
              [t(1414)]() {
                const e = t
                let i = this._scene[e(2325)],
                  s = Cesium.Texture[e(1814)]({ context: i, framebuffer: this[e(2397)] })
                s[e(1640)] = e(788)
                const n = Cesium[e(1643)](e(1673)),
                  o = {
                    lightDirection: new Cesium[e(310)](0, 0, 1),
                    rippleSize: this[e(343)][e(757)],
                    waterColor: this[e(343)][e(2096)],
                    waterAlpha: this._style.waterAlpha,
                    reflectivity: this._style[e(386)],
                    sunShiny: this[e(343)][e(254)],
                    distortionScale: this[e(343)].distortionScale,
                    waterColor: Cesium.Color[e(2008)](this[e(343)][e(1070)]),
                    normalTexture: n,
                    colorTexture: s,
                    fixedFrameToEastNorthUpTransform: Cesium.Matrix4.toArray(this[e(1869)]())
                  },
                  r = {}
                ;(r[e(1640)] = e(550)), (r[e(1285)] = o), (r[e(660)] = e(2195))
                let a = new Cesium[e(1637)]({
                  fabric: r,
                  translucent: !1,
                  minificationFilter: Cesium[e(485)][e(2523)],
                  magnificationFilter: Cesium[e(500)][e(2523)]
                })
                return (
                  Cesium[e(1283)]
                    [e(1030)](n)
                    [e(491)]()
                    [e(687)]((t) => {
                      const s = e
                      let n = new Cesium[s(1641)]({
                        context: i,
                        source: t,
                        sampler: new Cesium.Sampler({
                          wrapS: Cesium[s(917)][s(1563)],
                          wrapT: Cesium[s(917)][s(1563)],
                          minificationFilter: Cesium[s(485)][s(2523)],
                          magnificationFilter: Cesium[s(500)][s(2523)]
                        })
                      })
                      ;(n[s(1640)] = 'sampler2D'),
                        n[s(1130)](Cesium[s(1098)][s(1606)]),
                        (a.uniforms[s(1799)] = n)
                    }),
                  a
                )
              }
              [t(434)]() {
                const e = t
                this[e(2018)][e(1218)][e(704)](this[e(1218)], this), this[e(1515)]()
              }
              [t(1515)]() {
                const e = t
                this[e(840)] && this[e(840)][e(2039)](),
                  this[e(257)] && this[e(257)].destroy(),
                  this._colorFramebuffer && this[e(2397)][e(2039)](),
                  (this[e(840)] = void 0),
                  (this[e(257)] = void 0),
                  (this[e(2397)] = void 0)
              }
            })(e)
          case Q[i(653)]:
            return new (class extends Hi {
              constructor(e) {
                const i = t
                super(e),
                  (this[i(2336)] = Q.refractionWater),
                  (this[i(528)] = '折射水面'),
                  (this._style.color = Cesium[i(1960)](this[i(343)][i(1070)], 'rgb(112,144,201)')),
                  (this[i(343)][i(876)] = Cesium[i(1960)](this[i(343)].brightness, 2)),
                  (this[i(343)][i(2535)] = Cesium[i(1960)](this[i(343)][i(2535)], 3)),
                  (this[i(343)][i(730)] = Cesium[i(1960)](this._style[i(730)], 0.008))
              }
              [t(1218)](e) {
                const i = t
                if (!this[i(1144)]) return
                let s = e[i(1452)][i(1391)],
                  n = e[i(1694)],
                  o = e.globe[i(1482)],
                  r = e[i(2557)][i(255)]
                ;(this[i(1934)][i(1482)] = !1),
                  (e.shadowMap = void 0),
                  (e.globe[i(1482)] = !1),
                  (e[i(2557)].showSkirts = !1)
                const a = e[i(2325)],
                  h = a[i(2030)],
                  l = a[i(237)],
                  c = e[i(337)]
                this._createFramebuffer(a, h, l, c),
                  (function (t, e) {
                    const n = i
                    let o = t[n(1138)],
                      r = t[n(2325)],
                      a = r.uniformState,
                      h = t[n(1452)]
                    ;(t[n(1445)] = h),
                      t[n(2021)](),
                      (o[n(1809)][n(1914)] = !0),
                      (o.passes[n(1002)] = t.postProcessStages[n(389)]),
                      (o.tilesetPassState = Yi)
                    let l = Cesium[n(1960)](t[n(821)], Cesium[n(1154)][n(539)])
                    t[n(867)] &&
                      (((l = Cesium.Color.clone(l, Xi))[n(1626)] = Math[n(2394)](
                        l[n(1626)],
                        t[n(1488)]
                      )),
                      (l[n(1361)] = Math.pow(l.green, t[n(1488)])),
                      (l[n(2034)] = Math[n(2394)](l[n(2034)], t[n(1488)]))),
                      (o[n(821)] = l),
                      t[n(1927)][n(720)](o),
                      a[n(720)](o)
                    const c = t.shadowMap
                    Cesium[n(2330)](c) &&
                      c[s(147)] &&
                      (!Cesium[n(2330)](t[n(399)]) || t[n(399)] instanceof Cesium[s(137)]
                        ? Cesium[n(310)][n(875)](a[n(1160)], t[n(1560)][n(1657)])
                        : Cesium[n(310)][n(1902)](t.light[n(1657)], t[n(1560)][n(1657)]),
                      o[n(1747)][n(2553)](c)),
                      (t._computeCommandList[n(277)] = 0),
                      (t[n(1492)][n(277)] = 0)
                    const u = h[n(2067)]
                    ;(u.x = 0),
                      (u.y = 0),
                      (u[n(575)] = r[n(2030)]),
                      (u[n(2306)] = r.drawingBufferHeight)
                    const m = h.passState
                    ;(m[n(915)] = e),
                      (m.blendingEnabled = void 0),
                      (m.scissorTest = void 0),
                      (m[n(2067)] = Cesium[n(2111)][n(1902)](u, m[n(2067)])),
                      Cesium.defined(t[n(2557)]) && t.globe.beginFrame(o),
                      t[n(1273)](),
                      t[n(1645)](m, l),
                      t.resolveFramebuffers(m),
                      Cesium[n(2330)](t[n(2557)]) &&
                        (t[n(2557)].endFrame(o), !t[n(2557)][n(2508)] && (t[n(445)] = !0)),
                      r.endFrame()
                  })(e, this[i(2397)]),
                  this._primitive[i(1994)]
                const u = {}
                ;(u[i(2325)] = a), (u[i(915)] = this._colorFramebuffer)
                const m = Cesium[i(1641)].fromFramebuffer(u)
                ;(m.type = i(788)),
                  (this[i(1296)][i(1285)][i(1052)] = m),
                  (this._primitive.show = !0),
                  (e[i(1452)][i(1391)] = s),
                  (e[i(1694)] = n),
                  (e[i(2557)][i(1482)] = o),
                  (e[i(2557)][i(255)] = r)
              }
              [t(856)](e) {
                const i = t,
                  s = this[i(2018)][i(2325)]
                this[i(1494)](s, s[i(2030)], s[i(237)], this[i(2018)][i(337)]),
                  this[i(2018)][i(1218)][i(1973)](this[i(1218)], this)
              }
              [t(2231)]() {
                const e = t
                ;(this[e(1296)].uniforms[e(876)] = this[e(343)][e(876)]),
                  (this._material[e(1285)][e(2535)] = this._style[e(2535)]),
                  (this[e(1296)][e(1285)].waveStrength = this[e(343)][e(730)])
                let i = this[e(1496)] ? this[e(343)][e(506)] : this._style[e(1070)]
                ;(this[e(1296)][e(1285)].waterColor = Cesium[e(1154)].fromCssColorString(i)),
                  this[e(343)][e(2306)] != this[e(1516)] &&
                    ((this.height = this[e(343)][e(2306)]), this[e(1920)] && this[e(1920)]())
              }
              _createAppearence() {
                const e = t
                this[e(1296)] = this[e(1414)]()
                const i = {}
                ;(i[e(322)] = this[e(1296)]),
                  (i[e(2387)] =
                    '\n               in vec3 position3DHigh;\n               in vec3 position3DLow;\n               in vec3 normal;\n               in vec2 st;\n               in float batchId;\n               \n               out vec3 v_positionEC;\n               out vec3 v_normalEC;\n               out vec2 v_st;\n               \n               uniform mat4 reflectorProjectionMatrix;\n               uniform mat4 reflectorViewMatrix;\n               uniform mat4 reflectMatrix;\n               out vec4 v_worldPosition;  // 世界坐标\n               out vec4 v_uv;             // 纹理坐标\n               \n               \n               \n               void main()\n               {\n                   vec4 p = czm_computePosition();\n               \n                   v_positionEC = (czm_modelViewRelativeToEye * p).xyz;      // position in eye coordinates\n                   v_normalEC = czm_normal * normal;                         // normal in eye coordinates\n                   v_st = st;\n               \n               \n                   mat4 modelView = reflectorViewMatrix * reflectMatrix * czm_model;\n                   modelView[3][0] = 0.0;\n                   modelView[3][1] = 0.0;\n                   modelView[3][2] = 0.0;\n                   v_uv = reflectorProjectionMatrix * modelView * p;\n                   vec4 positionMC = vec4( position3DHigh + position3DLow, 1.0 );\n                   v_worldPosition = czm_model * positionMC;\n               \n                   gl_Position = czm_modelViewProjectionRelativeToEye * p;\n          }'),
                  (i[e(2381)] = !0)
                const s = new Cesium[e(985)](i)
                return (
                  (s[e(1285)] = {}),
                  (s.uniforms[e(361)] = Cesium.Matrix4[e(1789)](this[e(997)])),
                  (s[e(1285)].reflectorProjectionMatrix = Cesium.Matrix4[e(1789)](
                    this._reflectorProjectionMatrix
                  )),
                  (s[e(1285)].reflectorViewMatrix = Cesium[e(2066)][e(1789)](this[e(639)])),
                  s
                )
              }
              _createMaterial() {
                const e = t
                let i = this[e(2018)][e(2325)]
                new Cesium[e(1641)]({
                  context: this[e(2018)][e(2325)],
                  source: {
                    width: 1,
                    height: 1,
                    arrayBufferView: new Uint8Array([255, 0, 0, 255])
                  },
                  sampler: new Cesium[e(2203)]({
                    wrapS: Cesium[e(917)].REPEAT,
                    wrapT: Cesium.TextureWrap.REPEAT,
                    minificationFilter: Cesium[e(485)][e(2523)],
                    magnificationFilter: Cesium.TextureMinificationFilter[e(2523)]
                  })
                })[e(1640)] = e(788)
                const s = {}
                ;(s[e(2325)] = i), (s[e(915)] = this._colorFramebuffer)
                let n = Cesium[e(1641)].fromFramebuffer(s)
                n[e(1640)] = e(788)
                const o = this._style[e(2230)] || Cesium[e(1643)](e(1673)),
                  r = {
                    brightness: this._style.brightness,
                    speed: this[e(343)][e(2535)],
                    waveStrength: this[e(343)][e(730)],
                    waterColor: new Cesium[e(1154)][e(2008)](this[e(343)].color),
                    distortTexture: o,
                    colorTexture: n
                  },
                  a = {}
                ;(a.type = e(1348)),
                  (a[e(1285)] = r),
                  (a[e(660)] =
                    '\n                uniform vec4 waterColor;\n                uniform sampler2D colorTexture;\n                uniform sampler2D distortTexture;\n                uniform float brightness;\n                uniform float speed;\n                uniform float waveStrength;\n                \n                czm_material czm_getMaterial(czm_materialInput materialInput)\n                {\n                    czm_material material = czm_getDefaultMaterial(materialInput);\n        \n                    float time1 = czm_frameNumber/60. * speed;\n                    vec2 uv = gl_FragCoord.xy / czm_viewport.zw;\n                    vec2 vUv = fract(materialInput.st);\n        \n                    // float waveStrength = 0.1;\n                    // float waveSpeed = 0.03;\n        \n                    float waveSpeed = 0.03;\n        \n                    // simple distortion (ripple) via dudv map (see https://www.youtube.com/watch?v=6B7IF6GOu7s)\n        \n                    vec2 distortedUv = texture( distortTexture, vec2( vUv.x + time1 * waveSpeed, vUv.y ) ).rg * waveStrength;\n                    distortedUv = vUv.xy + vec2( distortedUv.x, distortedUv.y + time1 * waveSpeed );\n                    vec2 distortion = ( texture( distortTexture, distortedUv ).rg * 2.0 - 1.0 ) * waveStrength;\n                \n                    uv.xy += distortion;\n                    \n                    material.diffuse = texture(colorTexture, uv).rgb; \n                    material.diffuse = brightness * waterColor.rgb * material.diffuse;\n                    material.specular = 0.99;\n                    material.shininess = 20.0;\n                    material.alpha = waterColor.a;\n                    material.normal = materialInput.normalEC;\n        \n                    if ( materialInput.normalEC.z < 0.0 ) {\n                        material.diffuse.rgb = vec3(0, 0, 0); // back\n                        material.shininess = 0.0;\n                        material.specular = 0.0;\n                        material.alpha = waterColor.a;\n                    }  \n        \n                    return material;\n                }')
                let h = new Cesium[e(1637)]({
                  fabric: a,
                  translucent: !1,
                  minificationFilter: Cesium[e(485)][e(2523)],
                  magnificationFilter: Cesium[e(500)][e(2523)]
                })
                return (
                  Cesium[e(1283)]
                    [e(1030)](o)
                    [e(491)]()
                    [e(687)]((t) => {
                      const s = e
                      let n = new Cesium[s(1641)]({
                        context: i,
                        source: t,
                        sampler: new Cesium[s(2203)]({
                          wrapS: Cesium[s(917)][s(1563)],
                          wrapT: Cesium.TextureWrap[s(1563)],
                          minificationFilter: Cesium.TextureMinificationFilter.LINEAR,
                          magnificationFilter: Cesium[s(500)][s(2523)]
                        })
                      })
                      ;(n[s(1640)] = s(788)),
                        n.generateMipmap(Cesium.MipmapHint[s(1606)]),
                        (h[s(1285)][s(2230)] = n)
                    }),
                  h
                )
              }
              [t(434)]() {
                const e = t
                this[e(2018)][e(1218)][e(704)](this[e(1218)], this), this[e(1515)]()
              }
            })(e)
        }
      }
    }
    let Zi = {
        readFeature(e) {
          const i = t
          let s = e.geometry[i(2365)],
            n = e.properties
          return (n[i(2251)] = s), this[i(1951)](n)
        },
        create(e) {
          const i = t
          if (e[i(1158)] === Q[i(2040)])
            return new (class extends b {
              constructor(e = {}) {
                const i = t
                super(e),
                  (this[i(1912)] = e.position),
                  (this[i(528)] = i(1522)),
                  (this[i(314)] = i(1865)),
                  (this[i(1115)] = X),
                  (this._graphicType = Q[i(2040)]),
                  (this[i(350)] = 2),
                  (this[i(343)][i(1981)] = Cesium[i(1960)](e.style.radius, 10)),
                  (this[i(343)][i(1308)] = Cesium[i(1960)](e[i(1679)].visibleColor, i(552))),
                  (this._style[i(661)] = Cesium[i(1960)](e.style[i(506)], i(197))),
                  (this[i(343)][i(2345)] = Cesium.defaultValue(e.style[i(2345)], 0.8))
                let s = Cesium.defaultValue(e.position, [111, 28, 0])
                ;(this._cartesian3 = Cesium[i(310)][i(667)](s[0], s[1], s[2])),
                  (this[i(794)] = this[i(1912)] = s),
                  (this._boundingSphere = new Cesium.BoundingSphere(
                    this[i(2238)],
                    this[i(343)][i(1981)]
                  )),
                  (this._radius = this[i(343)][i(1981)])
              }
              get [t(2416)]() {
                return this[t(1934)]
              }
              [t(415)](e) {
                const i = t
                this[i(1934)] && (this[i(1934)][i(1482)] = e)
              }
              [t(2231)]() {
                const e = t
                ;(this[e(1024)] = new Cesium[e(1242)](this._cartesian3, this[e(343)][e(1981)])),
                  this._primitive[e(1994)] &&
                    ((this[e(1934)][e(1994)][e(322)][e(1285)][e(1743)] = Cesium[e(1154)][e(2008)](
                      this[e(343)][e(1308)]
                    )),
                    (this[e(1934)].appearance[e(322)][e(1285)][e(2308)] = Cesium[e(1154)][e(2008)](
                      this[e(343)][e(661)]
                    )),
                    (this._primitive[e(1994)][e(322)].uniforms[e(880)] = this[e(343)][e(1981)]),
                    this[e(920)] != this._style.radius &&
                      (this._addPrimitive(), (this[e(920)] = this[e(343)][e(1981)])))
              }
              [t(1253)](e) {
                const i = t
                ;(this[i(1912)] = e),
                  (this._cartesian3 = Cesium[i(310)][i(667)](e[0], e[1], e[2])),
                  (this[i(794)] = this[i(1912)]),
                  (this[i(1024)] = new Cesium[i(1242)](this[i(2238)], this._style[i(1981)])),
                  this[i(1365)] &&
                    ((this[i(1365)][i(2251)] = this[i(2238)]),
                    (this[i(1365)].up = Cesium.Cartesian3[i(379)](
                      this[i(2238)],
                      new Cesium[i(310)]()
                    )),
                    this[i(443)]())
              }
              [t(2515)](e) {
                const i = t
                ;(this[i(2469)] = e), (this[i(2469)][1][2] = this[i(2469)][0][2])
                let s = Cesium[i(310)][i(774)]([][i(1500)].apply([], e)),
                  n = Cesium[i(310)].distance(s[0], s[1])
                ;(this._style[i(1981)] = Number(n.toFixed(2)) || 1),
                  this[i(1253)](this._positions[0]),
                  this[i(1934)].appearance &&
                    ((this[i(1934)][i(1994)][i(322)][i(1285)][i(880)] = this[i(343)][i(1981)]),
                    (this[i(2444)][i(426)] = this[i(343)][i(1981)]))
              }
              [t(438)]() {
                const e = t
                return new Cesium[e(2186)]({
                  geometryInstances: new Cesium[e(2155)]({
                    geometry: new Cesium[e(1718)]({
                      vertexFormat: Cesium[e(2567)][e(1145)],
                      radius: this.style.radius
                    }),
                    modelMatrix: Cesium.Transforms[e(1648)](this._cartesian3)
                  }),
                  asynchronous: !1
                })
              }
              _createShadowMap(e) {
                const i = t
                ;(this[i(1365)] = new Cesium[i(2162)](e)),
                  (this[i(1365)][i(2251)] = this[i(2238)]),
                  (this[i(1365)].up = Cesium[i(310)][i(379)](
                    this[i(2238)],
                    new Cesium.Cartesian3()
                  )),
                  (this[i(2444)] = new Cesium[i(433)]({
                    lightCamera: this[i(1365)],
                    enable: !1,
                    darkness: 1,
                    isPointLight: !0,
                    isSpotLight: !1,
                    cascadesEnabled: !1,
                    context: e[i(2325)],
                    pointLightRadius: this[i(1679)][i(1981)],
                    fromLightSource: !1
                  }))
              }
              [t(1414)](e) {
                const i = t,
                  s = new Cesium[i(1637)]({
                    fabric: {
                      type: i(2245),
                      uniforms: {
                        u_visibleColor: Cesium[i(1154)][i(2008)](i(1263)),
                        u_hiddenColor: Cesium.Color[i(2008)](i(669)),
                        u_radius: this[i(1679)].radius
                      },
                      source: i(1300)
                    }
                  })
                ;(s[i(1685)][i(965)] = () =>
                  this[i(2444)][i(1247)] ? this[i(2444)][i(1247)] : e[i(2325)][i(1627)]),
                  (s[i(1685)][i(1097)] = () =>
                    Cesium[i(2141)][i(2434)](
                      1 / this._shadowMap[i(1188)].x,
                      1 / this[i(2444)][i(1188)].y,
                      this._shadowMap[i(2537)][i(1658)],
                      this._shadowMap[i(2537)][i(1379)]
                    )),
                  (s[i(1685)][i(1065)] = () => this[i(2444)][i(585)]),
                  (s._uniforms.u_lightPositionEC = () => this[i(2444)]._lightPositionEC),
                  (s._uniforms[i(2345)] = () => this[i(1679)].globalAlpha),
                  (this[i(1075)] = new Cesium[i(985)]({
                    flat: !0,
                    material: s,
                    vertexShaderSource: i(2191),
                    fragmentShaderSource: i(1442)
                  }))
              }
              [t(720)](e) {
                const i = t
                e[i(1747)].push(this[i(2444)]), this[i(1934)].update(e)
              }
              [t(2039)]() {}
              [t(443)]() {
                const e = t
                this[e(1934)] && this[e(2462)][e(245)][e(696)][e(1346)].remove(this[e(1934)]),
                  (this[e(1934)] = this._createPrimitive()),
                  (this[e(1934)][e(2122)] = this[e(2462)].id),
                  (this[e(1934)][e(2366)] = this[e(1570)]),
                  (this._primitive[e(1994)] = this._appearance)
              }
              [t(1545)](e) {
                const i = t
                ;(this._layer = e),
                  this[i(234)](e[i(245)][i(696)]),
                  this._createMaterial(e._viewer.scene),
                  this._addPrimitive(),
                  e[i(245)][i(696)][i(1346)][i(1861)](this)
              }
              [t(2389)](e) {
                const i = t
                e[i(245)][i(696)][i(1346)][i(1896)](this._primitive),
                  e[i(245)][i(696)][i(1346)][i(1896)](this),
                  (this._primitive = null)
              }
            })(e)
        }
      },
      Ki = {
        readFeatures(e) {
          const i = t
          let s = []
          return (
            e[i(1602)]((t) => {
              const e = i
              let n = this[e(1549)](t)
              n && s[e(2553)](n)
            }),
            s
          )
        },
        readFeature(e) {
          const i = t
          switch (e.properties[i(1059)]) {
            case X:
              return Zi[i(1549)](e)
            case j:
              return hi[i(1549)](e)
            case Y:
              return Qi[i(1549)](e)
            case q:
              return Fi[i(1549)](e)
            case A:
              return ue[i(1549)](e)
            case D:
              return Ee[i(1549)](e)
            case P:
              return Ft.readFeature(e)
            case R:
              return be[i(1549)](e)
            case k:
              return Se[i(1549)](e)
            case I:
              return Ze[i(1549)](e)
            case T:
              return pe[i(1549)](e)
            case M:
              return ce.readFeature(e)
            case O:
              return Je[i(1549)](e)
            case F:
              return ti.readFeature(e)
            case E:
              return fe[i(1549)](e)
            case L:
              return ve[i(1549)](e)
            case z:
              return ge[i(1549)](e)
            case B:
              return si[i(1549)](e)
            case N:
              return oi[i(1549)](e)
            case G:
              return ai[i(1549)](e)
            case W:
              return li[i(1549)](e)
            case U:
              return gi[i(1549)](e)
          }
        }
      }
    class Ji {
      constructor(e) {
        const i = t
        ;(this[i(245)] = e),
          (this[i(261)] = this[i(245)][i(696)][i(1346)].add(new Cesium.PointPrimitiveCollection()))
      }
      updatePosition(e) {
        const i = t,
          s = {}
        ;(s[i(2251)] = e),
          (s[i(225)] = 8),
          (s[i(1070)] = Cesium.Color[i(2147)]),
          (s.outlineWidth = 1),
          (s[i(725)] = Cesium[i(1154)][i(588)]),
          this[i(261)][i(1876)](),
          this[i(261)][i(1861)](s)
      }
      [t(1896)]() {
        const e = t
        this[e(245)][e(696)][e(1346)][e(1896)](this[e(261)])
      }
    }
    class $i {
      constructor(e) {
        const i = t
        ;(this[i(245)] = e), this[i(2271)](), (this._content = [])
      }
      [t(1520)](e) {
        const i = t
        ;(this[i(992)] = e), this[i(2161)](e)
      }
      [t(302)](e) {
        const i = t
        let s = this[i(992)][i(1500)](e)
        this[i(2161)](s)
      }
      _initContent(e) {
        const i = t
        this[i(2526)]()
        let s,
          n = Array[i(1108)](e) && e[i(277)] > 0
        ;(this._container[i(1679)].display = i(n ? 969 : 1772)),
          n &&
            e.forEach((t) => {
              const e = i
              ;((s = document[e(1945)](e(1120))).innerHTML = t), this[e(360)].appendChild(s)
            })
      }
      setPosition(e) {
        const i = t
        ;(this[i(248)] = e), this[i(731)]()
      }
      [t(731)]() {
        const e = t
        if (!this[e(360)] || !this._container[e(1679)]) return
        if (!this[e(248)]) return
        const i = this[e(360)][e(2127)] / 2
        ;(this[e(360)][e(1679)].top = this[e(248)].y - i + 'px'),
          (this[e(360)][e(1679)][e(1709)] = this[e(248)].x + 30 + 'px')
      }
      _createDom() {
        const e = t
        ;(this[e(360)] = document[e(1945)]('div')),
          this[e(360)][e(2453)][e(1861)](e(1327)),
          this[e(245)][e(436)][e(1392)].appendChild(this[e(360)])
      }
      [t(2526)]() {
        const e = t
        for (; this[e(360)][e(167)](); ) this[e(360)].removeChild(this[e(360)][e(535)])
      }
      [t(1896)]() {
        const e = t
        this._viewer[e(436)].container[e(2269)](this._container)
      }
    }
    function ts(e) {
      const i = t,
        s = Cesium[i(2285)][i(2579)](e)
      return [Cesium[i(475)].toDegrees(s[i(2106)]), Cesium.Math.toDegrees(s[i(199)])]
    }
    function es(e) {
      const i = t
      return Cesium[i(2285)].fromCartesian(e)[i(2306)]
    }
    function is(e) {
      const i = t,
        s = Cesium[i(2285)].fromCartesian(e)
      return [Cesium.Math[i(363)](s[i(2106)]), Cesium.Math.toDegrees(s[i(199)]), s.height]
    }
    function ss(e) {
      const i = t
      if (e[i(277)] < 2) return 0
      let s = 0
      for (let t = 0; t < e[i(277)] - 1; t++) s += Cesium.Cartesian3[i(1849)](e[t], e[t + 1])
      return s
    }
    function ns(e) {
      const i = t
      if (e[i(277)] < 3) return 0
      let s = []
      e.forEach((t) => {
        s[i(2553)](ts(t))
      }),
        s.push(s[0])
      let n = turf[i(2087)]([s])
      return turf[i(2346)](n)
    }
    class os extends v {
      constructor(e = {}) {
        const i = t
        if ((super(e), (this._viewer = e[i(395)]), !this[i(245)]))
          throw new Cesium[i(2360)]('viewer是必须的！')
        const s = { x: 7, y: 5 },
          n = { x: 0, y: -40 },
          o = { x: 0, y: 0, z: 0 },
          r = {}
        ;(r[i(1292)] = 0),
          (r.far = 1e4),
          (this[i(2571)] = Cesium[i(1960)](e[i(2540)], 2)),
          (this._measureEntities = []),
          (this[i(2469)] = []),
          (this[i(192)] = []),
          (this[i(437)] = new Cesium[i(1561)](this[i(245)].scene[i(1493)])),
          (this._labelOpts = {
            font: i(869),
            style: 2,
            scale: 0.5,
            showBackground: !0,
            backgroundColor: Cesium[i(1154)][i(2008)]('rgba(170, 170, 169, 0.72)'),
            backgroundPadding: s,
            pixelOffset: n,
            eyeOffset: o,
            horizontalOrigin: 0,
            verticalOrigin: 0,
            heightReference: 0,
            fillColor: Cesium.Color[i(2008)]('#FFFFFF'),
            outlineColor: Cesium[i(1154)][i(2008)](i(1593)),
            outlineWidth: 4,
            disableDepthTestDistance: 1e4,
            distanceDisplayCondition: r
          }),
          (this[i(387)] = {
            pixelSize: 10,
            color: Cesium.Color[i(2008)]('rgba(28,25,125,0.99)'),
            outline: !0,
            outlineWidth: 2,
            outlineColor: Cesium.Color[i(588)][i(329)](0.5),
            disableDepthTestDistance: 2e7
          }),
          (this._lineOpts = {
            width: 10,
            material: new Cesium[i(1847)]({ color: Cesium[i(1154)].YELLOW, glowPower: 0.1 })
          }),
          (this[i(1962)] = { material: Cesium[i(1154)][i(2474)][i(329)](0.5) }),
          (this[i(565)] = { material: Cesium[i(1154)][i(2474)][i(329)](0.3) }),
          (this._polygonOpts = Cesium.combine(Cesium[i(1960)](e[i(1598)], {}), this._polygonOpts)),
          (this._ellipseOpts = Cesium[i(938)](Cesium[i(1960)](e[i(1899)], {}), this._ellipseOpts)),
          (this._labelOpts = Cesium[i(938)](Cesium.defaultValue(e[i(2041)], {}), this._labelOpts)),
          (this[i(387)] = Cesium[i(938)](Cesium[i(1960)](e.pointOpts, {}), this[i(387)])),
          (this._lineOpts = Cesium[i(938)](Cesium.defaultValue(e[i(655)], {}), this[i(282)]))
      }
      start() {
        const e = t
        this[e(519)](),
          this._registerEvents(),
          (this[e(245)][e(471)] = !1),
          (this[e(245)]._element[e(1679)][e(402)] = e(1653)),
          (this[e(526)] = !0),
          (this[e(1795)] = 0),
          (this[e(960)] = new Ji(this._viewer)),
          (this[e(272)] = new $i(this._viewer)),
          this[e(188)](),
          this.fire(p.measureStart)
      }
      [t(188)]() {}
      [t(519)]() {
        const e = t
        this[e(526)] &&
          (this._unRegisterEvents(),
          (this[e(245)][e(374)][e(1679)][e(402)] = 'pointer'),
          (this[e(245)][e(471)] = !0),
          (this[e(526)] = !1),
          (this[e(2469)] = []),
          this[e(960)] &&
            (this[e(960)][e(1896)](),
            (this[e(960)] = null),
            this._mouseTip[e(1896)](),
            (this[e(272)] = null),
            this[e(240)]()))
      }
      [t(240)]() {}
      [t(2221)]() {
        const e = t
        this.stop(), this[e(450)](p[e(2070)], this._measureResult)
      }
      clear() {
        const e = t
        this[e(2243)][e(1602)]((t) => {
          const i = e
          this[i(245)].entities[i(1896)](t)
        }),
          (this[e(2243)] = [])
      }
      [t(2402)]() {
        const e = t
        this[e(1698)](), this._rightClickEvent(), this[e(1008)]()
      }
      [t(344)]() {
        const e = t
        this[e(437)][e(1638)](Cesium.ScreenSpaceEventType[e(1309)]),
          this._handler[e(1638)](Cesium[e(1827)][e(1282)]),
          this[e(437)][e(1638)](Cesium[e(1827)][e(1073)])
      }
      [t(2039)]() {
        const e = t
        this[e(1555)](), this[e(437)][e(2039)]()
      }
    }
    const rs = {}
    rs[t(922)] = 'Module'
    const as = Object[t(1731)](
      Object[t(1399)](
        {
          __proto__: null,
          CameraPicker: class extends v {
            constructor(e = {}) {
              const i = t
              super(e),
                (this[i(2137)] = new Date()[i(1566)]()),
                (this[i(245)] = e[i(395)]),
                (this[i(962)] = Cesium.defaultValue(e[i(675)], ',')),
                (this[i(2571)] = Cesium[i(1960)](e[i(2540)], 4)),
                this._createDom(),
                this[i(245)][i(436)][i(1392)][i(1621)](this[i(360)]),
                this._setPosition(),
                this[i(2077)](),
                this[i(831)]()
            }
            [t(2271)]() {
              const e = t
              ;(this[e(360)] = document[e(1945)](e(1120))), this[e(360)].classList.add(e(2083))
              let i = document[e(1945)](e(1120))
              i[e(2453)][e(1861)](e(1903)), i[e(2453)][e(1861)]('info-item-close')
              let s = document.createElement(e(1081))
              s.classList.add(e(945)),
                (s[e(176)] = 'x'),
                (s.title = '移除'),
                i.appendChild(s),
                this._container[e(1621)](i),
                (i = document.createElement(e(1120)))[e(2453)].add(e(1903)),
                (i[e(1679)][e(1811)] = e(2393))
              let n = document.createElement('span')
              ;(n[e(176)] = e(2437)), i[e(1621)](n)
              let o = document[e(1945)](e(1711))
              o[e(2453)][e(1861)](e(1714)),
                o[e(2453)][e(1861)](e(2255)),
                (o[e(922)] = this[e(2571)]),
                (o[e(1124)] = (t) => {
                  const i = e
                  ;(this._precision = o.value), this[i(2077)]()
                }),
                i[e(1621)](o),
                ((n = document[e(1945)]('span'))[e(176)] = '分隔符：'),
                (n[e(1679)][e(1811)] = e(889)),
                i[e(1621)](n)
              let r = document[e(1945)](e(1711))
              r[e(2453)][e(1861)](e(2255)),
                (r[e(922)] = this[e(962)]),
                i.appendChild(r),
                (r[e(1124)] = (t) => {
                  const i = e
                  ;(this[i(962)] = r[i(922)]), this[i(2077)]()
                }),
                this[e(360)][e(1621)](i),
                (i = document.createElement('div'))[e(2453)][e(1861)]('info-item'),
                (i = document[e(1945)](e(1120))).classList.add(e(1903)),
                (i[e(1679)][e(1811)] = e(648)),
                ((n = document[e(1945)](e(1081)))[e(176)] = e(483)),
                i.appendChild(n)
              let a = document.createElement(e(459))
              a[e(2453)].add(e(1714)),
                (a[e(1679)][e(1811)] = 'height: 58px; width: 184px;'),
                (this[e(1502)] = a),
                i[e(1621)](a),
                this[e(360)][e(1621)](i),
                (i = document[e(1945)](e(1120)))[e(2453)][e(1861)](e(1903)),
                (i[e(1679)][e(1811)] = e(648)),
                ((n = document[e(1945)]('span'))[e(176)] = e(908)),
                i.appendChild(n),
                (a = document.createElement(e(459)))[e(2453)][e(1861)](e(1714)),
                (a[e(1679)][e(1811)] = 'height: 58px; width: 140px;'),
                (this[e(488)] = a),
                i[e(1621)](a),
                this[e(360)][e(1621)](i)
              let h = document[e(1945)]('button')
              ;(h[e(176)] = '复制'),
                (h.title = e(169)),
                (h.onclick = (t) => {
                  const i = e
                  let s = this[i(1401)](),
                    n =
                      '\n                x: ' +
                      s[i(2543)].x +
                      i(984) +
                      s[i(2543)].y +
                      i(294) +
                      s[i(2543)].z +
                      ' , \n                heading: ' +
                      s[i(925)][i(1198)] +
                      ',\n                pitch: ' +
                      s.orientation.pitch +
                      i(1202) +
                      s[i(925)][i(1143)] +
                      ' '
                  navigator[i(1110)][i(728)](n)
                }),
                i.appendChild(h),
                this[e(360)][e(1621)](i)
            }
            _getCameraInfo() {
              const e = t
              let i = this[e(245)].scene
              return {
                destination: {
                  x: i[e(1391)][e(2251)].x[e(1268)](this._precision),
                  y: i[e(1391)][e(2251)].y[e(1268)](this[e(2571)]),
                  z: i.camera.position.z[e(1268)](this[e(2571)])
                },
                orientation: {
                  heading: i[e(1391)].heading[e(1268)](this._precision),
                  pitch: i[e(1391)][e(711)].toFixed(this[e(2571)]),
                  roll: i[e(1391)][e(1143)].toFixed(this._precision)
                }
              }
            }
            [t(1486)]() {
              const e = t
              if (!this[e(360)] || !this[e(360)][e(1679)]) return
              this[e(360)][e(1679)].position = e(2076)
              const i = this[e(360)][e(168)],
                s = this._container[e(2127)]
              ;(this[e(360)][e(1679)].bottom = e(1752) + s / 2 + 'px)'),
                (this._container[e(1679)][e(1709)] = 'calc(50% - ' + i / 2 + 'px)')
            }
            _setInputContent() {
              const e = t
              let i = this[e(1401)]()
              ;(this._inputP[e(922)] =
                i.destination.x +
                this[e(962)] +
                '\n' +
                i[e(2543)].y +
                this[e(962)] +
                '\n' +
                i[e(2543)].z),
                (this[e(488)][e(922)] =
                  i[e(925)][e(1198)] +
                  this[e(962)] +
                  '\n' +
                  i[e(925)][e(711)] +
                  this[e(962)] +
                  '\n' +
                  i[e(925)][e(1143)])
            }
            [t(831)]() {
              const e = t
              this._viewer[e(1391)][e(2552)][e(1973)](this[e(353)], this),
                (this[e(360)][e(846)]('info-close')[0][e(1803)] = (t) => {
                  this.remove()
                })
            }
            _cameraEvent(e) {
              this[t(2077)]()
            }
            [t(1896)]() {
              const e = t
              this[e(360)] &&
                (this[e(245)][e(1391)].moveEnd.removeEventListener(this[e(353)], this),
                this[e(245)][e(436)][e(1392)][e(2269)](this._container),
                (this._container = null))
            }
            [t(1712)]() {
              this.remove()
            }
          },
          DrawHandler: class extends v {
            constructor(e) {
              const i = t
              super(e),
                (this[i(245)] = e[i(395)]),
                (this[i(2469)] = []),
                (this[i(2290)] = []),
                (this[i(2055)] = !1),
                (this[i(1378)] = []),
                (this[i(264)] = []),
                (this[i(2490)] = new Cesium[i(1561)](this[i(245)][i(696)].canvas)),
                (this._polygonOpts = Cesium[i(1960)](e.polygonOpts, {})),
                (this[i(1995)] = Cesium[i(1960)](e.polylineOpts, {})),
                (this[i(960)] = new Ji(this._viewer)),
                (this[i(272)] = new $i(this[i(245)])),
                this._init()
            }
            _init() {
              const e = t
              this.handler[e(1068)]((t) => {
                const i = e
                if (!this[i(2055)]) return
                this[i(245)][i(374)][i(1679)].cursor = i(1653)
                let s = this[i(245)].scene.pickPosition(t[i(2251)])
                if (
                  (s ||
                    (s = this[i(245)][i(696)][i(1391)][i(1816)](
                      t[i(2251)],
                      this._viewer[i(696)][i(2557)].ellipsoid
                    )),
                  !s)
                )
                  return
                this[i(1378)][i(2553)](s), this._cartesian3ArrayTemp[i(2553)](s)
                let n = s,
                  o = Cesium[i(2285)].fromCartesian(n)
                ;(n = [
                  Cesium.Math[i(363)](o[i(2106)]),
                  Cesium.Math[i(363)](o[i(199)]),
                  o[i(2306)]
                ]),
                  this[i(2469)][i(2553)](n),
                  'point' == this._drawOpts[i(1640)] &&
                    (this[i(450)](p[i(182)], {
                      type: this[i(1485)][i(1640)],
                      position: this._positions[0]
                    }),
                    this.stop()),
                  this._positionsTemp[i(2553)](n),
                  this[i(2060)]()
              }, Cesium.ScreenSpaceEventType[e(1282)]),
                this[e(2490)].setInputAction((t) => {
                  const i = e
                  if (!this[i(2055)]) return
                  this[i(245)][i(374)][i(1679)].cursor = i(1653)
                  let s = this[i(245)][i(696)].pickPosition(t[i(954)])
                  if (
                    (s ||
                      (s = this[i(245)].scene[i(1391)][i(1816)](
                        t[i(813)],
                        this[i(245)][i(696)][i(2557)].ellipsoid
                      )),
                    !s)
                  )
                    return
                  this[i(960)].updatePosition(s),
                    this[i(272)][i(1253)](t[i(954)]),
                    (this[i(264)] = this._cartesian3Array[i(1500)]([s]))
                  let n = s,
                    o = Cesium.Cartographic[i(2579)](n)
                  ;(n = [
                    Cesium[i(475)][i(363)](o[i(2106)]),
                    Cesium[i(475)][i(363)](o[i(199)]),
                    o[i(2306)]
                  ]),
                    this._positionsTemp.concat([n])
                }, Cesium[e(1827)].MOUSE_MOVE),
                this[e(2490)][e(1068)]((t) => {
                  const i = e
                  if (this[i(2055)]) {
                    if (this._positions.length < 2) return this.stop(), void this.fire(p[i(289)])
                    if (i(2087) == this._drawOpts[i(1640)] && this._positions[i(277)] < 3)
                      return this.stop(), void this[i(450)](p[i(289)])
                    this[i(450)](p[i(182)], {
                      type: this._drawOpts[i(1640)],
                      positions: this[i(2469)]
                    }),
                      this[i(519)]()
                  }
                }, Cesium[e(1827)].RIGHT_CLICK)
            }
            [t(2060)]() {
              const e = t
              switch (this[e(1485)].type) {
                case e(811):
                  this[e(272)][e(1520)]([
                    e(1786),
                    '已有' +
                      this[e(2469)][e(277)] +
                      e(1898) +
                      (this._positions[e(277)] + 1) +
                      e(1471),
                    this[e(2469)][e(277)] < 2 ? e(2258) : e(921)
                  ])
                  break
                case e(2087):
                  this[e(272)][e(1520)]([
                    e(1456),
                    '已有' +
                      this[e(2469)][e(277)] +
                      e(1898) +
                      (this[e(2469)][e(277)] + 1) +
                      e(1471),
                    this[e(2469)][e(277)] < 3 ? '点击鼠标右键取消绘制' : e(921)
                  ])
              }
            }
            start(e) {
              const i = t
              switch (
                (this[i(2055)] && this[i(519)](),
                (this[i(1485)] = e),
                this._drawOpts[i(1598)] &&
                  Cesium[i(938)](this._polygonOpts, this[i(1485)].polygonOpts),
                this._drawOpts[i(1767)] && Cesium[i(938)](this[i(1995)], this[i(1485)][i(1767)]),
                (this[i(960)][i(261)][i(1482)] = !0),
                (this[i(2055)] = !0),
                e.type)
              ) {
                case i(1542):
                  this[i(272)][i(1520)]([i(2009), i(218)])
                  break
                case 'polyline':
                  this._mouseTip.setContent([i(1786), i(857), '点击鼠标右键取消绘制']),
                    this[i(2432)]()
                  break
                case i(2087):
                  this[i(272)][i(1520)]([i(1456), i(857), i(2258)]),
                    this[i(2432)](),
                    this._createPolygon()
              }
              this.fire(p[i(2279)])
            }
            [t(2432)]() {
              const e = t
              this[e(1967)] = this[e(245)][e(118)][e(1861)]({
                polyline: {
                  clampToGround: !0,
                  classificationType: Cesium[e(1597)].BOTH,
                  width: 2,
                  ...this[e(1995)],
                  positions: new Cesium.CallbackProperty(
                    (t) =>
                      e(811) == this[e(1485)][e(1640)]
                        ? this[e(264)]
                        : this._cartesian3ArrayTemp[e(277)] > 0
                          ? this._cartesian3ArrayTemp[e(1500)]([this._cartesian3ArrayTemp[0]])
                          : this[e(264)],
                    !1
                  )
                }
              })
            }
            [t(1706)]() {
              const e = t
              this[e(1967)] = this[e(245)].entities.add({
                polygon: {
                  classificationType: Cesium[e(1597)][e(1722)],
                  material: Cesium.Color[e(2474)][e(329)](0.2),
                  outline: !0,
                  outlineColor: Cesium[e(1154)][e(2474)],
                  outlineWidth: 1,
                  ...this[e(565)],
                  hierarchy: new Cesium[e(2569)]((t) => new Cesium[e(2229)](this[e(264)]), !1)
                }
              })
            }
            [t(2275)]() {
              const e = t
              this[e(1967)] &&
                (this._viewer[e(118)][e(1896)](this[e(1967)]), (this[e(1967)] = null))
            }
            [t(519)]() {
              const e = t
              this[e(2275)](),
                (this._positions = []),
                (this[e(2290)] = []),
                (this[e(1378)] = []),
                (this[e(264)] = []),
                (this[e(2055)] = !1),
                (this[e(960)][e(261)].show = !1),
                this[e(272)][e(1520)]()
            }
            [t(1712)]() {
              const e = t
              this[e(960)][e(1896)](),
                this[e(272)].remove(),
                this[e(519)](),
                this[e(2490)][e(2039)]()
            }
          },
          Measure: class extends v {
            constructor(e) {
              const i = t
              super(e), (this[i(245)] = e[i(395)]), (this[i(1932)] = null), (this[i(2243)] = [])
            }
            start(e) {
              const i = t
              switch (
                (this[i(519)](),
                (this[i(1656)] = Cesium[i(1960)](e[i(1640)], 'distance')),
                this[i(1656)])
              ) {
                case i(2346):
                  this[i(1932)] = new (class extends os {
                    constructor(t) {
                      super(t)
                    }
                    _start() {
                      const e = t
                      this[e(272)][e(1520)]([e(1352), '点击鼠标左键开始，点击鼠标右键取消'])
                    }
                    [t(240)]() {
                      const e = t
                      ;(this[e(192)] = []), (this[e(1182)] = null)
                    }
                    [t(514)]() {
                      const e = t
                      ;(this[e(715)] = this._viewer[e(118)].add({
                        polygon: {
                          ...this[e(565)],
                          perPositionHeight: !0,
                          hierarchy: new Cesium[e(2569)](
                            (t) => new Cesium[e(2229)](this[e(192)]),
                            !1
                          ),
                          extrudedHeight: new Cesium[e(2569)]((t) => this._getMaxHeight(), !1)
                        },
                        polyline: {
                          positions: new Cesium[e(2569)]((t) => this[e(446)](), !1),
                          ...this[e(282)]
                        }
                      })),
                        this._measureEntities.push(this[e(715)])
                    }
                    _createVertex() {
                      const e = t
                      let i = this[e(245)][e(118)][e(1861)]({
                        position: this._positions[this[e(2469)].length - 1],
                        type: e(2362),
                        point: { ...this[e(387)] }
                      })
                      this[e(2243)].push(i)
                    }
                    [t(2538)]() {
                      const e = t
                      ;(this[e(1182)] = this._viewer[e(118)][e(1861)]({
                        position: new Cesium[e(2569)]((t) => this[e(241)](), !1),
                        type: e(1310),
                        label: {
                          text: new Cesium.CallbackProperty(
                            (t) =>
                              '面积' +
                              ns(this._tempPositions).toFixed(this.fractionDigits) +
                              e(1212),
                            !1
                          ),
                          ...this._labelOpts
                        }
                      })),
                        this[e(2243)][e(2553)](this[e(1182)])
                    }
                    [t(446)]() {
                      const e = t
                      let i = this[e(1410)](),
                        s = []
                      if (this[e(192)][e(277)] < 1) return []
                      let n = [...this[e(192)]]
                      return (
                        n[e(2553)](n[0]),
                        n[e(1602)]((t) => {
                          const n = e,
                            o = Cesium[n(2285)][n(2579)](t),
                            r = Cesium[n(310)].fromRadians(o[n(2106)], o[n(199)], i)
                          s.push(r)
                        }),
                        s
                      )
                    }
                    [t(1410)]() {
                      const e = t
                      let i = -1e3
                      return (
                        this[e(192)][e(1602)]((t) => {
                          const s = e,
                            n = Cesium[s(2285)][s(2579)](t)
                          n[s(2306)] > i && (i = n[s(2306)])
                        }),
                        i
                      )
                    }
                    [t(241)]() {
                      const e = t
                      let i = -1e4,
                        s = new Cesium[e(310)]()
                      for (let t = 0; t < this[e(192)][e(277)]; t++) {
                        const n = this[e(192)][t],
                          o = Cesium.Cartographic[e(2579)](n)
                        o[e(2306)] > i && (i = o.height), Cesium[e(310)].add(s, n, s)
                      }
                      s = Cesium.Cartesian3[e(2317)](s, this[e(192)][e(277)], s)
                      const n = Cesium.Cartographic[e(2579)](s)
                      return (
                        (n.height = i),
                        Cesium[e(310)].fromRadians(n.longitude, n.latitude, n[e(2306)])
                      )
                    }
                    [t(1698)]() {
                      const e = t
                      this[e(437)].setInputAction((t) => {
                        const i = e
                        this[i(245)][i(374)].style[i(402)] = i(1653)
                        let s = this[i(245)][i(696)][i(1125)](t[i(2251)])
                        if (!s) {
                          const e = this[i(245)].scene[i(2557)].ellipsoid
                          s = this[i(245)][i(696)][i(1391)][i(1816)](t.position, e)
                        }
                        s &&
                          (this[i(2469)].push(s),
                          1 == this[i(2469)][i(277)] && this[i(514)](),
                          this[i(467)]())
                      }, Cesium[e(1827)][e(1282)])
                    }
                    [t(1008)]() {
                      const e = t
                      this[e(437)][e(1068)]((t) => {
                        const i = e
                        if (!this[i(526)]) return
                        this._viewer[i(374)].style[i(402)] = 'default'
                        let s = this[i(245)].scene.pickPosition(t.endPosition)
                        s ||
                          (s = this[i(245)][i(696)][i(1391)][i(1816)](
                            t[i(813)],
                            this[i(245)][i(696)][i(2557)].ellipsoid
                          )),
                          s &&
                            (this._mousePoint[i(1267)](s),
                            this[i(272)][i(1253)](t[i(954)]),
                            this[i(2469)][i(277)] < 3
                              ? this[i(272)].setContent([
                                  i(1352),
                                  '点击鼠标左键新增点，点击鼠标右键取消'
                                ])
                              : this[i(272)].setContent(['当前测量类型：水平面积', i(288)]),
                            this[i(2413)](s))
                      }, Cesium[e(1827)].MOUSE_MOVE)
                    }
                    [t(2413)](e) {
                      const i = t
                      this._positions[i(277)] < 1 ||
                        ((this._tempPositions = this[i(2469)][i(1500)](e)),
                        this[i(192)][i(277)] >= 3 && !this[i(1182)] && this[i(2538)]())
                    }
                    [t(284)]() {
                      const e = t
                      this[e(437)].setInputAction((t) => {
                        const i = e
                        if (!this[i(526)] || this[i(2469)][i(277)] < 3)
                          this[i(519)](), this[i(955)]()
                        else {
                          this[i(192)] = [...this[i(2469)]]
                          let t = this[i(446)]()
                          this[i(715)][i(811)].positions = new Cesium.CallbackProperty((e) => t, !1)
                          let e = new Cesium.PolygonHierarchy(this[i(192)])
                          ;(this[i(715)][i(2087)][i(603)] = new Cesium.CallbackProperty(
                            (t) => e,
                            !1
                          )),
                            (this[i(715)][i(2087)][i(398)] = this[i(1410)]()),
                            (this[i(1182)][i(2251)] = this[i(241)]()),
                            (this._measureResult = ns(this[i(2469)]).toFixed(this[i(632)])),
                            (this[i(1182)][i(2240)].text =
                              '总面积' + this._measureResult + i(1212)),
                            this._measureEnd()
                        }
                      }, Cesium[e(1827)][e(1309)])
                    }
                  })({ viewer: this[i(245)], ...e })
                  break
                case i(2306):
                  this[i(1932)] = new (class extends os {
                    constructor(t) {
                      super(t)
                    }
                    [t(188)]() {
                      const e = t
                      ;(this._circleRadius = 0.1),
                        (this[e(1795)] = 0),
                        this[e(272)][e(1520)]([e(847), e(1940)])
                    }
                    [t(1761)]() {
                      const e = t
                      ;(this[e(927)] = this._viewer[e(118)][e(1861)]({
                        polyline: {
                          ...this._lineOpts,
                          positions: new Cesium[e(2569)]((t) => this[e(2469)], !1)
                        }
                      })),
                        this[e(2243)][e(2553)](this[e(927)])
                    }
                    [t(1595)]() {
                      const e = t
                      ;(this[e(1634)] = this[e(245)][e(118)].add({
                        position: new Cesium[e(2569)](
                          (t) => this[e(2469)][this[e(2469)].length - 1],
                          !1
                        ),
                        label: { ...this._labelOpts }
                      })),
                        this[e(2243)][e(2553)](this._labelEntity)
                    }
                    _createVertex(e) {
                      const i = t
                      let s = this[i(245)].entities[i(1861)]({
                        position: this._positions[this[i(2469)][i(277)] - 1],
                        type: i(2362),
                        point: { ...this._pointOpts }
                      })
                      this[i(2243)].push(s), (this[i(923) + e] = s)
                    }
                    _createCircleEntitiy() {
                      const e = t
                      ;(this._circleEntity = this[e(245)][e(118)][e(1861)]({
                        position: new Cesium.CallbackProperty(
                          (t) => this[e(2469)][this[e(2469)].length - 1],
                          !1
                        ),
                        ellipse: {
                          numberOfVerticalLines: 36,
                          ...this[e(1962)],
                          height: new Cesium[e(2569)](
                            (t) => es(this[e(2469)][this[e(2469)].length - 1]),
                            !1
                          ),
                          semiMinorAxis: new Cesium[e(2569)]((t) => this[e(2312)], !1),
                          semiMajorAxis: new Cesium[e(2569)]((t) => this[e(2312)], !1)
                        }
                      })),
                        this._measureEntities[e(2553)](this[e(2170)])
                    }
                    _leftClickEvent() {
                      const e = t
                      this._handler.setInputAction((t) => {
                        const e = a0_0x3b79
                        this[e(245)]._element[e(1679)].cursor = e(1653)
                        let i = this[e(245)][e(696)][e(1125)](t[e(2251)])
                        if (!i) {
                          const s = this[e(245)].scene[e(2557)].ellipsoid
                          i = this[e(245)][e(696)].camera[e(1816)](t.position, s)
                        }
                        if (i)
                          if (0 == this[e(2469)].length)
                            this[e(2469)][e(2553)](i),
                              this[e(467)](0),
                              this[e(1761)](),
                              this[e(1771)](),
                              this[e(1595)]()
                          else {
                            let t = this[e(2469)][0],
                              i = this[e(2469)][1]
                            ;(this[e(1986)].position = t),
                              (this[e(861)][e(2251)] = i),
                              (this._lineEntity[e(811)][e(2333)] = new Cesium[e(2569)](
                                (e) => [t, i],
                                !1
                              )),
                              (this[e(1634)].position = this._positions[1]),
                              (this[e(2170)][e(2251)] = new Cesium[e(2569)]((e) => t, !1))
                            let s = this[e(2312)]
                            ;(this[e(2170)].ellipse.semiMinorAxis = s),
                              (this[e(2170)][e(1111)].semiMajorAxis = s),
                              (this._circleEntity[e(1111)][e(2306)] = es(
                                this[e(2469)][this[e(2469)][e(277)] - 1]
                              )),
                              this[e(272)][e(1520)](),
                              this._measureEnd()
                          }
                      }, Cesium.ScreenSpaceEventType[e(1282)])
                    }
                    [t(1008)]() {
                      const e = t
                      this[e(437)][e(1068)]((t) => {
                        const i = e
                        if (!this[i(526)]) return
                        this._viewer._element[i(1679)][i(402)] = i(1653)
                        let s = this[i(245)].scene[i(1125)](t[i(954)])
                        s ||
                          (s = this[i(245)][i(696)].camera[i(1816)](
                            t[i(813)],
                            this[i(245)].scene.globe[i(1565)]
                          )),
                          s &&
                            (this[i(272)][i(1253)](t[i(954)]),
                            1 == this._positions[i(277)] &&
                              this[i(272)][i(1520)]([i(847), i(1528)]),
                            this._mousePoint[i(1267)](s),
                            this[i(2413)](s))
                      }, Cesium.ScreenSpaceEventType[e(1073)])
                    }
                    [t(2413)](e) {
                      const i = t
                      if (this[i(2469)][i(277)] < 1) return
                      let s = is(this._positions[0]),
                        n = is(e)
                      const o = n[2] - s[2]
                      s[2] = n[2]
                      const r = Cesium.Cartesian3[i(667)](s[0], s[1], n[2])
                      var a, h
                      this[i(2469)][i(277)] < 2
                        ? (this[i(2469)][i(2553)](r), this._createVertex(1))
                        : ((this[i(2469)][1] = r),
                          (this[i(1795)] = o[i(1268)](this[i(2571)])),
                          (this[i(1634)].label[i(1982)] = '高度：' + this[i(1795)] + ' 米'),
                          (this[i(861)][i(2251)] = r)),
                        (this[i(2312)] =
                          ((a = this[i(2469)][0]),
                          (h = e),
                          (function (t, e) {
                            const i = a0_0x3b79
                            let s = (t[1] * Math.PI) / 180,
                              n = (e[1] * Math.PI) / 180,
                              o = s - n,
                              r = (t[0] * Math.PI) / 180 - (e[0] * Math.PI) / 180,
                              a =
                                2 *
                                Math.asin(
                                  Math.sqrt(
                                    Math[i(2394)](Math[i(884)](o / 2), 2) +
                                      Math[i(1272)](s) *
                                        Math[i(1272)](n) *
                                        Math.pow(Math[i(884)](r / 2), 2)
                                  )
                                )
                            return (a *= 6378.137), (a = Math[i(1647)](1e4 * a) / 10)
                          })(ts(a), ts(h))))
                    }
                    [t(284)]() {
                      const e = t
                      this[e(437)][e(1068)]((t) => {
                        const i = e
                        this[i(526)] && (this[i(519)](), this[i(1555)]())
                      }, Cesium[e(1827)].RIGHT_CLICK)
                    }
                  })({ viewer: this._viewer, ...e })
                  break
                default:
                  this[i(1932)] = new (class extends os {
                    constructor(t = {}) {
                      super(t)
                    }
                    _start() {
                      const e = t
                      this[e(272)][e(1520)]([e(141), e(1940)])
                    }
                    _stop() {
                      this[t(192)] = []
                    }
                    [t(1761)]() {
                      const e = t
                      ;(this[e(927)] = this[e(245)].entities.add({
                        polyline: {
                          positions: new Cesium.CallbackProperty((t) => this[e(192)], !1),
                          ...this[e(282)]
                        }
                      })),
                        this[e(2243)][e(2553)](this[e(927)])
                    }
                    [t(467)]() {
                      const e = t
                      let i = this[e(245)].entities.add({
                        position: this[e(2469)][this[e(2469)][e(277)] - 1],
                        type: e(1620),
                        label: {
                          text: ss(this._positions)[e(1268)](this[e(2571)]) + '米',
                          ...this._labelOpts
                        },
                        point: { ...this[e(387)] }
                      })
                      this[e(2243)].push(i), (this[e(1821)] = i)
                    }
                    [t(1825)]() {
                      const e = t
                      let i = this[e(245)].entities.add({
                        position: this[e(2469)][0],
                        type: e(1620),
                        label: { ...this[e(1088)], text: '起点' },
                        point: { ...this[e(387)] }
                      })
                      this[e(2243)][e(2553)](i)
                    }
                    [t(1808)]() {
                      const e = t
                      let i = this[e(1821)]
                      this[e(245)].entities.remove(i), this[e(245)][e(118)].remove(this[e(646)])
                      let s = this[e(245)][e(118)][e(1861)]({
                        position: this[e(2469)][this[e(2469)][e(277)] - 1],
                        type: e(1620),
                        label: {
                          ...this[e(1088)],
                          text: e(352) + ss(this._positions)[e(1268)](this._precision) + '米'
                        },
                        point: { ...this[e(387)] }
                      })
                      this._measureEntities.push(s)
                    }
                    [t(1698)]() {
                      const e = t
                      this[e(437)][e(1068)]((t) => {
                        const i = e
                        this[i(245)][i(374)].style[i(402)] = i(1653)
                        let s = this[i(245)][i(696)][i(1125)](t[i(2251)])
                        if (!s) {
                          const e = this[i(245)].scene.globe[i(1565)]
                          s = this[i(245)][i(696)][i(1391)][i(1816)](t[i(2251)], e)
                        }
                        if (s) {
                          if ((this[i(2469)][i(2553)](s), 1 == this[i(2469)][i(277)]))
                            return this[i(1761)](), void this[i(1825)]()
                          this[i(467)]()
                        }
                      }, Cesium[e(1827)][e(1282)])
                    }
                    [t(1008)]() {
                      const e = t
                      this[e(437)][e(1068)]((t) => {
                        const i = e
                        if (!this[i(526)]) return
                        this[i(272)][i(1253)](t[i(954)]),
                          this._positions.length < 2
                            ? this[i(272)][i(1520)]([i(141), i(1633)])
                            : this[i(272)].setContent([i(141), i(288)]),
                          (this[i(245)][i(374)][i(1679)].cursor = i(1653))
                        let s = this[i(245)][i(696)].pickPosition(t.endPosition)
                        s ||
                          (s = this[i(245)][i(696)].camera[i(1816)](
                            t[i(813)],
                            this[i(245)].scene[i(2557)][i(1565)]
                          )),
                          s && (this[i(960)][i(1267)](s), this[i(2413)](s))
                      }, Cesium[e(1827)][e(1073)])
                    }
                    [t(2413)](e) {
                      const i = t
                      this._positions[i(277)] < 1 ||
                        (this._tempPositions = this._positions[i(1500)](e))
                    }
                    [t(284)]() {
                      const e = t
                      this[e(437)][e(1068)]((t) => {
                        const i = e
                        if (!this[i(526)] || this[i(2469)][i(277)] < 2)
                          this[i(519)](), this[i(1555)]()
                        else {
                          this[i(1808)]()
                          let t = [...this[i(2469)]]
                          ;(this[i(927)].polyline[i(2333)] = new Cesium[i(2569)]((e) => t, !1)),
                            (this[i(1795)] = ss(this[i(2469)])[i(1268)](this._precision)),
                            this[i(2221)]()
                        }
                      }, Cesium.ScreenSpaceEventType[e(1309)])
                    }
                  })({ viewer: this[i(245)], ...e })
              }
              this[i(1932)][i(498)](),
                this[i(450)](p[i(244)]),
                this[i(1932)].on(p[i(2070)], (t) => {
                  const e = i
                  ;(this[e(2243)] = this[e(2243)][e(1500)](this[e(1932)]._measureEntities)),
                    this[e(450)](p.measureEnd, { result: Number(t), measureType: this[e(1656)] })
                })
            }
            [t(519)]() {
              const e = t
              this[e(1932)] && (this[e(1932)][e(519)](), (this._measureTool = null))
            }
            clear() {
              const e = t
              this[e(2243)].forEach((t) => {
                const i = e
                this[i(245)][i(118)].remove(t)
              }),
                (this._measureEntities = [])
            }
            [t(1712)]() {
              this[t(1555)](), this.stop()
            }
          },
          PointPicker: class extends v {
            constructor(e = {}) {
              const i = t
              super(e),
                (this[i(2137)] = new Date()[i(1566)]()),
                (this[i(245)] = e.viewer),
                (this._separator = Cesium[i(1960)](e.separator, ',')),
                (this._precision = Cesium[i(1960)](e[i(2540)], 4)),
                this._createDom(),
                this._viewer[i(436)][i(1392)][i(1621)](this[i(360)]),
                (this[i(1912)] = this[i(177)]()),
                (this[i(540)] = 0),
                this[i(1486)](this[i(1912)]),
                this._addPostRender(),
                this[i(2570)](),
                this._initEvents()
            }
            [t(2271)]() {
              const e = t
              ;(this[e(360)] = document.createElement(e(1120))),
                this[e(360)][e(2453)][e(1861)]('xt3d-point-picker-container')
              let i = document[e(1945)]('div')
              i[e(2453)][e(1861)](e(1903)), i.classList[e(1861)]('info-item-close')
              let s = document.createElement(e(1081))
              s[e(2453)].add(e(945)),
                (s[e(176)] = 'x'),
                (s[e(1038)] = '移除'),
                i.appendChild(s),
                this[e(360)][e(1621)](i),
                (i = document[e(1945)]('div'))[e(2453)].add(e(1903)),
                (i[e(1679)].cssText = 'margin:5px 1px;')
              let n = document[e(1945)](e(1081))
              ;(n[e(176)] = '小数位：'), i[e(1621)](n)
              let o = document[e(1945)]('input')
              o[e(2453)][e(1861)](e(1714)),
                o[e(2453)].add(e(2255)),
                (o[e(922)] = this._precision),
                (o.onchange = (t) => {
                  const i = e
                  ;(this[i(2571)] = o[i(922)]), this[i(1486)](this._position)
                }),
                i.appendChild(o),
                ((n = document[e(1945)](e(1081))).innerHTML = e(531)),
                (n[e(1679)][e(1811)] = e(889)),
                i[e(1621)](n)
              let r = document[e(1945)]('input')
              r.classList.add(e(2255)),
                (r[e(922)] = this[e(962)]),
                i[e(1621)](r),
                (r[e(1124)] = (t) => {
                  const i = e
                  ;(this[i(962)] = r.value), this[i(1486)](this[i(1912)])
                }),
                this[e(360)][e(1621)](i),
                (i = document[e(1945)](e(1120)))[e(2453)].add(e(1903)),
                ((n = document[e(1945)]('span'))[e(176)] = e(2228)),
                i[e(1621)](n)
              let a = document[e(1945)](e(1081))
              ;(a[e(176)] = '经纬度'),
                a[e(2453)][e(1861)]('value-type'),
                a.classList[e(1861)](e(850)),
                (a[e(1803)] = (t) => {
                  const i = e
                  h[i(2453)][i(1896)](i(850)),
                    l[i(2453)][i(1896)]('value-type-activate'),
                    a[i(2453)][i(1861)](i(850)),
                    (this[i(540)] = 0),
                    this[i(1486)](this[i(1912)])
                }),
                i[e(1621)](a)
              let h = document[e(1945)](e(1081))
              ;(h.innerHTML = e(1684)),
                (h[e(1679)][e(1811)] = e(2112)),
                h[e(2453)][e(1861)](e(1892)),
                (h[e(1803)] = (t) => {
                  const i = e
                  l[i(2453)][i(1896)](i(850)),
                    a[i(2453)][i(1896)](i(850)),
                    h[i(2453)].add(i(850)),
                    (this[i(540)] = 1),
                    this[i(1486)](this[i(1912)])
                }),
                i[e(1621)](h)
              let l = document.createElement('span')
              ;(l.innerHTML = '弧度'),
                (l[e(1679)].cssText = 'margin-left: 18px;'),
                l[e(2453)][e(1861)]('value-type'),
                (l[e(1803)] = (t) => {
                  const i = e
                  a.classList[i(1896)](i(850)),
                    h[i(2453)].remove(i(850)),
                    l[i(2453)][i(1861)]('value-type-activate'),
                    (this._valueType = 2),
                    this[i(1486)](this[i(1912)])
                }),
                i.appendChild(l),
                this[e(360)][e(1621)](i),
                (i = document[e(1945)](e(1120)))[e(2453)].add(e(1903)),
                (i[e(1679)][e(1811)] = e(648)),
                ((n = document.createElement('span')).innerHTML = e(483)),
                i.appendChild(n)
              let c = document[e(1945)](e(459))
              c[e(2453)][e(1861)](e(1714)),
                (c[e(1679)][e(1811)] = e(2588)),
                (this[e(1711)] = c),
                i[e(1621)](c)
              let u = document.createElement(e(2164))
              ;(u[e(176)] = '复制'),
                (u[e(1038)] = e(169)),
                (u[e(1803)] = (t) => {
                  const i = e
                  let s = this._getPoint(),
                    n = s.x + this[i(962)] + s.y + this[i(962)] + s.z
                  navigator[i(1110)][i(728)](n)
                }),
                i.appendChild(u),
                this[e(360)][e(1621)](i)
            }
            [t(1486)](e) {
              const i = t
              this[i(1912)] = e
              let s = this[i(296)]()
              return (
                (this[i(1711)].value =
                  s.x + this._separator + '\n' + s.y + this._separator + '\n' + s.z),
                s
              )
            }
            [t(296)]() {
              const e = t
              let i,
                s = Cesium[e(2285)].fromCartesian(this[e(1912)])
              return (i =
                0 == this[e(540)]
                  ? {
                      x: Cesium.Math[e(363)](s[e(2106)]).toFixed(this._precision),
                      y: Cesium[e(475)][e(363)](s.latitude)[e(1268)](this[e(2571)]),
                      z: s[e(2306)][e(1268)](this[e(2571)])
                    }
                  : 1 == this[e(540)]
                    ? {
                        x: this[e(1912)].x[e(1268)](this._precision),
                        y: this[e(1912)].y[e(1268)](this[e(2571)]),
                        z: this[e(1912)].z[e(1268)](this._precision)
                      }
                    : {
                        x: s[e(2106)].toFixed(this[e(2571)]),
                        y: s[e(199)].toFixed(this[e(2571)]),
                        z: s[e(2306)][e(1268)](this[e(2571)])
                      })
            }
            [t(177)]() {
              const e = t
              let i = new Cesium.Cartesian2(
                  this[e(245)].canvas[e(1481)] / 2,
                  this._viewer.canvas[e(281)] / 2
                ),
                s = this[e(245)][e(696)][e(1125)](i)
              if ((s || (s = this._viewer[e(1391)][e(1816)](i)), !s)) return
              let n = Cesium.Ellipsoid[e(276)][e(926)](s),
                o = (180 * n[e(2106)]) / Math.PI,
                r = (180 * n.latitude) / Math.PI
              return Cesium[e(310)][e(667)](o, r, n[e(2306)])
            }
            _addBillboard() {
              const e = t
              this[e(2386)] = this[e(245)].entities[e(1861)]({
                position: new Cesium[e(2569)]((t) => this[e(1912)], !1),
                type: e(2121),
                index: this[e(2137)],
                point: {
                  color: Cesium[e(1154)][e(971)],
                  pixelSize: 10,
                  outlineColor: Cesium[e(1154)][e(588)],
                  outlineWidth: 3,
                  disableDepthTestDistance: 1e3
                }
              })
            }
            [t(2185)]() {
              const e = t
              this[e(245)].scene[e(1271)][e(1973)](this[e(731)], this)
            }
            [t(731)]() {
              const e = t
              if (!this[e(360)] || !this[e(360)].style) return
              const i = this[e(245)][e(696)].canvas[e(2306)],
                s = new Cesium.Cartesian2()
              Cesium.SceneTransforms.worldToWindowCoordinates(
                this[e(245)][e(696)],
                this[e(1912)],
                s
              ),
                (this[e(360)][e(1679)][e(2251)] = 'absolute'),
                (this[e(360)][e(1679)][e(151)] = i - s.y + 70 + 'px')
              const n = this[e(360)].offsetWidth
              this[e(360)][e(1679)].left = s.x - n / 2 + 'px'
            }
            _initEvents() {
              const e = t
              ;(this[e(437)] = new Cesium[e(1561)](this[e(245)][e(696)][e(1493)])),
                this[e(2402)](),
                (this[e(360)].getElementsByClassName('info-close')[0][e(1803)] = (t) => {
                  this[e(1896)]()
                })
            }
            [t(2402)]() {
              const e = t
              this.initLeftClickEventHandler(), this[e(1117)](), this[e(1860)](), this[e(1925)]()
            }
            [t(344)]() {
              const e = t
              this[e(437)].removeInputAction(Cesium.ScreenSpaceEventType[e(1172)]),
                this[e(437)][e(1638)](Cesium[e(1827)].LEFT_UP),
                this[e(437)][e(1638)](Cesium[e(1827)][e(1073)]),
                this[e(437)][e(1638)](Cesium[e(1827)][e(1282)])
            }
            initLeftClickEventHandler() {
              const e = t
              this[e(437)][e(1068)]((t) => {
                const i = e
                let s = this[i(245)][i(696)][i(1125)](t[i(2251)])
                s && this._setPosition(s)
              }, Cesium[e(1827)].LEFT_CLICK)
            }
            [t(1117)]() {
              const e = t
              this[e(437)][e(1068)]((t) => {
                const i = e
                let s = this[i(245)][i(696)][i(1646)](t[i(2251)])
                if (s && s.id && s.id[i(1640)] && i(2121) == s.id[i(1640)]) {
                  if (s.id[i(2331)] != this._index) return
                  ;(this[i(1954)] = !0),
                    (this[i(245)].scene[i(691)].enableRotate = !1),
                    (this._viewer[i(471)] = !1),
                    (this._viewer._element[i(1679)][i(402)] = ''),
                    (document.body.style.cursor = i(1754))
                }
              }, Cesium[e(1827)][e(1172)])
            }
            [t(1925)]() {
              const e = t
              this[e(437)].setInputAction((t) => {
                const i = e
                ;(this[i(1954)] = !1),
                  (this[i(245)][i(471)] = !0),
                  (document[i(879)].style[i(402)] = i(1653)),
                  (this[i(245)][i(696)][i(691)][i(1551)] = !0)
              }, Cesium[e(1827)].LEFT_UP)
            }
            [t(1860)]() {
              const e = t
              this[e(437)][e(1068)]((t) => {
                const i = e
                let s = this[i(245)][i(696)].pickPosition(t.endPosition)
                s && this[i(1954)] && this[i(1486)](s)
              }, Cesium[e(1827)][e(1073)])
            }
            remove() {
              const e = t
              this[e(360)] &&
                (this._unRegisterEvents(),
                this[e(245)][e(118)][e(1896)](this[e(2386)]),
                this[e(245)][e(436)][e(1392)][e(2269)](this[e(360)]),
                this[e(245)][e(696)].postRender[e(704)](this[e(731)], this),
                (this[e(360)] = null))
            }
            [t(1712)]() {
              const e = t
              this.remove(), this[e(437)][e(2039)]()
            }
          }
        },
        Symbol[t(1540)],
        rs
      )
    )
    class hs {
      constructor(e) {
        const i = t
        ;(this[i(568)] = e), (this._type = '')
        let s = this[i(407)]()
        this[i(1259)] = this[i(1539)](s)
      }
      [t(1539)](t) {}
      get [t(1060)]() {
        return this[t(1259)]
      }
      set options(e) {
        ;(this[t(568)] = e), this._merge()
      }
      get type() {
        return this._type
      }
      [t(407)]() {
        const e = t
        let i = { ...this[e(568)] }
        if (i[e(1412)]) {
          const t = i[e(1412)]
          let s = new Cesium[e(613)](
            Cesium.Math.toRadians(t.west),
            Cesium[e(475)][e(1149)](t[e(1802)]),
            Cesium[e(475)][e(1149)](t[e(221)]),
            Cesium.Math[e(1149)](t[e(2149)])
          )
          i[e(1412)] = s
        }
        return i
      }
      [t(548)]() {
        const e = t
        let i = this[e(407)]()
        for (const t in i)
          if (Object[e(782)][e(1669)](i, t)) {
            const e = i[t]
            null != e && null != e && (this._provider[t] = e)
          }
      }
      get [t(1350)]() {
        return this[t(568)]
      }
      toJson() {
        const e = t,
          i = {}
        return (i.type = this[e(1640)]), (i.options = this[e(1350)]), i
      }
    }
    class ls extends hs {
      constructor(e) {
        const i = t
        super(e), (this[i(1436)] = i(2020))
      }
      [t(1539)](e) {
        return new Cesium[t(881)](e)
      }
    }
    class cs extends hs {
      constructor(e = {}) {
        const i = t
        super(e), (this._type = i(2194))
      }
      [t(1539)](e) {
        return new Cesium[t(641)](e)
      }
    }
    const us = {}
    ;(us[t(2303)] = ls), (us[t(307)] = cs)
    const ms = us,
      ps = 'none',
      ds = 'control',
      fs = t(1793),
      Cs = t(2448),
      vs = {}
    ;(vs[t(1793)] = t(1937)),
      (vs[t(125)] = t(679)),
      (vs.moveHeight = t(2043)),
      (vs.addMidPoint = 'rgba(245,49,232,0.99)')
    let _s = vs
    class gs extends v {
      constructor(e) {
        const i = t
        super(e),
          (this[i(245)] = e.graphicLayer[i(245)]),
          (this._graphicLayer = e[i(190)]),
          (this[i(2213)] = ps),
          (this[i(2153)] = []),
          (this[i(437)] = new Cesium[i(1561)](this[i(245)][i(696)][i(1493)])),
          this[i(1519)]()
      }
      get [t(2392)]() {
        const e = t
        return this[e(279)][e(2392)]
      }
      [t(2574)](e, i) {
        const s = t
        let n = this[s(245)][s(696)][s(1646)](e[s(2251)])
        n &&
          n.id &&
          'EditAnchor' == n.id[s(1640)] &&
          ((this._editAnchor = n.id),
          (this[s(2213)] = n.id.editAnchorType),
          (this._editAnchorIndex = n.id[s(2331)]),
          this.setCursor(s(1754)),
          this._writeLeftDow(),
          this[s(272)].setContent([s(2028)]),
          this[s(272)][s(1253)](e[s(2251)]),
          (this[s(245)][s(696)].screenSpaceCameraController[s(1551)] = !1),
          this[s(2213)] == Cs && this[s(1998)](s(2476)),
          this[s(772)](e, i))
      }
      [t(772)](t, e) {}
      [t(767)](e) {
        const i = t
        if ((this._mouseTip[i(1253)](e[i(954)]), this._editAnchorType == ps))
          return void this[i(2298)](e)
        let s = this[i(561)](e)
        if (s) {
          if (
            (this._editAnchorType == ds &&
              ((this._controlPositions[this._editAnchorIndex] = s),
              this[i(1806)](s, this._editAnchorIndex)),
            this[i(2213)] == fs && this[i(612)](s),
            this[i(2213)] == Cs)
          ) {
            let t = this[i(2082)](e)
            this[i(1704)](e, t)
          }
          this[i(2456)](e), this._initAnchorPositions(), this._writeLeftDow()
        }
      }
      [t(2456)](t) {}
      [t(816)]() {
        const e = t
        let i = this[e(1375)].position[e(1474)]()
        i &&
          ((this[e(317)] = i),
          (this[e(617)] = Cesium[e(2058)].eastNorthUpToFixedFrame(i)),
          (this._inverseMatrix4 = Cesium[e(2066)][e(814)](this[e(617)], new Cesium.Matrix4())))
      }
      [t(2082)](e) {
        const i = t
        let s = this[i(317)][i(1902)](),
          n = Cesium.Cartesian3[i(1676)](
            this[i(245)][i(696)][i(1391)][i(2251)],
            s,
            new Cesium.Cartesian3()
          )
        n = Cesium[i(310)][i(379)](n, n)
        const o = Cesium.Plane[i(1163)](s, n),
          r = this[i(245)].scene[i(1391)][i(2033)](e[i(954)])
        if (r && o) {
          const t = Cesium.IntersectionTests.rayPlane(r, o)
          return (
            this._viewer.scene[i(2557)][i(1565)].cartesianToCartographic(t)[i(2306)] -
            Cesium.Cartographic[i(2579)](s).height
          )
        }
        return 0
      }
      _mouseMoveALT(e) {
        const i = t
        if (this[i(2213)] == ps) return
        this[i(272)].setPosition(e[i(954)])
        let s = this[i(2082)](e),
          n = Cesium.Matrix4[i(1420)](new Cesium[i(310)](0, 0, s), new Cesium[i(2066)]())
        this[i(246)](e, n), this._writeLeftDow(), this[i(1998)](i(2476))
      }
      [t(246)](t, e) {}
      [t(2006)](e, i) {
        const s = t
        this[s(2153)][s(1602)]((t) => {
          t[s(1482)] = !0
        }),
          this[s(1998)](''),
          this[s(2213)] != ps &&
            ((this._editAnchorType = ps),
            (this[s(245)][s(696)][s(691)][s(1551)] = !0),
            this[s(2222)](e, i),
            this[s(970)](!0))
      }
      handleLeftUp(t, e) {}
      [t(2298)](e) {
        const i = t
        this._mouseTip[i(1520)]()
        let s = this[i(245)][i(696)].pick(e[i(954)])
        s &&
          s.id &&
          i(541) == s.id[i(1640)] &&
          (this[i(272)][i(1520)](s.id[i(155)]), this[i(576)](e))
      }
      [t(576)]() {}
      udpateControlPosition(e, i) {
        const s = t
        let n = Cesium[s(2285)][s(2579)](e)
        const o = [Cesium.Math[s(363)](n.longitude), Cesium[s(475)][s(363)](n[s(199)]), n[s(2306)]]
        let r = this[s(2392)][s(2333)]
        ;(r[i] = o), this[s(2392)][s(2515)](r)
      }
      updateMoveAllPosition(t) {}
      [t(1704)](t, e) {}
      [t(1862)](e, i) {
        const s = t
        return (
          (i = Cesium[s(2066)][s(2483)](this[s(2293)], i, new Cesium[s(310)]())),
          (i = Cesium[s(2066)].multiplyByPoint(e, i, new Cesium[s(310)]())),
          Cesium.Matrix4[s(2483)](this[s(617)], i, new Cesium[s(310)]())
        )
      }
      [t(561)](e) {
        const i = t
        let s = this[i(245)].scene[i(1125)](e[i(954)])
        return (
          s ||
            this._viewer.scene[i(1391)][i(1816)](e[i(954)], this._viewer[i(696)][i(2557)][i(1565)]),
          s
        )
      }
      [t(2402)]() {
        const e = t
        this[e(437)][e(1068)]((t) => {
          this._leftDown(t)
        }, Cesium[e(1827)][e(1172)]),
          this[e(437)][e(1068)]((t) => {
            this[e(767)](t)
          }, Cesium[e(1827)][e(1073)]),
          this[e(437)].setInputAction((t) => {
            this[e(2006)](t)
          }, Cesium[e(1827)][e(1338)]),
          this[e(437)][e(1068)](
            (t) => {
              this[e(2574)](t, !0)
            },
            Cesium.ScreenSpaceEventType[e(1172)],
            Cesium[e(1635)][e(542)]
          ),
          this[e(437)].setInputAction(
            (t) => {
              this._mouseMoveALT(t)
            },
            Cesium[e(1827)][e(1073)],
            Cesium[e(1635)][e(542)]
          ),
          this[e(437)][e(1068)](
            (t) => {
              this._leftUp(t, !0)
            },
            Cesium[e(1827)][e(1338)],
            Cesium[e(1635)][e(542)]
          )
      }
      [t(344)]() {
        const e = t
        this[e(437)][e(1638)](Cesium[e(1827)][e(1172)]),
          this[e(437)].removeInputAction(Cesium[e(1827)][e(1073)]),
          this[e(437)].removeInputAction(Cesium.ScreenSpaceEventType[e(1338)]),
          this[e(437)][e(1638)](Cesium[e(1827)][e(1172)], Cesium[e(1635)][e(542)]),
          this[e(437)][e(1638)](Cesium[e(1827)].MOUSE_MOVE, Cesium[e(1635)][e(542)]),
          this[e(437)][e(1638)](
            Cesium.ScreenSpaceEventType[e(1338)],
            Cesium.KeyboardEventModifier[e(542)]
          )
      }
      [t(1519)]() {
        const e = t
        this._activate || (this[e(2402)](), (this[e(2031)] = !0), this[e(546)]())
      }
      [t(546)]() {
        ;(this[t(272)] = new $i(this._viewer)),
          this._initAnchorPositions(!0),
          this._createEditAnchors(),
          this.activateHook()
      }
      [t(1972)]() {}
      [t(175)]() {
        const e = t
        let i = ''
        this._controlPositions[e(1602)]((t, s) => {
          const n = e
          let o = new Cesium.CallbackProperty((t) => this._controlPositions[s], !1)
          ;(i = this[n(1569)](ds)), this[n(1003)](s, ds, o, i)
        }),
          this[e(1441)][e(1602)]((t, s) => {
            const n = e
            let o = new Cesium[n(2569)]((t) => this[n(1441)][s], !1)
            ;(i = this[n(1569)](Cs)), this[n(1003)](s, Cs, o, i)
          })
        let s = new Cesium.CallbackProperty((t) => this[e(2109)], !1)
        ;(i = this[e(1569)](fs)), this[e(1003)](0, fs, s, i)
      }
      [t(970)](e) {
        const i = t
        ;(this[i(1094)] = []), (this[i(1441)] = []), (this._moveAllPosition = null)
        let s = new Cesium.Cartesian3()
        if (this[i(2392)][i(2333)] && this[i(2392)].positions[i(277)] > 0) {
          if (
            (this[i(2392)][i(2333)][i(1602)]((t) => {
              const e = i
              let n = Cesium[e(310)][e(667)](t[0], t[1], t[2])
              ;(s = Cesium[e(310)][e(1861)](s, n, new Cesium[e(310)]())), this[e(1094)][e(2553)](n)
            }),
            e && this[i(2392)][i(1679)][i(2080)])
          ) {
            let t = [],
              e = []
            this[i(1094)].forEach((s) => {
              const n = i,
                o = Cesium[n(2285)][n(2579)](s),
                r = this[n(245)][n(696)].sampleHeight(o, this[n(245)].entities.values)
              let a = Cesium[n(310)][n(1069)](o[n(2106)], o[n(199)], r)
              t[n(2553)](a)
              const h = [Cesium.Math[n(363)](o[n(2106)]), Cesium[n(475)][n(363)](o[n(199)]), r]
              e.push(h)
            }),
              (this[i(1094)] = t),
              this[i(2392)][i(2515)](e)
          }
          let t = this[i(2392)][i(2333)].length
          if (
            ((s = Cesium[i(310)].divideByScalar(s, t, new Cesium.Cartesian3())),
            (this[i(2109)] = s),
            e && this[i(2392)][i(1679)][i(2080)])
          ) {
            let t = Cesium.Cartographic.fromCartesian(s)
            const e = this[i(245)][i(696)].sampleHeight(t, this._viewer[i(118)][i(736)])
            this[i(2109)] = Cesium[i(310)][i(1069)](t[i(2106)], t[i(199)], e)
          }
        }
        this[i(567)](e)
      }
      [t(567)]() {}
      getTooltip(t) {}
      [t(701)]() {
        const e = t
        ;(this[e(245)][e(696)][e(691)][e(1551)] = !0),
          this.setCursor(''),
          this[e(344)](),
          (this._activate = !1),
          this[e(2459)]()
      }
      [t(2459)]() {
        const e = t
        this[e(1675)](),
          this[e(272)] && this._mouseTip[e(1896)](),
          (this[e(272)] = void 0),
          this[e(2016)]()
      }
      [t(2016)]() {}
      [t(1675)]() {
        const e = t
        this._editAnchors[e(1602)]((t) => {
          const i = e
          this[i(245)][i(118)][i(1896)](t)
        }),
          (this[e(2153)] = [])
      }
      [t(1003)](e, i, s, n) {
        const o = t
        let r = _s[i],
          a = this[o(245)].entities[o(1861)]({
            index: e,
            editAnchorType: i,
            type: o(541),
            position: s,
            tooltip: n,
            point: {
              pixelSize: 10,
              color: Cesium[o(1154)][o(2008)](r),
              outline: !0,
              outlineWidth: 2,
              outlineColor: Cesium[o(1154)].WHITE[o(329)](0.5),
              disableDepthTestDistance: 2e9,
              zIndex: 999
            }
          })
        this[o(2153)][o(2553)](a)
      }
      _beforDestroy() {
        const e = t
        this[e(701)](), this[e(437)][e(2039)]()
      }
      getZMoveMatrix() {}
      [t(1998)](e) {
        const i = t
        this[i(245)][i(374)][i(1134)][i(1679)][i(402)] = e
      }
    }
    class ys extends gs {
      constructor(t) {
        super(t)
      }
      initAnchorPositions(e) {
        const i = t
        ;(this._moveAllPosition = this[i(2392)].cartesian3),
          (this[i(1094)] = []),
          (this[i(1441)] = [])
      }
      [t(246)](e, i) {
        const s = t
        let n = this[s(2392)][s(2251)],
          o = this._leftDownPosition
        o = this[s(1862)](i, o)
        let r = Cesium.Cartographic[s(2579)](o)
        ;(n = [
          Cesium[s(475)][s(363)](r[s(2106)]),
          Cesium.Math[s(363)](r.latitude),
          r[s(2306)] < 0 ? 0 : r[s(2306)]
        ]),
          this.editGraphic[s(1253)](n),
          (this._moveAllPosition = this[s(2392)][s(1131)]),
          this[s(1839)]()
      }
      [t(1839)]() {}
      updateMoveAllPosition(e) {
        const i = t
        let s = Cesium.Cartographic[i(2579)](e)
        ;(e = [Cesium[i(475)].toDegrees(s.longitude), Cesium.Math[i(363)](s[i(199)]), s[i(2306)]]),
          this[i(2392)][i(1253)](e)
      }
      [t(1569)](e) {
        const i = t
        return [i(1266), i(1882)]
      }
    }
    class ws extends gs {
      constructor(t) {
        super(t)
      }
      handleLeftDown() {
        const e = t
        this[e(2392)].style.clampToGround &&
          this[e(2213)] == fs &&
          this._editAnchors[e(1602)]((t) => {
            t[e(1482)] = !1
          }),
          this[e(2213)] == fs && (this._editAnchor[e(1482)] = !1)
      }
      [t(567)]() {
        const e = t
        if (
          ((this._moveHeightPositions = []),
          !this[e(2392)][e(1679)][e(2080)] && 0 != this[e(2392)][e(1679)][e(398)])
        ) {
          let t = this.editGraphic._entity[e(2087)][e(398)]
          this[e(2392)][e(2333)][e(1602)]((i) => {
            const s = e,
              n = Cesium[s(310)][s(667)](i[0], i[1], t)
            this[s(1441)][s(2553)](n)
          })
        }
      }
      [t(246)](e, i) {
        const s = t
        if (
          !(
            this.editGraphic[s(1679)][s(2080)] ||
            this._editAnchorType == Cs ||
            (this[s(2213)] == ds && this[s(2392)][s(1679)][s(398)] > 0)
          )
        ) {
          if (this[s(2213)] == ds && 0 == this[s(2392)].style[s(398)]) {
            let t = this[s(1862)](i, this[s(317)])
            return (
              (this[s(1094)][this[s(930)]] = t), void this.udpateControlPosition(t, this[s(930)])
            )
          }
          this[s(562)](i), this[s(970)]()
        }
      }
      [t(561)](e) {
        const i = t
        let s
        return (
          this[i(2213)] == fs || (!this[i(2392)][i(1679)].clampToGround && this[i(2213)] == ds)
            ? (s = J(this[i(245)], e[i(954)], this[i(317)]))
            : (s = this[i(245)][i(696)][i(1125)](e.endPosition)) ||
              this[i(245)][i(696)][i(1391)][i(1816)](
                e.endPosition,
                this[i(245)][i(696)][i(2557)][i(1565)]
              ),
          s
        )
      }
      [t(1704)](e, i) {
        const s = t
        ;(this[s(2392)][s(1679)].extrudedHeight += i),
          this[s(2392)][s(2515)](this[s(2392)][s(2469)])
      }
      [t(612)](e) {
        const i = t
        ;(e = Cesium.Matrix4[i(2483)](this._inverseMatrix4, e, new Cesium[i(310)]())).z = 0
        let s = Cesium[i(2066)][i(1420)](e, new Cesium[i(2066)]())
        this[i(562)](s)
      }
      [t(562)](e) {
        const i = t
        let s = [],
          n = []
        this[i(2392)].positions[i(1602)]((t) => {
          const e = i,
            s = Cesium[e(310)][e(667)](t[0], t[1], t[2])
          n[e(2553)](s)
        })
        for (let t = 0; t < n.length; t++) {
          let o = n[t]
          o = this.getMoveByMatrix(e, o)
          let r = Cesium.Cartographic[i(2579)](o)
          this[i(1094)][t] = o
          const a = [
            Cesium[i(475)][i(363)](r.longitude),
            Cesium[i(475)][i(363)](r.latitude),
            r[i(2306)]
          ]
          s.push(a)
        }
        this.editGraphic[i(2515)](s)
      }
      [t(1569)](e) {
        const i = t
        return e == ds
          ? this[i(2392)][i(1679)][i(2080)] || 0 != this[i(2392)][i(1679)][i(398)]
            ? ['拖拽改变位置']
            : [i(1266), i(1882)]
          : e == fs
            ? this[i(2392)].style[i(2080)]
              ? ['拖拽改变位置']
              : [i(607), i(1800)]
            : e == Cs
              ? ['拖拽改变高度']
              : void 0
      }
      [t(1972)]() {
        const e = t
        this[e(2392)][e(1920)] = (t) => {
          const i = e
          this[i(1675)](),
            setTimeout(() => {
              const t = i
              this[t(970)](!0), this[t(175)]()
            }, 200)
        }
      }
    }
    class xs extends ws {
      constructor(t) {
        super(t)
      }
    }
    class bs extends ys {
      constructor(t) {
        super(t)
      }
      [t(561)](e) {
        const i = t
        return J(this[i(245)], e[i(954)], this[i(317)])
      }
    }
    class Ss extends gs {
      constructor(t) {
        super(t)
      }
      [t(246)](e, i) {
        const s = t
        if (!this[s(2392)].style.clampToGround && this[s(2213)] != Cs) {
          if (this._editAnchorType == ds) {
            let t = this[s(1862)](i, this[s(317)])
            return (
              (this[s(1094)][this[s(930)]] = t),
              this[s(1806)](t, this._editAnchorIndex),
              void this[s(970)]()
            )
          }
          this[s(562)](i), this[s(970)]()
        }
      }
      [t(567)]() {
        const e = t
        if (this[e(2392)][e(1158)] == Q.uprightLine) {
          ;(this[e(1094)] = []), (this[e(2109)] = this.editGraphic[e(1378)][0])
          let t = Cesium[e(2285)][e(2579)](this[e(2109)])
          this[e(1441)] = [
            Cesium[e(310)].fromRadians(
              t[e(2106)],
              t[e(199)],
              t[e(2306)] + this[e(2392)][e(1679)][e(2306)]
            )
          ]
        }
      }
      pickMousePosition(e) {
        const i = t
        let s
        return (
          this._editAnchorType == fs
            ? (s = J(this[i(245)], e.endPosition, this[i(317)]))
            : (s = this[i(245)][i(696)].pickPosition(e[i(954)])) ||
              this[i(245)][i(696)][i(1391)].pickEllipsoid(
                e.endPosition,
                this[i(245)].scene.globe[i(1565)]
              ),
          s
        )
      }
      [t(612)](e) {
        const i = t
        ;(e = Cesium[i(2066)].multiplyByPoint(this._inverseMatrix4, e, new Cesium.Cartesian3())).z =
          0
        let s = Cesium.Matrix4[i(1420)](e, new Cesium.Matrix4())
        this[i(562)](s)
      }
      [t(562)](e) {
        const i = t
        let s = [],
          n = []
        this[i(2392)][i(2333)][i(1602)]((t) => {
          const e = i,
            s = Cesium.Cartesian3[e(667)](t[0], t[1], t[2])
          n[e(2553)](s)
        })
        for (let t = 0; t < n[i(277)]; t++) {
          let o = n[t]
          o = this[i(1862)](e, o)
          let r = Cesium.Cartographic[i(2579)](o)
          this[i(1094)][t] = o
          const a = [
            Cesium.Math[i(363)](r[i(2106)]),
            Cesium[i(475)][i(363)](r.latitude),
            r[i(2306)]
          ]
          s[i(2553)](a)
        }
        this.editGraphic[i(2515)](s)
      }
      [t(1569)](e) {
        const i = t
        let s = []
        return (
          e == ds &&
            ((s = [i(1266)]),
            this.editGraphic[i(343)].clampToGround || s[i(2553)]('按住alt拖拽改变高度')),
          e == fs && ((s = [i(607)]), this[i(2392)][i(343)][i(2080)] || s[i(2553)](i(1800))),
          e == Cs && (s = [i(1874)]),
          s
        )
      }
      [t(1704)](e, i) {
        const s = t
        ;(this[s(2392)][s(1679)][s(2306)] += i), this.editGraphic[s(1253)](this[s(2392)][s(1912)])
      }
      [t(1972)]() {
        const e = t
        this[e(2392)][e(1920)] = (t) => {
          this[e(970)]()
        }
      }
    }
    class Ps extends gs {
      constructor(t) {
        super(t)
      }
      [t(246)](e, i) {
        const s = t
        if (this[s(2213)] == ds) {
          let t = this[s(1862)](i, this._leftDownPosition)
          return (
            (this._controlPositions[this[s(930)]] = t),
            this[s(1806)](t, this[s(930)]),
            void this._initAnchorPositions()
          )
        }
        this.updateAllPositionByMatrix(i), this[s(970)]()
      }
      [t(561)](e) {
        const i = t
        return J(this[i(245)], e.endPosition, this[i(317)])
      }
      [t(612)](e) {
        const i = t
        ;(e = Cesium.Matrix4[i(2483)](this[i(2293)], e, new Cesium[i(310)]())).z = 0
        let s = Cesium[i(2066)][i(1420)](e, new Cesium[i(2066)]())
        this.updateAllPositionByMatrix(s)
      }
      [t(562)](e) {
        const i = t
        let s = [],
          n = []
        this[i(2392)].positions[i(1602)]((t) => {
          const e = i,
            s = Cesium[e(310)][e(667)](t[0], t[1], t[2])
          n[e(2553)](s)
        })
        for (let t = 0; t < n[i(277)]; t++) {
          let o = n[t]
          o = this[i(1862)](e, o)
          let r = Cesium[i(2285)][i(2579)](o)
          this[i(1094)][t] = o
          const a = [
            Cesium[i(475)][i(363)](r[i(2106)]),
            Cesium[i(475)][i(363)](r[i(199)]),
            r[i(2306)]
          ]
          s[i(2553)](a)
        }
        this[i(2392)][i(2515)](s)
      }
      [t(1569)](e) {
        const i = t
        let s = []
        return e == ds && (s = [i(121), i(1882)]), e == fs && (s = [i(607), i(1800)]), s
      }
    }
    let Ms = t(1754),
      As = t(494),
      Ts = t(2313),
      Es = t(2022),
      zs = 'rotate'
    class Ds extends gs {
      constructor(t) {
        super(t)
      }
      [t(2574)](t, e) {
        this._svgEditAnchors.forEach((e) => {
          const i = a0_0x3b79,
            s = e[i(2132)]
          let n = t[i(2251)]
          Math.abs(n.x - (s[i(805)] + 6)) < 4 &&
            Math[i(636)](n.y - (s[i(1129)] + 6)) < 4 &&
            (this[i(2450)](s),
            (this[i(1403)] = s[i(1548)][i(789)]),
            (this[i(1605)] = !0),
            (this[i(245)][i(696)].screenSpaceCameraController.enableRotate = !1))
        })
      }
      _mouseMove(e) {
        const i = t
        if ((this[i(272)][i(1253)](e[i(954)]), !this[i(1605)])) return void this[i(2298)](e)
        if (this._editType == Ms) return void this[i(601)](e)
        let s = Cesium[i(2176)][i(708)](
          this._viewer[i(696)],
          this[i(2392)][i(2238)],
          new Cesium.Cartesian2()
        )
        if (!s) return
        let n = Cesium.Cartesian2[i(1849)](s, e[i(954)]) - Cesium[i(194)][i(1849)](s, e[i(813)]),
          o = this.editGraphic[i(1275)]()
        if ((this[i(1403)] == zs && this[i(1774)](e, s), this[i(1403)] == As)) {
          let t = ((o[i(2306)] + n) / this[i(2392)][i(1462)])[i(1268)](3),
            e = ((o.width + n) / this[i(2392)][i(1262)]).toFixed(3)
          this[i(2392)][i(544)](t, e)
        }
        if (this._editType == Ts) {
          let t = ((o.height + n) / this[i(2392)][i(1462)])[i(1268)](3)
          this[i(2392)].setScale(t, null)
        }
        if (this[i(1403)] == Es) {
          let t = ((o[i(575)] + n) / this[i(2392)][i(1262)]).toFixed(3)
          this.editGraphic[i(544)](null, t)
        }
      }
      [t(601)](e) {
        const i = t
        let s = this[i(245)][i(696)][i(1125)](e.endPosition)
        if (!s) return
        let n = Cesium[i(2285)][i(2579)](s)
        n[i(2306)] < 0 && (n[i(2306)] = 0),
          this[i(2392)][i(1253)]([
            Cesium[i(475)][i(363)](n.longitude),
            Cesium[i(475)][i(363)](n[i(199)]),
            n.height
          ])
      }
      _setRotate(e, i) {
        const s = t
        let n = Cesium.Cartesian2.subtract(e[s(954)], i, new Cesium[s(194)]())
        n.y = -n.y
        let o = Cesium[s(194)][s(379)](n, new Cesium[s(194)]()),
          r = new Cesium[s(194)](0, 1),
          a = Cesium.Cartesian2.dot(r, o),
          h = Math[s(196)](a),
          l = Cesium.Math.toDegrees(h)
        n.x > 0 && (l = 360 - l), this[s(2392)].setRotate(parseInt(l))
      }
      [t(2298)](e) {
        const i = t
        let s = null,
          n = e[i(954)]
        this[i(1848)][i(1602)]((t) => {
          const e = i,
            o = t.dom
          if (
            Math[e(636)](n.x - (o.offsetLeft + 6)) < 4 &&
            Math[e(636)](n.y - (o.offsetTop + 6)) < 4
          ) {
            let t = o[e(1548)][e(789)]
            switch ((this[e(2450)](o), (this[e(2233)] = o), t)) {
              case Ms:
                ;(s = [e(502)]), (this._viewer[e(374)].style[e(402)] = e(1754))
                break
              case As:
                ;(s = [e(1562)]), (this._viewer._element.style.cursor = o[e(1548)][e(402)])
                break
              case Ts:
                ;(s = [e(405)]), (this[e(245)]._element[e(1679)][e(402)] = 'ns-resize')
                break
              case Es:
                ;(s = [e(452)]), (this[e(245)][e(374)].style[e(402)] = e(371))
                break
              case zs:
                ;(s = [e(1043)]), (this[e(245)][e(374)][e(1679)].cursor = e(1754))
            }
          }
        }),
          s ||
            (this[i(2233)] && (this[i(1372)](this._highlightDom), (this[i(2233)] = null)),
            (this[i(245)]._element[i(1679)][i(402)] = '')),
          this[i(272)][i(1520)](s)
      }
      [t(2006)](e) {
        const i = t
        ;(this[i(245)][i(374)][i(1679)][i(402)] = ''),
          (this._editType = null),
          (this[i(1605)] = !1),
          (this[i(245)][i(696)][i(691)][i(1551)] = !0),
          this[i(1848)][i(1602)]((t) => {
            const e = i
            this[e(1372)](t[e(2132)])
          })
      }
      [t(1972)]() {
        const e = t
        this[e(245)].scene[e(1271)][e(1973)](this[e(731)], this), this._createSvgEditAnchors()
      }
      [t(731)]() {
        const e = t
        this._svgEditAnchors &&
          this[e(1848)].forEach((t) => {
            t.update()
          })
      }
      [t(818)]() {
        const e = t
        ;(this[e(1848)] = []),
          this._addMoveAnchor(),
          this[e(1196)](),
          this[e(1686)](),
          this[e(1817)]()
      }
      [t(1196)]() {
        const e = t
        let i = document[e(1945)](e(1120))
        ;(i.attributes.eType = zs), this[e(1372)](i), this._viewer[e(436)][e(1392)][e(1621)](i)
        const s = {}
        s.dom = i
        let n = s
        ;(n.update = (t) => {
          const i = e
          let s = this[i(2392)].getSize()[i(2306)] / 2 + 20,
            o = Cesium[i(475)].toRadians(this[i(2392)][i(1679)][i(1837)])
          this._updateAnchor(n, 0, s, o)
        }),
          this._svgEditAnchors[e(2553)](n)
      }
      [t(2253)]() {
        const e = t
        let i = document[e(1945)](e(1120))
        ;(i[e(1548)][e(789)] = Ms), this[e(1372)](i), this[e(245)][e(436)][e(1392)][e(1621)](i)
        const s = {}
        s.dom = i
        let n = s
        ;(n.update = (t) => {
          this._updateAnchor(n, 0, 0, 0)
        }),
          this[e(1848)].push(n)
      }
      _addHWAnchors() {
        const e = t
        this[e(2399)](Ts, 0, 1),
          this[e(2399)](Es, -1, 0),
          this.addHWAnchor(Ts, 0, -1),
          this[e(2399)](Es, 1, 0)
      }
      [t(2399)](e, i, s) {
        const n = t
        let o = document[n(1945)](n(1120))
        ;(o[n(1548)][n(789)] = e), this[n(1372)](o), this._viewer.cesiumWidget.container[n(1621)](o)
        const r = {}
        r.dom = o
        let a = r
        ;(a[n(720)] = (t) => {
          const e = n
          let o = this.editGraphic[e(1275)](),
            r = (i * o[e(575)]) / 2,
            h = (s * o[e(2306)]) / 2,
            l = Cesium[e(475)].toRadians(this[e(2392)][e(1679)][e(1837)])
          this[e(1984)](a, r, h, l)
        }),
          this[n(1848)][n(2553)](a)
      }
      [t(1817)]() {
        const e = t
        this[e(1935)](1, 1, e(2454)),
          this[e(1935)](-1, 1, e(2514)),
          this[e(1935)](-1, -1, e(2454)),
          this[e(1935)](1, -1, e(2514))
      }
      [t(1935)](e, i, s) {
        const n = t
        let o = document[n(1945)](n(1120))
        ;(o.attributes[n(789)] = As),
          (o[n(1548)][n(402)] = s),
          this._resetAnchorStyle(o),
          this[n(245)][n(436)][n(1392)].appendChild(o)
        const r = {}
        r[n(2132)] = o
        let a = r
        ;(a[n(720)] = (t) => {
          const s = n
          let o = this[s(2392)].getSize(),
            r = (e * o.width) / 2,
            h = (i * o[s(2306)]) / 2,
            l = Cesium[s(475)][s(1149)](this[s(2392)][s(1679)][s(1837)])
          this[s(1984)](a, r, h, l)
        }),
          this._svgEditAnchors[n(2553)](a)
      }
      [t(1372)](e) {
        const i = t
        switch (e[i(1548)].eType) {
          case zs:
            e[i(1679)].cssText = i(2244)
            break
          case Ms:
            e[i(1679)][i(1811)] = i(2460)
            break
          default:
            e[i(1679)].cssText = i(1305)
        }
      }
      _highlightAnchor(e) {
        const i = t
        ;(e.style[i(1277)] = i(320)), (e.style[i(1900)] = '1px solid yellow')
      }
      [t(1984)](e, i, s, n) {
        const o = t
        let r = i * Math[o(1272)](n) - s * Math[o(884)](n),
          a = s * Math[o(1272)](n) + i * Math[o(884)](n),
          h = e[o(2114)] ? e.windowPosition : new Cesium[o(194)]()
        e.windowPosition = h
        const l = this._viewer[o(696)][o(1493)][o(2306)]
        Cesium.SceneTransforms[o(708)](this[o(245)][o(696)], this[o(2392)]._cartesian3, h),
          (e.dom[o(1679)][o(2251)] = o(2076)),
          (e.dom.style[o(1709)] = h.x + r - 5 + 'px'),
          (e[o(2132)][o(1679)][o(151)] = l - h.y + a - 5 + 'px')
      }
      _removeSvgEditAnchors() {
        const e = t
        this[e(1848)] &&
          this[e(1848)].forEach((t) => {
            const i = e
            this[i(245)][i(436)].container[i(2269)](t[i(2132)])
          }),
          (this[e(1848)] = void 0)
      }
      [t(2016)]() {
        const e = t
        this._viewer[e(696)][e(1271)].removeEventListener(this._postRender, this), this[e(1425)]()
      }
    }
    let Is = class extends gs {
      constructor(t) {
        super(t)
      }
      [t(970)](e) {
        const i = t,
          s = this.editGraphic[i(1679)]
        let n = this[i(2392)][i(1131)][i(1902)]()
        if (((this[i(2109)] = n), e && s[i(2080)])) {
          const t = Cesium.Cartographic[i(2579)](n),
            e = this[i(245)][i(696)].sampleHeight(t, this[i(245)][i(118)][i(736)])
          this[i(2109)] = Cesium.Cartesian3[i(1069)](t[i(2106)], t[i(199)], e)
        }
        let o = this[i(2109)][i(1902)](),
          r = Cesium[i(2058)][i(1648)](o),
          a = Cesium[i(2066)].multiplyByPoint(
            r,
            new Cesium[i(310)](0, s[i(1981)], 0),
            new Cesium[i(310)]()
          )
        if (((this[i(1094)] = [a]), e && s.clampToGround)) {
          o = Cesium[i(2285)].fromCartesian(a)
          const t = this[i(245)][i(696)].sampleHeight(o, this[i(245)][i(118)][i(736)])
          ;(o = Cesium[i(310)][i(1069)](o.longitude, o[i(199)], t)), (this[i(1094)] = [o])
        }
        this._moveHeightPositions = []
        let h = this[i(2109)][i(1902)]()
        r = Cesium[i(2058)][i(1648)](h)
        let l = Cesium.Matrix4[i(2483)](
          r,
          new Cesium.Cartesian3(s[i(1981)], 0, 0),
          new Cesium[i(310)]()
        )
        s[i(2080)] ||
          0 == this.editGraphic.style[i(398)] ||
          (((h = Cesium[i(2285)][i(2579)](l))[i(2306)] += s[i(398)]),
          (l = Cesium[i(310)][i(1069)](h[i(2106)], h[i(199)], h[i(2306)])),
          (this._moveHeightPositions = [l]))
      }
      handleMouseMoveALT(e, i) {
        const s = t
        if (this[s(2213)] != fs) return
        if (this[s(2392)][s(1679)][s(2080)]) return
        let n = this[s(2392)][s(2251)],
          o = this[s(317)]
        o = this[s(1862)](i, o)
        let r = Cesium[s(2285)][s(2579)](o)
        ;(n = [
          Cesium.Math[s(363)](r[s(2106)]),
          Cesium.Math[s(363)](r[s(199)]),
          r[s(2306)] < 0 ? 0 : r[s(2306)]
        ]),
          (this[s(317)] = o),
          this[s(2392)].setPosition(n)
        let a = n[2]
        this[s(272)][s(302)](['高度' + a[s(1268)](3) + ' 米']), this[s(970)]()
      }
      [t(1806)](e, i) {
        const s = t
        let n = Cesium.Cartographic[s(2579)](e)
        const o = [
          Cesium[s(475)][s(363)](n.longitude),
          Cesium[s(475)][s(363)](n[s(199)]),
          n[s(2306)]
        ]
        let r = [this[s(2392)][s(2251)], o]
        this.editGraphic.setPositions(r)
      }
      [t(561)](e) {
        const i = t
        let s
        return (
          this.editGraphic[i(1679)][i(2080)]
            ? (s = this[i(245)][i(696)].pickPosition(e[i(954)])) &&
              (s = J(this[i(245)], e[i(954)], s))
            : (s = J(this[i(245)], e[i(954)], this.editGraphic[i(1131)])),
          s
        )
      }
      [t(612)](e) {
        const i = t
        let s = Cesium[i(2285)].fromCartesian(e)
        ;(e = [Cesium[i(475)][i(363)](s[i(2106)]), Cesium[i(475)][i(363)](s[i(199)]), s[i(2306)]]),
          this[i(2392)][i(1253)](e)
      }
      getTooltip(e) {
        const i = t
        if (e == ds) return ['拖拽改变半径']
        if (e == fs) {
          let t = [i(1266)]
          return this[i(2392)][i(1679)][i(2080)] || t.push(i(1882)), t
        }
        return e == Cs ? ['上下拖拽改变拉伸高度'] : void 0
      }
      [t(772)]() {
        const e = t
        this[e(2153)][e(1602)]((t) => {
          const i = e
          t[i(1682)] == ds && (t[i(1482)] = !1)
        })
      }
      [t(2456)]() {
        const e = t
        if (this[e(1375)].editAnchorType == ds) {
          let t = this[e(2392)][e(1679)][e(1981)]
          this[e(272)][e(302)](['半径' + t + ' 米'])
        }
      }
      [t(1704)](e) {
        const i = t,
          s = K(this[i(245)], e.startPosition),
          n = K(this[i(245)], e[i(954)])
        let o = Z(this[i(245)], this[i(617)], s, n),
          r = Cesium.Matrix4.getTranslation(o, new Cesium[i(310)]())
        ;(this.editGraphic.style[i(398)] += r.z), this[i(2392)][i(1253)](this.editGraphic[i(1912)])
      }
      [t(1972)]() {
        const e = t
        this[e(2392)][e(1920)] = (t) => {
          const i = e
          this[i(1675)](),
            setTimeout(() => {
              const t = i
              this[t(970)](!0), this[t(175)]()
            }, 200)
        }
      }
    }
    class ks extends gs {
      constructor(t) {
        super(t)
      }
      [t(567)](e) {
        const i = t
        let s = this[i(2392)][i(1131)]
        this[i(2109)] = s
        let n = this._moveAllPosition,
          o = Cesium.Transforms[i(1648)](n),
          r = Cesium[i(2066)][i(2483)](
            o,
            new Cesium.Cartesian3(0, this[i(2392)].style[i(1981)], 0),
            new Cesium[i(310)]()
          )
        ;(this[i(1094)] = [r]), (this[i(1441)] = [])
      }
      [t(246)](e, i) {
        const s = t
        if (this[s(2213)] != fs) return
        let n = this[s(2392)][s(2251)],
          o = this[s(317)]
        o = this[s(1862)](i, o)
        let r = Cesium[s(2285)][s(2579)](o)
        ;(n = [
          Cesium[s(475)].toDegrees(r[s(2106)]),
          Cesium[s(475)].toDegrees(r[s(199)]),
          r.height
        ]),
          (this[s(317)] = o),
          this[s(2392)].setPosition(n),
          this[s(567)]()
      }
      [t(1806)](e, i) {
        const s = t
        let n = Cesium.Cartographic[s(2579)](e)
        const o = [
          Cesium[s(475)][s(363)](n.longitude),
          Cesium[s(475)][s(363)](n.latitude),
          n[s(2306)]
        ]
        let r = [this[s(2392)][s(2251)], o]
        this.editGraphic[s(2515)](r)
      }
      pickMousePosition(e) {
        const i = t
        return J(this._viewer, e[i(954)], this.editGraphic[i(1131)])
      }
      [t(612)](e) {
        const i = t
        let s = Cesium[i(2285)][i(2579)](e)
        ;(e = [Cesium[i(475)][i(363)](s.longitude), Cesium[i(475)][i(363)](s[i(199)]), s[i(2306)]]),
          this[i(2392)][i(1253)](e)
      }
      [t(1569)](e) {
        const i = t
        return e == ds ? [i(2053)] : e == fs ? ['拖拽改变位置', i(1882)] : void 0
      }
      [t(2456)]() {
        const e = t
        if (this[e(1375)][e(1682)] == ds) {
          let t = this.editGraphic[e(1679)][e(1981)]
          this._mouseTip[e(302)](['半径' + t + ' 米'])
        }
      }
    }
    class Fs extends gs {
      constructor(t) {
        super(t)
      }
      initAnchorPositions(e) {
        const i = t
        if (this[i(2392)][i(2336)] != Q[i(738)]) {
          let t = this[i(2392)][i(2251)]
          this[i(2109)] = Cesium[i(310)][i(667)](t[0], t[1], t[2])
          let e = this._moveAllPosition,
            s = Cesium[i(2058)][i(1648)](e),
            n = Cesium[i(2066)].multiplyByPoint(
              s,
              new Cesium[i(310)](0, this[i(2392)].style[i(1981)], 0),
              new Cesium[i(310)]()
            )
          ;(this._controlPositions = [n]),
            (n = Cesium[i(2066)][i(2483)](
              s,
              new Cesium[i(310)](0, this[i(2392)].style[i(1981)], this[i(2392)][i(1679)].height),
              new Cesium[i(310)]()
            )),
            (this._moveHeightPositions = [n])
        } else {
          this[i(1441)] = []
          let t = this[i(2392)][i(1679)].height
          this[i(2392)][i(2333)][i(1602)]((e) => {
            const s = i,
              n = Cesium[s(310)][s(667)](e[0], e[1], e[2] + t)
            this[s(1441)][s(2553)](n)
          })
        }
      }
      [t(246)](e, i) {
        const s = t
        if (this[s(2213)] == ds) {
          let t = this[s(1862)](i, this[s(317)])
          return (
            (this[s(317)] = t),
            (this[s(1094)][this._editAnchorIndex] = t),
            this.udpateControlPosition(t, this[s(930)]),
            void this._initAnchorPositions()
          )
        }
        this[s(2213)] != Cs && (this.updateAllPositionByMatrix(i), this._initAnchorPositions(!0))
      }
      [t(1806)](e, i) {
        const s = t
        let n = Cesium[s(2285)].fromCartesian(e)
        const o = [
          Cesium[s(475)][s(363)](n.longitude),
          Cesium[s(475)][s(363)](n[s(199)]),
          n[s(2306)]
        ]
        let r = this[s(2392)][s(2333)]
        if (this[s(2392)]._graphicType != Q[s(738)])
          return (r = [this[s(2392)][s(2251)], o]), void this[s(2392)][s(2515)](r)
        ;(r[i] = o), this[s(2392)][s(2515)](r)
      }
      [t(561)](e) {
        const i = t
        let s
        return (
          this[i(2213)] == fs || (this[i(2213)] == ds && this[i(2392)][i(2336)] != Q[i(738)])
            ? (s = J(this[i(245)], e[i(954)], this[i(317)]))
            : (s = this._viewer[i(696)][i(1125)](e[i(954)])) ||
              this._viewer[i(696)][i(1391)][i(1816)](
                e[i(954)],
                this[i(245)].scene[i(2557)][i(1565)]
              ),
          s
        )
      }
      [t(1704)](e, i) {
        const s = t,
          n = K(this[s(245)], e.startPosition),
          o = K(this[s(245)], e[s(954)])
        let r = Z(this[s(245)], this._matrix4, n, o),
          a = Cesium[s(2066)][s(309)](r, new Cesium[s(310)]())
        ;(this[s(2392)][s(1679)].height += a.z),
          (this[s(2392)][s(1679)][s(2306)] = Number(this[s(2392)].style[s(2306)][s(1268)](2))),
          this[s(2392)][s(2231)]()
      }
      [t(612)](e) {
        const i = t
        ;(e = Cesium[i(2066)][i(2483)](this[i(2293)], e, new Cesium[i(310)]())).z = 0
        let s = Cesium[i(2066)][i(1420)](e, new Cesium[i(2066)]())
        this[i(562)](s)
      }
      [t(562)](e) {
        const i = t
        if (this[i(2392)][i(2336)] != Q.simpleWall) {
          let t = this[i(1862)](e, this._leftDownPosition),
            s = Cesium.Cartographic[i(2579)](t)
          this[i(317)] = t
          const n = [
            Cesium.Math.toDegrees(s[i(2106)]),
            Cesium[i(475)][i(363)](s[i(199)]),
            s[i(2306)]
          ]
          return void this[i(2392)][i(1253)](n)
        }
        let s = [],
          n = []
        this.editGraphic.positions[i(1602)]((t) => {
          const e = i,
            s = Cesium[e(310)][e(667)](t[0], t[1], t[2])
          n[e(2553)](s)
        })
        for (let t = 0; t < n[i(277)]; t++) {
          let o = n[t]
          o = this[i(1862)](e, o)
          let r = Cesium[i(2285)].fromCartesian(o)
          this[i(1094)][t] = o
          const a = [
            Cesium.Math.toDegrees(r.longitude),
            Cesium[i(475)].toDegrees(r[i(199)]),
            r[i(2306)]
          ]
          s.push(a)
        }
        this.editGraphic.setPositions(s)
      }
      [t(1569)](e) {
        const i = t
        return e == ds
          ? [i(1266), i(1882)]
          : e == fs
            ? [i(607), '按住alt拖拽改变整体高度']
            : e == Cs
              ? ['拖拽改变高度']
              : void 0
      }
      [t(1972)]() {
        const e = t
        this.editGraphic[e(1920)] = (t) => {
          this._initAnchorPositions(!0)
        }
      }
    }
    class Rs extends ys {
      constructor(t) {
        super(t)
      }
      [t(561)](e) {
        const i = t
        return J(this[i(245)], e[i(954)], this[i(317)])
      }
    }
    class Ls extends ys {
      constructor(t) {
        super(t)
      }
      [t(561)](e) {
        const i = t
        return J(this[i(245)], e[i(954)], this._leftDownPosition)
      }
    }
    class Os extends gs {
      constructor(t) {
        super(t)
      }
      [t(246)](e, i) {
        const s = t
        if (!this[s(2392)][s(1679)].clampToGround) {
          if (this[s(2213)] == ds) {
            let t = this[s(1862)](i, this[s(317)])
            return (
              (this[s(317)] = t),
              (this[s(1094)][this[s(930)]] = t),
              this[s(1806)](t, this[s(930)]),
              void this[s(970)]()
            )
          }
          this[s(562)](i), this[s(970)]()
        }
      }
      [t(1806)](e, i) {
        const s = t
        let n = Cesium[s(2285)].fromCartesian(e)
        const o = [
          Cesium[s(475)][s(363)](n[s(2106)]),
          Cesium.Math.toDegrees(n.latitude),
          n[s(2306)]
        ]
        let r = this[s(2392)].positions
        ;(r[i] = o), this[s(2392)][s(2515)](r)
      }
      [t(561)](e) {
        const i = t
        return J(this[i(245)], e[i(954)], this[i(317)])
      }
      [t(612)](e) {
        const i = t
        ;(e = Cesium[i(2066)][i(2483)](this[i(2293)], e, new Cesium[i(310)]())).z = 0
        let s = Cesium[i(2066)][i(1420)](e, new Cesium.Matrix4())
        this[i(562)](s)
      }
      updateAllPositionByMatrix(e) {
        const i = t
        let s = [],
          n = []
        this.editGraphic.positions[i(1602)]((t) => {
          const e = i,
            s = Cesium[e(310)][e(667)](t[0], t[1], t[2])
          n[e(2553)](s)
        })
        for (let t = 0; t < n.length; t++) {
          let o = n[t]
          o = this[i(1862)](e, o)
          let r = Cesium[i(2285)][i(2579)](o)
          this[i(1094)][t] = o
          const a = [
            Cesium[i(475)][i(363)](r[i(2106)]),
            Cesium[i(475)][i(363)](r[i(199)]),
            r[i(2306)]
          ]
          s[i(2553)](a)
        }
        this[i(2392)][i(2515)](s)
      }
      getTooltip(e) {
        const i = t
        return e == ds ? [i(2410), '按住alt拖拽改变高度'] : e == fs ? [i(607), , i(1800)] : void 0
      }
    }
    class Bs extends ys {
      constructor(t) {
        super(t)
      }
      [t(561)](e) {
        const i = t
        return J(this[i(245)], e[i(954)], this[i(317)])
      }
    }
    class Vs extends ys {
      constructor(t) {
        super(t)
      }
      [t(561)](e) {
        const i = t
        return J(this[i(245)], e[i(954)], this[i(317)])
      }
    }
    class Ns extends ys {
      constructor(t) {
        super(t)
      }
      pickMousePosition(e) {
        const i = t
        return J(this[i(245)], e.endPosition, this[i(317)])
      }
    }
    class Hs extends ys {
      constructor(t) {
        super(t)
      }
      [t(561)](e) {
        const i = t
        return J(this[i(245)], e[i(954)], this._leftDownPosition)
      }
    }
    class Gs extends ys {
      constructor(t) {
        super(t)
      }
      [t(561)](e) {
        const i = t
        return J(this[i(245)], e[i(954)], this[i(317)])
      }
    }
    class Ws extends ys {
      constructor(t) {
        super(t)
      }
      pickMousePosition(e) {
        const i = t
        return J(this[i(245)], e[i(954)], this[i(317)])
      }
      handleMouseMove() {
        const e = t
        this[e(245)].camera[e(1499)](1e-4)
      }
      [t(1839)]() {
        const e = t
        this._viewer[e(1391)][e(1499)](1e-4)
      }
      [t(772)]() {}
      [t(2222)]() {
        const e = t
        this[e(245)].camera[e(1499)](1e-4)
      }
    }
    class Us extends gs {
      constructor(t) {
        super(t)
      }
      [t(772)]() {}
      [t(2222)]() {}
      [t(246)](e, i) {
        const s = t
        if (
          this[s(2213)] == ds &&
          this[s(2392)].graphicType == Q[s(1213)] &&
          null == this[s(2392)].style[s(2306)]
        ) {
          let t = this[s(1862)](i, this[s(317)])
          return (this[s(1094)][this[s(930)]] = t), void this[s(1806)](t, this[s(930)])
        }
        if (
          this.editGraphic[s(1158)] == Q.waterPrimitive &&
          null == this[s(2392)][s(1679)][s(2306)]
        )
          this.updateAllPositionByMatrix(i)
        else {
          let t = this[s(1862)](i, this[s(317)])
          this.editGraphic[s(2306)] = Cesium[s(2285)][s(2579)](t)[s(2306)]
        }
        this._initAnchorPositions()
      }
      [t(561)](e) {
        const i = t
        return J(this._viewer, e[i(954)], this[i(317)])
      }
      updateMoveAllPosition(e) {
        const i = t
        ;(e = Cesium[i(2066)][i(2483)](this._inverseMatrix4, e, new Cesium.Cartesian3())).z = 0
        let s = Cesium[i(2066)].fromTranslation(e, new Cesium[i(2066)]())
        this.updateAllPositionByMatrix(s)
      }
      [t(562)](e) {
        const i = t
        let s = [],
          n = []
        this[i(2392)][i(2333)][i(1602)]((t) => {
          const e = i,
            s = Cesium[e(310)][e(667)](t[0], t[1], t[2])
          n.push(s)
        })
        for (let t = 0; t < n[i(277)]; t++) {
          let o = n[t]
          o = this[i(1862)](e, o)
          let r = Cesium[i(2285)][i(2579)](o)
          this[i(1094)][t] = o
          const a = [
            Cesium[i(475)].toDegrees(r[i(2106)]),
            Cesium[i(475)][i(363)](r[i(199)]),
            r[i(2306)]
          ]
          s.push(a)
        }
        this[i(2392)][i(2515)](s)
      }
      getTooltip(e) {
        const i = t
        return e == ds
          ? this[i(2392)][i(1158)] == Q.waterPrimitive && null == this[i(2392)][i(1679)][i(2306)]
            ? [i(1266), '按住alt拖拽改变高度']
            : ['拖拽改变位置']
          : e == fs
            ? ['拖拽整体平移', i(1800)]
            : void 0
      }
      [t(1972)]() {
        const e = t
        this[e(2392)][e(1920)] = (t) => {
          const i = e
          this[i(1675)](),
            setTimeout(() => {
              const t = i
              this[t(970)](!0), this[t(175)]()
            }, 200)
        }
      }
      deactivateHook() {
        const e = t
        this[e(2392)][e(1920)] = null
      }
    }
    class js extends v {
      constructor(e) {
        const i = t
        super(e),
          (this[i(245)] = e[i(190)]._viewer),
          (this._graphicLayer = e[i(190)]),
          (this[i(1252)] = e.graphicOpts),
          (this[i(2026)] = !0),
          (this[i(2031)] = !1),
          (this._positions = []),
          (this[i(1214)] = null),
          this[i(245)][i(436)][i(321)][i(1638)](Cesium[i(1827)][i(421)]),
          (this[i(2490)] = new Cesium.ScreenSpaceEventHandler(this[i(245)][i(696)][i(1493)]))
      }
      [t(2539)](t) {}
      _handleMouseMove(t) {}
      _handleRightClick(t) {}
      _registerEvents() {
        const e = t
        this[e(2490)].setInputAction((t) => {
          this[e(2539)](t)
        }, Cesium[e(1827)][e(1282)]),
          this.handler[e(1068)]((t) => {
            this[e(673)](t)
          }, Cesium[e(1827)].MOUSE_MOVE),
          this.handler[e(1068)]((t) => {
            this[e(146)](t)
          }, Cesium.ScreenSpaceEventType[e(1309)])
      }
      [t(344)]() {
        const e = t
        this.handler.removeInputAction(Cesium[e(1827)][e(1282)]),
          this[e(2490)][e(1638)](Cesium.ScreenSpaceEventType[e(1073)]),
          this[e(2490)][e(1638)](Cesium.ScreenSpaceEventType.RIGHT_CLICK)
      }
      [t(1519)]() {
        const e = t
        this[e(2031)] ||
          (this._registerEvents(),
          (this._activate = !0),
          (this[e(1214)] = this[e(279)][e(2300)]),
          (this[e(279)].selectedEnable = !1),
          (this[e(855)] = this._graphicLayer[e(778)]),
          (this[e(279)][e(778)] = !1),
          this[e(546)]())
      }
      [t(701)]() {
        const e = t
        this[e(344)](),
          (this[e(2031)] = !1),
          (this._positions = []),
          (this._graphicLayer[e(2300)] = this[e(1214)]),
          (this[e(279)][e(778)] = this._hasEdit),
          this[e(2459)]()
      }
      [t(546)]() {}
      [t(2459)]() {}
      [t(1712)]() {
        const e = t
        this[e(701)](), this[e(2490)].destroy()
      }
    }
    class qs extends js {
      constructor(t) {
        super(t)
      }
      [t(2539)](e) {
        const i = t
        let s = this._pickPosition(e[i(2251)])
        if (!s) return
        let o = Cesium[i(2285)][i(2579)](s)
        const r = [
          Cesium.Math.toDegrees(o[i(2106)]),
          Cesium[i(475)][i(363)](o[i(199)]),
          o[i(2306)] > 0 ? o[i(2306)] : 0
        ]
        if (0 == this[i(2469)][i(277)]) this[i(2469)][i(2553)](r), this[i(279)].add(this[i(759)])
        else {
          if (De(r, this[i(2469)][this[i(2469)].length - 1]) < n) return
          this[i(2469)][i(2553)](r), this[i(759)][i(2515)](this[i(2469)])
        }
        this._setTipContent(),
          this[i(759)][i(2368)] == this._positions[i(277)] &&
            ((this[i(759)].editMode = !1),
            (this[i(759)]._drawEnd = !0),
            this[i(279)][i(450)](p.drawEnd, this[i(759)]),
            (this[i(759)] = null),
            this.deactivate(),
            this[i(2026)] || this[i(1519)]())
      }
      [t(1608)](e) {
        const i = t
        let s
        return (
          (s = this[i(245)][i(696)].pickPosition(e)) ||
            (s = this[i(245)][i(696)][i(1391)][i(1816)](e, this[i(245)][i(696)].globe[i(1565)])),
          s
        )
      }
      [t(673)](e) {
        const i = t
        let s = this[i(1608)](e[i(954)])
        if (!s) return
        if ((this[i(272)][i(1253)](e[i(954)]), this[i(960)][i(1267)](s), this[i(2469)][i(277)] < 1))
          return
        let o = Cesium[i(2285)][i(2579)](s)
        const r = [
          Cesium[i(475)][i(363)](o.longitude),
          Cesium[i(475)][i(363)](o[i(199)]),
          o[i(2306)] > 0 ? o[i(2306)] : 0
        ]
        if (De(r, this[i(2469)][this._positions[i(277)] - 1]) < n) return
        const a = this._positions[i(1500)]([r])
        this[i(759)].setPositions(a), this[i(2456)]()
      }
      [t(2456)]() {}
      [t(146)]() {
        const e = t
        return 0 == this[e(2469)][e(277)] ||
          this[e(759)][e(2368)] ||
          (this[e(759)][e(293)] && this[e(2469)][e(277)] < this._graphic[e(293)])
          ? (this[e(279)][e(1896)](this[e(759)]),
            this[e(279)].fire(p[e(289)]),
            (this._graphic = null),
            void this[e(701)]())
          : (this[e(759)][e(2515)](this._positions),
            (this._graphic[e(1804)] = !1),
            (this._graphic._drawEnd = !0),
            this[e(279)][e(450)](p[e(182)], this._graphic),
            (this[e(759)] = null),
            this[e(701)](),
            void (this._singleEnd || this[e(1519)]()))
      }
      [t(546)]() {
        const e = t
        ;(this[e(759)] = this[e(1341)]()),
          (this[e(759)][e(1804)] = !0),
          (this._graphic._drawEnd = !1),
          (this[e(272)] = new $i(this[e(245)])),
          (this[e(960)] = new Ji(this[e(245)])),
          this[e(1841)](),
          this[e(279)][e(450)](p[e(2279)], this[e(759)])
      }
      [t(1341)]() {}
      _initTipContent() {
        const e = t,
          i = this[e(759)][e(695)],
          s = e(1614)
        let n
        ;(n = this._graphic[e(293)]
          ? [
              '绘制图元名称：' + i + e(527) + this[e(759)][e(293)] + '个点',
              '按下鼠标左键确定第1个点的位置',
              s
            ]
          : 1 == this[e(759)][e(2368)]
            ? [e(1405) + i + e(1517) + this[e(759)][e(2368)] + '个点', '按下鼠标左键确定位置', s]
            : [
                e(1405) + i + '，需要' + this[e(759)].fixPointCount + '个点',
                '按下鼠标左键确定第1个点的位置',
                s
              ]),
          this._mouseTip[e(1520)](n)
      }
      _setTipContent() {
        const e = t,
          i = this[e(759)][e(695)]
        let s
        ;(s = this[e(759)][e(293)]
          ? [
              e(1405) + i + e(527) + this[e(759)][e(293)] + '个点',
              '已有' +
                this[e(2469)][e(277)] +
                '个点，按下鼠标左键确定第' +
                (this._positions.length + 1) +
                '个点',
              ,
              this[e(2469)].length < this._graphic[e(293)] ? '按下鼠标右键取消绘制' : e(1938)
            ]
          : [
              e(1405) + i + '，需要' + this._graphic[e(2368)] + '个点',
              '已有' + this._positions[e(277)] + e(114) + (this[e(2469)].length + 1) + '个点',
              e(1614)
            ]),
          this[e(272)].setContent(s)
      }
      _handleDeActivate() {
        const e = t
        this[e(759)] && this[e(279)].remove(this[e(759)]),
          (this[e(759)] = null),
          this._mouseTip && this[e(272)][e(1896)](),
          this[e(960)] && this[e(960)][e(1896)](),
          (this[e(272)] = void 0),
          (this._mousePoint = void 0)
      }
    }
    class Ys extends qs {
      constructor(t) {
        super(t)
      }
      [t(1608)](e) {
        const i = t
        let s
        if (
          !this._graphic[i(1679)][i(2080)] &&
          this[i(759)][i(1679)][i(398)] > 0 &&
          this[i(2469)][i(277)] > 0
        ) {
          const t = this._positions[0],
            n = Cesium[i(310)].fromDegrees(t[0], t[1], t[2])
          s = J(this._viewer, e, n)
        } else
          (s = this[i(245)][i(696)].pickPosition(e)) ||
            (s = this._viewer[i(696)][i(1391)][i(1816)](e, this[i(245)][i(696)].globe[i(1565)]))
        return s
      }
      _createGraphic() {
        return new me(this[t(1252)])
      }
    }
    class Xs extends Ys {
      constructor(t) {
        super(t)
      }
      [t(1341)]() {
        const e = t
        return Ze.create(this[e(1252)])
      }
    }
    class Qs extends qs {
      constructor(t) {
        super(t)
      }
      _createGraphic() {
        const e = t
        return ce[e(1951)](this[e(1252)])
      }
    }
    class Zs extends qs {
      constructor(t) {
        super(t)
      }
      [t(1341)]() {
        const e = t
        return ue[e(1951)](this[e(1252)])
      }
    }
    class Ks extends js {
      constructor(t) {
        super(t)
      }
      _handleLeftClick(e) {
        const i = t
        ;(this[i(759)][i(1804)] = !1),
          (this[i(759)][i(2371)] = !0),
          this[i(279)][i(450)](p[i(182)], this[i(759)]),
          (this[i(759)] = null),
          this.deactivate(),
          this._singleEnd || this[i(1519)]()
      }
      [t(1608)](e) {
        const i = t
        let s = this._viewer.scene[i(1125)](e)
        return (
          s ||
            (s = this[i(245)][i(696)][i(1391)][i(1816)](e, this[i(245)][i(696)][i(2557)][i(1565)])),
          s
        )
      }
      [t(673)](e) {
        const i = t
        let s = this[i(1608)](e[i(954)])
        if (!s) return
        this._mouseTip[i(1253)](e[i(954)]), this._mousePoint[i(1267)](s)
        let n = Cesium[i(2285)][i(2579)](s)
        const o = [
          Cesium[i(475)][i(363)](n[i(2106)]),
          Cesium.Math[i(363)](n[i(199)]),
          n.height > 0 ? n[i(2306)] : 0
        ]
        0 == this[i(2469)][i(277)] &&
          (this[i(2469)][i(2553)](o), this._graphicLayer.add(this._graphic)),
          this[i(759)][i(1253)](o)
      }
      [t(146)]() {
        const e = t
        this._graphicLayer[e(450)](p[e(289)]), this.deactivate()
      }
      [t(546)]() {
        const e = t
        ;(this[e(759)] = this[e(1341)]()),
          (this[e(759)].editMode = !0),
          (this[e(759)]._drawEnd = !1),
          (this[e(272)] = new $i(this[e(245)])),
          this._initTipContent(),
          (this[e(960)] = new Ji(this[e(245)])),
          this._graphicLayer[e(450)](p[e(2279)], this[e(759)])
      }
      [t(1341)]() {}
      [t(1841)]() {
        const e = t
        let i
        ;(i = [
          e(842) + this[e(759)][e(695)] + e(1517) + this[e(759)].fixPointCount + '个点',
          e(1269),
          e(1614)
        ]),
          this[e(272)].setContent(i)
      }
      [t(2459)]() {
        const e = t
        this[e(759)] && this[e(279)][e(1896)](this[e(759)]),
          (this._graphic = null),
          this._mouseTip && this._mouseTip[e(1896)](),
          this[e(960)] && this[e(960)][e(1896)](),
          (this._mouseTip = void 0),
          (this[e(960)] = void 0)
      }
    }
    class Js extends Ks {
      constructor(t) {
        super(t)
      }
      [t(1341)]() {
        const e = t
        return ce.create(this[e(1252)])
      }
    }
    class $s extends Ks {
      constructor(t) {
        super(t)
      }
      _createGraphic() {
        const e = t
        return Ee[e(1951)](this[e(1252)])
      }
    }
    class tn extends Ks {
      constructor(t) {
        super(t)
      }
      [t(1341)]() {
        const e = t
        return Ft[e(1951)](this[e(1252)])
      }
    }
    class en extends qs {
      constructor(t) {
        super(t)
      }
      [t(1608)](e) {
        const i = t
        let s
        if (this[i(2469)][i(277)] > 0) {
          const t = this[i(2469)][0],
            n = Cesium.Cartesian3[i(667)](t[0], t[1], t[2])
          s = J(this._viewer, e, n)
        } else {
          if (
            ((s = this[i(245)][i(696)].pickPosition(e)) ||
              (s = this[i(245)].scene.camera.pickEllipsoid(
                e,
                this[i(245)][i(696)][i(2557)][i(1565)]
              )),
            !s)
          )
            return
          const t = Cesium[i(2285)].fromCartesian(s)
          t[i(2306)] < 0 && (t.height = 0),
            (s = Cesium.Cartesian3[i(1069)](t[i(2106)], t[i(199)], t[i(2306)]))
        }
        return s
      }
      [t(2456)]() {
        const e = t
        let i = this._graphic[e(1679)][e(1981)]
        this[e(272)][e(302)](['当前半径' + i + ' 米'])
      }
      _createGraphic() {
        return new de(this._graphicOpts)
      }
    }
    class sn extends qs {
      constructor(t) {
        super(t)
      }
      [t(1608)](e) {
        const i = t
        let s
        if (this[i(2469)][i(277)] > 0) {
          const t = this._positions[0],
            n = Cesium[i(310)][i(667)](t[0], t[1], t[2])
          s = J(this[i(245)], e, n)
        } else
          (s = this[i(245)][i(696)][i(1125)](e)) ||
            (s = this[i(245)][i(696)][i(1391)].pickEllipsoid(
              e,
              this._viewer[i(696)].globe[i(1565)]
            ))
        return s
      }
      [t(2456)]() {
        const e = t
        let i = this[e(759)][e(1679)].radius
        this[e(272)][e(302)]([e(1284) + i + ' 米'])
      }
      _createGraphic() {
        return new Ce(this._graphicOpts)
      }
    }
    class nn extends Ks {
      constructor(t) {
        super(t)
      }
      [t(1341)]() {
        const e = t
        return (this[e(1252)][e(2251)] = [0, 0, 0]), Se[e(1951)](this[e(1252)])
      }
    }
    class on extends qs {
      constructor(t) {
        super(t)
      }
      [t(1608)](e) {
        const i = t
        let s = this[i(245)][i(696)].pickPosition(e)
        if (
          (s ||
            (s = this[i(245)][i(696)][i(1391)][i(1816)](e, this[i(245)][i(696)][i(2557)][i(1565)])),
          0 == this[i(2469)].length)
        ) {
          const t = Cesium.Cartographic[i(2579)](s)
          ;(t[i(2306)] += 0.2), (s = Cesium.Cartesian3[i(1069)](t[i(2106)], t[i(199)], t[i(2306)]))
        }
        if (this._graphic[i(1699)]) {
          const t = Cesium[i(2285)][i(2579)](s)
          ;(t[i(2306)] = this[i(245)][i(696)][i(2501)](t, [this._graphic[i(1699)]])),
            t[i(2306)] < 0 && (t[i(2306)] = 0),
            (s = Cesium[i(310)][i(1069)](t[i(2106)], t[i(199)], t[i(2306)]))
        }
        return s
      }
      [t(1341)]() {
        const e = t
        return (this[e(1252)][e(2251)] = [0, 0, 0]), be.create(this._graphicOpts)
      }
    }
    class rn extends qs {
      constructor(t) {
        super(t)
      }
      _createGraphic() {
        const e = t
        return ge.create(this[e(1252)])
      }
      [t(1608)](e) {
        const i = t
        let s
        if (this[i(2469)].length > 0 && this._graphic[i(2336)] != Q[i(738)]) {
          const t = this[i(2469)][0],
            n = Cesium[i(310)][i(667)](t[0], t[1], t[2])
          s = J(this[i(245)], e, n)
        } else
          (s = this._viewer[i(696)].pickPosition(e)) ||
            (s = this[i(245)][i(696)][i(1391)][i(1816)](e, this._viewer[i(696)][i(2557)][i(1565)]))
        return s
      }
    }
    class an extends Ks {
      constructor(t) {
        super(t)
      }
      [t(1341)]() {
        const e = t
        return Je[e(1951)](this[e(1252)])
      }
      [t(1608)](e) {
        const i = t
        let s = this[i(245)].scene[i(1125)](e)
        if (
          (s ||
            (s = this._viewer[i(696)][i(1391)][i(1816)](e, this[i(245)].scene[i(2557)][i(1565)])),
          !s)
        )
          return
        let n = Cesium.Cartographic[i(2579)](s)
        return (s = Cesium.Cartesian3[i(1069)](n.longitude, n[i(199)], n[i(2306)] + 20))
      }
    }
    class hn extends Ks {
      constructor(t) {
        super(t)
      }
      [t(1341)]() {
        const e = t
        return (this[e(1252)][e(2251)] = [0, 0, 0]), si[e(1951)](this[e(1252)])
      }
    }
    let ln = {
      readFeature(e) {
        const i = t
        let s = e[i(138)][i(2365)],
          n = e[i(1004)]
        return (n.position = s), this[i(1951)](n)
      },
      create(e) {
        const i = t
        if (e[i(1158)] === Q[i(441)])
          return new (class extends b {
            constructor(e) {
              const i = t
              super(e),
                (this[i(1115)] = V),
                (this[i(2336)] = Q.model),
                (this[i(528)] = '模型图元'),
                (this[i(314)] = i(1865)),
                (this[i(350)] = 1),
                (this[i(1487)] = Cesium[i(1960)](e.uri, '')),
                (this[i(1912)] = Cesium.defaultValue(e.position, [111, 28, 0])),
                (this[i(343)][i(1198)] = Cesium[i(1960)](this[i(343)].heading, 0)),
                (this[i(343)].pitch = Cesium[i(1960)](this[i(343)].pitch, 0)),
                (this._style[i(1143)] = Cesium.defaultValue(this[i(343)][i(1143)], 0)),
                (this[i(343)][i(2341)] = Cesium.defaultValue(this[i(343)][i(2341)], 1)),
                (this._entity = this[i(2547)]()),
                this[i(1253)](this[i(1912)]),
                this[i(2231)]()
            }
            set [t(634)](e) {
              const i = t
              ;(this._uri = e), (this._entity[i(441)].uri = e)
            }
            get [t(634)]() {
              return this[t(1487)]
            }
            [t(415)](e) {
              const i = t
              this[i(987)][i(1482)] = e
            }
            [t(2547)]() {
              const e = t
              return new Cesium.Entity({
                graphicId: this[e(1570)],
                position: this[e(2238)],
                model: { uri: this[e(1487)], ...this[e(343)], silhouetteColor: Cesium.Color.YELLOW }
              })
            }
            [t(1253)](e) {
              const i = t
              ;(this[i(1912)] = e),
                (this._cartesian3 = Cesium[i(310)][i(667)](e[0], e[1], e[2])),
                (this[i(794)] = this._position),
                (this[i(987)][i(2251)] = this[i(2238)]),
                (this[i(987)][i(925)] = this[i(2120)]()),
                (this[i(1024)] = new Cesium[i(1242)](this._cartesian3, 3))
            }
            [t(2231)]() {
              const e = t
              for (const t in this[e(343)])
                if (Object[e(782)][e(1669)](this._style, t)) {
                  const e = this._style[t]
                  this._entity.model[t] = e
                }
              this._entity.orientation = this[e(2120)]()
              let i = this[e(1496)] ? 1 : 0
              this._entity[e(441)][e(351)] = i
            }
            [t(2120)]() {
              const e = t,
                i = Cesium[e(475)][e(1149)](this[e(343)][e(1198)]),
                s = Cesium[e(475)][e(1149)](this[e(343)].pitch),
                n = Cesium[e(475)].toRadians(this._style[e(1143)]),
                o = new Cesium[e(183)](i, s, n)
              return Cesium[e(2058)][e(950)](this[e(2238)], o)
            }
            [t(415)](e) {
              const i = t
              this[i(987)][i(1482)] = e
            }
            [t(1545)](e) {
              const i = t
              ;(this[i(987)][i(2122)] = e.id),
                (this[i(2462)] = e),
                e._viewer.entities.add(this[i(987)])
            }
            [t(2389)](e) {
              const i = t
              e[i(245)][i(118)][i(1896)](this[i(987)])
            }
          })(e)
      }
    }
    class cn extends Ks {
      constructor(t) {
        super(t)
      }
      [t(1341)]() {
        const e = t
        return (this[e(1252)][e(2251)] = [0, 0, 0]), ln[e(1951)](this._graphicOpts)
      }
      _pickPosition(e) {
        const i = t
        let s = this[i(245)][i(696)][i(1125)](e)
        s ||
          (s = this._viewer[i(696)][i(1391)].pickEllipsoid(
            e,
            this._viewer[i(696)][i(2557)][i(1565)]
          ))
        const n = Cesium[i(2285)].fromCartesian(s)
        return (
          (n[i(2306)] = this[i(245)].scene[i(2501)](n, [this[i(759)][i(987)]])),
          n[i(2306)] < 0 && (n[i(2306)] = 0),
          Cesium[i(310)][i(1069)](n[i(2106)], n[i(199)], n[i(2306)])
        )
      }
    }
    class un extends qs {
      constructor(t) {
        super(t)
      }
      _createGraphic() {
        const e = t
        return oi[e(1951)](this[e(1252)])
      }
    }
    class mn extends Ks {
      constructor(t) {
        super(t)
      }
      [t(1341)]() {
        const e = t
        return (this._graphicOpts.position = [0, 0, 0]), ai.create(this[e(1252)])
      }
    }
    class pn {
      constructor(e = {}) {
        const i = t
        ;(this[i(1165)] = '箭头轴图元'),
          (this[i(343)] = Cesium[i(1960)](e[i(1679)], {})),
          (this[i(343)][i(1070)] = Cesium[i(1960)](e[i(1679)][i(1070)], i(1655))),
          (this[i(343)][i(506)] = Cesium.defaultValue(e.style[i(506)], 'rgba(255,255,0,1)')),
          (this[i(343)][i(359)] = Cesium[i(1960)](e[i(1679)].lineLength, 2)),
          (this._style[i(2113)] = Cesium.defaultValue(e[i(1679)][i(2113)], 0.1)),
          (this[i(343)][i(1299)] = Cesium[i(1960)](e.style[i(1299)], 1)),
          (this._style.headRadius = Cesium[i(1960)](e.style[i(1411)], 0.2)),
          (this[i(343)][i(1592)] = Cesium[i(1960)](e[i(1679)][i(1592)], Cesium[i(1980)].Z)),
          (this._position = e[i(2251)]),
          (this[i(1934)] = this[i(438)]()),
          (this._isSelected = !1),
          this[i(2231)]()
      }
      get [t(468)]() {
        return this._isSelected
      }
      set [t(468)](e) {
        const i = t
        ;(this[i(1496)] = e), this[i(2231)]()
      }
      get primitive() {
        return this._primitive
      }
      _setStyle() {
        const e = t,
          i = this._isSelected ? this[e(343)][e(506)] : this[e(343)][e(1070)]
        this[e(1934)][e(1994)].material[e(1285)][e(1070)] = Cesium[e(1154)][e(2008)](i)
      }
      [t(438)]() {
        const e = t
        let i = this[e(2214)](),
          s = this[e(1464)](this[e(1912)])
        const n = {}
        n[e(138)] = i[0]
        const o = {}
        return (
          (o[e(138)] = i[1]),
          new Cesium.Primitive({
            modelMatrix: s,
            geometryInstances: [new Cesium[e(2155)](n), new Cesium[e(2155)](o)],
            appearance: new Cesium[e(985)]({ material: Cesium[e(1637)][e(1082)](e(1154)) }),
            asynchronous: !1
          })
        )
      }
      _createGeometry() {
        const e = t,
          i = {}
        ;(i[e(277)] = this[e(343)].lineLength),
          (i[e(523)] = this[e(343)].lineRadius),
          (i[e(2395)] = this[e(343)][e(2113)])
        const s = {}
        ;(s[e(277)] = this[e(343)][e(1299)]), (s[e(523)] = 0), (s[e(2395)] = this[e(343)][e(1411)])
        const n = Cesium[e(1355)][e(933)](new Cesium[e(1355)](i)),
          o = Cesium[e(1355)].createGeometry(new Cesium[e(1355)](s))
        let r = (this[e(343)].lineLength + this[e(343)].headLength) / 2,
          a = new Cesium[e(310)](0, 0, r)
        return this[e(1755)](o, a), [n, o]
      }
      [t(1464)]() {
        const e = t
        let i = Cesium[e(819)][e(712)](Cesium.Math[e(1149)](90))
        switch (this._style[e(1592)]) {
          case Cesium[e(1980)].X:
            i = Cesium[e(819)][e(1693)](Cesium.Math[e(1149)](90))
            break
          case Cesium[e(1980)].Y:
            i = Cesium[e(819)].fromRotationX(Cesium[e(475)].toRadians(-90))
        }
        let s = Cesium[e(2066)][e(942)](i),
          n = Cesium[e(2058)][e(1648)](this[e(1912)])
        return Cesium.Matrix4.multiply(n, s, n), n
      }
      [t(1253)](e) {
        const i = t
        this[i(1912)] = e
        let s = this[i(1464)]()
        this._primitive[i(223)] = s
      }
      [t(1755)](e, i) {
        const s = t
        for (let t = 0; t < e[s(1548)].position[s(736)][s(277)]; t += 3)
          (e[s(1548)][s(2251)][s(736)][t] += i.x),
            (e.attributes[s(2251)][s(736)][t + 1] += i.y),
            (e[s(1548)][s(2251)][s(736)][t + 2] += i.z)
      }
    }
    let dn = {
      readFeature(e) {
        const i = t
        let s = e.geometry[i(2365)],
          n = e[i(1004)]
        return (n.position = s), this[i(1951)](n)
      },
      create(e) {
        const i = t
        if (e[i(1158)] === Q[i(1369)])
          return new (class extends b {
            constructor(e = {}) {
              const i = t
              super(e),
                (this[i(1912)] = e[i(2251)]),
                (this[i(528)] = i(1728)),
                (this[i(314)] = i(1865)),
                (this[i(1115)] = H),
                (this[i(2336)] = Q.arrowAxis),
                (this[i(343)][i(1070)] = Cesium[i(1960)](e[i(1679)].color, i(1655))),
                (this._style[i(506)] = Cesium[i(1960)](e[i(1679)][i(506)], i(615))),
                (this[i(343)][i(359)] = Cesium.defaultValue(e[i(1679)][i(359)], 2)),
                (this[i(343)][i(2113)] = Cesium.defaultValue(e[i(1679)][i(2113)], 0.1)),
                (this[i(343)][i(1299)] = Cesium[i(1960)](e[i(1679)][i(1299)], 1)),
                (this[i(343)].headRadius = Cesium[i(1960)](e[i(1679)][i(1411)], 0.2)),
                (this[i(343)].axis = Cesium[i(1960)](e[i(1679)][i(1592)], Cesium[i(1980)].Z))
              let s = Cesium[i(1960)](e[i(2251)], [111, 28, 0])
              ;(this[i(2238)] = Cesium.Cartesian3[i(667)](s[0], s[1], s[2])),
                (this[i(794)] = this[i(1912)] = s),
                (this._boundingSphere = new Cesium[i(1242)](this[i(2238)], this._style[i(359)])),
                (this[i(1934)] = this[i(438)]()),
                (this[i(2549)] = !1)
            }
            get [t(2416)]() {
              return this[t(1934)]
            }
            [t(415)](e) {
              const i = t
              this[i(1934)][i(1482)] = e
            }
            [t(2231)]() {
              const e = t
              let i = !1
              this[e(2549)] &&
                (this[e(2389)](this[e(2462)]), (i = !0), (this._primitive = this[e(438)]()))
              const s = this[e(1496)] ? this[e(343)][e(506)] : this._style.color
              ;(this._primitive[e(1994)].material.uniforms[e(1070)] = Cesium[e(1154)][e(2008)](s)),
                (this[e(1024)] = new Cesium[e(1242)](this[e(2238)], this[e(343)].lineLength)),
                i && this[e(1545)](this[e(2462)])
            }
            [t(1253)](e) {
              const i = t
              ;(this[i(1912)] = e),
                (this._cartesian3 = Cesium[i(310)][i(667)](e[0], e[1], e[2])),
                (this[i(794)] = this[i(1912)]),
                this[i(1137)][i(1253)](this[i(2238)]),
                (this._boundingSphere = new Cesium.BoundingSphere(
                  this[i(2238)],
                  this._style[i(359)]
                ))
            }
            [t(438)]() {
              const e = t
              return (
                (this._arrowPrimitive = new pn({ style: this[e(343)], position: this[e(2238)] })),
                this[e(1137)].primitive
              )
            }
            [t(1545)](e) {
              const i = t
              ;(this._primitive.layerId = e.id),
                (this[i(1934)][i(2366)] = this._id),
                (this[i(2462)] = e),
                (this[i(2549)] = !0),
                e[i(245)][i(696)].primitives[i(1861)](this[i(1934)])
            }
            [t(2389)](e) {
              const i = t
              e[i(245)].scene[i(1346)][i(1896)](this._primitive),
                (this[i(2549)] = !1),
                (this._arrowPrimitive = null)
            }
          })(e)
      }
    }
    class fn extends Ks {
      constructor(t) {
        super(t)
      }
      [t(1341)]() {
        const e = t
        return dn[e(1951)](this[e(1252)])
      }
      [t(1608)](e) {
        const i = t
        let s = this[i(245)].scene[i(1125)](e)
        s || (s = this[i(245)][i(696)].camera[i(1816)](e, this[i(245)][i(696)][i(2557)].ellipsoid))
        const n = Cesium[i(2285)][i(2579)](s)
        return (
          (n[i(2306)] = this[i(245)][i(696)][i(2501)](n, [this._graphic[i(1934)]])),
          n[i(2306)] < 0 && (n[i(2306)] = 0),
          Cesium[i(310)].fromRadians(n.longitude, n[i(199)], n[i(2306)])
        )
      }
    }
    class Cn extends Ks {
      constructor(t) {
        super(t)
      }
      [t(1341)]() {
        const e = t
        return li[e(1951)](this[e(1252)])
      }
      [t(1608)](e) {
        const i = t
        let s = this[i(245)].scene.pickPosition(e)
        const n = Cesium.Cartographic[i(2579)](s)
        return (
          n[i(2306)] < 0 && (n[i(2306)] = 0),
          Cesium[i(310)][i(1069)](n[i(2106)], n.latitude, n[i(2306)])
        )
      }
    }
    class vn extends Ks {
      constructor(t) {
        super(t)
      }
      _createGraphic() {
        const e = t
        return gi[e(1951)](this[e(1252)])
      }
      [t(1608)](e) {
        const i = t
        let s = this[i(245)][i(696)].pickPosition(e)
        s || (s = this[i(245)][i(696)].camera[i(1816)](e, this[i(245)][i(696)][i(2557)][i(1565)]))
        const n = Cesium[i(2285)][i(2579)](s)
        return Cesium.Cartesian3[i(1069)](n.longitude, n[i(199)], n[i(2306)])
      }
    }
    class _n extends Ks {
      constructor(t) {
        super(t)
      }
      [t(1341)]() {
        const e = t
        return hi.create(this[e(1252)])
      }
      _pickPosition(e) {
        const i = t
        let s = this[i(245)][i(696)][i(1125)](e)
        s || (s = this[i(245)].scene.camera[i(1816)](e, this[i(245)][i(696)].globe.ellipsoid))
        const n = Cesium[i(2285)].fromCartesian(s)
        return (
          (n[i(2306)] = this[i(245)][i(696)][i(2501)](n, [this._graphic[i(987)]])),
          n[i(2306)] < 0 && (n[i(2306)] = 0),
          Cesium[i(310)][i(1069)](n.longitude, n[i(199)], n.height)
        )
      }
    }
    class gn extends Ks {
      constructor(t) {
        super(t)
      }
      _createGraphic() {
        const e = t
        return Fi[e(1951)](this[e(1252)])
      }
      [t(1608)](e) {
        const i = t
        let s = this[i(245)].scene[i(1125)](e)
        s ||
          (s = this[i(245)][i(696)][i(1391)][i(1816)](e, this._viewer[i(696)][i(2557)].ellipsoid))
        const n = Cesium[i(2285)][i(2579)](s)
        return Cesium[i(310)][i(1069)](n[i(2106)], n[i(199)], n[i(2306)])
      }
    }
    class yn extends qs {
      constructor(t) {
        super(t)
      }
      _pickPosition(e) {
        const i = t
        let s
        if (this._positions[i(277)] > 0) {
          const t = this[i(2469)][0],
            n = Cesium.Cartesian3[i(667)](t[0], t[1], t[2])
          s = J(this[i(245)], e, n)
        } else {
          if (
            ((s = this[i(245)][i(696)][i(1125)](e)) ||
              (s = this[i(245)][i(696)][i(1391)].pickEllipsoid(
                e,
                this[i(245)][i(696)][i(2557)][i(1565)]
              )),
            !s)
          )
            return
          const t = Cesium[i(2285)].fromCartesian(s)
          t.height < 0 && (t.height = 0),
            (s = Cesium[i(310)].fromRadians(t[i(2106)], t[i(199)], t.height))
        }
        return s
      }
      [t(2456)]() {
        const e = t
        let i = this[e(759)][e(1679)][e(1981)]
        this[e(272)][e(302)](['当前半径' + i + ' 米'])
      }
      [t(1341)]() {
        return Zi[t(1951)](this._graphicOpts)
      }
    }
    class wn extends v {
      constructor(e = {}) {
        const i = t
        if ((super(e), (this._viewer = e[i(395)]), !this[i(245)]))
          throw new Cesium[i(2360)](i(1928))
        ;(this._id = Cesium[i(1960)](e.id, Cesium[i(721)]())),
          (this[i(1165)] = Cesium.defaultValue(e[i(1916)], '')),
          (this._show = Cesium.defaultValue(e[i(1482)], !0)),
          (this._remarks = Cesium[i(1960)](e[i(755)], ''))
      }
      get id() {
        return this[t(1570)]
      }
      set [t(1482)](t) {
        ;(this._show = t), this._setVisible(t)
      }
      get [t(1482)]() {
        return this[t(1144)]
      }
      get type() {
        return this._type
      }
      get [t(1916)]() {
        return this._name
      }
      set [t(1916)](e) {
        this[t(1165)] = e
      }
      get [t(755)]() {
        return this[t(1478)]
      }
      set remarks(t) {
        this._remarks = t
      }
      _setVisible(t) {}
      async [t(2322)](e) {
        const i = t
        return fetch(e)
          .then((t) => t[i(2241)]())
          .then((t) => this[i(949)](t))
      }
      async [t(949)](t) {}
      toGeoJson() {}
      async loadFromLocalFile() {
        const e = t
        return _(e(2241))[e(687)]((t) => ((t = JSON[e(201)](t)), this.loadFromGeoJson(t)))
      }
      downloadToLocalFile(t) {
        g(this.toGeoJson(), t)
      }
      [t(1712)]() {}
    }
    const xn = t(444),
      bn = t(1538),
      Sn = t(1759)
    class Pn {
      constructor(e = {}) {
        const i = t
        ;(this._options = e), (this[i(1436)] = i(1070))
      }
      get [t(1640)]() {
        return this[t(1436)]
      }
      get [t(1350)]() {
        return this._options
      }
      set [t(1350)](e) {
        const i = t
        ;(this._options = e), this[i(2378)]()
      }
      [t(2378)]() {}
      clear() {}
      [t(1931)](e) {
        const i = t
        ;(this[i(966)] = e), this[i(2378)]()
      }
      [t(1930)]() {}
    }
    const Mn = {
      create(e, i) {
        const s = t
        switch (e) {
          case s(1070):
          default:
            return new (class extends Pn {
              constructor(e) {
                const i = t
                super(e),
                  (this._type = 'color'),
                  (this[i(568)][i(1070)] = Cesium[i(1960)](e[i(1070)], i(287)))
              }
              [t(2378)]() {
                const e = t
                if (!this._tileset) return
                let i = this._options.color
                const s = {}
                s[e(1556)] = [[e(2430), e(1745) + i + '")']]
                const n = {}
                ;(n[e(1070)] = s), (this[e(966)].style = new Cesium[e(1753)](n))
              }
              clear() {
                const e = t
                this._tileset && (this._tileset[e(1679)] = null)
              }
            })(i)
          case s(762):
            return new (class extends Pn {
              constructor(e) {
                const i = t
                super(e),
                  (this._type = i(762)),
                  (this[i(568)][i(1070)] = Cesium[i(1960)](e.color, i(287))),
                  (this[i(568)][i(732)] = Cesium[i(1960)](e[i(732)], 0)),
                  (this[i(568)][i(1971)] = Cesium.defaultValue(e[i(1971)], 0)),
                  (this[i(568)][i(1312)] = Cesium.defaultValue(e[i(1312)], 0)),
                  (this[i(568)][i(1289)] = Cesium[i(1960)](e.glowEnable, !1)),
                  (this[i(568)][i(1832)] = Cesium[i(1960)](e[i(1832)], 100)),
                  (this[i(568)][i(666)] = Cesium[i(1960)](e.baseHeight, 100)),
                  (this._options[i(2261)] = Cesium.defaultValue(e[i(2261)], 'z'))
              }
              [t(2378)]() {
                const e = t
                if (!this._tileset) return
                let i = this[e(568)][e(1070)]
                i = Cesium[e(1154)][e(2008)](i)
                let s = Cesium[e(1960)](this[e(568)][e(762)], this[e(602)]())
                const n = {}
                ;(n.lightingModel = this[e(568)][e(1971)]),
                  (n[e(1312)] = this[e(568)].mode),
                  (n[e(732)] = this[e(568)][e(732)]),
                  (n[e(388)] = s)
                const o = new Cesium[e(1468)](n)
                this[e(966)][e(1639)] = o
              }
              [t(602)]() {
                const e = t
                let i = this[e(568)][e(1070)]
                return (
                  (i = Cesium.Color[e(2008)](i)),
                  e(2468) +
                    i.red.toFixed(1) +
                    ',' +
                    i[e(1361)].toFixed(1) +
                    ', ' +
                    i[e(2034)][e(1268)](1) +
                    ', ' +
                    i.alpha[e(1268)](1) +
                    e(1237) +
                    this[e(568)].upAxis +
                    e(897) +
                    this[e(568)][e(666)] +
                    e(2002) +
                    this[e(568)].glowEnable +
                    e(1045) +
                    this[e(568)][e(2261)] +
                    e(897) +
                    this._options[e(1832)] +
                    '.0, 0.0, 1.0) - time));\n                color.rgb += color.rgb * (1.0 - diff);\n               \n                material.diffuse =color.rgb;\n              } '
                )
              }
              [t(1555)]() {
                const e = t
                this[e(966)] && (this[e(966)].customShader = null)
              }
            })(i)
          case 'texture':
            return new (class extends Pn {
              constructor(e) {
                const i = t
                super(e),
                  (this[i(1436)] = i(1660)),
                  (this[i(568)][i(259)] = Cesium[i(1960)](e[i(259)], ''))
              }
              [t(2378)]() {
                const e = t
                if (!this[e(966)]) return
                if (!this[e(568)][e(259)]) return void this.clear()
                const i = {}
                i[e(2171)] = this[e(568)][e(259)]
                const s = new Cesium[e(1468)]({
                  lightingModel: Cesium.LightingModel[e(1194)],
                  mode: Cesium[e(327)][e(2411)],
                  lightingModel: Cesium[e(2555)][e(215)],
                  varyings: { v_normalMC: Cesium[e(1187)][e(785)] },
                  uniforms: {
                    u_texture: { value: new Cesium[e(764)](i), type: Cesium.UniformType[e(1853)] }
                  },
                  vertexShaderText: e(1977),
                  fragmentShaderText:
                    '\n            void fragmentMain(FragmentInput fsInput, inout czm_modelMaterial material) {\n              vec3 positionMC = fsInput.attributes.positionMC;\n              if (dot(vec3(0.0, 0.0, 1.0), v_normalMC) > 0.95) {\n                //处理楼顶:统一处理成深色。\n                material.diffuse = vec3(0.079, 0.107, 0.111);\n              } else {\n                //处理四个侧面: 贴一样的图\n                float width = 100.0;\n                float height = 50.0;\n                float textureX = 0.0;\n                float dotXAxis = dot(vec3(0.0, 1.0, 0.0), v_normalMC);\n                if (dotXAxis > 0.52 || dotXAxis < -0.52) {\n                  textureX = mod(positionMC.x, width) / width;\n                } else {\n                  textureX = mod(positionMC.y, width) / width; //positionMC.z\n                }\n                float textureY = mod(positionMC.z, height) / height; //positionMC.y\n                material.diffuse = texture(u_texture, vec2(textureX, textureY)).rgb;\n              }\n            }'
                })
                this[e(966)].customShader = s
              }
              [t(1555)]() {
                this[t(966)] && (this._tileset.customShader = null)
              }
            })(i)
          case s(1222):
            return new (class extends Pn {
              constructor(e) {
                const i = t
                super(e),
                  (this[i(1436)] = i(1222)),
                  (this[i(568)][i(1556)] = Cesium.defaultValue(e[i(1556)], [
                    [i(2430), "color('#FFFFFF')"]
                  ]))
              }
              [t(2378)]() {
                const e = t
                this[e(966)] &&
                  (this._tileset.style = new Cesium[e(1753)]({
                    color: { conditions: this[e(568)][e(1556)] }
                  }))
              }
              [t(1555)]() {
                const e = t
                this._tileset && (this[e(966)][e(1679)] = null)
              }
            })(i)
        }
      }
    }
    const An = {
        create(e) {
          const i = t
          switch (e[i(1640)]) {
            case i(2194):
              return new cs(e[i(1350)])
            case i(2020):
              return new ls(e[i(1350)])
          }
        }
      },
      Tn = {
        GraphicLayer: class extends wn {
          constructor(e = {}) {
            const i = t
            super(e),
              (this._type = xn),
              (this[i(528)] = i(1470)),
              (this[i(1214)] = Cesium[i(1960)](e.selectedEnable, !0)),
              (this[i(855)] = Cesium.defaultValue(e.hasEdit, !1)),
              (this[i(2295)] = null),
              (this._mousePosition = null),
              (this[i(2013)] = 0),
              (this[i(1032)] = new Cesium[i(1561)](this[i(245)][i(696)][i(1493)])),
              (this._graphics = {})
            let s = this[i(245)][i(903)] || {}
            ;(this[i(903)] = s), (this[i(245)][i(903)] = s)
            let n = this[i(245)][i(1483)] || {}
            ;(this[i(1483)] = n),
              (this[i(245)][i(1483)] = n),
              this[i(245)][i(696)].preRender.addEventListener(this[i(1779)], this),
              (this._viewer._drawing = !1),
              this.on(p.drawEnd, (t) => {
                const e = i
                ;(this[e(245)][e(2055)] = !1),
                  this[e(855)] && this._autoEditing && this._handlePickObject(t)
              }),
              this.on(p[i(289)], (t) => {
                const e = i
                this[e(245)][e(2055)] = !1
              }),
              (this[i(1091)] = new (class {
                constructor(e) {
                  const i = t
                  ;(this[i(245)] = e[i(245)]), (this[i(279)] = e), (this[i(457)] = null)
                }
                startDraw(e) {
                  const i = t
                  this[i(457)] = e
                  const s = e[i(1158)]
                  if (!s || !Q[s]) throw new Cesium[i(2360)](i(2532) + s)
                  this[i(2486)](), (this._graphicDraw = this[i(1697)]()), this[i(1941)][i(1519)]()
                }
                [t(1697)]() {
                  const e = t
                  let i
                  switch (this[e(457)].graphicType) {
                    case Q[e(2040)]:
                      const t = {}
                      ;(t[e(190)] = this[e(279)]), (t.graphicOpts = this[e(457)]), (i = new yn(t))
                      break
                    case Q[e(2513)]:
                    case Q[e(1581)]:
                    case Q[e(2580)]:
                    case Q[e(1910)]:
                    case Q[e(482)]:
                    case Q[e(1629)]:
                    case Q[e(1820)]:
                    case Q[e(2065)]:
                    case Q[e(1956)]:
                    case Q[e(1147)]:
                    case Q[e(2175)]:
                    case Q[e(2014)]:
                    case Q[e(1553)]:
                    case Q[e(2498)]:
                    case Q[e(1155)]:
                    case Q[e(1301)]:
                    case Q[e(651)]:
                    case Q.shaderEffetMagicBox:
                    case Q[e(1171)]:
                    case Q[e(797)]:
                    case Q.shaderEffetNeonPoint:
                    case Q[e(243)]:
                    case Q[e(2167)]:
                    case Q.shaderEffetRadarScan:
                    case Q[e(1322)]:
                    case Q[e(2188)]:
                    case Q[e(928)]:
                    case Q[e(1183)]:
                    case Q.shaderEffetVirus:
                    case Q.shaderEffetWavePetals:
                    case Q.shaderEffetGlowPyramid:
                      const s = {}
                      ;(s[e(190)] = this[e(279)]), (s.graphicOpts = this[e(457)]), (i = new gn(s))
                      break
                    case Q[e(1458)]:
                    case Q[e(203)]:
                      const n = {}
                      ;(n[e(190)] = this[e(279)]), (n[e(563)] = this[e(457)]), (i = new vn(n))
                      break
                    case Q.plane:
                      const o = {}
                      ;(o.graphicLayer = this[e(279)]), (o[e(563)] = this[e(457)]), (i = new _n(o))
                      break
                    case Q[e(649)]:
                    case Q[e(1625)]:
                    case Q.rectSensor:
                    case Q[e(427)]:
                    case Q.probeRadar:
                    case Q[e(1650)]:
                      const r = {}
                      ;(r[e(190)] = this[e(279)]), (r[e(563)] = this[e(457)]), (i = new Cn(r))
                      break
                    case Q[e(2087)]:
                      const a = {}
                      ;(a[e(190)] = this._graphicLayer),
                        (a[e(563)] = this._drawGraphicOpts),
                        (i = new Ys(a))
                      break
                    case Q.doubleArrow:
                    case Q[e(941)]:
                    case Q[e(2286)]:
                    case Q[e(1181)]:
                    case Q[e(2568)]:
                    case Q[e(181)]:
                    case Q[e(2197)]:
                    case Q[e(1668)]:
                    case Q[e(1394)]:
                    case Q.tailedFineArrow:
                      const h = {}
                      ;(h[e(190)] = this[e(279)]),
                        (h[e(563)] = this._drawGraphicOpts),
                        (i = new Xs(h))
                      break
                    case Q[e(1107)]:
                    case Q[e(2052)]:
                    case Q[e(1810)]:
                    case Q[e(163)]:
                    case Q[e(2352)]:
                    case Q[e(835)]:
                    case Q[e(1944)]:
                    case Q[e(665)]:
                    case Q.uprightLabel:
                    case Q.customTemplate:
                    case Q[e(1193)]:
                    case Q.divIndicator:
                      const l = {}
                      ;(l[e(190)] = this[e(279)]),
                        (l[e(563)] = this._drawGraphicOpts),
                        (i = new $s(l))
                      break
                    case Q.simpleMarker:
                    case Q[e(784)]:
                    case Q.alertMarker:
                    case Q[e(143)]:
                    case Q[e(1582)]:
                    case Q[e(1554)]:
                    case Q.textMarker:
                      const c = {}
                      ;(c.graphicLayer = this[e(279)]),
                        (c.graphicOpts = this._drawGraphicOpts),
                        (i = new tn(c))
                      break
                    case Q.simpleLine:
                    case Q[e(975)]:
                      const u = {}
                      ;(u[e(190)] = this[e(279)]), (u[e(563)] = this[e(457)]), (i = new Qs(u))
                      break
                    case Q.uprightLine:
                      const m = {}
                      ;(m[e(190)] = this[e(279)]), (m[e(563)] = this[e(457)]), (i = new Js(m))
                      break
                    case Q[e(1121)]:
                      const p = {}
                      ;(p.graphicLayer = this[e(279)]), (p[e(563)] = this[e(457)]), (i = new Zs(p))
                      break
                    case Q.svg:
                      const d = {}
                      ;(d[e(190)] = this._graphicLayer),
                        (d[e(563)] = this._drawGraphicOpts),
                        (i = new nn(d))
                      break
                    case Q[e(2098)]:
                    case Q[e(608)]:
                      const f = {}
                      ;(f[e(190)] = this[e(279)]), (f[e(563)] = this[e(457)]), (i = new on(f))
                      break
                    case Q[e(2584)]:
                      const C = {}
                      ;(C[e(190)] = this[e(279)]), (C[e(563)] = this[e(457)]), (i = new en(C))
                      break
                    case Q.sphere:
                      const v = {}
                      ;(v[e(190)] = this._graphicLayer), (v[e(563)] = this[e(457)]), (i = new sn(v))
                      break
                    case Q[e(738)]:
                    case Q[e(1921)]:
                    case Q[e(549)]:
                      const _ = {}
                      ;(_[e(190)] = this[e(279)]),
                        (_[e(563)] = this._drawGraphicOpts),
                        (i = new rn(_))
                      break
                    case Q[e(1421)]:
                    case Q[e(1989)]:
                      const g = {}
                      ;(g[e(190)] = this._graphicLayer), (g[e(563)] = this[e(457)]), (i = new an(g))
                      break
                    case Q[e(2316)]:
                      const y = {}
                      ;(y[e(190)] = this[e(279)]),
                        (y.graphicOpts = this._drawGraphicOpts),
                        (i = new hn(y))
                      break
                    case Q.model:
                      const w = {}
                      ;(w[e(190)] = this[e(279)]), (w[e(563)] = this[e(457)]), (i = new cn(w))
                      break
                    case Q[e(265)]:
                      const x = {}
                      ;(x[e(190)] = this[e(279)]), (x[e(563)] = this[e(457)]), (i = new un(x))
                      break
                    case Q[e(1369)]:
                    case Q[e(2500)]:
                      const b = {}
                      ;(b[e(190)] = this[e(279)]), (b.graphicOpts = this[e(457)]), (i = new fn(b))
                      break
                    case Q[e(1946)]:
                    case Q[e(726)]:
                      const S = {}
                      ;(S[e(190)] = this._graphicLayer), (S[e(563)] = this[e(457)]), (i = new mn(S))
                  }
                  return i
                }
                [t(2486)]() {
                  const e = t
                  this[e(1941)] && (this._graphicDraw[e(2039)](), (this[e(1941)] = null))
                }
              })(this)),
              (this[i(372)] = new (class {
                constructor(e) {
                  const i = t
                  ;(this[i(245)] = e._viewer),
                    (this[i(279)] = e),
                    (this[i(2068)] = null),
                    (this[i(372)] = null),
                    (this._editAnchorType = null),
                    this[i(279)].on(p[i(308)], (t) => {
                      const e = i
                      ;(this[e(2068)] = t),
                        (this[e(2068)][e(1804)] = !0),
                        (this[e(372)] = this[e(1408)]())
                    }),
                    this[i(279)].on(p[i(957)], (t) => {
                      const e = i
                      ;(this[e(2068)][e(1804)] = !1), this[e(1438)]()
                    })
                }
                _createGraphicEditor() {
                  const e = t,
                    i = this[e(2068)][e(1059)],
                    s = {}
                  s[e(190)] = this._graphicLayer
                  let n,
                    o = s
                  switch (i) {
                    case j:
                      n = new Gs(o)
                      break
                    case q:
                      n = new Ws(o)
                      break
                    case U:
                      n = new Hs(o)
                      break
                    case W:
                      n = new Ns(o)
                      break
                    case I:
                      n = new xs(o)
                      break
                    case T:
                      n = new ws(o)
                      break
                    case D:
                    case P:
                    case B:
                      n = new ys(o)
                      break
                    case R:
                      n = new bs(o)
                      break
                    case k:
                      n = new Ds(o)
                      break
                    case M:
                      n = new Ss(o)
                      break
                    case A:
                      n = new Ps(o)
                      break
                    case E:
                      n = new Is(o)
                      break
                    case L:
                      n = new ks(o)
                      break
                    case z:
                      n = new Fs(o)
                      break
                    case O:
                      n = new Rs(o)
                      break
                    case V:
                      n = new Ls(o)
                      break
                    case N:
                      n = new Os(o)
                      break
                    case G:
                      n = new Bs(o)
                      break
                    case H:
                      n = new Vs(o)
                      break
                    case Y:
                      n = new Us(o)
                  }
                  return n
                }
                [t(1438)]() {
                  const e = t
                  this[e(372)] && (this[e(372)][e(2039)](), (this[e(372)] = null))
                }
              })(this)),
              this[i(2342)]()
          }
          get [t(1979)]() {
            const e = t
            return Object.values(this[e(716)])
          }
          get [t(778)]() {
            return this[t(855)]
          }
          set hasEdit(e) {
            const i = t
            ;(this._hasEdit = e), this[i(503)](), e || this[i(147)]()
          }
          get [t(2300)]() {
            return this[t(1214)]
          }
          set [t(2300)](e) {
            const i = t
            ;(this._selectedEnable = e), e || this[i(503)]()
          }
          get [t(1356)]() {
            return this[t(2295)]
          }
          set [t(1356)](t) {
            this._selectedGraphic = t
          }
          get [t(2392)]() {
            return this._editGraphic
          }
          set [t(2392)](e) {
            const i = t
            this[i(855)] &&
              (this._clearEditObject(),
              (this._editGraphic = e),
              this[i(450)](p[i(308)], this[i(2068)]))
          }
          [t(415)](e) {
            const i = t
            Object[i(2576)](this._graphics)[i(1602)]((t) => {
              const s = i
              this._graphics[t][s(1482)] = e
            }),
              e || (this[i(503)](), this._clearEditObject())
          }
          [t(2342)]() {
            const e = t
            this[e(1032)][e(1068)]((t) => {
              const i = e
              if (!this._selectedEnable && !this._hasEdit) return
              if (this._viewer._drawing) return
              let s = this[i(245)][i(696)][i(1646)](t[i(2251)])
              if (!s) return void this[i(1965)](null)
              if (!s.id && s[i(2416)] && s.collection) {
                let t = s.collection[i(2129)]
                if (!t) return
                let e = Object[i(736)](this._particleGraphics)
                for (let s = 0; s < e[i(277)]; s++) {
                  const n = e[s]
                  n._layer.id == this[i(1570)] &&
                    t == n[i(2496)]._billboardCollection[i(2129)] &&
                    this[i(1965)](n)
                }
                return
              }
              if (s[i(2416)] && s[i(2416)][i(2122)] == this[i(1570)]) {
                let t = s[i(2416)].graphicId,
                  e = this[i(127)](t)
                return void this[i(1965)](e)
              }
              let n = s.id
              if (n)
                if (n[i(2122)] != this.id) this[i(1965)](null)
                else {
                  const t = this[i(127)](n.graphicId)
                  this[i(1965)](t)
                }
              else this[i(1965)](null)
            }, Cesium[e(1827)].LEFT_CLICK)
          }
          [t(503)]() {
            const e = t
            this._selectedGraphic &&
              ((this[e(2295)][e(468)] = !1),
              this[e(450)](p[e(786)], this[e(2295)]),
              (this._selectedGraphic = null))
          }
          [t(147)]() {
            const e = t
            this._editGraphic &&
              (this[e(450)](p[e(957)], this._editGraphic),
              (this[e(2068)][e(468)] = !1),
              (this[e(2068)] = null))
          }
          [t(1965)](e) {
            const i = t
            if (this[i(855)]) {
              if (!e) return void this[i(147)]()
              if (this[i(2068)] && this[i(2068)].id == e.id) return
              return this[i(147)](), (this[i(2068)] = e), void this.fire(p[i(308)], this[i(2068)])
            }
            e
              ? (this[i(2295)] && this[i(2295)].id == e.id) ||
                (this._clearSelectedObject(),
                (this[i(2295)] = e),
                (this[i(2295)][i(468)] = !0),
                this.fire(p.selectedGraphic, this[i(2295)]))
              : this[i(503)]()
          }
          [t(2224)](e) {
            e ? this[t(2497)](e) : this._removeTooltip()
          }
          _addTooltip(e) {
            const i = t
            ;(this[i(674)] && this[i(674)].id == e.id) ||
              ((this[i(674)] = e),
              this._tooltipContainer && this[i(2247)](),
              (this[i(2477)] = document.createElement('div')),
              (this._tooltipContainer[i(1679)][i(1811)] = i(1913)),
              this[i(245)][i(1392)][i(1621)](this[i(2477)]))
          }
          _removeTooltip() {
            const e = t
            this[e(2477)] &&
              (this[e(2477)][e(1896)](), (this[e(2477)] = null), (this[e(674)] = null))
          }
          [t(1779)]() {
            const e = t
            if (!this[e(1482)]) return
            if (!this[e(1233)]) return void (this[e(1233)] = Date[e(1710)]())
            if (Date[e(1710)]() - this[e(1233)] < 10) return
            this[e(1233)] = Date[e(1710)]()
            let i = Number[e(1664)],
              s = Number.MIN_VALUE
            Object[e(2576)](this[e(903)])[e(1602)]((t) => {
              const n = e,
                o = this[n(903)][t]
              if (o[n(1482)]) {
                const t = this._viewer[n(1391)].position,
                  e = o[n(1131)]
                ;(o[n(779)] = Math[n(2492)](Cesium[n(310)][n(1849)](t, e))),
                  i > o[n(779)] && (i = o[n(779)]),
                  s < o[n(779)] && (s = o._cameraDistance)
              }
            }),
              Object.keys(this[e(903)]).forEach((t) => {
                const i = e,
                  n = this._htmlGraphics[t]
                n[i(1482)] && n[i(1550)](s)
              })
          }
          [t(1861)](e) {
            const i = t
            if (this._graphics[e.id]) throw new Cesium.DeveloperError('id ' + e.id + ' 已存在！')
            ;(this[i(716)][e.id] = e),
              (e.index = ++this[i(2013)]),
              (e[i(1482)] = this[i(1144)]),
              e[i(1545)](this),
              this[i(450)](p[i(1362)], e),
              e._graphicClassType == D &&
                ((this._htmlGraphics[e.id] = e),
                (e.onDomClick = (t) => {
                  const s = i
                  if (!this._selectedEnable) return
                  let n = []
                  Object[s(736)](this[s(903)])[s(1602)]((t) => {
                    const i = s
                    t[i(2462)].id != e[i(2462)].id &&
                      n.indexOf(t._layer.id) < 0 &&
                      (n[i(2553)](t._layer.id), t[i(2462)][i(1965)](null))
                  }),
                    this._handlePickObject(t)
                })),
              e[i(1115)] == G && (this[i(1483)][e.id] = e)
          }
          [t(1896)](e) {
            const i = t
            this[i(2295)] && this[i(2295)].id == e.id && this[i(503)](),
              this._editGraphic && this[i(2068)].id == e.id && this[i(147)](),
              this[i(450)](p[i(2051)], e),
              e._removeHook(this),
              delete this[i(716)][e.id],
              this[i(903)][e.id] && delete this._htmlGraphics[e.id],
              this[i(1483)][e.id] && delete this[i(1483)][e.id]
          }
          [t(474)](e) {
            const i = t,
              s = this[i(716)][e]
            if (!s) throw new Cesium[i(2360)](i(132))
            this[i(1896)](s)
          }
          [t(1876)]() {
            const e = t
            this[e(503)](),
              this[e(147)](),
              Object[e(2576)](this._graphics).forEach((t) => {
                const i = e,
                  s = this[i(716)][t]
                this[i(1896)](s)
              })
          }
          [t(127)](e) {
            return this[t(716)][e]
          }
          [t(2138)](e, i) {
            const s = t
            if (e) {
              let t = e[s(1611)]
              return void this[s(245)][s(1391)][s(1644)](t, i)
            }
            let n = null
            for (let t = 0; t < this[s(1979)][s(277)]; t++) {
              let e = this[s(1979)][t][s(1611)]
              n ? Cesium[s(1242)][s(1974)](e, n, n) : (n = e[s(1902)]())
            }
            n && this[s(245)][s(1391)][s(1644)](n, i)
          }
          [t(2124)](e, i) {
            const s = t,
              n = this[s(716)][e]
            if (!n) throw new Cesium[s(2360)](s(132))
            this[s(2138)](n, i)
          }
          loadFromGeoJson(e) {
            const i = t
            return (
              this[i(1876)](),
              new Promise((t, s) => {
                const n = i
                if (!e) return void s(n(2376))
                let o = e[n(1004)]
                if (o && o[n(1640)] == xn)
                  if (
                    ((this._id = o.id),
                    (this._name = o[n(1916)]),
                    (this[n(1478)] = o[n(755)]),
                    (this[n(1144)] = o[n(1482)]),
                    (this[n(855)] = o[n(778)]),
                    (this[n(1214)] = o[n(2300)]),
                    'FeatureCollection' == e[n(1640)] && e.features)
                  ) {
                    const i = Ki.readFeatures(e[n(2211)])
                    i[n(1602)]((t) => {
                      this[n(1861)](t)
                    }),
                      t(i)
                  } else s('数据格式有问题！')
                else s('导入失败！类型错误，请确保是通过toGeoJson（）方法导出的数据')
              })
            )
          }
          async [t(702)](e) {
            const i = t
            return fetch(e)
              [i(687)]((t) => t[i(2241)]())
              [i(687)]((t) => this[i(1349)](t))
          }
          async [t(1349)](t) {
            return new Promise((e, i) => {
              const s = a0_0x3b79
              if (!t) return void i('失败，数据为空！')
              let n = t.properties
              if (n && n[s(1640)] == xn)
                if (s(2449) == t[s(1640)] && t[s(2211)]) {
                  const i = Ki[s(1256)](t.features)
                  i.forEach((t) => {
                    const e = s
                    this.getById(t.id) && (t[e(1570)] = Cesium[e(721)]()), this[e(1861)](t)
                  }),
                    e(i)
                } else i(s(161))
              else i(s(2090))
            })
          }
          async appendFromLocalFile() {
            const e = t
            return _(e(2241))[e(687)]((t) => ((t = JSON[e(201)](t)), this[e(1349)](t)))
          }
          [t(1930)]() {
            const e = t
            let i = []
            return (
              Object[e(2576)](this[e(716)])[e(1602)]((t) => {
                const s = e,
                  n = this[s(716)][t]
                i[s(2553)](n[s(1930)]())
              }),
              {
                type: e(2449),
                features: i,
                properties: {
                  id: this._id,
                  name: this[e(1165)],
                  remarks: this[e(1478)],
                  type: this._type,
                  show: this[e(1144)],
                  hasEdit: this[e(855)],
                  selectedEnable: this[e(1214)]
                }
              }
            )
          }
          [t(204)](e) {
            const i = t
            this._clearSelectedObject(), this[i(1091)].startDraw(e), (this[i(245)][i(2055)] = !0)
          }
          [t(2486)]() {
            const e = t
            this[e(1091)][e(2486)](), (this[e(245)][e(2055)] = !1)
          }
          [t(1712)]() {
            const e = t
            this[e(245)][e(696)].preRender[e(704)](this[e(1779)], this),
              this._screenSpaceEventHandler.destroy(),
              this[e(2486)](),
              this[e(1876)]()
          }
        },
        TilesetLayer: class extends wn {
          constructor(e = {}) {
            const i = t
            super(e),
              (this[i(1934)] = null),
              (this[i(1436)] = bn),
              (this[i(528)] = i(2439)),
              (this[i(136)] = Cesium[i(1960)](e[i(522)], {})),
              (this[i(1912)] = e.position),
              (this._rotation = Cesium[i(1960)](e[i(2307)], { x: 0, y: 0, z: 0 })),
              (this[i(2035)] = Cesium[i(1960)](e[i(2341)], 1)),
              (this[i(1885)] = e.buildingStyle)
          }
          _setVisible(e) {
            const i = t
            this[i(1934)] && (this[i(1934)][i(1482)] = e)
          }
          get [t(2416)]() {
            return this._primitive
          }
          get [t(522)]() {
            return this[t(136)]
          }
          get [t(2038)]() {
            return this[t(1885)]
          }
          set buildingStyle(e) {
            const i = t
            e
              ? ((this[i(1885)] = e), this[i(1885)][i(1931)](this[i(1934)]))
              : this[i(1885)] && (this._buildingStyle[i(1555)](), (this[i(1885)] = null))
          }
          get url() {
            return this[t(1844)]
          }
          get position() {
            return this[t(1912)]
          }
          get [t(2307)]() {
            return this[t(219)]
          }
          set [t(2307)](e) {
            const i = t
            ;(this[i(219)] = e), this[i(1112)]()
          }
          get [t(2341)]() {
            return this._scale
          }
          set scale(e) {
            const i = t
            ;(this[i(2035)] = e), this[i(1112)]()
          }
          setPosition(e) {
            ;(this[t(1912)] = e), this._setMatrix()
          }
          async [t(2488)](e, i) {
            const s = t
            return (
              (this[s(1844)] = e),
              (this[s(136)] = Cesium[s(938)](this[s(136)], i)),
              Cesium.Cesium3DTileset[s(2488)](e, i)[s(687)]((t) => {
                const e = s
                this[e(1934)] && this[e(245)][e(696)][e(1346)][e(1896)](this[e(1934)]),
                  (this._primitive = this[e(245)][e(696)][e(1346)].add(t)),
                  (this[e(1934)][e(1482)] = this._show),
                  (this[e(2038)] = this._buildingStyle)
                let i = this[e(1934)][e(249)].transform,
                  n = Cesium.Matrix4[e(309)](i, new Cesium.Cartesian3()),
                  o = Cesium[e(2285)][e(2579)](n),
                  r = Cesium[e(475)].toDegrees(o[e(2106)]),
                  a = Cesium[e(475)][e(363)](o.latitude),
                  h = o[e(2306)]
                return (
                  this[e(1912)] ? this[e(1112)]() : (this[e(1912)] = [r, a, h]),
                  Promise[e(1672)](this)
                )
              })
            )
          }
          [t(949)](e) {
            const i = t
            if (!e) return Promise[i(1139)](i(2376))
            let s = e[i(1004)]
            return s && s[i(1640)] == bn
              ? ((this[i(1570)] = s.id),
                (this[i(1165)] = s[i(1916)]),
                (this[i(136)] = s[i(522)]),
                (this._show = s.show),
                (this._url = s[i(2171)]),
                (this[i(219)] = s[i(2307)]),
                (this[i(2035)] = s[i(2341)]),
                (this[i(1912)] = e[i(138)][i(2365)]),
                '' != s[i(238)] && (this[i(1885)] = Mn[i(1951)](s[i(238)], s.styleOpts)),
                this[i(2488)](this[i(1844)], this[i(136)]))
              : Promise[i(1139)](i(2409))
          }
          [t(1930)]() {
            const e = t
            let i = '',
              s = {}
            return (
              this[e(1885)] && ((i = this[e(1885)][e(1640)]), (s = this[e(1885)][e(1350)])),
              {
                type: 'feature',
                geometry: { type: e(1865), coordinates: this[e(2251)] },
                properties: {
                  id: this[e(1570)],
                  type: this[e(1436)],
                  name: this[e(1165)],
                  remarks: this[e(1478)],
                  url: this[e(1844)],
                  tilesetOpts: this[e(522)],
                  show: this[e(1144)],
                  rotation: this[e(2307)],
                  scale: this[e(2341)],
                  styleType: i,
                  styleOpts: s
                }
              }
            )
          }
          async [t(2138)](e) {
            const i = t
            if (this[i(1934)]) return this._viewer[i(2138)](this[i(1934)], e)
          }
          [t(1112)]() {
            const e = t
            let i = Cesium[e(310)].fromDegrees(
                this[e(1912)][0],
                this._position[1],
                this[e(1912)][2]
              ),
              s = Cesium[e(819)][e(2377)](Cesium[e(475)][e(1149)](this._rotation.x)),
              n = Cesium[e(819)][e(1693)](Cesium[e(475)].toRadians(this[e(219)].y)),
              o = Cesium.Matrix3.fromRotationZ(Cesium[e(475)][e(1149)](this[e(219)].z)),
              r = Cesium[e(2066)][e(942)](s),
              a = Cesium[e(2066)].fromRotationTranslation(n),
              h = Cesium[e(2066)].fromRotationTranslation(o),
              l = Cesium[e(2058)][e(1648)](i)
            Cesium.Matrix4[e(999)](l, r, l),
              Cesium[e(2066)].multiply(l, a, l),
              Cesium[e(2066)][e(999)](l, h, l)
            let c = Cesium[e(2066)][e(605)](this[e(2035)])
            Cesium[e(2066)][e(999)](l, c, l), (this[e(1934)][e(249)][e(1863)] = l)
          }
          remove() {
            const e = t
            this[e(1934)] &&
              (this[e(245)][e(696)][e(1346)][e(1896)](this[e(1934)]), (this[e(1934)] = null))
          }
          _beforDestroy() {
            this.remove()
          }
        },
        ImageryLayer: class extends wn {
          constructor(e = {}) {
            const i = t
            super(e),
              (this[i(1436)] = Sn),
              (this[i(528)] = i(1087)),
              (this[i(1796)] = null),
              (this[i(1385)] = e[i(1063)]),
              (this[i(1387)] = Cesium[i(1960)](e.layerOpts, {})),
              (this[i(1387)].show = this[i(1144)]),
              (this[i(745)] = Cesium.defaultValue(e.filter, {})),
              (this[i(745)].enabled = Cesium[i(1960)](this[i(745)][i(299)], !1)),
              (this[i(745)][i(1070)] = Cesium.defaultValue(this._filter[i(1070)], '#000000')),
              (this[i(745)].alpha = Cesium[i(1960)](this._filter[i(2315)], 0.5)),
              (this[i(745)][i(113)] = Cesium.defaultValue(this._filter[i(113)], !0)),
              (this[i(745)][i(458)] = Cesium.defaultValue(this[i(745)][i(458)], !0)),
              this[i(1385)] && this._imageryProvider[i(1259)] && this._createLayer(),
              (this[i(714)] = this[i(745)])
          }
          [t(415)](e) {
            const i = t
            this._imageryLayer && (this[i(1796)].show = e), (this[i(1387)][i(1482)] = e)
          }
          get [t(2510)]() {
            return this._imageryLayer
          }
          get [t(714)]() {
            return this[t(745)]
          }
          set filter(e) {
            const i = t
            ;(this._filter = e),
              this[i(745)][i(299)]
                ? this._layerFilter
                  ? this[i(2025)][i(720)](this[i(745)])
                  : (this._layerFilter = new (class {
                      constructor(e) {
                        const i = t
                        ;(this[i(245)] = e.viewer),
                          (this._imageryLayer = e[i(2510)]),
                          (this[i(343)] = e.style)
                        let s = this[i(343)]
                        const n = {}
                        ;(n[i(913)] = () => s[i(113)]),
                          (n[i(2470)] = () => s[i(2315)]),
                          (n[i(1195)] = () => s[i(2208)])
                        const o = {}
                        ;(o[i(791)] = i(341)),
                          (o[i(1725)] = i(2534)),
                          (o.uniformMap = n),
                          (this[i(343)][i(2208)] = Cesium[i(1154)][i(2008)](this[i(343)][i(1070)])),
                          (this[i(2384)] = o),
                          (this[i(743)] = this[i(1796)]._createTextureWebGL)
                        let r = this
                        ;(this[i(1796)][i(743)] = function (t, e) {
                          const s = i
                          let n = r[s(743)][s(1505)](this)(t, e)
                          return r[s(1238)](t, n) || n
                        }),
                          this[i(720)](s)
                      }
                      [t(720)](e) {
                        const i = t
                        ;(this[i(343)] = e),
                          (this[i(343)][i(2208)] = Cesium.Color[i(2008)](this._style.color)),
                          this[i(245)].imageryLayers[i(1896)](this[i(1796)], !1),
                          this[i(245)][i(872)].add(this[i(1796)])
                      }
                      [t(1238)](e, i) {
                        const s = t
                        let n = this[s(2384)],
                          o = new Cesium[s(1641)]({
                            context: e,
                            pixelFormat: Cesium[s(2277)][s(810)],
                            pixelDatatype: Cesium[s(547)][s(1382)],
                            sampler: i.sampler,
                            width: i[s(575)],
                            height: i[s(2306)],
                            flipY: i[s(681)],
                            target: Cesium[s(1576)][s(859)],
                            preMultiplyAlpha: this._style[s(458)]
                          }),
                          r = new Cesium[s(401)]({
                            context: e,
                            colorTextures: [o],
                            destroyAttachments: !1
                          })
                        const a = {}
                        a[s(299)] = !1
                        const h = { x: 0, y: 0 }
                        ;(h[s(575)] = i[s(575)]), (h[s(2306)] = i[s(2306)])
                        const l = Cesium[s(515)][s(2407)]({
                          depthTest: a,
                          blending: this[s(343)][s(458)]
                            ? Cesium[s(1870)][s(1465)]
                            : Cesium[s(1870)].ALPHA_BLEND,
                          viewport: h
                        })
                        return (
                          (n[s(1720)][s(710)] = function () {
                            return i
                          }),
                          !Cesium[s(1154)][s(1525)](
                            this[s(343)].cesiumColor,
                            Cesium[s(1154)][s(466)]
                          ) &&
                            e[s(2010)](n.bgFS, {
                              framebuffer: r,
                              renderState: l,
                              uniformMap: n[s(1720)]
                            }).execute(e),
                          e
                            .createViewportQuadCommand(n[s(1725)], {
                              framebuffer: r,
                              renderState: l,
                              uniformMap: n[s(1720)]
                            })
                            [s(154)](e),
                          r[s(2039)](),
                          i.destroy(),
                          o
                        )
                      }
                      [t(1896)]() {
                        const e = t
                        ;(this[e(1796)]._createTextureWebGL = this[e(743)]),
                          this._viewer[e(872)].remove(this[e(1796)], !1),
                          this._viewer.imageryLayers.add(this[e(1796)])
                      }
                    })({
                      viewer: this._viewer,
                      imageryLayer: this.imageryLayer,
                      style: this[i(745)]
                    }))
                : this._layerFilter && (this[i(2025)][i(1896)](), (this[i(2025)] = null))
          }
          [t(336)]() {
            const e = t
            this[e(1387)].index = this._viewer[e(872)].length
            const i = new Cesium[e(1457)](this[e(1385)][e(1259)], this._layerOpts)
            this[e(245)][e(872)][e(1861)](i), (this._imageryLayer = i), this._mergeOpts()
          }
          async [t(949)](e) {
            const i = t
            return (
              this[i(1896)](),
              new Promise((t, s) => {
                const n = i
                if (!e) return void s(n(2376))
                let o = e[n(1004)]
                if (!o || o.type != Sn) return void s(n(2409))
                ;(this._id = o.id),
                  (this[n(1165)] = o.name),
                  (this[n(1144)] = o[n(1482)]),
                  (this._layerOpts = o[n(753)]),
                  (this[n(1387)][n(1482)] = this[n(1144)]),
                  (this[n(745)] = o[n(714)])
                let r = o[n(1063)]
                ;(this[n(1385)] = An[n(1951)](r)),
                  this._createLayer(),
                  (this[n(714)] = this._filter),
                  t(this)
              })
            )
          }
          get layerOpts() {
            return this._layerOpts
          }
          set [t(753)](e) {
            const i = t
            ;(this[i(1387)] = e), this[i(1078)]()
          }
          [t(1078)]() {
            const e = t
            for (const t in this[e(1387)])
              if (Object.hasOwnProperty[e(1669)](this[e(1387)], t)) {
                const i = this[e(1387)][t]
                null != i && null != i && (this._imageryLayer[t] = i)
              }
          }
          [t(1930)]() {
            const e = t,
              i = {}
            return (
              (i[e(1640)] = 'Point'),
              (i.coordinates = null),
              {
                type: e(543),
                geometry: i,
                properties: {
                  id: this[e(1570)],
                  type: this[e(1436)],
                  name: this[e(1165)],
                  show: this[e(1144)],
                  filter: this.filter,
                  layerOpts: this.layerOpts,
                  imageryProvider: this[e(1385)][e(1901)]()
                }
              }
            )
          }
          [t(895)]() {
            const e = t
            this._imageryLayer && this[e(245)][e(872)][e(895)](this[e(1796)])
          }
          [t(1896)]() {
            const e = t
            this[e(1796)] && (this[e(245)][e(872)].remove(this[e(1796)]), (this[e(1796)] = null)),
              (this[e(2025)] = null)
          }
          [t(1712)]() {
            this[t(1896)]()
          }
        },
        ClusterLayer: class extends wn {
          constructor(e = {}) {
            const i = t
            super(e),
              (this[i(528)] = i(583)),
              (this[i(1436)] = 'cluster'),
              (this[i(1077)] = e.labelField),
              (this[i(1766)] = Cesium[i(1960)](e[i(1589)], {})),
              (this._clusterEnable = Cesium.defaultValue(e.clusterEnable, !0)),
              (this[i(509)] = {}),
              (this[i(343)] = Cesium[i(1960)](e.style, {})),
              (this[i(343)].colorMap = Cesium.defaultValue(this[i(343)].colorMap, [])),
              (this[i(343)][i(566)] = Cesium.defaultValue(this[i(343)][i(566)], 1)),
              (this[i(343)][i(2452)] = Cesium[i(1960)](this._style.minimumClusterSize, 2)),
              (this[i(343)][i(454)] = Cesium[i(1960)](this[i(343)].pixelRange, 10)),
              (this[i(2119)] = new Dt(e[i(2240)])),
              (this[i(556)] = new It(e.billboard)),
              this[i(1201)]().then((t) => {
                this[i(2138)]()
              })
          }
          get [t(2116)]() {
            return this[t(2069)]
          }
          set clusterEnable(e) {
            const i = t
            ;(this[i(2069)] = e), this[i(2278)] && (this[i(2278)][i(490)][i(299)] = e)
          }
          get [t(1679)]() {
            return this._style
          }
          set [t(1679)](e) {
            const i = t
            ;(this[i(343)] = e), this[i(2231)]()
          }
          [t(415)](e) {
            const i = t
            this[i(2278)] && (this[i(2278)][i(1482)] = e)
          }
          set [t(1589)](e) {
            const i = t
            if (((this[i(1766)] = e), 'string' == typeof e)) {
              if (e.indexOf(i(756)) < 0) return
              this[i(2485)]()
            }
            if (i(1784) == typeof e) {
              if (i(2449) != e[i(1640)]) return
              this[i(2485)]()
            }
          }
          get [t(1589)]() {
            return this[t(1766)]
          }
          [t(1201)]() {
            const e = t
            let i = this[e(1766)]
            if ('string' == typeof i) {
              if (i[e(877)]('.json') < 0) return
              return this[e(2485)]()
            }
            if ('object' == typeof i) {
              if (e(2449) != i[e(1640)]) return
              return this._addData()
            }
          }
          [t(2485)]() {
            const e = t
            return (
              this[e(1876)](),
              Cesium[e(709)]
                [e(2338)](this[e(1766)])
                .then(
                  (t) => (
                    (t[e(1482)] = this[e(1144)]),
                    this[e(245)][e(1652)].add(t),
                    (this._geoJsonDataSource = t),
                    (t[e(490)][e(299)] = this._clusterEnable),
                    (t[e(490)][e(454)] = this[e(343)][e(454)]),
                    (t.clustering[e(2452)] = this._style.minimumClusterSize),
                    this[e(1176)](),
                    this[e(2231)](),
                    Promise.resolve(t)
                  )
                )
            )
          }
          [t(1176)]() {
            const e = t
            this[e(2278)][e(490)][e(117)][e(1973)]((t, i) => {
              const s = e
              ;(i[s(909)][s(1482)] = !0),
                (i.label[s(1482)] = !1),
                (i[s(909)].id = i[s(2240)].id),
                (i.billboard[s(1990)] = Cesium.VerticalOrigin[s(1240)]),
                (i[s(909)].image = this._getImage(t[s(277)]))
            }, this)
          }
          [t(626)](e) {
            const i = t
            if (this[i(509)][e]) return this._imageCache[e]
            const s = this[i(343)][i(566)],
              n = this[i(907)](e)
            let o
            switch (s) {
              case 1:
                o = this._createImage_1(e, n)
                break
              case 2:
                o = this[i(469)](e, n)
                break
              case 3:
                o = this[i(1904)](e, n)
                break
              case 4:
                o = this[i(590)](e, n)
            }
            return (this[i(509)][e] = o), o
          }
          [t(2516)](e, i) {
            const s = t
            let n = document.createElement(s(1493))
            const o = 12 * (e + '')[s(277)] + 50
            n.width = n.height = o
            let r = n[s(1019)]('2d')
            r[s(1735)](),
              (r[s(2345)] = 0.5),
              (r[s(537)] = i),
              r[s(1855)](o / 2, o / 2, o / 2 - 5, 0, 2 * Math.PI),
              r[s(2327)](),
              r[s(1735)](),
              (r[s(2345)] = 0.8),
              (r[s(537)] = i),
              r.arc(o / 2, o / 2, o / 2 - 10, 0, 2 * Math.PI),
              r[s(2327)](),
              (r.font = s(1480)),
              (r[s(2345)] = 1),
              (r[s(537)] = s(1416))
            const a = o / 2 - (12 * e[s(174)]()[s(277)]) / 2
            return r[s(1104)](e, a, o / 2 + 10), n
          }
          [t(469)](e, i) {
            const s = t
            let n = document[s(1945)](s(1493)),
              o = (n.height = n.width = 12 * (e + '')[s(277)] + 50),
              r = n[s(1019)]('2d')
            r[s(557)](),
              (r[s(537)] = i),
              r.arc(o / 2, o / 2, o / 4, 0, 2 * Math.PI),
              r.fill(),
              r.closePath(),
              (r[s(839)] = 4)
            for (let t = 0; t < 1; t++)
              r[s(1735)](),
                r[s(1855)](o / 2, o / 2, o / 4 + 4, 0, 2 * Math.PI, !1),
                (r.strokeStyle = Cesium[s(1154)].fromCssColorString(i)[s(329)](0.4)[s(1886)]()),
                r[s(189)](),
                r[s(1855)](o / 2, o / 2, o / 4 + 10, 0, 2 * Math.PI, !1),
                (r.strokeStyle = Cesium[s(1154)][s(2008)](i)[s(329)](0.2)[s(1886)]()),
                r[s(189)](),
                r[s(418)]()
            ;(r.font = s(1480)), (r[s(2345)] = 1), (r[s(537)] = s(1416))
            const a = o / 2 - (12 * e[s(174)]()[s(277)]) / 2
            return r[s(1104)](e, a, o / 2 + 8), r[s(555)](), n
          }
          [t(1904)](e, i) {
            const s = t
            let n = document[s(1945)](s(1493)),
              o = (n[s(2306)] = n[s(575)] = 12 * (e + '')[s(277)] + 50),
              r = n.getContext('2d')
            r[s(557)](),
              (r[s(537)] = i),
              r.arc(o / 2, o / 2, o / 4, 0, 2 * Math.PI),
              r[s(2327)](),
              r.closePath()
            let a = 0,
              h = (90 * Math.PI) / 180,
              l = (30 * Math.PI) / 180
            r[s(839)] = 4
            for (let t = 0; t < 3; t++)
              r[s(1735)](),
                r[s(1855)](o / 2, o / 2, o / 4 + 4, a, a + h, !1),
                (r[s(1354)] = Cesium[s(1154)][s(2008)](i)[s(329)](0.4)[s(1886)]()),
                r[s(189)](),
                r[s(1855)](o / 2, o / 2, o / 4 + 10, a, a + h, !1),
                (r.strokeStyle = Cesium[s(1154)][s(2008)](i)[s(329)](0.2).toCssColorString()),
                r[s(189)](),
                (a = a + h + l),
                r[s(418)]()
            ;(r.font = s(1480)), (r[s(2345)] = 1), (r[s(537)] = s(1416))
            const c = o / 2 - (12 * e[s(174)]()[s(277)]) / 2
            return r[s(1104)](e, c, o / 2 + 8), r[s(555)](), n
          }
          [t(590)](e, i) {
            const s = t
            let n = document[s(1945)](s(1493)),
              o = (n.height = n[s(575)] = 12 * (e + '')[s(277)] + 50),
              r = n[s(1019)]('2d')
            r[s(557)](),
              (r[s(537)] = i),
              r[s(1855)](o / 2, o / 2, o / 4, 0, 2 * Math.PI),
              r.fill(),
              r[s(418)]()
            let a = 0,
              h = (80 * Math.PI) / 180,
              l = (40 * Math.PI) / 180
            r[s(839)] = o / 4
            for (let t = 0; t < 3; t++)
              r[s(1735)](),
                r.arc(o / 2, o / 2, o / 4 + 4, a, a + h, !1),
                (r[s(1354)] = Cesium[s(1154)]
                  .fromCssColorString(i)
                  [s(329)](0.4)
                  .toCssColorString()),
                r.stroke(),
                (a = a + h + l),
                r[s(418)]()
            ;(r.font = s(1480)), (r[s(2345)] = 1), (r[s(537)] = s(1416))
            const c = o / 2 - (12 * e[s(174)]().length) / 2
            return r[s(1104)](e, c, o / 2 + 8), r[s(555)](), n
          }
          [t(907)](e) {
            const i = t
            for (let t = 0; t < this[i(343)][i(2130)][i(277)]; t++) {
              const s = this._style.colorMap[t]
              if (e >= s[i(922)]) return s[i(1070)]
            }
          }
          [t(2231)]() {
            const e = t
            this._geoJsonDataSource[e(118)][e(736)][e(1602)]((t) => {
              const i = e
              this[i(556)][i(1567)](t[i(909)]),
                this[i(1077)] &&
                  ((t[i(2240)] = {}),
                  (this[i(2119)].text = t[i(1004)][this[i(1077)]]._value),
                  this[i(2119)][i(1567)](t[i(2240)]))
            })
          }
          [t(2138)](e) {
            const i = t
            if (this[i(2278)]) return this[i(245)][i(2138)](this[i(2278)][i(118)][i(736)], e)
          }
          [t(949)](e) {
            const i = t
            return (
              this[i(1876)](),
              new Promise((t, s) => {
                const n = i
                if (!e) return void s(n(2376))
                let o = e[n(1004)]
                o && 'cluster' == o[n(1640)]
                  ? ((this[n(1570)] = o.id),
                    (this._name = o[n(1916)]),
                    (this._remarks = o[n(755)]),
                    (this[n(1144)] = o[n(1482)]),
                    (this[n(343)] = o[n(1679)]),
                    (this[n(2119)] = new Dt(o[n(2240)])),
                    (this._billboard = new It(o[n(909)])),
                    (this[n(1077)] = o.labelField),
                    (this[n(1766)] = o[n(1589)]),
                    this[n(1201)]()
                      [n(687)]((e) => {
                        t(e)
                      })
                      [n(472)]((t) => {
                        s()
                      }))
                  : s(n(2409))
              })
            )
          }
          [t(1930)]() {
            const e = t,
              i = {}
            ;(i.id = this[e(1570)]),
              (i[e(1916)] = this[e(1165)]),
              (i[e(755)] = this._remarks),
              (i.type = this[e(1436)]),
              (i[e(1482)] = this[e(1144)]),
              (i.clusterEnable = this[e(2069)]),
              (i[e(1679)] = this.style),
              (i[e(2240)] = this._label),
              (i[e(909)] = this[e(556)]),
              (i[e(2466)] = this[e(1077)]),
              (i[e(1589)] = this[e(1589)])
            const s = {}
            return (s[e(1640)] = e(1829)), (s.properties = i), s
          }
          [t(1876)]() {
            const e = t
            ;(this._imageCache = {}),
              this[e(2278)] &&
                (this[e(245)].dataSources.remove(this[e(2278)]), (this[e(2278)] = null))
          }
          [t(1712)]() {
            const e = t
            this.removeAll(), Cesium[e(776)](this[e(437)])
          }
        }
      },
      En = { satellite: Et, Parser: Ki },
      zn = {}
    ;(zn[t(1653)] = t(1653)), (zn[t(845)] = 'custom')
    const Dn = {}
    ;(Dn[t(2166)] = h),
      (Dn.Resource = m),
      (Dn.Event = d),
      (Dn[t(2302)] = p),
      (Dn[t(1531)] = C),
      (Dn[t(356)] = zn),
      (Dn[t(444)] = En),
      (Dn[t(322)] = re),
      (Dn[t(1622)] = Tn),
      (Dn[t(1538)] = {}),
      (Dn[t(362)] = as),
      (Dn[t(1063)] = ms),
      (Dn.fs = x)
    let In = window.xt3d || {},
      kn = Dn
    for (const e in kn)
      if (Object[t(782)].call(kn, e)) {
        const t = kn[e]
        In[e] = t
      }
    return (
      (In[t(1135)] = t(2210)),
      (In[t(1628)] = t(1884)),
      (In[t(158)] = t(2548)),
      (In[t(538)] = t(2495)),
      (In[t(627)] = t(2524)),
      In
    )
  })
