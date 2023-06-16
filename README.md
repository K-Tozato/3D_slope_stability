# 3D_slope_stability (Hovland法を用いた斜面安定解析)

- 3次元極限平衡法に基づく斜面安定解析を行うためのFortranコード
- Hovland法を用いた斜面安定解析
- すべり面に楕円体を仮定した円弧すべり計算


## Hovland法による安全率の算出

<img src="./imgs/LEM_image.png" width="60%">

ある崩壊面（すべり面）を仮定したときに，その土塊を土柱に分割して作用する力を算定．  
滑動力と抵抗力の比として安全率$`F`$を以下のように表現する．

```math
F = \frac{\sum_i \sum_j \{cA + (N_{ij}-u_{ij}A)\tan \phi \}}{\sum_i \sum_j T_{ij}}
```

$i, j$：それぞれ土柱の$x,y$方向の番号  
$T_{ij}$：土柱のせん断力  
$N_{ij}$：すべり面上の垂直力  
$W_{ij}$：土柱重量  
$u_{ij}$：間隙水圧  
$c$：粘着力  
$A$：土柱のすべり面の面積   
$\phi$：内部摩擦角

間隙水圧$`\mathbf{u}_ij`$については，Green-Amptモデルによる浸透解析の結果を使用．
浸透解析のコードについては，[Infiltration_sflow](https://github.com/K-Tozato/infiltration_sflow)のページを参照．

Hovland法では，土柱側面に働く力をゼロとして仮定して安全率を算出する．  
垂直力方向の力のつり合いとすべり体全体でのすべり方向に垂直な方向$\bm{v}$のモーメントのつり合いを考えれば以下のようになる．

```math
( \mathbf{T}_{ij} + \mathbf{N}_{ij} + \mathbf{W}_{ij}) \cdot \mathbf{n}_{ij} = 0  
```

```math
\sum_i \sum_j ( \mathbf{r}_b\times\mathbf{T}_{ij} + \mathbf{r}_b\times\mathbf{N}_{ij} + \mathbf{r}_g\times\mathbf{W}_{ij} )\cdot \mathbf{v} = 0  
```


$\\mathbf{u}$：すべり方向を表現する単位ベクトル  
$\\mathbf{w}$：地表面に垂直で地中方向を向く単位ベクトル  
$\mathbf{v}$：$`\mathbf{v}=\mathbf{u}\times\mathbf{w}`$で求められる方向（単位ベクトル）  
$`\mathbf{T}_{ij}=T_{ij}\mathbf{t}_{ij}, \mathbf{N}_{ij}=N_{ij}\mathbf{n}_{ij}, \mathbf{W}_{ij}=W_{ij}\mathbf{g}`$：せん断力，垂直力，重力のベクトル表現  
$`\mathbf{t}_{ij}, \mathbf{n}_{ij}, \mathbf{g}`$：：それぞれの力の方向を表す単位ベクトル   
$`\mathbf{r}_b, \mathbf{r}_g`$：楕円体の回転中心から土柱底面中心，土柱重心までの位置ベクトル  

上記3式から，安全率は以下の式のように書き換えられる．
```math 
F= \frac{\sum_i\sum_j (\mathbf{t}_{ij}\times\mathbf{r}_b)\cdot\mathbf{v} [cA - \{W_{ij}(\mathbf{g}\cdot\mathbf{n}_{ij}) + u_{ij}A)\}\tan\phi]}{\sum_i\sum_j W_{ij}\{ -(\mathbf{g}\cdot\mathbf{n}_{ij})(\mathbf{r}_b\times\mathbf{n}_{ij})\cdot\mathbf{v}+(\mathbf{r}_g\times\mathbf{g})\cdot\mathbf{v} \}}
```



## 解析コードの説明

### 入力ファイル

- **coordinate.txt**  
  座標データ．左の行からx座標，y座標，z座標
- **num_node.txt**  
  対象領域のx方向，y方向の節点数

上記の2つについては，[Infiltration_sflow](https://github.com/K-Tozato/infiltration_sflow)と同じ．   


- **gwater_case.txt**  
  どの浸透解析結果を読み込むかのファイル．  
  浸透解析結果(water_depth_level_××××)の××××の数字の範囲を指定．  
- **infsupdip.dat**  
  安全率の算出の対象とする斜面角の範囲
- **parameter.dat**  
  斜面安定解析のパラメータ   
  左から(初期の)粘着力，内部摩擦角，単位体積重量，(飽和時の)粘着力，内部摩擦角，単位体積重量    

- **斜面形状に関する入力ファイル**




### 出力ファイル







