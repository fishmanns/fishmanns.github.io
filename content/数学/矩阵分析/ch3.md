# Chapter3 内积空间与正交分解

## 3.1 欧几里德空间与酉空间

### 3.1.1 内积的引入

线性空间上只定义了**加法**与**数乘**两种线性运算。仅由这两种运算出发，我们只能讨论线性组合、线性相关与线性无关、基与坐标等**线性结构**问题，却无法回答下面这些**度量性**问题：

- 一个向量有多长？
- 两个向量之间的夹角是多少？
- 两个向量是否"垂直"？
- 一个向量与一个子空间"贴近"到什么程度？

例如在 $\mathbb{R}^2$ 中，若只知道线性结构，就无法判断向量 $(1,0)$ 与 $(0,1)$ 是否垂直。因此需要引入一种新的运算——**内积**，用数量化的方式同时刻画"长度"与"夹角"。

**定义**: 设 $V$ 是实数域 $\mathbb{R}$ 上的线性空间，若对 $V$ 中任意一对有序向量 $\alpha,\beta$，都对应一个唯一确定的实数 $(\alpha,\beta)$，且满足：

1. **对称性**：$(\alpha,\beta)=(\beta,\alpha)$；
2. **线性性**：$(\alpha+\beta,\gamma)=(\alpha,\gamma)+(\beta,\gamma),\quad (k\alpha,\beta)=k(\alpha,\beta)$；
3. **正定性**：$(\alpha,\alpha)\ge 0$，且 $(\alpha,\alpha)=0\iff \alpha=0$；

则称 $(\alpha,\beta)$ 为 $\alpha$ 与 $\beta$ 的**内积**，称定义了内积的实线性空间 $V$ 为**欧几里德空间**（Euclidean space），简称为**欧氏空间**。

---

**性质**:

由对称性 1 与线性性 2 立刻可得,实内积对**两个变元都是线性的**,即它是**双线性**的:

$$
\begin{align}
(\alpha,\beta+\gamma)&=(\beta+\gamma,\alpha)=(\beta,\alpha)+(\gamma,\alpha)=(\alpha,\beta)+(\alpha,\gamma)\\
(\alpha,k\beta)&=(k\beta,\alpha)=k(\beta,\alpha)=k(\alpha,\beta)\\
(\alpha,0)&=(0,\alpha)=0
\end{align}
$$

推广到有限和的形式:

$$
\left(\sum_{i=1}^{m}k_i\alpha_i,\ \sum_{j=1}^{n}l_j\beta_j\right)=\sum_{i=1}^{m}\sum_{j=1}^{n}k_il_j(\alpha_i,\beta_j)
$$

---

*eg1: $\mathbb{R}^n$ 上的标准内积*

对 $\mathbb{R}^n$ 中任意两个向量

$$
x=(x_1,x_2,\cdots,x_n)^T,\qquad y=(y_1,y_2,\cdots,y_n)^T
$$

令

$$
(x,y)=\sum_{i=1}^{n}x_iy_i=x^Ty
$$

容易验证它满足内积的三条公理,称为 $\mathbb{R}^n$ 上的**标准内积**。以后若无特别说明,$\mathbb{R}^n$ 均指带标准内积的欧氏空间。

---

*eg2: 加权内积*

设 $w_1,w_2,\cdots,w_n$ 为给定的**正数**,令

$$
(x,y)=\sum_{i=1}^{n}w_ix_iy_i
$$

它同样满足内积公理。这说明**同一个线性空间上可以定义不同的内积**,从而得到不同的欧氏空间。

---

*eg3: 连续函数空间*

在闭区间 $[a,b]$ 上全体实连续函数构成的线性空间 $C[a,b]$ 中,令

$$
(f,g)=\int_{a}^{b}f(t)g(t)\,dt
$$

由定积分的性质可得它满足内积公理;但需注意,若 $(f,f)=\int_a^bf^2(t)dt=0$,由连续函数的性质必有 $f\equiv 0$,这正是正定性所要求的。

---

*eg4: 矩阵空间*

在实矩阵空间 $\mathbb{R}^{m\times n}$ 中,令

$$
(A,B)=tr(A^TB)=\sum_{i=1}^{m}\sum_{j=1}^{n}a_{ij}b_{ij}
$$

则 $(\cdot,\cdot)$ 是一个内积,称为矩阵的 **Frobenius 内积**。

---

**注意**:

- 内积是**依赖于数域**的运算:同一集合在不同的数域上可以有不同的内积。
- 内积是一个**数**,不是向量;$(k\alpha,\beta)=k(\alpha,\beta)$ 与 $(k\alpha,k\beta)=k^2(\alpha,\beta)$ 不可混淆。
- 欧氏空间 = **实内积空间**;它是"线性空间 + 内积"这样一个整体结构,而不只是线性空间本身。

### 3.1.2 复向量空间中的内积与酉空间

复线性空间上没有"大小顺序"可言,但复数的**模**与**共轭**可以承担类似实数的角色。为使得 $(\alpha,\alpha)$ 是实数并可用于定义长度,复内积必须引入**共轭**。

**定义**: 设 $V$ 是复数域 $\mathbb{C}$ 上的线性空间,若对 $V$ 中任意一对有序向量 $\alpha,\beta$,都对应一个唯一确定的复数 $(\alpha,\beta)$,且满足：

1. **共轭对称性**（Hermite 性）：$(\alpha,\beta)=\overline{(\beta,\alpha)}$；
2. **对第二变元线性**：
$$
(\alpha,\beta+\gamma)=(\alpha,\beta)+(\alpha,\gamma),\qquad (\alpha,k\beta)=k(\alpha,\beta)
$$
3. **对第一变元共轭线性**：
$$
(\alpha+\beta,\gamma)=(\alpha,\gamma)+(\beta,\gamma),\qquad (k\alpha,\beta)=\overline{k}(\alpha,\beta)
$$
4. **正定性**：$(\alpha,\alpha)\ge 0$,且 $(\alpha,\alpha)=0\iff\alpha=0$；

则称 $(\alpha,\beta)$ 为 $\alpha$ 与 $\beta$ 的**内积**,称 $V$ 为**酉空间**（unitary space）,有时也称为**复欧氏空间**。

---

**注意**:

- 条件 1 说明 $(\alpha,\alpha)=\overline{(\alpha,\alpha)}$,因此 $(\alpha,\alpha)$ 必为**实数**,这样条件 4 中的不等式才有意义。
- 条件 2、3 合起来称为**半双线性**（一个变元线性、另一个变元共轭线性）。共轭线性也常称为**反线性**。
- 本课程统一采用"**对第二变元线性、对第一变元共轭线性**"的约定。此时坐标表示自然地写成

$$
(x,y)=x^HGy
$$

并与 "$A$ 的列向量组的 Gram 矩阵为 $A^HA$" 这一常用形式保持一致。若采用相反的约定(对第一变元线性),只需把下面所有公式中的矩阵换成其共轭转置即可,内容完全等价。

---

**性质**:

复内积由共轭对称性 1 与条件 2、3 可得:

$$
\begin{align}
(\alpha,0)=0,\qquad (0,\alpha)&=0\\
(\alpha,k\beta)&=k(\alpha,\beta)\\
(k\alpha,\beta)&=\overline{k}(\alpha,\beta)\\
(k\alpha,l\beta)&=k\overline{l}(\alpha,\beta)
\end{align}
$$

---

*eg1: $\mathbb{C}^n$ 上的标准内积*

对 $\mathbb{C}^n$ 中任意向量 $x=(x_1,\cdots,x_n)^T,\ y=(y_1,\cdots,y_n)^T$,令

$$
(x,y)=\sum_{i=1}^{n}\overline{x_i}y_i=x^Hy
$$

它满足复内积的四条公理,称为 $\mathbb{C}^n$ 上的**标准内积**。

---

*eg2: 复函数空间*

在复值连续函数空间 $C[a,b]$ 中,令

$$
(f,g)=\int_{a}^{b}\overline{f(t)}g(t)\,dt
$$

则 $(f,g)$ 满足复内积的四条公理,是一个内积。这里**共轭加在第一个函数上**,正与本课程"对第二变元线性"的约定一致。

---

*eg3: 复矩阵空间*

在 $\mathbb{C}^{m\times n}$ 中,令

$$
(A,B)=tr(A^HB)=\sum_{i=1}^{m}\sum_{j=1}^{n}\overline{a_{ij}}b_{ij}
$$

则 $(A,B)$ 是一个内积。

---

**注意**: 实内积与复内积的区别

| 比较项 | 实内积（欧氏空间） | 复内积（酉空间） |
| :--- | :--- | :--- |
| 对称性 | 对称 $(\alpha,\beta)=(\beta,\alpha)$ | 共轭对称 $(\alpha,\beta)=\overline{(\beta,\alpha)}$ |
| 线性性 | 对两个变元都线性（双线性） | 对第二变元线性,对第一变元共轭线性 |
| 数乘公式 | $(k\alpha,k\beta)=k^2(\alpha,\beta)$ | $(k\alpha,k\beta)=|k|^2(\alpha,\beta)$ |
| $(\alpha,\alpha)$ | 非负实数 | 非负实数（由共轭对称性保证） |
| 标准空间 | $\mathbb{R}^n$, $(x,y)=x^Ty$ | $\mathbb{C}^n$, $(x,y)=x^Hy$ |
| 正交 | $(\alpha,\beta)=0$ | $(\alpha,\beta)=0$ |

### 3.1.3 向量长度（范数）

有了内积,就可以在空间中定义"长度"。由于 $(\alpha,\alpha)\ge 0$,下面的定义是合理的。

**定义**: 设 $V$ 是内积空间,$\alpha\in V$,称非负实数

$$
\|\alpha\|=\sqrt{(\alpha,\alpha)}
$$

为向量 $\alpha$ 的**长度**（或**范数**）。若 $\|\alpha\|=1$,则称 $\alpha$ 为**单位向量**。

---

**性质**:

1. **非负性**：$\|\alpha\|\ge 0$,且 $\|\alpha\|=0\iff\alpha=0$；
2. **齐次性**：$\|k\alpha\|=|k|\cdot\|\alpha\|$；
3. **单位化**：若 $\alpha\neq 0$,则 $\dfrac{\alpha}{\|\alpha\|}$ 是单位向量。

---

*eg*

在 $\mathbb{R}^3$ 中,设 $x=(1,2,-1)^T,\ y=(2,0,3)^T$,则

$$
(x,y)=1\cdot2+2\cdot0+(-1)\cdot3=-1
$$

$$
\|x\|=\sqrt{1^2+2^2+(-1)^2}=\sqrt6,\qquad
\|y\|=\sqrt{2^2+0^2+3^2}=\sqrt{13}
$$

其单位化向量分别为

$$
\frac{x}{\|x\|}=\frac{1}{\sqrt6}(1,2,-1)^T,\qquad
\frac{y}{\|y\|}=\frac{1}{\sqrt{13}}(2,0,3)^T
$$

### 3.1.4 向量之间的夹角与正交

**定义**（实内积空间）: 设 $\alpha,\beta$ 是欧氏空间中的两个**非零**向量,称

$$
\theta=\arccos\frac{(\alpha,\beta)}{\|\alpha\|\,\|\beta\|},\qquad \theta\in[0,\pi]
$$

为 $\alpha$ 与 $\beta$ 的**夹角**。

要使该定义有意义,必须保证

$$
-1\le \frac{(\alpha,\beta)}{\|\alpha\|\,\|\beta\|}\le 1
$$

这正是下一小节要证明的 **Cauchy–Schwarz 不等式**。

---

**定义**: 设 $\alpha,\beta$ 是内积空间中的两个向量,若

$$
(\alpha,\beta)=0
$$

则称 $\alpha$ 与 $\beta$ **正交**,记作 $\alpha\perp\beta$。

---

**注意**:

- **零向量与任意向量正交**,但零向量不参与夹角与单位化的讨论。
- 在**复内积空间**中,由于 $(\alpha,\beta)$ 一般是复数,通常**不定义夹角**,但"正交"的定义仍完全有效。
- 由对称性,正交关系具有对称性:$\alpha\perp\beta\iff\beta\perp\alpha$。

### 3.1.5 Cauchy–Schwarz 不等式

**定理**（Cauchy–Schwarz 不等式）:

设 $V$ 为内积空间,则对任意 $\alpha,\beta\in V$ 有

$$
|(\alpha,\beta)|\le\|\alpha\|\,\|\beta\|
$$

等号成立当且仅当 $\alpha$ 与 $\beta$ **线性相关**。

---

**证明**:

若 $\beta=0$,则两边均为 $0$,结论显然成立。

设 $\beta\neq 0$。作**正交分解**

$$
\gamma=\alpha-\frac{(\alpha,\beta)}{(\beta,\beta)}\beta
$$

由内积对第二变元的线性性,

$$
(\gamma,\beta)=(\alpha,\beta)-\frac{(\alpha,\beta)}{(\beta,\beta)}(\beta,\beta)=0
$$

即 $\gamma\perp\beta$。于是由勾股定理(见 3.7)

$$
\|\alpha\|^2=\left\|\gamma+\frac{(\alpha,\beta)}{(\beta,\beta)}\beta\right\|^2
=\|\gamma\|^2+\left\|\frac{(\alpha,\beta)}{(\beta,\beta)}\beta\right\|^2
=\|\gamma\|^2+\frac{|(\alpha,\beta)|^2}{\|\beta\|^2}
\ge \frac{|(\alpha,\beta)|^2}{\|\beta\|^2}
$$

两边同乘 $\|\beta\|^2$ 再开方,即得

$$
|(\alpha,\beta)|\le\|\alpha\|\,\|\beta\|
$$

且等号成立 $\iff \gamma=0\iff \alpha=\dfrac{(\alpha,\beta)}{(\beta,\beta)}\beta\iff\alpha,\beta$ 线性相关。

---

**推论**:

在 $\mathbb{R}^n$ 中取标准内积,Cauchy–Schwarz 不等式即为

$$
\left(\sum_{i=1}^{n}a_ib_i\right)^2\le\left(\sum_{i=1}^{n}a_i^2\right)\left(\sum_{i=1}^{n}b_i^2\right)
$$

在 $C[a,b]$ 中取标准内积,即为

$$
\left(\int_a^bf(t)g(t)\,dt\right)^2\le\int_a^bf^2(t)\,dt\int_a^bg^2(t)\,dt
$$

---

**注意**: Cauchy–Schwarz 不等式是"夹角"定义的合法性依据,也是三角不等式、正交投影误差估计等结论的基础。

### 3.1.6 三角不等式

**定理**（三角不等式）:

设 $V$ 为内积空间,则对任意 $\alpha,\beta\in V$ 有

$$
\|\alpha+\beta\|\le\|\alpha\|+\|\beta\|
$$

---

**证明**:

由内积的半双线性与共轭对称性,

$$
\|\alpha+\beta\|^2=(\alpha+\beta,\alpha+\beta)=\|\alpha\|^2+(\alpha,\beta)+(\beta,\alpha)+\|\beta\|^2
$$

而 $(\alpha,\beta)+(\beta,\alpha)=2\operatorname{Re}(\alpha,\beta)\le 2|(\alpha,\beta)|$,于是

$$
\|\alpha+\beta\|^2\le\|\alpha\|^2+2|(\alpha,\beta)|+\|\beta\|^2
\le\|\alpha\|^2+2\|\alpha\|\,\|\beta\|+\|\beta\|^2
=\big(\|\alpha\|+\|\beta\|\big)^2
$$

两边开方即得结论。

---

**推论**:

$$
\big|\,\|\alpha\|-\|\beta\|\,\big|\le\|\alpha-\beta\|
$$

---

**性质**:

由长度与内积还可以得到下面两个恒等式,它们表明"长度"中已经包含了内积的全部信息。

1. **平行四边形法则**：

$$
\|\alpha+\beta\|^2+\|\alpha-\beta\|^2=2\|\alpha\|^2+2\|\beta\|^2
$$

2. **极化恒等式**（实内积空间）：

$$
(\alpha,\beta)=\frac14\Big(\|\alpha+\beta\|^2-\|\alpha-\beta\|^2\Big)
$$

3. **勾股定理**：若 $\alpha\perp\beta$,则

$$
\|\alpha+\beta\|^2=\|\alpha\|^2+\|\beta\|^2
$$

---

**注意**: 极化恒等式说明:在实内积空间中,**保持范数**与**保持内积**是等价的。这一事实将在 3.6 讨论酉变换时反复使用。

### 3.1.7 例题

*eg1: 判断一个二元函数是否为内积*

在 $\mathbb{C}^2$ 中,令

$$
(x,y)=2\overline{x_1}y_1+\overline{x_1}y_2+\overline{x_2}y_1+3\overline{x_2}y_2
$$

即由矩阵

$$
G=\begin{pmatrix}2&1\\1&3\end{pmatrix}
$$

给出 $(x,y)=x^HGy$。判断它是否为内积。

**解**: 容易验证 1、2、3 三条公理;关键看正定性。对任意 $x=(x_1,x_2)^T\neq0$,

$$
x^HGy=2|x_1|^2+2\operatorname{Re}(\overline{x_1}x_2)+3|x_2|^2
\ge 2|x_1|^2-2|x_1||x_2|+3|x_2|^2
$$

配方得

$$
=2\left(|x_1|-\frac12|x_2|\right)^2+\frac52|x_2|^2\ge 0
$$

且等号成立只能 $|x_2|=0$ 且 $|x_1|=0$,即 $x=0$。所以它是 $\mathbb{C}^2$ 上的一个内积。

---

*eg2: 异于标准内积的正交*

仍取 eg1 中的内积,设 $x=(1,0)^T,\ y=(1,-2)^T$,则

$$
(x,y)=2\cdot1+1\cdot(-2)+0+0=0
$$

即在此内积下 $x\perp y$;而在标准内积下 $(x,y)=x^Hy=1\neq 0$。这说明**正交性与内积的选取有关**。

---

*eg3: 用 Cauchy–Schwarz 不等式证明*

证明:对任意正实数 $a,b$ 有

$$
\sqrt{(a+b)^2}\ \le\ \sqrt{a^2+b^2}\cdot\sqrt{2}
$$

**证明**: 在 $\mathbb{R}^2$ 中取 $x=(a,b)^T,\ y=(1,1)^T$,由 Cauchy–Schwarz 不等式

$$
(a+b)^2=(x,y)^2\le\|x\|^2\|y\|^2=(a^2+b^2)\cdot2
$$

两端开方即得。

---

**方法总结**:

1. **验证内积**:逐条核对对称（共轭对称）性、线性（半双线性）、正定性,其中正定性往往需要配方或利用二次型的判别法。
2. **计算范数、夹角**:先算 $(x,y)$,再算 $\|x\|,\|y\|$,最后代入定义。
3. **证明不等式**:把待证式改写成某个内积空间中已知不等式的形式（Cauchy–Schwarz、三角不等式、勾股定理),即"构造合适的内积"是关键技巧。

## 3.2 Gramian矩阵与内积度量矩阵

### 3.2.1 向量组的 Gram 矩阵

**定义**: 设 $V$ 是内积空间,$\alpha_1,\alpha_2,\cdots,\alpha_k\in V$,称 $k$ 阶方阵

$$
G=
\begin{pmatrix}
(\alpha_1,\alpha_1)&(\alpha_1,\alpha_2)&\cdots&(\alpha_1,\alpha_k)\\
(\alpha_2,\alpha_1)&(\alpha_2,\alpha_2)&\cdots&(\alpha_2,\alpha_k)\\
\vdots&\vdots&\ddots&\vdots\\
(\alpha_k,\alpha_1)&(\alpha_k,\alpha_2)&\cdots&(\alpha_k,\alpha_k)
\end{pmatrix}
$$

为向量组 $\alpha_1,\alpha_2,\cdots,\alpha_k$ 的 **Gram 矩阵**（Gramian matrix）,称其行列式 $\det G$ 为 **Gram 行列式**。

即

$$
G=\big[(\alpha_i,\alpha_j)\big]_{k\times k}
$$

---

*eg*

在 $\mathbb{R}^3$ 中取标准内积,$\alpha_1=(1,1,0)^T,\ \alpha_2=(1,0,1)^T,\ \alpha_3=(0,1,1)^T$,则 $(\alpha_i,\alpha_j)$ 恰为 $\alpha_i,\alpha_j$ 的内积,于是

$$
G=
\begin{pmatrix}
2&1&1\\
1&2&1\\
1&1&2
\end{pmatrix}
$$

### 3.2.2 Gram 矩阵与内积的关系

**定理**: Gram 矩阵必为 **Hermite 矩阵**（实空间为**对称矩阵**）。

---

**证明**:

由共轭对称性,

$$
g_{ji}=(\alpha_j,\alpha_i)=\overline{(\alpha_i,\alpha_j)}=\overline{g_{ij}}
$$

即 $G^H=G$。

---

**定理**（Gram 矩阵与线性组合的内积）:

设 $\alpha_1,\cdots,\alpha_k\in V$,$x=(x_1,\cdots,x_k)^T,\ y=(y_1,\cdots,y_k)^T$ 为任意坐标向量,则

$$
\left(\sum_{i=1}^{k}x_i\alpha_i,\ \sum_{j=1}^{k}y_j\alpha_j\right)=x^HGy
$$

特别地,

$$
\left\|\sum_{i=1}^{k}x_i\alpha_i\right\|^2=x^HGx
$$

---

**证明**:

由内积对第二变元线性、对第一变元共轭线性,

$$
\left(\sum_{i=1}^{k}x_i\alpha_i,\ \sum_{j=1}^{k}y_j\alpha_j\right)
=\sum_{i=1}^{k}\sum_{j=1}^{k}\overline{x_i}y_j(\alpha_i,\alpha_j)
=\sum_{i=1}^{k}\sum_{j=1}^{k}\overline{x_i}g_{ij}y_j
=x^HGy
$$

第二式取 $y=x$ 即得。

---

**定理**（Gram 矩阵的半正定性）:

Gram 矩阵 $G$ 是 **Hermite 半正定矩阵**,即对任意 $x\in\mathbb{C}^k$ 有

$$
x^HGx\ge 0
$$

---

**证明**:

由上一定理,

$$
x^HGx=\left\|\sum_{i=1}^{k}x_i\alpha_i\right\|^2\ge 0
$$

再由 $G^H=G$ 知 $G$ 为 Hermite 半正定矩阵。

---

**定理**（Gram 矩阵正定的判别）:

Gram 矩阵 $G$ **正定**的充分必要条件是 $\alpha_1,\alpha_2,\cdots,\alpha_k$ **线性无关**。

---

**证明**:

**必要性**:若 $\alpha_1,\cdots,\alpha_k$ 线性相关,则存在非零向量 $x=(x_1,\cdots,x_k)^T$ 使

$$
\sum_{i=1}^{k}x_i\alpha_i=0
$$

于是

$$
x^HGx=\left\|\sum_{i=1}^{k}x_i\alpha_i\right\|^2=0
$$

与非零 $x$ 下 $x^HGx>0$ 矛盾。

**充分性**:若 $G$ 不正定,则存在 $x\neq 0$ 使 $x^HGx=0$,从而

$$
\left\|\sum_{i=1}^{k}x_i\alpha_i\right\|^2=0
\ \Longrightarrow\ \sum_{i=1}^{k}x_i\alpha_i=0
$$

与线性无关矛盾。

---

**推论**:

1. $rank(G)=\dim\operatorname{span}\{\alpha_1,\alpha_2,\cdots,\alpha_k\}$,即 Gram 矩阵的秩等于向量组的秩。
2. Gram 行列式 $\det G\ge 0$;且 $\det G>0\iff\alpha_1,\cdots,\alpha_k$ 线性无关。
3. 当 $k=n$ 且 $\alpha_1,\cdots,\alpha_n$ 构成基时,$G$ 必为正定矩阵。

---

**注意**:

- 结论 2 是 Cauchy–Schwarz 不等式的高维推广:当 $k=2$ 时,

$$
\det G=(\alpha_1,\alpha_1)(\alpha_2,\alpha_2)-|(\alpha_1,\alpha_2)|^2\ge 0
$$

即为 $|(\alpha_1,\alpha_2)|\le\|\alpha_1\|\,\|\alpha_2\|$。

- Gram 矩阵只能给出"是否线性相关",不能给出"秩是多少"以外的更多线性信息;但它完整刻画了向量组的内积结构。

### 3.2.3 内积的矩阵表示:度量矩阵

**定义**: 设 $\varepsilon_1,\varepsilon_2,\cdots,\varepsilon_n$ 是 $n$ 维内积空间 $V$ 的一组**基**,向量 $\alpha,\beta$ 在该基下的坐标分别为

$$
x=(x_1,\cdots,x_n)^T,\qquad y=(y_1,\cdots,y_n)^T
$$

即

$$
\alpha=\sum_{i=1}^{n}x_i\varepsilon_i,\qquad \beta=\sum_{j=1}^{n}y_j\varepsilon_j
$$

令

$$
G=\big[(\varepsilon_i,\varepsilon_j)\big]_{n\times n}
$$

则

$$
(\alpha,\beta)=x^HGy
$$

称 $G$ 为基 $\varepsilon_1,\cdots,\varepsilon_n$ 下的**度量矩阵**（也称**内积的矩阵表示**）。它就是基向量组的 Gram 矩阵。

---

**定理**:

$n$ 维内积空间的任意一组基的度量矩阵都是 **Hermite 正定矩阵**;反之,任意给定的 $n$ 阶 Hermite 正定矩阵 $G$,都可以通过 $(\alpha,\beta)=x^HGy$ 在该空间中定义一个内积。

---

**证明**:

由 3.2.2,$G$ 是 Hermite 矩阵,且

$$
x^HGx=\left\|\sum_{i=1}^{n}x_i\varepsilon_i\right\|^2
$$

因为 $\varepsilon_1,\cdots,\varepsilon_n$ 是基,坐标 $x\neq 0$ 当且仅当 $\sum x_i\varepsilon_i\neq 0$,故 $x^HGx>0$ 对一切 $x\neq0$ 成立,即 $G$ 正定。

反之,若 $G$ 为 Hermite 正定矩阵,直接验证 $(\alpha,\beta)=x^HGy$ 满足内积四条公理即可(正定性由 $G$ 的正定性保证),细节留作练习。

---

**定理**（不同基下度量矩阵的关系——共轭合同）:

设 $\varepsilon_1,\cdots,\varepsilon_n$ 与 $\eta_1,\cdots,\eta_n$ 是同一内积空间的两组基,且

$$
(\eta_1,\eta_2,\cdots,\eta_n)=(\varepsilon_1,\varepsilon_2,\cdots,\varepsilon_n)P
$$

即 $P$ 为从基 $\varepsilon$ 到基 $\eta$ 的**过渡矩阵**。若两组基下的度量矩阵分别为 $G_\varepsilon,\ G_\eta$,则

$$
G_\eta=P^HG_\varepsilon P
$$

---

**证明**:

设向量 $\alpha$ 在两组基下的坐标分别为 $x,x'$,由 $\eta=\varepsilon P$ 得

$$
\alpha=\varepsilon x=\eta x'=\varepsilon Px'
\quad\Longrightarrow\quad x=Px'
$$

于是

$$
(\alpha,\beta)=x^HG_\varepsilon y=(Px')^HG_\varepsilon(Py')=x'^H\big(P^HG_\varepsilon P\big)y'
$$

由度量矩阵的唯一性即得 $G_\eta=P^HG_\varepsilon P$。

---

**注意**:

- $G_\eta=P^HG_\varepsilon P$ 称为**共轭合同**（复空间）或**合同**（实空间）。合同关系保持 Hermite 性与（半）正定性,这正是不同基下度量矩阵都正定的原因。
- 过渡矩阵 $P$ 必可逆,因此 $rank(G_\eta)=rank(G_\varepsilon)=n$。
- 特别地,若 $\varepsilon_1,\cdots,\varepsilon_n$ 是**标准正交基**(见 3.4),则

$$
g_{ij}=(\varepsilon_i,\varepsilon_j)=\delta_{ij}
\quad\Longrightarrow\quad G=I
$$

此时内积就是坐标的**标准形式**

$$
(\alpha,\beta)=\sum_{i=1}^{n}\overline{x_i}y_i=x^Hy
$$

这正是"标准正交基使内积形式最简"的含义。

### 3.2.4 实空间与复空间中的对应形式

| 概念 | 实内积空间（欧氏空间） | 复内积空间（酉空间） |
| :--- | :--- | :--- |
| Gram 矩阵元素 | $g_{ij}=(\alpha_i,\alpha_j)$ | $g_{ij}=(\alpha_i,\alpha_j)$ |
| 矩阵性质 | $G^T=G$（对称） | $G^H=G$（Hermite） |
| 内积表示 | $(\alpha,\beta)=x^TGy$ | $(\alpha,\beta)=x^HGy$ |
| 半正定 | $x^TGx\ge 0$ | $x^HGx\ge 0$ |
| 正定判别 | 线性无关 | 线性无关 |
| 换基公式 | $G_\eta=P^TG_\varepsilon P$ | $G_\eta=P^HG_\varepsilon P$ |
| 标准正交基 | $G=I$ | $G=I$ |
| 列向量组 $A$ 的 Gram 矩阵 | $A^TA$ | $A^HA$ |

### 3.2.5 例题

*eg1: 用 Gram 矩阵判断线性相关性*

判断

$$
\alpha_1=(1,1,0)^T,\quad \alpha_2=(1,0,1)^T,\quad \alpha_3=(0,1,1)^T
$$

的线性相关性。

**解**: 其 Gram 矩阵为

$$
G=
\begin{pmatrix}
2&1&1\\
1&2&1\\
1&1&2
\end{pmatrix}
$$

计算各阶顺序主子式:

$$
D_1=2>0,\qquad D_2=\begin{vmatrix}2&1\\1&2\end{vmatrix}=3>0,\qquad
D_3=\det G=4>0
$$

故 $G$ 正定,于是 $\alpha_1,\alpha_2,\alpha_3$ 线性无关。

---

*eg2: 由内积求度量矩阵*

在 $\mathbb{R}^2$ 中取标准内积,基为

$$
\varepsilon_1=(1,0)^T,\qquad \varepsilon_2=(1,2)^T
$$

求度量矩阵,并用它计算 $\alpha=(1,0)^T$ 与 $\beta=(1,2)^T$ 的内积。

**解**:

$$
G=
\begin{pmatrix}
(\varepsilon_1,\varepsilon_1)&(\varepsilon_1,\varepsilon_2)\\
(\varepsilon_2,\varepsilon_1)&(\varepsilon_2,\varepsilon_2)
\end{pmatrix}
=
\begin{pmatrix}
1&1\\
1&5
\end{pmatrix}
$$

$\alpha$ 在基下的坐标为 $x=(1,0)^T$,$\beta$ 的坐标为 $y=(0,1)^T$,于是

$$
(\alpha,\beta)=x^TGy=(1,0)\begin{pmatrix}1&1\\1&5\end{pmatrix}\begin{pmatrix}0\\1\end{pmatrix}=1
$$

与直接计算一致。

---

*eg3: 非标准内积下的度量矩阵*

在 $\mathbb{R}^2$ 中取标准基 $\varepsilon_1=(1,0)^T,\ \varepsilon_2=(0,1)^T$,定义

$$
(x,y)=2x_1y_1+x_1y_2+x_2y_1+3x_2y_2
$$

则

$$
G=\begin{pmatrix}2&1\\1&3\end{pmatrix}
$$

因为 $2>0,\ \det G=5>0$,所以 $G$ 正定,它确实是一个内积;并且在此内积下,标准基不再是标准正交基。

---

**方法总结**:

1. **求 Gram 矩阵**:逐个计算所有内积 $(\alpha_i,\alpha_j)$,写成矩阵;注意 $G^H=G$,只需算上三角。
2. **判别线性相关性**:计算 $\det G$ 或判断 $G$ 是否正定（顺序主子式、特征值、配方均可）。$\det G>0$ 时线性无关。
3. **内积的矩阵表示**:先确定基,再按 $g_{ij}=(\varepsilon_i,\varepsilon_j)$ 写出度量矩阵,最后用 $(\alpha,\beta)=x^HGy$ 计算。
4. **换基**:过渡矩阵满足 $\eta=\varepsilon P$ 时,$G_\eta=P^HG_\varepsilon P$;反之若已知 $G_\eta$,则 $G_\varepsilon=(P^H)^{-1}G_\eta P^{-1}$。

## 3.3 内积空间的度量

### 3.3.1 范数及其公理化刻画

由 3.1.3,内积可以诱导出长度 $\|\alpha\|=\sqrt{(\alpha,\alpha)}$。下面说明它满足一般的**范数公理**,从而"长度"确实是一个合理的度量。

**定义**: 设 $V$ 是数域 $F$ 上的线性空间,若对任意 $\alpha\in V$ 都对应一个非负实数 $\|\alpha\|$,且满足：

1. **正定性**：$\|\alpha\|\ge 0$,且 $\|\alpha\|=0\iff\alpha=0$；
2. **齐次性**：$\|k\alpha\|=|k|\,\|\alpha\|$；
3. **三角不等式**：$\|\alpha+\beta\|\le\|\alpha\|+\|\beta\|$；

则称 $\|\cdot\|$ 为 $V$ 上的一个**范数**,称 $(V,\|\cdot\|)$ 为**赋范线性空间**。

---

**定理**:

内积空间上由 $\|\alpha\|=\sqrt{(\alpha,\alpha)}$ 定义的 $\|\cdot\|$ 必为范数。

---

**证明**:

正定性由内积的正定性给出;齐次性由

$$
\|k\alpha\|^2=(k\alpha,k\alpha)=k\overline{k}(\alpha,\alpha)=|k|^2\|\alpha\|^2
$$

得到;三角不等式即 3.1.6 的结论。

---

**注意**:

- 反过来,**并非每个范数都可由内积诱导**。由内积诱导的范数必须满足**平行四边形法则**

$$
\|\alpha+\beta\|^2+\|\alpha-\beta\|^2=2\|\alpha\|^2+2\|\beta\|^2
$$

例如 $\mathbb{R}^2$ 上的 $\|x\|_\infty=\max\{|x_1|,|x_2|\}$ 就不满足平行四边形法则,因此它不能由任何内积诱导。

- 在 $\mathbb{C}^n$ 上由标准内积诱导的范数

$$
\|x\|_2=\left(\sum_{i=1}^{n}|x_i|^2\right)^{1/2}
$$

称为 **2-范数**（Euclid 范数）;矩阵的 Frobenius 范数

$$
\|A\|_F=\sqrt{tr(A^HA)}=\left(\sum_{i,j}|a_{ij}|^2\right)^{1/2}
$$

正是矩阵空间在 Frobenius 内积下诱导的范数。

### 3.3.2 距离

**定义**: 设 $V$ 为内积空间,$\alpha,\beta\in V$,称

$$
d(\alpha,\beta)=\|\alpha-\beta\|
$$

为 $\alpha$ 与 $\beta$ 的**距离**。

---

**性质**（距离三公理）:

1. **非负性与正定性**：$d(\alpha,\beta)\ge 0$,且 $d(\alpha,\beta)=0\iff\alpha=\beta$；
2. **对称性**：$d(\alpha,\beta)=d(\beta,\alpha)$；
3. **三角不等式**：$d(\alpha,\gamma)\le d(\alpha,\beta)+d(\beta,\gamma)$。

---

**证明**:

1、2 由范数的正定性与齐次性直接得到。对 3,由三角不等式

$$
d(\alpha,\gamma)=\|\alpha-\gamma\|=\|(\alpha-\beta)+(\beta-\gamma)\|
\le\|\alpha-\beta\|+\|\beta-\gamma\|=d(\alpha,\beta)+d(\beta,\gamma)
$$

---

**注意**:

- 距离具有**平移不变性**：$d(\alpha+\gamma,\beta+\gamma)=d(\alpha,\beta)$。
- 距离具有**齐次性**：$d(k\alpha,k\beta)=|k|\,d(\alpha,\beta)$。
- 正交投影的"最佳逼近"性质（见 3.7）正是用距离来描述的。

### 3.3.3 正交与线性无关

**定理**:

内积空间中**两两正交的非零向量组必线性无关**。

---

**证明**:

设 $\alpha_1,\alpha_2,\cdots,\alpha_k$ 两两正交且均非零,若有

$$
k_1\alpha_1+k_2\alpha_2+\cdots+k_k\alpha_k=0
$$

两端与 $\alpha_i$ 作内积。由正交性 $(\alpha_j,\alpha_i)=0\ (j\neq i)$,得

$$
k_i(\alpha_i,\alpha_i)=0
$$

又 $\alpha_i\neq0$,故 $(\alpha_i,\alpha_i)=\|\alpha_i\|^2>0$,于是 $k_i=0\ (i=1,\cdots,k)$,即向量组线性无关。

---

**定理**（勾股定理的推广）:

设 $\alpha_1,\alpha_2,\cdots,\alpha_k$ 两两正交,则

$$
\|\alpha_1+\alpha_2+\cdots+\alpha_k\|^2=\|\alpha_1\|^2+\|\alpha_2\|^2+\cdots+\|\alpha_k\|^2
$$

---

**证明**:

由内积的半双线性,

$$
\left\|\sum_{i=1}^{k}\alpha_i\right\|^2
=\sum_{i=1}^{k}\sum_{j=1}^{k}(\alpha_i,\alpha_j)
=\sum_{i=1}^{k}(\alpha_i,\alpha_i)
=\sum_{i=1}^{k}\|\alpha_i\|^2
$$

其中交叉项 $(\alpha_i,\alpha_j)=0\ (i\neq j)$。

---

**注意**: 正交性与线性无关的关系是"单向"的:线性无关的向量组未必两两正交,但两两正交的非零向量组一定线性无关。把线性无关组改造为两两正交组,就是 3.8 的 Schmidt 正交化。

### 3.3.4 度量中的基本不等式

**定理**（Cauchy–Schwarz 不等式的矩阵形式）:

对 $A\in\mathbb{C}^{m\times n}$, $x,y\in\mathbb{C}^n$ 有

$$
|x^HA^HAy|\le\|Ax\|_2\|Ay\|_2
$$

特别地,对任意 $x,y\in\mathbb{C}^n$,

$$
|x^Hy|^2\le\big(x^Hx\big)\big(y^Hy\big)=\|x\|_2^2\|y\|_2^2
$$

---

**证明**:

在 $\mathbb{C}^n$ 中对向量 $Ax,Ay$ 应用 Cauchy–Schwarz 不等式,

$$
|(Ax,Ay)|=|x^HA^HAy|\le\|Ax\|\,\|Ay\|
$$

取 $A=I$ 即得第二式。

---

**定理**（Bessel 不等式）:

设 $\varepsilon_1,\varepsilon_2,\cdots,\varepsilon_k$ 是内积空间 $V$ 中的**标准正交组**,则对任意 $\alpha\in V$ 有

$$
\sum_{i=1}^{k}|(\varepsilon_i,\alpha)|^2\le\|\alpha\|^2
$$

---

**证明**:

令

$$
\gamma=\alpha-\sum_{i=1}^{k}(\varepsilon_i,\alpha)\varepsilon_i
$$

则对每个 $j$,

$$
(\varepsilon_j,\gamma)=(\varepsilon_j,\alpha)-\sum_{i=1}^{k}(\varepsilon_i,\alpha)(\varepsilon_j,\varepsilon_i)
=(\varepsilon_j,\alpha)-(\varepsilon_j,\alpha)=0
$$

即 $\gamma$ 与每个 $\varepsilon_j$ 正交。由勾股定理,

$$
\|\alpha\|^2=\left\|\gamma+\sum_{i=1}^{k}(\varepsilon_i,\alpha)\varepsilon_i\right\|^2
=\|\gamma\|^2+\sum_{i=1}^{k}|(\varepsilon_i,\alpha)|^2
\ge\sum_{i=1}^{k}|(\varepsilon_i,\alpha)|^2
$$

---

**注意**: 当 $\varepsilon_1,\cdots,\varepsilon_k$ 是**标准正交基**时,$\gamma=0$,Bessel 不等式成为**等式**

$$
\|\alpha\|^2=\sum_{i=1}^{n}|(\varepsilon_i,\alpha)|^2
$$

这就是 **Parseval 等式**。

### 3.3.5 例题

*eg1: 计算距离与夹角*

在 $\mathbb{R}^3$ 中取标准内积,$\alpha=(1,1,1)^T,\ \beta=(1,2,0)^T$,求 $\|\alpha\|,\ \|\beta\|,\ d(\alpha,\beta)$ 及夹角 $\theta$。

**解**:

$$
\|\alpha\|=\sqrt3,\qquad\|\beta\|=\sqrt5
$$

$$
\alpha-\beta=(0,-1,1)^T,\qquad d(\alpha,\beta)=\|\alpha-\beta\|=\sqrt2
$$

$$
(\alpha,\beta)=1+2+0=3,\qquad
\cos\theta=\frac{3}{\sqrt3\cdot\sqrt5}=\frac{3}{\sqrt{15}}
$$

$$
\theta=\arccos\frac{3}{\sqrt{15}}
$$

---

*eg2: 正交向量组的线性无关性*

设 $\alpha_1=(1,1,1)^T,\ \alpha_2=(1,-1,0)^T,\ \alpha_3=(1,1,-2)^T$,验证它们两两正交,并由此断言它们线性无关。

**解**: 计算内积

$$
(\alpha_1,\alpha_2)=1-1+0=0,\quad
(\alpha_1,\alpha_3)=1+1-2=0,\quad
(\alpha_2,\alpha_3)=1-1+0=0
$$

即它们两两正交且均非零,故线性无关。

---

*eg3: 用勾股定理证明恒等式*

设 $\alpha\perp\beta$,证明

$$
\|\alpha+\beta\|^2+\|\alpha-\beta\|^2=2\|\alpha\|^2+2\|\beta\|^2
$$

**证明**: 由 $\alpha\perp\beta$ 得 $(\alpha,\beta)=0$,于是 $(\alpha,\beta)+(\beta,\alpha)=2\operatorname{Re}(\alpha,\beta)=0$,所以

$$
\|\alpha+\beta\|^2=\|\alpha\|^2+\|\beta\|^2,\qquad
\|\alpha-\beta\|^2=\|\alpha\|^2+\|\beta\|^2
$$

两式相加即得结论。

---

**方法总结**:

1. **算距离**:先作差 $\alpha-\beta$,再求其范数。
2. **判正交**:直接计算 $(\alpha,\beta)$ 是否为零;注意内积的选取会影响正交性。
3. **证明线性无关**:若能验证两两正交且非零,即可直接断言线性无关,无需解方程组。
4. **证明不等式**:优先套用 Cauchy–Schwarz、三角不等式、Bessel 不等式与勾股定理,而不是展开坐标。

## 3.4 标准正交基与酉矩阵

### 3.4.1 正交向量组与标准正交基

**定义**: 设 $V$ 为 $n$ 维内积空间,$\varepsilon_1,\varepsilon_2,\cdots,\varepsilon_n$ 是 $V$ 的一组基。

- 若 $(\varepsilon_i,\varepsilon_j)=0\ (i\neq j)$,则称它为一组**正交基**；
- 若进一步 $\|\varepsilon_i\|=1\ (\forall i)$,即

$$
(\varepsilon_i,\varepsilon_j)=\delta_{ij}=
\begin{cases}
1,&i=j\\
0,&i\neq j
\end{cases}
$$

则称它为一组**标准正交基**（规范正交基,复空间中也称为**酉基**）。

---

**定理**:

$n$ 维内积空间必存在标准正交基。

---

**注意**: 存在性由 Schmidt 正交化方法给出,具体构造见 3.8。这里先承认结论并讨论其性质。

---

**定理**（标准正交基下的坐标公式与内积公式）:

设 $\varepsilon_1,\cdots,\varepsilon_n$ 是标准正交基,向量

$$
\alpha=\sum_{i=1}^{n}x_i\varepsilon_i,\qquad \beta=\sum_{j=1}^{n}y_j\varepsilon_j
$$

则

$$
x_i=(\varepsilon_i,\alpha),\qquad
(\alpha,\beta)=\sum_{i=1}^{n}\overline{x_i}y_i=x^Hy,\qquad
\|\alpha\|^2=\sum_{i=1}^{n}|x_i|^2
$$

---

**证明**:

由内积对第二变元线性,

$$
(\varepsilon_i,\alpha)=\left(\varepsilon_i,\sum_{j=1}^{n}x_j\varepsilon_j\right)
=\sum_{j=1}^{n}x_j(\varepsilon_i,\varepsilon_j)=x_i
$$

再由半双线性,

$$
(\alpha,\beta)=\sum_{i=1}^{n}\sum_{j=1}^{n}\overline{x_i}y_j(\varepsilon_i,\varepsilon_j)
=\sum_{i=1}^{n}\overline{x_i}y_i=x^Hy
$$

取 $\beta=\alpha$ 得第三式。

---

**推论**:

1. 标准正交基下的**度量矩阵为单位矩阵** $G=I$,即内积化为坐标的标准内积。
2. 在标准正交基下,向量的坐标可以通过内积 $x_i=(\varepsilon_i,\alpha)$ 直接读出,无需解线性方程组。
3. 标准正交基使**向量的范数等于坐标的 2-范数**:$\|\alpha\|=\|x\|_2$。

---

*eg*

在 $\mathbb{R}^3$ 中,

$$
\varepsilon_1=\frac{1}{\sqrt2}(1,1,0)^T,\quad
\varepsilon_2=\frac{1}{\sqrt2}(1,-1,0)^T,\quad
\varepsilon_3=(0,0,1)^T
$$

两两内积为 $0$,且范数均为 $1$,故为标准正交基。向量 $\alpha=(3,1,2)^T$ 在此基下的坐标为

$$
x_1=(\varepsilon_1,\alpha)=\frac{4}{\sqrt2}=2\sqrt2,\quad
x_2=(\varepsilon_2,\alpha)=\frac{2}{\sqrt2}=\sqrt2,\quad
x_3=(\varepsilon_3,\alpha)=2
$$

于是 $\|\alpha\|^2=8+2+4=14$,与 $\|\alpha\|^2=9+1+4=14$ 一致。

### 3.4.2 酉矩阵与正交矩阵

**定义**: 设 $U\in\mathbb{C}^{n\times n}$,若

$$
U^HU=I
$$

则称 $U$ 为**酉矩阵**（unitary matrix）。若 $Q\in\mathbb{R}^{n\times n}$ 满足

$$
Q^TQ=I
$$

则称 $Q$ 为**正交矩阵**（orthogonal matrix）。

---

**定理**（酉矩阵的等价刻画）:

设 $U\in\mathbb{C}^{n\times n}$,则下列命题等价：

1. $U^HU=I$；
2. $U$ 可逆,且 $U^{-1}=U^H$；
3. $UU^H=I$；
4. $U$ 的 $n$ 个**列向量**构成 $\mathbb{C}^n$ 的标准正交组；
5. $U$ 的 $n$ 个**行向量**构成 $\mathbb{C}^n$ 的标准正交组。

---

**证明**:

记 $U=(u_1,u_2,\cdots,u_n)$,$u_i$ 为列向量。

**(1)$\iff$(4)**: $U^HU$ 的第 $(i,j)$ 元素为

$$
(U^HU)_{ij}=u_i^Hu_j=(u_i,u_j)
$$

故 $U^HU=I\iff(u_i,u_j)=\delta_{ij}$。

**(1)$\Rightarrow$(2)**: 由 $U^HU=I$ 知 $U$ 可逆(满秩),两边左乘 $U^{-1}$ 得 $U^{-1}=U^H$。

**(2)$\Rightarrow$(3)**: 由 $UU^{-1}=I$ 及 $U^{-1}=U^H$ 得 $UU^H=I$。

**(3)$\iff$(5)**: 对 $UU^H=I$ 转置共轭,或直接计算 $UU^H$ 的第 $(i,j)$ 元素为 $U$ 的第 $i$ 行与第 $j$ 行的内积。

**(3)$\Rightarrow$(1)**: 由 $UU^H=I$ 知 $U$ 可逆,于是 $U^H=(U)^{-1}$,从而 $U^HU=I$。

---

**性质**:

设 $U,V$ 为 $n$ 阶酉矩阵,则

1. $|\det U|=1$；
2. $U^{-1}=U^H$ 仍是酉矩阵,且 $U^H$ 也是酉矩阵；
3. $UV$ 仍是酉矩阵；
4. $U$ 的任意特征值 $\lambda$ 都满足 $|\lambda|=1$；
5. $U$ 保持内积与范数:对任意 $x,y\in\mathbb{C}^n$,

$$
(Ux,Uy)=(x,y),\qquad \|Ux\|=\|x\|
$$

---

**证明**:

1. 由 $\det(U^HU)=\overline{\det U}\det U=|\det U|^2=\det I=1$。
2. $(U^{-1})^HU^{-1}=(U^H)^HU^H=UU^H=I$。
3. $(UV)^H(UV)=V^HU^HUV=V^HV=I$。
4. 设 $Ux=\lambda x,\ x\neq0$,则

$$
x^Hx=x^HU^HUx=(Ux)^H(Ux)=\overline{\lambda}\lambda\,x^Hx=|\lambda|^2x^Hx
$$

故 $|\lambda|=1$。

5. $(Ux,Uy)=(Ux)^H(Uy)=x^HU^HUy=x^Hy=(x,y)$;取 $y=x$ 得范数保持。

---

**定理**（标准正交基之间的过渡矩阵是酉矩阵）:

设 $\varepsilon_1,\cdots,\varepsilon_n$ 与 $\eta_1,\cdots,\eta_n$ 是同一内积空间的两组标准正交基,且

$$
(\eta_1,\cdots,\eta_n)=(\varepsilon_1,\cdots,\varepsilon_n)P
$$

则过渡矩阵 $P$ 是酉矩阵。

---

**证明**:

由 3.2.3 的换基公式 $G_\eta=P^HG_\varepsilon P$。因为两组基都是标准正交基,故 $G_\varepsilon=G_\eta=I$,于是

$$
P^HP=I
$$

即 $P$ 为酉矩阵。

---

**注意**:

- 该定理说明:**标准正交基的"坐标变换"保持内积**,这是标准正交基在数值计算中特别稳定的根本原因。
- 实空间中的正交矩阵满足 $Q^{-1}=Q^T$,且 $\det Q=\pm1$。$\det Q=1$ 对应**旋转**,$\det Q=-1$ 对应**反射**。

---

*eg1: 旋转矩阵*

对任意 $\theta\in\mathbb{R}$,

$$
R_\theta=
\begin{pmatrix}
\cos\theta&-\sin\theta\\
\sin\theta&\cos\theta
\end{pmatrix}
$$

则

$$
R_\theta^TR_\theta=
\begin{pmatrix}
\cos^2\theta+\sin^2\theta&0\\
0&\sin^2\theta+\cos^2\theta
\end{pmatrix}
=I
$$

故 $R_\theta$ 是正交矩阵,且 $\det R_\theta=1$,它是平面上的**旋转**。

---

*eg2: 反射矩阵与 Hadamard 矩阵*

$$
H=\frac{1}{\sqrt2}\begin{pmatrix}1&1\\1&-1\end{pmatrix}
$$

满足 $H^T=H$ 且

$$
H^TH=\frac12\begin{pmatrix}1&1\\1&-1\end{pmatrix}\begin{pmatrix}1&1\\1&-1\end{pmatrix}=I
$$

故 $H$ 是正交矩阵,且 $\det H=-1$,对应**反射**。

---

*eg3: 判断酉矩阵*

设

$$
U=\frac{1}{\sqrt2}\begin{pmatrix}1&i\\ i&1\end{pmatrix}
$$

计算

$$
U^HU=\frac12\begin{pmatrix}1&-i\\-i&1\end{pmatrix}\begin{pmatrix}1&i\\ i&1\end{pmatrix}
=\frac12\begin{pmatrix}1+1&i-i\\ i-i&1+1\end{pmatrix}=I
$$

故 $U$ 是酉矩阵。

---

**方法总结**:

1. **验证标准正交基**:逐一计算内积 $(\varepsilon_i,\varepsilon_j)$,核对是否为 $\delta_{ij}$。
2. **验证酉矩阵**:计算 $U^HU$ 看是否等于 $I$;也可只看列向量是否标准正交。
3. **求标准正交基下的坐标**:直接用 $x_i=(\varepsilon_i,\alpha)$,不必解方程组。
4. **判断旋转与反射**:先用 $Q^TQ=I$ 判断正交,再看 $\det Q=1$ 还是 $-1$。

## 3.5 正交子空间

### 3.5.1 正交补的定义

**定义**: 设 $V$ 是内积空间,$W$ 是 $V$ 的子空间,称集合

$$
W^\perp=\{\alpha\in V\ :\ (\alpha,\beta)=0,\ \forall\beta\in W\}
$$

为 $W$ 的**正交补**（orthogonal complement）。

即 $W^\perp$ 由 $V$ 中所有与 $W$ 中**每个**向量都正交的向量组成。

---

**注意**: 正交是对**整个子空间**提出的要求:$\alpha\in W^\perp$ 必须与 $W$ 中所有向量正交,而不是只与 $W$ 中某几个向量正交。

### 3.5.2 正交补的基本性质

**性质**:

1. $W^\perp$ 是 $V$ 的**子空间**；
2. $\{0\}^\perp=V$,$V^\perp=\{0\}$；
3. $W\cap W^\perp=\{0\}$；
4. 若 $W\subseteq U$,则 $U^\perp\subseteq W^\perp$（反单调性）；
5. $\dim W+\dim W^\perp=\dim V$；
6. $(W^\perp)^\perp=W$（有限维内积空间）。

---

**证明**:

**1** 设 $\alpha_1,\alpha_2\in W^\perp$,$k\in\mathbb{C}$。对任意 $\beta\in W$,

$$
(k\alpha_1+\alpha_2,\beta)=\overline{k}(\alpha_1,\beta)+(\alpha_2,\beta)=0
$$

故 $k\alpha_1+\alpha_2\in W^\perp$,即 $W^\perp$ 对加法和数乘封闭,是子空间。

**2** 零向量与任意向量正交,故 $\{0\}^\perp=V$;若 $\alpha\in V^\perp$,则特别地 $(\alpha,\alpha)=0$,故 $\alpha=0$,即 $V^\perp=\{0\}$。

**3** 若 $\alpha\in W\cap W^\perp$,则 $(\alpha,\alpha)=0$,由正定性得 $\alpha=0$。

**4** 设 $W\subseteq U$,$\alpha\in U^\perp$。对任意 $\beta\in W\subseteq U$ 有 $(\alpha,\beta)=0$,故 $\alpha\in W^\perp$,即 $U^\perp\subseteq W^\perp$。

**5** 取 $W$ 的一组标准正交基 $\varepsilon_1,\cdots,\varepsilon_k$,由 3.4.1 的结论可将其扩充为 $V$ 的一组标准正交基

$$
\varepsilon_1,\cdots,\varepsilon_k,\varepsilon_{k+1},\cdots,\varepsilon_n
$$

下面证明

$$
W^\perp=\operatorname{span}\{\varepsilon_{k+1},\cdots,\varepsilon_n\}
$$

一方面,当 $j>k$ 时 $(\varepsilon_j,\varepsilon_i)=0\ (i\le k)$,故 $\varepsilon_j$ 与 $W$ 的标准正交基正交,从而与 $W$ 中每个向量正交,即 $\varepsilon_j\in W^\perp$。

另一方面,设 $\alpha=\sum_{i=1}^{n}x_i\varepsilon_i\in W^\perp$。由坐标公式 $x_i=(\varepsilon_i,\alpha)$,当 $i\le k$ 时 $\varepsilon_i\in W$,故 $x_i=0$,于是 $\alpha$ 只含 $\varepsilon_{k+1},\cdots,\varepsilon_n$ 的项。

因此 $\dim W^\perp=n-k$,即 $\dim W+\dim W^\perp=n$。

**6** 由 $W\subseteq(W^\perp)^\perp$ 直接可得一边包含。再由 5,

$$
\dim (W^\perp)^\perp=n-\dim W^\perp=n-(n-\dim W)=\dim W
$$

故 $(W^\perp)^\perp=W$。

---

**定理**（正交直和分解）:

设 $W$ 是 $n$ 维内积空间 $V$ 的子空间,则 $V$ 可以分解为 $W$ 与 $W^\perp$ 的**正交直和**:

$$
V=W\oplus W^\perp
$$

即每个 $\alpha\in V$ 都可以**唯一**地表示为

$$
\alpha=\beta+\gamma,\qquad \beta\in W,\ \gamma\in W^\perp
$$

且其中的 $\beta$ 就是 $\alpha$ 在 $W$ 上的**正交投影**（见 3.7）。

---

**证明**:

取 $W$ 的标准正交基 $\varepsilon_1,\cdots,\varepsilon_k$,令

$$
\beta=\sum_{i=1}^{k}(\varepsilon_i,\alpha)\varepsilon_i\in W,\qquad \gamma=\alpha-\beta
$$

则对每个 $j\le k$,

$$
(\varepsilon_j,\gamma)=(\varepsilon_j,\alpha)-\sum_{i=1}^{k}(\varepsilon_i,\alpha)(\varepsilon_j,\varepsilon_i)
=(\varepsilon_j,\alpha)-(\varepsilon_j,\alpha)=0
$$

故 $\gamma\in W^\perp$,即 $\alpha=\beta+\gamma$ 是 $W+W^\perp$ 中的分解。

唯一性由 $W\cap W^\perp=\{0\}$ 得到:$W+W^\perp$ 是直和。

---

**性质**（正交补与子空间运算）:

设 $W_1,W_2$ 是 $V$ 的子空间,则

$$
(W_1+W_2)^\perp=W_1^\perp\cap W_2^\perp,\qquad
(W_1\cap W_2)^\perp=W_1^\perp+W_2^\perp
$$

---

**证明**:

第一式:$\alpha\in(W_1+W_2)^\perp\iff\alpha$ 与 $W_1+W_2$ 中所有向量正交 $\iff\alpha$ 分别与 $W_1$ 和 $W_2$ 中所有向量正交 $\iff\alpha\in W_1^\perp\cap W_2^\perp$。

第二式:在第一式中把 $W_1,W_2$ 换成 $W_1^\perp,W_2^\perp$,并利用性质 6,

$$
(W_1^\perp+W_2^\perp)^\perp=(W_1^\perp)^\perp\cap(W_2^\perp)^\perp=W_1\cap W_2
$$

两端再取正交补即得。

### 3.5.3 四个基本子空间的正交关系

正交补把矩阵的四个基本子空间整齐地联系起来。

**定理**:

设 $A\in\mathbb{C}^{m\times n}$,则

$$
\operatorname{range}(A)^\perp=\operatorname{null}(A^H),\qquad
\operatorname{null}(A)^\perp=\operatorname{range}(A^H)
$$

从而

$$
\mathbb{C}^m=\operatorname{range}(A)\oplus\operatorname{null}(A^H),\qquad
\mathbb{C}^n=\operatorname{range}(A^H)\oplus\operatorname{null}(A)
$$

---

**证明**:

设 $y\in\mathbb{C}^m$。则

$$
y\in\operatorname{range}(A)^\perp
\iff y^H(Ax)=0,\quad \forall x\in\mathbb{C}^n
\iff (A^Hy)^Hx=0,\quad \forall x
\iff A^Hy=0
\iff y\in\operatorname{null}(A^H)
$$

即第一式成立。对矩阵 $A^H$ 应用第一式得

$$
\operatorname{range}(A^H)^\perp=\operatorname{null}\big((A^H)^H\big)=\operatorname{null}(A)
$$

两端取正交补,并利用 $(W^\perp)^\perp=W$,得第二式。最后两式由正交直和分解给出。

---

**推论**:

1. $rank(A)=rank(A^H)$；
2. **秩–零度定理**：$rank(A)+\dim\operatorname{null}(A)=n$；
3. 方程组 $Ax=b$ 有解 $\iff b\perp\operatorname{null}(A^H)$,即 $b$ 与 $A^H$ 的零空间正交。

---

**注意**: 上述结论给出了**最小二乘问题**的理论基础:$Ax=b$ 无解时,可在 $\operatorname{range}(A)$ 中寻找与 $b$ 距离最近的向量,而 $b$ 到 $\operatorname{range}(A)$ 的正交投影正是最小二乘解所对应的像。

### 3.5.4 例题

*eg1: 求平面子空间的正交补*

在 $\mathbb{R}^3$ 中取标准内积,设

$$
W=\{(x_1,x_2,x_3)^T\ :\ x_1+x_2+x_3=0\}
$$

求 $W^\perp$。

**解**: $W$ 是由向量 $(1,1,1)^T$ 的方程定义的一张过原点的平面,而

$$
W=\{(1,-1,0)^T,\ (1,0,-1)^T\}^\perp
$$

事实上,对 $\alpha=(a_1,a_2,a_3)^T$,

$$
\alpha\in W^\perp
\iff (\alpha,\beta)=0,\ \forall\beta\in W
$$

取 $\beta_1=(1,-1,0)^T,\ \beta_2=(1,0,-1)^T$（它们张成 $W$）,得

$$
a_1-a_2=0,\qquad a_1-a_3=0
\ \Longrightarrow\ a_1=a_2=a_3
$$

故

$$
W^\perp=\operatorname{span}\{(1,1,1)^T\}
$$

这与 $\dim W+\dim W^\perp=2+1=3$ 一致。

---

*eg2: 正交直和分解*

在 $\mathbb{R}^3$ 中,设 $W=\operatorname{span}\{(1,1,1)^T\}$,$\alpha=(1,2,3)^T$,将 $\alpha$ 分解为 $W$ 与 $W^\perp$ 中向量之和。

**解**: 令 $\varepsilon=\dfrac{1}{\sqrt3}(1,1,1)^T$ 为 $W$ 的单位向量,则

$$
\beta=(\varepsilon,\alpha)\varepsilon=\frac{6}{3}(1,1,1)^T=(2,2,2)^T\in W
$$

$$
\gamma=\alpha-\beta=(1,2,3)^T-(2,2,2)^T=(-1,0,1)^T\in W^\perp
$$

验证:$(\gamma,(1,1,1)^T)=-1+0+1=0$,确实正交。

---

*eg3: 四个基本子空间*

设

$$
A=\begin{pmatrix}1&0&1\\0&1&1\end{pmatrix}
$$

求 $\operatorname{range}(A)$、$\operatorname{null}(A)$、$\operatorname{range}(A^H)$、$\operatorname{null}(A^H)$,并验证正交关系。

**解**:

$$
\operatorname{range}(A)=\mathbb{C}^2,\qquad \operatorname{null}(A^H)=\{0\}
$$

解 $Ax=0$:

$$
\begin{cases}
x_1+x_3=0\\
x_2+x_3=0
\end{cases}
\ \Longrightarrow\
\operatorname{null}(A)=\operatorname{span}\{(-1,-1,1)^T\}
$$

$$
\operatorname{range}(A^H)=\operatorname{span}\{(1,0,1)^T,\ (0,1,1)^T\}
$$

验证:$\operatorname{null}(A)^\perp=\{(x_1,x_2,x_3)^T:-x_1-x_2+x_3=0\}=\operatorname{span}\{(1,0,1)^T,(0,1,1)^T\}=\operatorname{range}(A^H)$。

---

**方法总结**:

1. **求正交补**:把 $W$ 的基作为行（或列）写出齐次方程组 $(\alpha,\beta_i)=0$,解出即得 $W^\perp$;再由维数公式核对 $\dim W^\perp=n-\dim W$。
2. **求正交分解**:先由正交投影公式求 $\beta$,再令 $\gamma=\alpha-\beta$,最后用 $(\gamma,\beta)=0$ 验证。
3. **验证正交关系**:利用 $\operatorname{range}(A)^\perp=\operatorname{null}(A^H)$ 与 $\operatorname{null}(A)^\perp=\operatorname{range}(A^H)$,以及维数公式。

## 3.6 酉变换

### 3.6.1 酉变换的定义

**定义**: 设 $V$ 是内积空间,$T:V\to V$ 是**线性变换**。若对任意 $\alpha,\beta\in V$ 都有

$$
(T\alpha,T\beta)=(\alpha,\beta)
$$

则称 $T$ 为 $V$ 上的**酉变换**（复内积空间）或**正交变换**（实内积空间）。二者统称为**保内积变换**。

---

**注意**:

- 酉变换的定义中包含**线性**这一前提:只是在集合上保持内积的映射未必是酉变换。
- 由定义直接可得,酉变换**保持向量的长度**:在定义中取 $\beta=\alpha$ 得 $\|T\alpha\|=\|\alpha\|$。

### 3.6.2 酉变换的等价刻画

**定理**:

设 $T:V\to V$ 是有限维内积空间上的线性变换,则下列命题等价：

1. $T$ 保持内积:$ (T\alpha,T\beta)=(\alpha,\beta)$；
2. $T$ 保持范数:$\|T\alpha\|=\|\alpha\|$；
3. $T$ 把**某一组**标准正交基映为标准正交基；
4. $T$ 把**任意一组**标准正交基映为标准正交基；
5. $T$ 在**任意**标准正交基下的矩阵是**酉矩阵**；
6. $T$ 可逆,且其**伴随变换**满足 $T^*T=TT^*=I$,即 $T^*=T^{-1}$。

---

**证明**:

**(1)$\Rightarrow$(2)**:取 $\beta=\alpha$。

**(2)$\Rightarrow$(1)**:由范数保持,对任意 $\alpha,\beta$,

$$
\|T(\alpha+\beta)\|^2=\|T\alpha\|^2+\|T\beta\|^2+2\operatorname{Re}(T\alpha,T\beta)
$$

又 $\|T(\alpha+\beta)\|^2=\|\alpha+\beta\|^2=\|\alpha\|^2+\|\beta\|^2+2\operatorname{Re}(\alpha,\beta)$,于是

$$
\operatorname{Re}(T\alpha,T\beta)=\operatorname{Re}(\alpha,\beta)
$$

把 $\beta$ 换成 $i\beta$,注意 $T(i\beta)=iT\beta$,得

$$
\operatorname{Re}\big(i(T\alpha,T\beta)\big)=\operatorname{Re}\big(i(\alpha,\beta)\big)
\ \Longrightarrow\
-\operatorname{Im}(T\alpha,T\beta)=-\operatorname{Im}(\alpha,\beta)
$$

故 $\operatorname{Im}(T\alpha,T\beta)=\operatorname{Im}(\alpha,\beta)$,实部虚部都相等,即 $(T\alpha,T\beta)=(\alpha,\beta)$。

**(1)$\Rightarrow$(4)**:设 $\varepsilon_1,\cdots,\varepsilon_n$ 是标准正交基,则

$$
(T\varepsilon_i,T\varepsilon_j)=(\varepsilon_i,\varepsilon_j)=\delta_{ij}
$$

即 $T\varepsilon_1,\cdots,T\varepsilon_n$ 仍是标准正交组,且个数为 $n$,故为标准正交基。

**(4)$\Rightarrow$(3)**:显然。

**(3)$\Rightarrow$(1)**:设 $\varepsilon_1,\cdots,\varepsilon_n$ 是标准正交基且 $T\varepsilon_1,\cdots,T\varepsilon_n$ 也是标准正交基。对 $\alpha=\sum x_i\varepsilon_i,\ \beta=\sum y_j\varepsilon_j$,由 $T$ 的线性,

$$
(T\alpha,T\beta)=\left(\sum_i x_iT\varepsilon_i,\ \sum_j y_jT\varepsilon_j\right)
=\sum_{i,j}\overline{x_i}y_j\delta_{ij}
=\sum_i\overline{x_i}y_i=(\alpha,\beta)
$$

**(4)$\iff$(5)**:设 $T$ 在标准正交基 $\varepsilon$ 下的矩阵为 $U$,即 $T\varepsilon_j=\sum_i u_{ij}\varepsilon_i$。由 3.4.1 的坐标公式,

$$
u_{ij}=(\varepsilon_i,T\varepsilon_j)
$$

于是 $T\varepsilon_1,\cdots,T\varepsilon_n$ 为标准正交基 $\iff$ 矩阵 $U=(u_{ij})$ 的列向量标准正交 $\iff U^HU=I\iff U$ 为酉矩阵。

**(5)$\Rightarrow$(6)**:矩阵酉意味着可逆,且 $T^{-1}$ 在标准正交基下的矩阵为 $U^H$,故 $T^*=T^{-1}$。

**(6)$\Rightarrow$(1)**:$(T\alpha,T\beta)=(\alpha,T^*T\beta)=(\alpha,\beta)$。

---

**注意**: 由该定理,**判断一个线性变换是否为酉变换,最方便的做法是看它在标准正交基下的矩阵是否为酉矩阵**。

### 3.6.3 酉变换的性质

**性质**:

设 $T,S$ 是 $n$ 维内积空间上的酉变换,则

1. $T$ 是可逆变换,且 $T^{-1}$ 也是酉变换；
2. $TS$ 是酉变换；
3. $T$ 保持距离:$d(T\alpha,T\beta)=d(\alpha,\beta)$；
4. $T$ 的特征值的模都为 $1$；
5. 在复空间中 $|\det T|=1$;在实空间中 $\det T=\pm1$；
6. 恒等变换 $I$ 与 $-I$ 都是酉变换（$-I$ 在实空间是中心对称）。

---

**证明**:

1. 由 $T^*T=I$ 知 $T$ 可逆;又 $(T^{-1})^*T^{-1}=(T^*)^{-1}T^{-1}=(TT^*)^{-1}=I$,故 $T^{-1}$ 是酉变换。
2. $(TS)^*(TS)=S^*T^*TS=S^*S=I$。
3. $d(T\alpha,T\beta)=\|T\alpha-T\beta\|=\|T(\alpha-\beta)\|=\|\alpha-\beta\|=d(\alpha,\beta)$。
4. 设 $T\alpha=\lambda\alpha,\ \alpha\neq0$,则

$$
\|\alpha\|=\|T\alpha\|=\|\lambda\alpha\|=|\lambda|\,\|\alpha\|
\ \Longrightarrow\ |\lambda|=1
$$

5. 由 $T^*T=I$ 取行列式,$|\det T|^2=1$;实变换的行列式为实数,故 $\det T=\pm1$。

---

**注意**:

- 酉变换是**等距同构**,它把内积空间"原样"地搬到自身,因此保持一切由内积导出的几何量:长度、距离、夹角、正交性。
- 实空间中 $\det T=1$ 的酉变换称为**第一类**（旋转）,$\det T=-1$ 的称为**第二类**（反射或旋转加反射）。

### 3.6.4 实正交变换的几何分类

**定理**（二维正交变换的标准形）:

设 $Q$ 是 $2$ 阶实正交矩阵,则

$$
Q=\begin{pmatrix}\cos\theta&-\sin\theta\\ \sin\theta&\cos\theta\end{pmatrix}
\quad(\det Q=1)
$$

或

$$
Q=\begin{pmatrix}\cos\theta&\sin\theta\\ \sin\theta&-\cos\theta\end{pmatrix}
\quad(\det Q=-1)
$$

前者是绕原点的**旋转**,后者是关于过原点某条直线的**反射**。

---

**证明**: 设

$$
Q=\begin{pmatrix}a&b\\ c&d\end{pmatrix}
$$

由 $Q^TQ=I$ 得 $a^2+c^2=1,\ b^2+d^2=1,\ ab+cd=0$。由前两式可令

$$
a=\cos\theta,\ c=\sin\theta
$$

再由 $ab+cd=0$ 与 $b^2+d^2=1$ 解得 $(b,d)=(-\sin\theta,\cos\theta)$ 或 $(\sin\theta,-\cos\theta)$。前者对应 $\det Q=1$,后者对应 $\det Q=-1$。

---

**注意**: 在 $n$ 维实内积空间中,任何正交变换都可以在适当的正交基下化为若干二维旋转块（以及可能的 $-1$ 块）的直和,即

$$
\begin{pmatrix}
R_{\theta_1}&&&\\
&\ddots&&\\
&&R_{\theta_k}&\\
&&&\pm I
\end{pmatrix}
$$

这是实正交矩阵在正交相似意义下的标准形。

### 3.6.5 例题

*eg1: 判断酉变换*

在 $\mathbb{C}^2$ 上定义

$$
T\binom{x_1}{x_2}=\frac{1}{\sqrt2}\binom{x_1+ix_2}{ix_1+x_2}
$$

判断 $T$ 是否为酉变换。

**解**: $T$ 在标准正交基下的矩阵为

$$
U=\frac{1}{\sqrt2}\begin{pmatrix}1&i\\ i&1\end{pmatrix}
$$

而由 3.4.2 的 eg3,$U^HU=I$,故 $U$ 是酉矩阵,从而 $T$ 是酉变换。事实上,任取 $x$,直接计算可得 $\|Tx\|=\|x\|$。

---

*eg2: 旋转是正交变换*

在 $\mathbb{R}^2$ 上定义 $T(x,y)=(x\cos\theta-y\sin\theta,\ x\sin\theta+y\cos\theta)$。由 3.4.2 的 eg1,$R_\theta$ 是正交矩阵,故 $T$ 是正交变换。它保持长度、距离与夹角,且 $\det R_\theta=1$,是**旋转**。

---

*eg3: Householder 变换*

设 $u\in\mathbb{C}^n$ 且 $\|u\|=1$,定义

$$
H=I-2uu^H
$$

证明 $H$ 是酉矩阵,且 $H^2=I$、$\det H=-1$。

**证明**: 首先

$$
H^H=(I-2uu^H)^H=I-2uu^H=H
$$

其次,

$$
H^2=I-4uu^H+4u(u^Hu)u^H=I-4uu^H+4uu^H=I
$$

于是 $H^HH=H^2=I$,即 $H$ 是酉矩阵。又 $H$ 的特征值只有 $\pm1$,且 $Hu=-u$,其余特征值为 $1$(在 $u^\perp$ 上 $H=I$),故 $\det H=-1$。从几何上看,$H$ 是关于超平面 $u^\perp$ 的**反射**。

---

**方法总结**:

1. **判断酉变换**:写出变换在标准正交基下的矩阵,验证 $U^HU=I$;或验证 $\|T\alpha\|=\|\alpha\|$ 对所有 $\alpha$ 成立。
2. **利用伴随变换**:$T$ 为酉变换 $\iff T^*=T^{-1}$,这一形式在证明中往往最简洁。
3. **几何解释**:实正交变换按 $\det$ 分为旋转与反射;二维情形可直接写出 $R_\theta$ 或反射矩阵。
4. **特征值方法**:酉变换的特征值模为 $1$,可用于构造例子或反证。

## 3.7 幂等矩阵、正交投影与勾股定理

### 3.7.1 幂等矩阵

**定义**: 设 $P\in\mathbb{C}^{n\times n}$,若

$$
P^2=P
$$

则称 $P$ 为**幂等矩阵**（idempotent matrix）。

---

**性质**:

1. $P$ 的特征值只能是 $0$ 或 $1$；
2. $I-P$ 也是幂等矩阵；
3. $P^H$ 也是幂等矩阵；
4. $\operatorname{range}(P)=\{x:Px=x\}=\operatorname{null}(I-P)$；
5. $\mathbb{C}^n=\operatorname{range}(P)\oplus\operatorname{null}(P)$（一般不是正交直和）；
6. $rank(P)=tr(P)$；
7. 幂等矩阵一定**可对角化**,且相似于

$$
\begin{pmatrix}
I_r&0\\
0&0
\end{pmatrix},\qquad r=rank(P)
$$

---

**证明**:

**1** 设 $Px=\lambda x,\ x\neq0$,则

$$
\lambda x=Px=P^2x=\lambda Px=\lambda^2x
\ \Longrightarrow\ (\lambda^2-\lambda)x=0
\ \Longrightarrow\ \lambda\in\{0,1\}
$$

**2** $(I-P)^2=I-2P+P^2=I-2P+P=I-P$。

**3** $(P^H)^2=(P^2)^H=P^H$。

**4** 若 $x\in\operatorname{range}(P)$,则 $x=Py$,于是 $Px=P^2y=Py=x$,故 $x\in\{x:Px=x\}$;反之若 $Px=x$,则显然 $x\in\operatorname{range}(P)$。又 $Px=x\iff(I-P)x=0\iff x\in\operatorname{null}(I-P)$。

**5** 对任意 $x$,

$$
x=Px+(I-P)x,\qquad Px\in\operatorname{range}(P),\quad (I-P)x\in\operatorname{null}(P)
$$

故 $\mathbb{C}^n=\operatorname{range}(P)+\operatorname{null}(P)$;若 $x\in\operatorname{range}(P)\cap\operatorname{null}(P)$,则 $x=Px=0$,故为直和。

**6** 由 7 知 $P$ 相似于 $\operatorname{diag}(I_r,0)$,相似矩阵的迹相同,故

$$
tr(P)=tr\begin{pmatrix}I_r&0\\0&0\end{pmatrix}=r=rank(P)
$$

**7** $P$ 的最小多项式整除 $\lambda(\lambda-1)$,它没有重根,故 $P$ 可对角化;对角元只能是 $0$ 或 $1$,重排后即为 $\operatorname{diag}(I_r,0)$。

---

**注意**: 幂等矩阵不一定是正交投影。例如

$$
P=\begin{pmatrix}1&1\\0&0\end{pmatrix}
$$

满足 $P^2=P$,但它不是 Hermite 矩阵。此时

$$
\operatorname{range}(P)=\operatorname{span}\{(1,0)^T\},\qquad
\operatorname{null}(P)=\operatorname{span}\{(1,-1)^T\}
$$

二者并不正交。$P$ 是沿 $(1,-1)^T$ 方向到 $x$ 轴的**斜投影**。

### 3.7.2 正交投影

**定义**: 设 $P\in\mathbb{C}^{n\times n}$（或 $\mathbb{R}^{n\times n}$）,若

$$
P^2=P,\qquad P^H=P
$$

则称 $P$ 为**正交投影矩阵**（orthogonal projection matrix）。实空间中条件为 $P^2=P,\ P^T=P$。

---

**定理**（正交投影的判别）:

设 $P$ 是幂等矩阵,则 $P$ 是正交投影矩阵的充分必要条件是

$$
\operatorname{range}(P)\perp\operatorname{null}(P)
$$

---

**证明**:

**必要性**:设 $x\in\operatorname{range}(P),\ y\in\operatorname{null}(P)$。由 $P^H=P$,

$$
(x,y)=(Px,y)=(x,P^Hy)=(x,Py)=(x,0)=0
$$

故 $\operatorname{range}(P)\perp\operatorname{null}(P)$。

**充分性**:设 $P^2=P$ 且 $\operatorname{range}(P)\perp\operatorname{null}(P)$。把任意 $x,y$ 按 $\mathbb{C}^n=\operatorname{range}(P)\oplus\operatorname{null}(P)$ 分解为

$$
x=x_1+x_2,\qquad y=y_1+y_2,\qquad x_1,y_1\in\operatorname{range}(P),\ x_2,y_2\in\operatorname{null}(P)
$$

则

$$
(Px,y)=(x_1,y_1+y_2)=(x_1,y_1),\qquad
(x,Py)=(x_1+x_2,y_1)=(x_1,y_1)
$$

故 $(Px,y)=(x,Py)$ 对一切 $x,y$ 成立,即 $P^H=P$。

---

**定理**（正交投影的存在性与唯一性）:

设 $S$ 是 $\mathbb{C}^n$ 的子空间,则存在唯一的正交投影矩阵 $P$,使得

$$
\operatorname{range}(P)=S,\qquad \operatorname{null}(P)=S^\perp
$$

称 $P$ 为**在 $S$ 上的正交投影**。

---

**定理**（正交投影的两种标准构造）:

1. 若 $S$ 的**标准正交基**构成矩阵 $Q$ 的列（即 $Q^HQ=I$）,则

$$
P=QQ^H
$$

2. 若 $A$ 的列向量张成 $S$ 且 $A$ **列满秩**,则

$$
P=A(A^HA)^{-1}A^H
$$

---

**证明**:

1. 由 $Q^HQ=I$ 得 $(QQ^H)^2=QQ^HQQ^H=QQ^H$,且 $(QQ^H)^H=QQ^H$,故 $P=QQ^H$ 是正交投影矩阵。又对任意 $x$,

$$
Px=Q(Q^Hx)\in\operatorname{range}(Q)=S
$$

且若 $Q^Hx=0$,则 $Px=0$;反之 $Px=0\Rightarrow Q^HQ(Q^Hx)=Q^Hx=0$。故 $\operatorname{range}(P)=S$。

2. 当 $A$ 列满秩时 $A^HA$ 正定,故可逆。令 $P=A(A^HA)^{-1}A^H$,则

$$
P^2=A(A^HA)^{-1}\underline{A^HA(A^HA)^{-1}}A^H=A(A^HA)^{-1}A^H=P
$$

$$
P^H=A(A^HA)^{-H}A^H=A(A^HA)^{-1}A^H=P
$$

(利用了 $A^HA$ 为 Hermite 矩阵)。又 $x\in\operatorname{range}(P)\iff x=Az$,且

$$
Px=0\iff A(A^HA)^{-1}A^Hx=0\iff A^Hx=0
$$

故 $\operatorname{null}(P)=\operatorname{null}(A^H)=\operatorname{range}(A)^\perp=S^\perp$。

---

**性质**:

设 $P$ 是到 $S$ 的正交投影,则

1. $I-P$ 是到 $S^\perp$ 的正交投影；
2. $\|Px\|\le\|x\|$,且等号成立 $\iff x\in S$；
3. $tr(P)=\dim S$；
4. **分解公式**：

$$
x=Px+(I-P)x,\qquad Px\in S,\ (I-P)x\in S^\perp
$$

5. 若 $P_1,P_2$ 分别是到 $S_1\perp S_2$ 的正交投影,则 $P_1P_2=P_2P_1=0$;若 $S_1\subseteq S_2$,则 $P_1P_2=P_2P_1=P_1$。

---

**证明**:

1. $(I-P)^2=I-P$,$(I-P)^H=I-P$;且 $\operatorname{range}(I-P)=\operatorname{null}(P)=S^\perp$。
2. 由正交分解 $x=Px+(I-P)x$ 与勾股定理,

$$
\|x\|^2=\|Px\|^2+\|(I-P)x\|^2\ge\|Px\|^2
$$

等号成立 $\iff(I-P)x=0\iff x\in\operatorname{range}(P)=S$。

3. 由 3.7.1 性质 6,$tr(P)=rank(P)=\dim\operatorname{range}(P)=\dim S$。

---

### 3.7.3 勾股定理与最佳逼近

**定理**（勾股定理）:

设 $x\perp y$,则

$$
\|x+y\|^2=\|x\|^2+\|y\|^2
$$

更一般地,若 $x_1,x_2,\cdots,x_k$ 两两正交,则

$$
\|x_1+x_2+\cdots+x_k\|^2=\|x_1\|^2+\|x_2\|^2+\cdots+\|x_k\|^2
$$

---

**定理**（正交投影的最佳逼近性质）:

设 $S$ 是内积空间 $V$ 的子空间,$P$ 是到 $S$ 的正交投影,则对任意 $x\in V$ 与任意 $y\in S$,

$$
\|x-Px\|\le\|x-y\|
$$

且等号成立当且仅当 $y=Px$。即 $Px$ 是 $S$ 中**距离 $x$ 最近**的向量。

---

**证明**:

由 $P$ 为到 $S$ 的正交投影,

$$
x-Px\in S^\perp,\qquad Px-y\in S
$$

二者正交,由勾股定理

$$
\|x-y\|^2=\big\|(x-Px)+(Px-y)\big\|^2=\|x-Px\|^2+\|Px-y\|^2\ge\|x-Px\|^2
$$

等号成立 $\iff\|Px-y\|=0\iff y=Px$。

---

**注意**: 该定理是**最小二乘方法**的几何本质:当 $Ax=b$ 无解时,最小二乘解就是使 $Ax$ 为 $b$ 在 $\operatorname{range}(A)$ 上正交投影的 $x$,即

$$
A^HAx=A^Hb
$$

这也解释了 3.8 中 QR 分解在求解最小二乘问题中的作用。

### 3.7.4 例题

*eg1: 用标准正交基构造投影矩阵*

在 $\mathbb{R}^3$ 中,$S=\operatorname{span}\{(1,1,1)^T\}$,求到 $S$ 的正交投影矩阵,并把 $x=(1,2,3)^T$ 作正交分解。

**解**: $S$ 的单位向量为

$$
q=\frac{1}{\sqrt3}(1,1,1)^T
$$

故

$$
P=qq^T=\frac13\begin{pmatrix}1&1&1\\1&1&1\\1&1&1\end{pmatrix}
$$

$$
Px=\frac13\begin{pmatrix}6\\6\\6\end{pmatrix}=(2,2,2)^T,\qquad
x-Px=(-1,0,1)^T
$$

验证:$(2,2,2)^T\in S$,$(-1,0,1)^T\perp(1,1,1)^T$,且

$$
\|x\|^2=14=\|Px\|^2+\|x-Px\|^2=12+2
$$

勾股定理成立。

---

*eg2: 用一般公式求投影矩阵*

求 $\mathbb{R}^3$ 中到平面

$$
S=\{x=(x_1,x_2,x_3)^T:x_1+x_2+x_3=0\}
$$

的正交投影矩阵。

**解**: 取 $S$ 的一组基为列构成

$$
A=\begin{pmatrix}1&1\\-1&0\\0&-1\end{pmatrix}
$$

则

$$
A^HA=\begin{pmatrix}2&1\\1&2\end{pmatrix},\qquad
(A^HA)^{-1}=\frac13\begin{pmatrix}2&-1\\-1&2\end{pmatrix}
$$

$$
P=A(A^HA)^{-1}A^H
=\frac13\begin{pmatrix}
2&-1&-1\\
-1&2&-1\\
-1&-1&2
\end{pmatrix}
=I-\frac13J
$$

其中 $J$ 为全 $1$ 矩阵。这与 $S^\perp=\operatorname{span}\{(1,1,1)^T\}$、$P=I-qq^T$ 一致。

---

*eg3: 幂等但非正交的投影*

验证 $P=\begin{pmatrix}1&1\\0&0\end{pmatrix}$ 是幂等矩阵但不是正交投影矩阵,并求出它的投影方向。

**解**: $P^2=P$ 直接计算可得,故 $P$ 幂等。但 $P^T\neq P$,故它不是正交投影。

$$
\operatorname{range}(P)=\operatorname{span}\{(1,0)^T\},\qquad
\operatorname{null}(P)=\operatorname{span}\{(1,-1)^T\}
$$

二者内积为 $1\neq0$,不正交。$P$ 把向量沿 $(1,-1)^T$ 方向投影到 $x$ 轴上,是**斜投影**。

---

**方法总结**:

1. **验证正交投影**:核对 $P^2=P$ 与 $P^H=P$ 两条;只满足第一条是斜投影。
2. **求投影矩阵**:
   - 若 $S$ 的标准正交基容易得到,用 $P=QQ^H$;
   - 否则取 $S$ 的一组基作列构成 $A$,用 $P=A(A^HA)^{-1}A^H$。
3. **求正交分解与距离**:先算 $Px$,再算 $x-Px$;距离即 $\|x-Px\|$。
4. **验证结果**:用 $\|x\|^2=\|Px\|^2+\|x-Px\|^2$（勾股定理）检验,并确认 $Px\in S$、$x-Px\in S^\perp$。

## 3.8 Schmidt正交化方法与QR分解

### 3.8.1 Schmidt 正交化方法

**问题**: 给定内积空间中一组**线性无关**的向量 $\alpha_1,\alpha_2,\cdots,\alpha_k$,如何构造一组与之"等价"的标准正交向量组?

**方法**（Schmidt 正交化 / Gram–Schmidt 过程）:

**第一步**（正交化）:依次令

$$
\begin{align}
\beta_1&=\alpha_1\\
\beta_2&=\alpha_2-\frac{(\beta_1,\alpha_2)}{(\beta_1,\beta_1)}\beta_1\\
\beta_3&=\alpha_3-\frac{(\beta_1,\alpha_3)}{(\beta_1,\beta_1)}\beta_1-\frac{(\beta_2,\alpha_3)}{(\beta_2,\beta_2)}\beta_2\\
&\ \vdots\\
\beta_k&=\alpha_k-\sum_{i=1}^{k-1}\frac{(\beta_i,\alpha_k)}{(\beta_i,\beta_i)}\beta_i
\end{align}
$$

**第二步**（单位化）:令

$$
\varepsilon_i=\frac{\beta_i}{\|\beta_i\|},\qquad i=1,2,\cdots,k
$$

则 $\varepsilon_1,\varepsilon_2,\cdots,\varepsilon_k$ 是一组**标准正交向量组**。

---

**定理**:

设 $\alpha_1,\cdots,\alpha_k$ 线性无关,按上述方法得到的 $\beta_1,\cdots,\beta_k$ 均非零,且对每个 $j$,

$$
\operatorname{span}\{\beta_1,\cdots,\beta_j\}=\operatorname{span}\{\alpha_1,\cdots,\alpha_j\}
$$

$$
\operatorname{span}\{\varepsilon_1,\cdots,\varepsilon_j\}=\operatorname{span}\{\alpha_1,\cdots,\alpha_j\}
$$

从而 $\varepsilon_1,\cdots,\varepsilon_k$ 是 $\operatorname{span}\{\alpha_1,\cdots,\alpha_k\}$ 的一组标准正交基。

---

**证明**:

对 $j$ 归纳。$j=1$ 时 $\beta_1=\alpha_1\neq0$。

设 $\beta_1,\cdots,\beta_{j-1}$ 两两正交且非零,且 $\operatorname{span}\{\beta_1,\cdots,\beta_{j-1}\}=\operatorname{span}\{\alpha_1,\cdots,\alpha_{j-1}\}$。由构造,

$$
\beta_j=\alpha_j-\sum_{i=1}^{j-1}\frac{(\beta_i,\alpha_j)}{(\beta_i,\beta_i)}\beta_i
$$

对每个 $l<j$,

$$
(\beta_l,\beta_j)=(\beta_l,\alpha_j)-\sum_{i=1}^{j-1}\frac{(\beta_i,\alpha_j)}{(\beta_i,\beta_i)}(\beta_l,\beta_i)
=(\beta_l,\alpha_j)-\frac{(\beta_l,\alpha_j)}{(\beta_l,\beta_l)}(\beta_l,\beta_l)=0
$$

故 $\beta_j$ 与 $\beta_1,\cdots,\beta_{j-1}$ 都正交。

若 $\beta_j=0$,则 $\alpha_j\in\operatorname{span}\{\beta_1,\cdots,\beta_{j-1}\}=\operatorname{span}\{\alpha_1,\cdots,\alpha_{j-1}\}$,与 $\alpha_1,\cdots,\alpha_k$ 线性无关矛盾。因此 $\beta_j\neq0$,可以单位化。

最后,由于 $\beta_j$ 与 $\beta_1,\cdots,\beta_{j-1}$ 的系数关系是上三角的,两个张成关系同时成立。

---

**注意**:

- Schmidt 正交化的**几何意义**是:从 $\alpha_j$ 中减去它在前面各方向上的投影,剩下的部分 $\beta_j$ 自然与前面所有方向正交。
- 若 $\alpha_j$ 与前面的向量线性相关,则 $\beta_j=0$,此时只能得到正交组而不能得到标准正交组;因此该过程要求输入**线性无关**。
- **经典 Gram–Schmidt** 与**修正 Gram–Schmidt**（逐次减去投影）在数学上等价,但后者数值稳定性更好。

---

*eg*

在 $\mathbb{R}^3$ 中取标准内积,对线性无关组

$$
\alpha_1=(1,1,0)^T,\qquad \alpha_2=(1,0,1)^T,\qquad \alpha_3=(0,1,1)^T
$$

作 Schmidt 正交化。

**解**:

$$
\beta_1=\alpha_1=(1,1,0)^T,\qquad \|\beta_1\|^2=2
$$

$$
(\beta_1,\alpha_2)=1,\qquad
\beta_2=\alpha_2-\frac12\beta_1=\left(\frac12,-\frac12,1\right)^T,\qquad \|\beta_2\|^2=\frac32
$$

$$
(\beta_1,\alpha_3)=1,\qquad(\beta_2,\alpha_3)=\frac12
$$

$$
\beta_3=\alpha_3-\frac12\beta_1-\frac{1/2}{3/2}\beta_2
=\left(0,1,1\right)^T-\frac12(1,1,0)^T-\frac13\left(\frac12,-\frac12,1\right)^T
=\left(-\frac23,\frac23,\frac23\right)^T
$$

单位化得

$$
\varepsilon_1=\frac{1}{\sqrt2}(1,1,0)^T,\qquad
\varepsilon_2=\frac{1}{\sqrt6}(1,-1,2)^T,\qquad
\varepsilon_3=\frac{1}{\sqrt3}(-1,1,1)^T
$$

可以验证它们两两正交且范数均为 $1$。

### 3.8.2 QR 分解

**定理**（QR 分解）:

设 $A\in\mathbb{C}^{n\times k}$ **列满秩**（$k\le n$）,则存在分解

$$
A=QR
$$

其中 $Q\in\mathbb{C}^{n\times k}$ 的列向量构成**标准正交组**（$Q^HQ=I$）,$R\in\mathbb{C}^{k\times k}$ 是**上三角矩阵**且对角元为**正实数**。这样的分解是**唯一**的。

---

**证明**（存在性）:

设 $A=(a_1,a_2,\cdots,a_k)$,$a_1,\cdots,a_k$ 线性无关。对它们作 Schmidt 正交化,得到标准正交组 $q_1,\cdots,q_k$。由正交化过程,

$$
a_j=\sum_{i=1}^{j-1}(q_i,a_j)q_i+\|b_j\|q_j,\qquad j=1,\cdots,k
$$

其中 $b_j$ 为正交化过程中的中间向量。令

$$
r_{ij}=(q_i,a_j)\ (i<j),\qquad r_{jj}=\|b_j\|>0,\qquad r_{ij}=0\ (i>j)
$$

则 $R=(r_{ij})$ 为上三角矩阵且对角元为正,并且

$$
A=QR
$$

---

**证明**（唯一性）:

设

$$
A=Q_1R_1=Q_2R_2
$$

其中 $Q_i^HQ_i=I$,$R_i$ 为上三角且对角元为正。由列满秩知 $R_1,R_2$ 可逆,于是

$$
Q_2^HQ_1=R_2R_1^{-1}
$$

左端为酉矩阵（列标准正交且为方阵）,右端为上三角矩阵。上三角的酉矩阵必为对角矩阵,且对角元模为 $1$;而 $R_2R_1^{-1}$ 的对角元为正实数之商,仍为正实数,故只能为 $1$。因此 $R_2R_1^{-1}=I$,即 $R_1=R_2$,进而 $Q_1=Q_2$。

---

**定理**（QR 分解与 Cholesky 分解的关系）:

设 $A=QR$ 为列满秩矩阵 $A$ 的 QR 分解,则

$$
A^HA=R^HR
$$

即 $R$ 恰为 Hermite 正定矩阵 $A^HA$ 的 **Cholesky 因子**（上三角、对角元为正）。

---

**证明**:

$$
A^HA=(QR)^H(QR)=R^HQ^HQR=R^HR
$$

由 $R$ 上三角且对角元为正知它正是 $A^HA$ 的 Cholesky 分解因子;而正定矩阵的 Cholesky 分解是唯一的,这也从另一角度说明了 QR 分解的唯一性。

---

**推论**:

若 $A\in\mathbb{C}^{n\times n}$ 可逆且 $A=QR$,则

$$
|\det A|=\prod_{i=1}^{n}r_{ii}
$$

即 $|\det A|$ 等于 $R$ 的对角元之积。几何上,这是 $A$ 的列向量张成的平行体的体积。

---

*eg*

求

$$
A=\begin{pmatrix}1&1&0\\1&0&1\\0&1&1\end{pmatrix}
$$

的 QR 分解。

**解**: 由 3.8.1 的 eg 已得

$$
q_1=\frac{1}{\sqrt2}\begin{pmatrix}1\\1\\0\end{pmatrix},\quad
q_2=\frac{1}{\sqrt6}\begin{pmatrix}1\\-1\\2\end{pmatrix},\quad
q_3=\frac{1}{\sqrt3}\begin{pmatrix}-1\\1\\1\end{pmatrix}
$$

逐项计算 $r_{ij}=(q_i,a_j)$:

$$
r_{11}=(q_1,a_1)=\sqrt2,\qquad
r_{22}=(q_2,a_2)=\frac{3}{\sqrt6}=\sqrt{\frac32},\qquad
r_{33}=(q_3,a_3)=\frac{2}{\sqrt3}
$$

$$
r_{12}=(q_1,a_2)=\frac{1}{\sqrt2},\qquad
r_{13}=(q_1,a_3)=\frac{1}{\sqrt2},\qquad
r_{23}=(q_2,a_3)=\frac{1}{\sqrt6}
$$

于是

$$
Q=\begin{pmatrix}
\frac{1}{\sqrt2}&\frac{1}{\sqrt6}&-\frac{1}{\sqrt3}\\[4pt]
\frac{1}{\sqrt2}&-\frac{1}{\sqrt6}&\frac{1}{\sqrt3}\\[4pt]
0&\frac{2}{\sqrt6}&\frac{1}{\sqrt3}
\end{pmatrix},
\qquad
R=\begin{pmatrix}
\sqrt2&\frac{1}{\sqrt2}&\frac{1}{\sqrt2}\\[4pt]
0&\sqrt{\frac32}&\frac{1}{\sqrt6}\\[4pt]
0&0&\frac{2}{\sqrt3}
\end{pmatrix}
$$

验证 $A=QR$:例如第三列

$$
\frac{1}{\sqrt2}q_1+\frac{1}{\sqrt6}q_2+\frac{2}{\sqrt3}q_3
=\frac12\begin{pmatrix}1\\1\\0\end{pmatrix}+\frac16\begin{pmatrix}1\\-1\\2\end{pmatrix}+\frac23\begin{pmatrix}-1\\1\\1\end{pmatrix}
=\begin{pmatrix}0\\1\\1\end{pmatrix}=a_3
$$

再由 $A^HA=R^HR$ 可核对上三角部分的正确性。

---

### 3.8.3 QR 分解的应用

**应用 1:求解最小二乘问题**

当 $Ax=b$ 无解时,最小二乘问题为

$$
\min_{x}\|Ax-b\|_2
$$

由 3.7.3 的结论,其正规方程（法方程）为

$$
A^HAx=A^Hb
$$

把 $A=QR$ 代入,由 $A^HA=R^HR$ 得

$$
R^HRx=R^HQ^Hb
\ \Longrightarrow\
Rx=Q^Hb
$$

由于 $R$ 是可逆上三角矩阵,可直接用**回代**求解,既避免了显式构造 $A^HA$（会放大条件数）,又提高了数值稳定性。

---

**应用 2:求解线性方程组**

若 $A$ 可逆,则 $Ax=b\iff Rx=Q^Hb$,回代即可。

---

**应用 3:计算行列式与正交基**

- $|\det A|=\prod_i r_{ii}$;
- $Q$ 的列给出了 $\operatorname{range}(A)$ 的一组标准正交基;
- 由 $A=QR$ 可对矩阵进行正交三角化,是特征值算法（如 QR 迭代）的起点。

---

**注意**:

- QR 分解要求 $A$ **列满秩**以保证 $R$ 的对角元为正且分解唯一;若 $A$ 不满秩,则需借助列主元技术得到更一般的分解。
- 实际计算中常用 **Householder 变换**（见 3.6.5 eg3）或 **Givens 旋转** 来构造 QR 分解,其数值稳定性优于经典 Gram–Schmidt;修正 Gram–Schmidt 也常被使用。

---

**方法总结**:

1. **Schmidt 正交化**:先取 $\beta_1=\alpha_1$;再用 $\beta_j=\alpha_j-\sum_{i<j}\dfrac{(\beta_i,\alpha_j)}{(\beta_i,\beta_i)}\beta_i$;最后单位化。每一步都要核对 $(\beta_i,\beta_j)=0$。
2. **求 QR 分解**:对 $A$ 的列向量作 Schmidt 正交化得 $Q$,再用 $r_{ij}=(q_i,a_j)\ (i\le j)$ 写出 $R$。
3. **检验 QR 分解**:验证 $Q^HQ=I$（列标准正交）、$R$ 为上三角且对角元为正,并核对 $QR=A$ 与 $A^HA=R^HR$。
4. **求解最小二乘**:优先用 $Rx=Q^Hb$ 回代,而不是直接求 $(A^HA)^{-1}$。
