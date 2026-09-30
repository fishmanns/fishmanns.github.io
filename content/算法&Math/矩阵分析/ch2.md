# Chapter2 $\lambda$矩阵与 Jordan 标准型

## 2.1 $\lambda$矩阵及标准型

### 2.1.1 l矩阵

**定义**: 设 $a_{ij}(\lambda)$ 为数域 $F(\mathbb{R}实数域/\mathbb{C}复数域/\mathbb{Q}有理数域)$ 上的多项式,称

$$
A(\lambda)=
\begin{bmatrix}
a_{11}(\lambda) & a_{12}(\lambda) & \cdots  & a_{1n}(\lambda)\\
a_{21}(\lambda) & a_{22}(\lambda) & \cdots & a_{2n}(\lambda)\\
\vdots & \vdots & \ddots & \vdots \\
a_{m1}(\lambda) & a_{m2}(\lambda) & \cdots & a_{mn}(\lambda)
\end{bmatrix}
$$

为多项式矩阵/l矩阵.**不为0的多项式 $a_{ij}(\lambda)$ 中最高的次数为 $A(\lambda) 的次数$**

---

*eg1: 非零数字矩阵为 0 次 l矩阵*

*eg2: 当 $A\in\mathbb{C}^{n\times n},\lambda I-A$ 是 1 次l矩阵*


### 2.1.2 子式

**定义**: 在 $m \times n$ 矩阵 $A$ 中，任取 $k$ 行与 $k$ 列（$1 \le k \le \min\{m,n\}$），位于这些行与列交叉处的 $k^2$ 个元素，按原来相对位置组成的 $k$ 阶行列式，称为矩阵 $A$ 的一个 $k$ 阶子式。

---

*eg1: 子式计算例子*

例：求矩阵 $A$ 的所有 2 阶子式

设

$$
A=
\begin{pmatrix}
1 & 2 & 3 \\
4 & 5 & 6
\end{pmatrix}
$$

$A$ 是 $2 \times 3$ 矩阵，所以 $k$ 最大为 $2$。

2 阶子式需要取 2 行 2 列。因为 $A$ 只有 2 行，所以行只能全取；列从 3 列中任取 2 列，共有

$$
\binom{3}{2}=3
$$

个 2 阶子式。

取第 1、2 列：

$$
\begin{vmatrix}
1 & 2 \\
4 & 5
\end{vmatrix}
=1\cdot 5-2\cdot 4=5-8=-3
$$

取第 1、3 列：

$$
\begin{vmatrix}
1 & 3 \\
4 & 6
\end{vmatrix}
=1\cdot 6-3\cdot 4=6-12=-6
$$

取第 2、3 列：

$$
\begin{vmatrix}
2 & 3 \\
5 & 6
\end{vmatrix}
=2\cdot 6-3\cdot 5=12-15=-3
$$

因为存在不为 $0$ 的 2 阶子式，例如

$$
\begin{vmatrix}
1 & 2 \\
4 & 5
\end{vmatrix}
=-3 \neq 0
$$

且没有 3 阶子式（行数不够），所以

$$
rank(A)=2
$$

### 2.1.3 正规秩

**定义**: 如果l矩阵 $A(\lambda)$ 中**有一个** r 阶 (r>=1) 子式不为 0,而所有 r+1 阶子式(如果有)全为 0,称 $A(\lambda)$ 的(正规)秩为 r ,记作:

$$
rank{A(\lambda)}=r
$$

**(r阶子式!=0只需要一个,而r+1阶子式=0需要所有)**

---

*eg1: 零矩阵的秩为 0*

*eg2: 当 $A\in\mathbb{C}^{n\times n},\lambda I-A$ 的秩为n*

### 2.1.4 l矩阵的逆矩阵

**定义**：设 $A(\lambda)$ 是一个 $n$ 阶 $\lambda$-矩阵，如果存在一个 $n$ 阶 $\lambda$-矩阵 $B(\lambda)$，使得

$$
A(\lambda)B(\lambda)=B(\lambda)A(\lambda)=I
$$

其中 $I$ 是 $n$ 阶单位矩阵，则称 $A(\lambda)$ 是**可逆的**，并称 $B(\lambda)$ 为 $A(\lambda)$ 的**逆矩阵**，记作

$$
B(\lambda)=A^{-1}(\lambda)
$$

同时可逆的l矩阵又称为单位模矩阵(幺模矩阵)

---

**定理**：$n$ 阶 $\lambda$-矩阵 $A(\lambda)$ 是一个单位模矩阵/可逆的**充分必要条件**是：$A(\lambda)$ 的行列式

$$
\det A(\lambda)=|A(\lambda)|\neq0
$$

是一个**非零常数**,或

$$
\forall\lambda,rank(A(\lambda))\equiv n
$$

---

*eg1*

$$
A(\lambda)=
\begin{pmatrix}
\lambda & 1 \\
0 & \lambda
\end{pmatrix}
$$

则

$$
|A(\lambda)|=
\begin{vmatrix}
\lambda & 1 \\
0 & \lambda
\end{vmatrix}
=\lambda^2
$$

不是非零常数，所以 $A(\lambda)$ 不可逆。

---

**注意**: 

$\lambda$-矩阵可逆与数字矩阵可逆有本质区别：

- 数字矩阵 $A$ 可逆 $\iff \det A\neq 0$
- $\lambda$-矩阵 $A(\lambda)$ 可逆 $\iff \det A(\lambda)$ 是**非零常数**

也就是说，$\det A(\lambda)$ 仅仅“不为零”还不够，必须是**与 $\lambda$ 无关的非零常数**。

### 2.1.5 初等变换

**定义**：对 $\lambda$-矩阵 $A(\lambda)$ 施行以下三种变换，称为 $\lambda$-矩阵的**初等变换**：

1. **换行（列）变换**：交换 $A(\lambda)$ 的第 $i$ 行与第 $j$ 行（或第 $i$ 列与第 $j$ 列），记作

$$
r_i \leftrightarrow r_j \qquad (\text{或 } c_i \leftrightarrow c_j)
$$

2. **倍行（列）变换**：用非零常数 $k$ 乘 $A(\lambda)$ 的第 $i$ 行（或第 $i$ 列），记作

$$
k r_i \qquad (\text{或 } k c_i)
$$

3. **消行（列）变换**：把 $A(\lambda)$ 的第 $j$ 行（或第 $j$ 列）的 $\varphi(\lambda)$ 倍加到第 $i$ 行（或第 $i$ 列）上，其中 $\varphi(\lambda)$ 是一个多项式，记作

$$
r_i + \varphi(\lambda) r_j \qquad (\text{或 } c_i + \varphi(\lambda) c_j)
$$

对单位矩阵 $I$ 施行一次初等变换所得到的矩阵，称为**初等矩阵**。

三种初等变换分别对应三种初等矩阵：

1. 交换 $I$ 的第 $i$ 行与第 $j$ 行，得

$$
P(i,j)=
\begin{pmatrix}
1 & & & & \\
& \ddots & & & \\
& & 0 & \cdots & 1 \\
& & \vdots & \ddots & \vdots \\
& & 1 & \cdots & 0 \\
& & & & \ddots & \\
& & & & & 1
\end{pmatrix}
$$

2. 用非零常数 $k$ 乘 $I$ 的第 $i$ 行，得

$$
P(i(k))=
\begin{pmatrix}
1 & & & \\
& \ddots & & \\
& & k & \\
& & & \ddots & \\
& & & & 1
\end{pmatrix}
$$

3. 把 $I$ 的第 $j$ 行的 $\varphi(\lambda)$ 倍加到第 $i$ 行，得

$$
P(i,j(\varphi))=
\begin{pmatrix}
1 & & & & \\
& \ddots & & & \\
& & 1 & \cdots & \varphi(\lambda) \\
& & & \ddots & \vdots \\
& & & & 1 \\
& & & & & \ddots & \\
& & & & & & 1
\end{pmatrix}
$$

---

**注意**：

- 初等变换需满足可逆变换,初等矩阵为单位模矩阵/可逆矩阵
- 第三种变换不能自乘 $\varphi(\lambda)$,因为初等变换需可逆.当 $r_i\rightarrow r_i+\varphi(\lambda)r_j$,可以通过 $r_i+\varphi(\lambda)r_j-\varphi(\lambda)r_j\rightarrow r_i$ 完成可逆变换.当 $r_i\rightarrow \big(1+\varphi(\lambda)\big)r_i$,无法通过 $\big(1+\varphi(\lambda)\big)r_i\rightarrow r_i$,**因为 $\frac{1}{1+\varphi(\lambda)}$ 并不是一个多项式**

### 2.1.6 等价

**定义**：对一个 $m\times n$ 的 $\lambda$-矩阵 $A(\lambda)$ 相当于用 $m$-阶初等矩阵左乘 $A(\lambda)$

对 $A(\lambda)$ 进行初等列变换相当于用 $n$-阶初等矩阵右乘 $A(\lambda)$. 

如果 $\lambda$-矩阵 $A(\lambda)$ 经过**有限次**初等变换后变成 $B(\lambda)$，则称 $A(\lambda)$ 与 $B(\lambda)$ **等价**，记作

$$
A(\lambda) \simeq B(\lambda)
$$

---

**性质**：

1. **反身性**：$A(\lambda) \simeq A(\lambda)$。

2. **对称性**：若 $A(\lambda) \simeq B(\lambda)$，则 $B(\lambda) \simeq A(\lambda)$。

3. **传递性**：若 $A(\lambda) \simeq B(\lambda)$，$B(\lambda) \simeq C(\lambda)$，则 $A(\lambda) \simeq C(\lambda)$。

---

**注意**: 等价描述的是**有限次**的初等变换

### 2.1.7 $\lambda$-矩阵的 Smith 标准型

**定义**：**任意一个**非零的 $m \times n$ 的 $\lambda$-矩阵 $A(\lambda)$ 都**等价**于一个如下形式的对角 $\lambda$-矩阵：

$$
A(\lambda)\simeq
\begin{pmatrix}
d_1(\lambda) & & & & \\
& d_2(\lambda) & & & \\
& & \ddots & & \\
& & & d_r(\lambda) & \\
& & & & 0 \\
& & & & & \ddots \\
& & & & & & 0
\end{pmatrix}
$$

其中 $r=rank A(\lambda)$，且 $d_i(\lambda)$ 都是首项系数为 $1$ 的多项式，满足

$$
d_i(\lambda) \mid d_{i+1}(\lambda) \qquad (i=1,2,\dots,r-1)\\
d_{i+1}(\lambda)=d_i(\lambda)q(\lambda)
$$

这个矩阵称为 $A(\lambda)$ 的 **Smith 标准型**

非零对角元 $d_1(\lambda), d_2(\lambda), \dots, d_r(\lambda)$ 称为 $A(\lambda)$ 的**不变因子**。

---

**注意**：

- 1 也算做不变因子中
- $d_i(\lambda) \mid d_{i+1}(\lambda)$ 表示 $d_i(\lambda)$ 整除 $d_{i+1}(\lambda)$。
- **首 1 多项式不等于左上角一定为 1!!!**,$d_1(\lambda)=\lambda^2+2\lambda$ 也是合法的

---

**定理**: 设 $A(\lambda)$ 是一个非零 $\lambda$-矩阵，$a_{11}(\lambda) \neq 0$.如果 $A(\lambda)$ 中至少有一个元素不能被 $a_{11}(\lambda)$ 整除,那么一定可以通过一系列初等变换,把 $A(\lambda)$ 化为 $B(\lambda)$,使得 $B(\lambda)$ 的左上角元素 $b_{11}(\lambda) \neq 0$,且

$$
\deg b_{11}(\lambda) < \deg a_{11}(\lambda)
$$

即**左上角元素的次数严格降低**

---

**注意**: 不能被 $a_{11}$ 整除的元素是任意的

---

**证明**: 

情况 1：同一行或同一列中有元素不能被整除

设第 $i$ 行第 $1$ 列的元素 $a_{i1}(\lambda)$ 不能被 $a_{11}(\lambda)$ 整除

做带余除法

$$
a_{i1}(\lambda) = q(\lambda) a_{11}(\lambda) + r(\lambda)
$$

其中

$$
r(\lambda) \neq 0, \qquad \deg r(\lambda) < \deg a_{11}(\lambda)
$$

于是做初等行变换

$$
r_i - q(\lambda) r_1
$$

第 $i$ 行第 $1$ 列就变成 $r(\lambda)$

再交换第 $1$ 行与第 $i$ 行，左上角元素就变成了 $r(\lambda)$，次数严格降低

情况 2：只有不同行不同列的元素不能被整除

设 $a_{ij}(\lambda)$ 不能被 $a_{11}(\lambda)$ 整除（$i \neq 1$，$j \neq 1$）。

先把第 $i$ 行加到第 $1$ 行：

$$
r_1 + r_i
$$

此时第 $1$ 行第 $j$ 列变成 $a_{1j}(\lambda) + a_{ij}(\lambda)$。

如果这个新元素仍不能被 $a_{11}(\lambda)$ 整除，就回到情况 1。

如果它能被整除，就再把第 $j$ 列加到第 $1$ 列，继续调整，最终一定能使第 $1$ 行或第 $1$ 列中出现不能被 $a_{11}(\lambda)$ 整除的元素，从而回到情况 1。

---

**性质**:

- 通过引理,每次能将左上角元素次数严格降低
- 多项式有限次下,不能无限降低次数,一定能达到**左上角的元素能整除所有元素**
- 整除时,可以轻易将第一行第一列所有元素变为 0

$$
\begin{pmatrix}
d_1(\lambda) & 0 & \cdots & 0 \\
0 & & & \\
\vdots & & A_1(\lambda) & \\
0 & & &
\end{pmatrix}
$$

- 然后对右下角的子矩阵 $A_1(\lambda)$ 重复同样的步骤,变换成为 Smith 矩阵

---

*eg1*

$$
A(\lambda)=
\begin{pmatrix}
\lambda+1 & \lambda & 0 \\
\lambda & \lambda+1 & 1 \\
0 & 1 & \lambda+1
\end{pmatrix}
$$

$R_2\leftarrow R_2-R_1,C_2\leftarrow C_2-C_1$

$$
\simeq
\begin{pmatrix}
\lambda+1 & -1 & 0 \\
-1 & 2 & 1 \\
0 & 1 & \lambda+1
\end{pmatrix}
$$

交换 $R_1,R_2$，交换 $C_1,C_2$：

$$
\simeq
\begin{pmatrix}
-1 & -1 & 2 \\
0 & \lambda+1 & -1 \\
1 & 0 & 1
\end{pmatrix}
$$

$
R_1\leftarrow-R_1,
R_3\leftarrow R_3+R_1,
$

$
C_2\leftarrow C_2-C_1,
C_3\leftarrow C_3-2C_1,
$

$$
\simeq
\begin{pmatrix}
1 & 0 & 0 \\
0 & \lambda+1 & -1 \\
0 & 1 & 1
\end{pmatrix}
$$

$
R_2\leftrightarrow R_3,
C_2\leftrightarrow C_3
$

$$
\simeq
\begin{pmatrix}
1 & 0 & 0 \\
0 & 1 & 1 \\
0 & -1 & \lambda+1
\end{pmatrix}
$$

$
C_3\leftarrow C_3-C_2,
R_3\leftarrow R_3+R_2
$

$$
\simeq
\begin{pmatrix}
1 & 0 & 0 \\
0 & 1 & 0 \\
0 & 0 & \lambda+2
\end{pmatrix}
$$

$$
\boxed{
D(\lambda)=
\begin{pmatrix}
1 & 0 & 0 \\
0 & 1 & 0 \\
0 & 0 & \lambda+2
\end{pmatrix}
}
$$

## 2.2 Smith 标准型的唯一性

### 2.2.1 行列式因子

**定义**: $A(\lambda)$ 为一个 $\lambda$-矩阵且 $rank(A(\lambda))=r$.对于任意正整数 $k,1\leq k\leq r,A(\lambda)$ 必存在非零 $k$ 阶子式, $A(\lambda)$ 的全部 $k$ 阶子式的首 1 (降幂排列)最高公因式 $D_k(\lambda)$ 称为 $A(\lambda)$ 的 **$k$ 阶行列式因子**,对于

$$
A(\lambda)\simeq
\begin{pmatrix}
d_1(\lambda) & & & & \\
& d_2(\lambda) & & & \\
& & \ddots & & \\
& & & d_r(\lambda) & \\
& & & & 0 \\
& & & & & \ddots \\
& & & & & & 0
\end{pmatrix}
$$

其各阶行列式因子

$$
\begin{align}
    D_1(\lambda)&=d_1(\lambda)\\
    D_2(\lambda)&=d_1(\lambda)d_2(\lambda)\\
    &\vdots\\
    D_r(\lambda)&=d_1(\lambda)d_2(\lambda)\cdots d_r(\lambda)
\end{align}
$$

或者:
$$
\begin{align}
d_1(\lambda)=&D_1(\lambda)\\
d_2(\lambda)=&\frac{D_2(\lambda)}{D_1(\lambda)}\\
\vdots&\\
d_r(\lambda)=&\frac{D_r(\lambda)}{D_{r-1}(\lambda)}
\end{align}
$$

---

**注意**: 

- 行列式因子一共有 $r$ 个
- 对于一般来说,流程为:先求所有 $k$ 阶子式,将各个 $k$ 阶子式化为首 1 且降幂排列,对所有 $k$ 阶子式求最大公因式
- 对于 Smith 标准型而言,$k$ 阶子式只有一种,所以就等于 $k$ 阶子式

---

**证明**:

由 Smith 标准型定义,假设 $c_i(\lambda)$ 为整除后的余项

$$
c_i(\lambda) =\frac{d_{i+1}(\lambda)}{d_{i}(\lambda)}
$$

则对于 Smith 的所有二阶子式而言

$$
d_1d_2=c_1d_1^2\\
d_1d_3=c_2c_1d_1^2\\
\cdots\\
d_1d_r=c_{r-1}\cdots c_1 d_1^2\\
$$

同理还有

$$
d_2d_3 = c_2c_1^2d_1^2\\
d_id_j = \prod_{t=1}^{i-1}c_t\prod_{q=1}^{j-1}c_qd_1^2
$$

由于要选最大公因式,所以选择 $c_1d_1^2=d_2d_1$ 作为行列式因子而不是 $d_1^2$
其他因子可以类比得到

---

**定理**: 初等变换不改变 $\lambda$-矩阵的行列式因子,因此等价矩阵拥有相同行列式因子,因此具有相同的秩(**秩=行列式因子个数=不变因子个数)**

### 2.2.2 Smith 标准型的唯一性

**定理**: $A(\lambda)$ 的 Smith 标准型是唯一的

---

**证明**: 

假设 $A(\lambda)$ 具有 Smith 标准型 $S_1(\lambda), S_2(\lambda)$

由 Smith 标准型定义,$A(\lambda)\simeq S_1(\lambda), A(\lambda)\simeq S_2(\lambda)$

初等变换不改变行列式因子,所以 $S_1,S_2$ 行列式因子一致

不变因子可由行列式因子导出,所以 $S_1,S_2$ 不变因子一致

$S_1=S_2$, Smith 标准型唯一

---

**定理**: 

$m\times n$ 维 $\lambda$ 矩阵 $A(\lambda)\simeq B(\lambda)$:

$$
\Leftrightarrow 1.对于任意的 k,其 k 阶行列式因子一致
\Leftrightarrow 2.
$$