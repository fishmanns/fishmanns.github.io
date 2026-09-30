# 5 近端策略优化算法(PPO)

> callback  
> 同策略:交互与学习的策略一致  
> 异策略:交互与学习的策略不同

回忆策略梯度,使用 $n$ 次采样,然后利用这 $n$ 次采样丢进网络里对 $\theta$ 进行改进,这是典型的同策略算法,而对于两个策略下,一个进行采样,一个进行迭代,需要对策略梯度进行重要性采样改进

## 5.1 重要性采样

回顾策略梯度表示

$$
\nabla \bar{R}_\theta=\mathbb{E}_{\tau\sim p_\theta(\tau)}[R(\tau)\nabla\log p_\theta(\tau)]
$$

**问题: 假设使用同一个策略,当 $\theta\rightarrow\theta'$ 后,$p_\theta(\tau)\rightarrow p_{\theta'}(\tau)$,也就是说概率模型改变,第一次采样完全失效了**,对于传统的策略梯度,往往需要大批量的采样才能实现

假设 $x\sim p$,从分布 $p$ 中得到 $f(x)$ 的期望

$$
\mathbb{E}_{x\sim p}[f(x)]\approx\frac{1}{N}\sum_{i=1}^{N}f(x_i)
$$

现在假设使用另一个分布 $q$ 采样,而仍对原式进行计算

$$
\begin{align}
\mathbb{E}_{x\sim p}[f(x)]&=\int f(x)p(x)dx\\
&=\int f(x)\frac{p(x)}{q(x)}q(x)dx\\
&=\mathbb{E}_{x\sim q}[f(x)\frac{p(x)}{q(x)}]
\end{align}
$$

其中 $q$ 可以是任意分布,因为根据定义,这两者在期望上是严格相等的,其中 $\frac{p(x)}{q(x)} $ 称为 **重要性权重**,由此,我们使用采样分布 $q$ 成功估计了 $p$ 下的期望

尽管 $q$ 不会影响期望,但根据方差公式

$$
\mathrm{Var}[X]=\mathbb{E}[X^2]-(\mathbb{E}[X])^2
$$

可以将重要性采样方差表示为

$$
\begin{align}
\mathrm{Var}_{x\sim p}[f(x)]&=\mathbb{E}_{x\sim p}[f(x)^2]-(\mathbb{E}_{x\sim p}[f(x)])^2\\
\mathrm{Var}_{x\sim q}[f(x)\frac{p(x)}{q(x)}]&=\mathbb{E}_{x\sim q}[\big(f(x)\frac{p(x)}{q(x)}\big)^2]-\big(\mathbb{E}_{x\sim q}[f(x)\frac{p(x)}{q(x)}]\big)^2\\
&=\int f(x)^2\frac{p(x)^2}{q(x)^2}q(x)dx-\big(\int f(x)\frac{p(x)}{q(x)}q(x)dx \big)^2\\
&=\int f(x)^2\frac{p(x)}{q(x)}p(x)dx-\big(\int f(x)p(x)dx \big)^2\\
&=\mathbb{E}_{x\sim p}[f(x)^2\frac{p(x)}{q(x)}]-(\mathbb{E}_{x\sim p}[f(x)])^2\\
\end{align}
$$

简单来说,即使 $q$ 不会影响期望,但是会影响方差中的第一项,假如 $\frac{p(x)}{q(x)}$ 太大,或者说**用来采样的分布p和估计的分布q**相差太大,那么期望上升,方差变大,由于采样不可能无限多,进而会影响估计结果

## 5.2 PPO

有了重要性采样,便能将原策略梯度转换为采样下的策略梯度

$$
\begin{align}
\nabla\bar{R}_\theta=\mathbb{E}_{\tau\sim p_\theta(\tau)}[R(\tau)\nabla\log p_\theta(\tau)]\\
\Downarrow\\
\nabla\bar{R}_\theta=\mathbb{E}_{\tau\sim p_{\theta'}(\tau)}[R(\tau)\frac{p_\theta(\tau)}{p_{\theta'}(\tau)}\nabla\log p_{\theta}(\tau)]
\end{align}
$$

即保持原式不变,乘以重要性权重

这样的好处是,$\theta'$ 用于交互,$\theta'$ 与 $\theta$ 本身并没有任何关系,可以先用 $\theta'$ 进行一大批采样,然后不断训练 $\theta$,当训练到一定程度,再进行下一次采样

引入critic网络替代全局回报

$$
\begin{align}
\nabla\bar{R}_\theta=\mathbb{E}_{(s_t, a_t)\sim \pi_{\theta'}}[\frac{p_\theta(s_t,a_t)}{p_{\theta'}(s_t,a_t)}A^{\theta'}(s_t,a_t)\nabla\log p_\theta(a_t|s_t)]\\
(\nabla\log是数学推导,并不是直接的代换,详见上一章)\\
\Downarrow\\
\mathbb{E}_{(s_t, a_t)\sim \pi_{\theta'}}[\frac{p_\theta(a_t|s_t)p_\theta(s_t)}{p_{\theta'}(a_t|s_t)p_{\theta'}(s_t)}A^{\theta'}(s_t,a_t)\nabla\log p_\theta(a_t|s_t)]\\
\Downarrow\\
\mathbb{E}_{(s_t, a_t)\sim \pi_{\theta'}}[\frac{p_\theta(a_t|s_t)}{p_{\theta'}(a_t|s_t)}A^{\theta'}(s_t,a_t)\nabla\log p_\theta(a_t|s_t)]
\end{align}
$$

关键的疑问在于 $p_\theta(s_t)=p_{\theta'}(s_t)$ 是如何成立的,即怎么能直接将其约分

简单的回答就是: 算不出来,$p_\theta(s_t)$ 代表着在策略 $\pi_\theta$ 下,$s_t$ 在所有状态中出现的占比,显然想要计算这样数据是几乎无法做到的,只能*委婉地*将其忽略

通过梯度能够获取原来的优化目标

$$
\begin{align}
\nabla f(x)=f(x)\nabla\log f(x)\\
\Downarrow\\
J^{\theta'}(\theta)=\sum \frac{p_\theta(a_t|s_t)}{p_{\theta'}(a_t|s_t)}A^{\theta'}(s_t,a_t) p_{\theta'}(a_t|s_t)\\
(因为原式就是在 \theta' 下的期望)\\
\Downarrow\\
J^{\theta'}(\theta)=\mathbb{E}_{(s_t,a_t)\sim \pi_{\theta'}}[\frac{p_\theta(a_t|s_t)}{p_{\theta'}(a_t|s_t)}A^{\theta'}(s_t,a_t)]
\end{align}
$$

## 5.3 PPO 优化

在重要性采样中提到,如果分布 $p$ 与分布 $q$ 相差太大,则会产生较大方差,因此除了原本的优化目标外,还需添加一项对分布约束进行限制

$$
\begin{align}
J^{\theta'}_{PPO}(\theta)=J^{\theta'}(\theta)-\beta\mathrm{KL}(\theta,\theta')\\
J^{\theta'}(\theta)=\mathbb{E}_{(s_t,a_t)\sim \pi_{\theta'}}[\frac{p_\theta(a_t|s_t)}{p_{\theta'}(a_t|s_t)}A^{\theta'}(s_t,a_t)]
\end{align}
$$

其中 $\mathrm{KL(\theta,\theta')}$ 可认为是关于 $(\theta,\theta')$ 的函数,能够**衡量两个参数在行为上的距离而不是参数本身的差距**,也就是在动作概率分布上的差距,实际上求参数的差距并没有意义,因为每个参数值对最后的结果影响比重都是不同的

### 5.3.1 PPO-penalty

即

$$
J^{\theta'}_{PPO}(\theta)=J^{\theta'}(\theta)-\beta\mathrm{KL}(\theta,\theta')
$$

KL 项提供约束,当采样和优化策略差距过大,给一个负梯度惩罚, $\beta$ 项用于调节惩罚的权重大小

### 5.3.2 PPO-clip

PPO-penalty 的问题在于,$\beta$ 并不是一个好调的参数,且 KL 散度的计算也不简单,引入裁剪操作,省略 KL 计算过程

给出裁剪函数

$$
\mathrm{clip}(\frac{p_\theta(a_t|s_t)}{p_{\theta'}(a_t|s_t)},1-\epsilon,1+\epsilon)
$$

其意义为:

$$
\begin{align}
\mathrm{if}\hspace{1em}\frac{p_\theta(a_t|s_t)}{p_{\theta'}(a_t|s_t)}\geq1+\epsilon,\frac{p_\theta(a_t|s_t)}{p_{\theta'}(a_t|s_t)}=1+\epsilon\\
\mathrm{if}\hspace{1em}\frac{p_\theta(a_t|s_t)}{p_{\theta'}(a_t|s_t)}\leq1-\epsilon,\frac{p_\theta(a_t|s_t)}{p_{\theta'}(a_t|s_t)}=1-\epsilon
\end{align}
$$

通常取 $\epsilon=0.1/0.2$,即

$$
\mathrm{clip}(\frac{p_\theta(a_t|s_t)}{p_{\theta'}(a_t|s_t)},0.8,1.2)
$$

同时将优化目标改写为

$$
\begin{align}
J^{\theta'}_{PPO-clip}(\theta)=\mathbb{E}_{(s_t,a_t)\sim \pi_{\theta'}}[\min\Big(\\
\frac{p_\theta(a_t|s_t)}{p_{\theta'}(a_t|s_t)}A^{\theta'}(s_t,a_t),\\
\mathrm{clip}(\frac{p_\theta(a_t|s_t)}{p_{\theta'}(a_t|s_t)},1-\epsilon,1+\epsilon)A^{\theta'}(s_t,a_t)\\
\Big)]
\end{align}
$$

通俗而言

$$
\boxed{PPO-Clip将复杂的KL散度和\beta调节转换为简易的对重要性采样的截断约束,进而将两个策略捆在一起}
$$

## 5.4 关于 On-policy 和 Off-policy

从 $PPO-clip$ 中可以看出**近端策略优化**的核心表示: **迭代策略理应只在采样策略附近游动**

对于典型的 $Off-policy: Q-learning$ 而言,采样所用的策略可以天马行空,可以是随机策略,只要保证采样数够大,让每个状态和动作都有访问的机会

而对于 $PPO$ 而言,我们必须加上强硬的限制: 用于采样的策略必须和优化的策略相捆绑,所以,我们通常认为其是 **On-policy**

