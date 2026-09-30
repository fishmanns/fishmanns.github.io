# 6 深度 Q 网络(DQN)

为了解决表格型办法无法处理的连续状态和动作空间的问题,通常采用神经网络对值函数进行拟合,即值函数近似

$$
Q_\theta(s,a)\approx Q_\pi(s,a)
$$

将使用网络估计值函数的方法称为**评论员**

## 6.1 估计状态价值函数

两种办法:

1.基于蒙特卡洛,让

$$
V_\theta(s_t)\rightarrow G_t
$$

2.基于时序差分,让

$$
V_\theta(s_t)\rightarrow r_{t+1}+\gamma V_\theta(s_{t+1})
$$

当然和表格型方法一样,由估计状态价值再计算 $Q(s_t,a_t)$,不如直接估计 $Q(s_t,a_t)$

## 6.2 估计动作价值函数

网络构造上有不同方法  

1. $(s_t,a_t)\rightarrow Q_\theta\rightarrow Q_\theta(s_t,a_t)$,输入状态和动作,直接给出对应 $Q$ 值,**适用于连续状态空间和连续动作空间**
2. $s_t\rightarrow Q_\theta\rightarrow Q_\theta(s_t,a_1)/Q_\theta(s_t,a_2)\cdots$,输入状态,输出是各个动作下的 $Q$,**只适用于离散的动作空间,因为输出没法做成连续的**

将谁作为更新目标也有两种方案

1. 类 SARSA: $Q_\theta(s_t,a_t)\rightarrow r_{t+1}+\gamma Q_\theta(s_{t+1},a_{t+1})$,使用采样时的 $a_{t+1}$
2. 类 Q-learning: $Q_\theta(s_t,a_t)\rightarrow r_{t+1}+\gamma\max\limits_{a_{t+1}}Q_\theta(s_{t+1},a_{t+1})$,使用贪心策略中的 $a_{t+1}$

在更新一定步骤后,拥有优化策略:

$$
\pi^*(a_t|s_t) = 
\begin{cases}
    1,a_t=\operatorname*{arg\,max}\limits_{a^*_t}Q_\theta(s_t,a^*_t)
    \\
    0, a_t\neq\operatorname*{arg\,max}\limits_{a^*_t}Q_\theta(s_t,a^*_t)
\end{cases}
$$

当然,在探索时,为了保留探索性,必须加上 $\epsilon$-贪心策略

$$
a=
\begin{cases}
    \operatorname*{arg\,max}\limits_{a_t}Q_\theta(s_t,a_t), 有 1-\epsilon 的概率  
    \\
    \mathrm{random}, \mathrm{else}
\end{cases}
$$

**$\epsilon$-贪心只用于探索和训练,假如训练完成,进行实验应当使用纯贪心**

一般情况下,认为 DQN 使用类 Q-learning 进行迭代,作为标准的 off-policy 方法

## 6.3 策略迭代

在表格型方法中,不需要显式地策略迭代: $Q(s_t,a_t)$ 直接通过 TD-target 更新,策略通过 $\operatorname*{arg\,max} + \epsilon-贪心$ 直接更改

然而在 DQN 中,我们必须显式迭代,原因在于我们无法直接接触 $Q_\theta(s_t,a_t)$,我们只能通过修正网络参数 $\theta$,让 $Q_\theta(s_t,a_t)$ 趋近于 TD-target

一般使用 MSE 计算 Loss

$$
\boxed{L(\theta) = \frac{1}{N}\sum_i^N[r_{i+1}+\gamma\max\limits_{a_{i+1}}Q_\theta(s_{i+1},a_{i+1})-Q_\theta(s_i,a_i)]^2
}
$$

通过反响传播更新 $\theta$

$$
\boxed{\theta\leftarrow\theta-\alpha\nabla_\theta L(\theta)}
$$

## 6.4 经验回放

经验回放构建一个回放缓冲区(replay buffer),由于异策略,所以使用采样时来自各个策略的经验 $(s_t,a_t,r_{t+1},s_{t+1})$ 是完全可行的

在训练时,从 replay buffer 中挑选出一个批量(batch)出来,进行更新

## 6.5 目标网络(target net)

假设 TD-target

$$
Q_\theta(s_t,a_t)\rightarrow r_{t+1}+\gamma\max\limits_{a_{t+1}}Q_\theta(s_{t+1},a_{t+1})
$$

那么发生 $Q_\theta(s_t,a_t)趋近\rightarrow Q_\theta(s_{t+1},a_{t+1})\rightarrow\theta 改变\rightarrow Q_\theta(s_{t+1},a_{t+1}) 改变\rightarrow Q_\theta(s_t,a_t)趋近\rightarrow\cdots$

这样自锁的迭代必然不太稳定,一般做法是引入**目标网络**将迭代与采样隔离开

$$
Q_\theta: 在线网络,用于梯度下降迭代策略\\
Q_{\theta^-}: 目标网络,用于计算 TD-target
$$

初始化时

$$
\theta=\theta^-
$$

在接下来的迭代中,相当于固定住目标网络,让在线网络趋近于目标网络

$$
Q_\theta(s_t,a_t)\rightarrow r_{t+1}+\gamma\max\limits_{a_{t+1}}Q_{\theta^-}(s_{t+1},a_{t+1})
$$

当迭代一定次数后,将在线网络 copy 至目标网络,防止二者差距太大

$$
\mathrm{when\space iter==n}:\theta^-=\theta
$$

## 6.6 DQN 算法

给出 suedocode 如下

$$
\begin{aligned}
&\textbf{Initialize } Q_\theta,\; Q_{\theta^-},\; \mathcal D,
\qquad \theta^- \leftarrow \theta;\\
&\textbf{for each episode do}\\
&\quad \text{Initialize }s_0;\\
&\quad \textbf{for }t=0,1,\ldots,T\textbf{ do}\\
&\qquad
a_t\sim\epsilon\text{-greedy}\left(Q_\theta(s_t,\cdot)\right);\\
&\qquad
\text{Execute }a_t,\text{ obtain }r_t,s_{t+1},d_t;\\
&\qquad
\mathcal D\leftarrow
\mathcal D\cup(s_t,a_t,r_t,s_{t+1},d_t);\\
&\qquad
\text{Sample }(s_i,a_i,r_i,s'_i,d_i)_{i=1}^{N}
\sim\mathcal D;\\
&\qquad
y_i=
r_i+\gamma(1-d_i)
\max_{a'}Q_{\theta^-}(s'_i,a');\\
&\qquad
L(\theta)=
\frac{1}{N}
\sum_{i=1}^{N}
\left(y_i-Q_\theta(s_i,a_i)\right)^2;\\
&\qquad
\theta\leftarrow
\theta-\alpha\nabla_\theta L(\theta);\\
&\qquad
\textbf{if }t\bmod C=0:
\quad\theta^-\leftarrow\theta;\\
&\quad \textbf{end for}\\
&\textbf{end for}
\end{aligned}
$$

## 6.7 优化

### 6.7.1 双深度 Q 网络(DDQN)

DQN 存在一个问题: $Q_\theta(s_t,a_t)$ 一般情况下会被高估

其原因在于,每一个迭代 

$$
Q_\theta(s_t,a_t)\rightarrow r_{t+1}+\gamma\max\limits_{a_{t+1}}Q_{\theta^-}(s_{t+1},a_{t+1})
$$

而网络预测必然存在误差,假设存在三个动作

$$
Q_{\theta^-}(s_2,a_1)=10\\
Q_{\theta^-}(s_2,a_2)=10\\
Q_{\theta^-}(s_2,a_3)=10
$$

他们 $Q$ 值相等, 在 $Q_\theta(s_1,a)$ 的迭代中应当有同等效力,但考虑误差

$$
Q_{\theta^-}(s_2,a_1)=10.001\\
Q_{\theta^-}(s_2,a_2)=9\\
Q_{\theta^-}(s_2,a_3)=9.99
$$

即使一点偏差也会使得 $(s_1,a)\rightarrow(s_2,a_1)$ 更具有依赖性

常见方法是**将选择动作和计算 Q 值的网络隔离**,即一个选择动作,另一个计算 Q 值

$$
Q_\theta(s_t,a_t)\rightarrow r_{t+1}+\gamma\max\limits_{a_{t+1}}Q_{\theta^-}(s_{t+1},a_{t+1})\\
\Downarrow\\
a_{t+1}=\operatorname*{arg\,max}\limits_{a}Q_\theta(s_{t+1},a)\\
Q_\theta(s_t,a_t)\rightarrow r_{t+1}+\gamma Q_{\theta^-}(s_{t+1},a_{t+1})\\
$$

### 6.7.2 竞争 Q 网络(Dueling DQN)

引入优势函数 $A(s_t,a_t)$ 概念,将原先的网络由 

$$
s_t\rightarrow \theta\rightarrow Q(s_t,a_t)
$$ 

改造为 

$$
s_t\rightarrow\theta\rightarrow V(s_t)+A(s_t,a_t)
$$

其中

$$
A(s_t,a_t)+V(s_t)=Q(s_t,a_t)
$$

最后合成 $Q$

主要体现思想为,既然某个比较好的状态 $s_1$ 得分比较高,那么对于其类似的状态 $s_t$,我们有理由给他一个较好的平均分,而 $A(s_t,a_t)$ 用于表示相对于这个平均分,哪个动作更好 

### 6.7.3 优先级经验回放(PER)

为每个经验赋予权重,让那些不方便采样的样本能够大概率被选上

### 6.7.4 噪声网络(Noise Net)

在 $\epsilon$-贪心中,我们通过随机概率完成自主探索,同样也有一种为参数添加噪声的方法

假定网络参数 $\theta$,每回合开始时,为其添加噪声 $\theta+\varepsilon\rightarrow\tilde{\theta}$

然后在**该回合中,保持 $\tilde{\theta}$ 不变**进行游戏,每次采样选择

$$
a=\operatorname*{arg\,max}\limits_{a_t}\tilde{Q}_{\tilde{\theta}}(s_t,a_t)
$$

这样做的好处是,对于 $\epsilon$-贪心而言,所谓的探索只是概率游戏,我们有概率选择最好的,或者有概率随机乱走,假设另一条采样又来到这个状态,其动作不能确定.而对于噪声网络而言,能保证当前回合,参数不变的情况下,每次来到这个状态,都能复现上一次的结果,然后下个回合我们再添加不同的噪声,达到不同的结果

### 6.7.5 彩虹(rainbow)

将以上多个模块进行耦合拼凑,称为 $DQN_{promax}$ 或者 $rainbow$