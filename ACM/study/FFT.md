# FFT

## 简介

FFT( ~~Fast Fast TLE~~ Fast Fourier Transformation) 是一种用来计算向量叉积 / 卷积的工具。在**高精度乘法**和图像处理应用广泛

## 前置芝士

### 卷积(Convolution)

卷积的定义式

$$
g(x) = \int_{-\infty}^{+\infty} f(x - t)\cdot g(t) \mathop{}\!\mathrm{d} t
$$

在一个数组中，我们这样计算

```cpp
vector Convolution(vector a, vector b) {
    vector c
    for(int i = 0; i < n; ++i)
        for(int j = ; j < m; ++j)
            c[i + j] += a[i] * b[i]
    return c
}
```

容易看出这是一个 $O(n^2)$ 的过程

### 多项式

我们把一个 $n$ 维向量与一个 $n - 1$ 次多项式对应

$$
f(x) = a_0 + a_1x + a_2x^2 + \cdots + a_{n-1}x^{n-1}
$$

显然，$a\times b$ ($\times$ 表示卷积) 的结果就是 $(f\cdot g)(x)$ 再转换成向量

这里提出多项式正是为了利用多项式的一些性质来计算卷积

#### 多项式的性质

##### 性质1

$\partial\degree(fg) = \partial\degree(f) + \partial\degree(g)$，其中 $\partial\degree(f)$ 表示多项式的次数

##### 性质2

$n$ 个不同的点可以确定一个 $n - 1$ 次的多项式

### 总体思路

有了多项式，计算卷积的过程就变成了:

根据性质 $1$ 知道 $f\cdot g$ 的次数为 $(n - 1) + (m - 1) = n + m - 2$，所以我们可以通过找到 $n + m - 1$ 个不同的点，确定多项式，进而得到卷积

所以我们取 $x_0, x_1, \cdots, x_{n + m - 2}$，分别计算出 $f(x_0)\cdots f(x_{n + m - 2})$ 和 $g(x_0)\cdots g(x_{n + m - 2})$。进而计算出 $fg(x_0)\cdots fg(x_{n + m - 2})$，即找到了 $fg(x)$ 上的 $n  +m- 1$ 个点

显然，如果随意取 $n + m - 1$ 个点，那么时间复杂度还是会变成 $O(n^2)$，所以我们需要找到一些具有特殊性质的点帮助我们提高效率。这是 FFT 算法最屌的地方，这里寻找的点是在复平面的单位圆上的。

### 复数

记 $\omega_n^k = e^{\small{i\frac{2k\pi}{n}}}$，如果 $n=2^t, t\in \N_+$，复数有以下性质:

$$
\begin{aligned}
&(\omega_n^k)^t = \omega_n^{kt}\\\\
&\omega_n^{k+n} = \omega_n^k\\\\
&\omega_n^{k + \frac{n}{2}} = -\omega_n^k\\\\
&\omega_{dn}^{dk} = \omega_n^k
\end{aligned}
$$

看起来证明都不难，留给读者了

## 正式过程

观前提醒，由于上述的复数性质提到: $n=2^t, t\in \N_+$，所以以下的多项式 / 数组都是由 $2^t$ 项组成

### 将系数转化成点

对于一个 $n - 1$ 次多项式 $f(x)$，需要找 $n$ 个点 ($n = 2^t, t \in \N_+$)。取 $x_k = w_n^k, k = 0, 1, 2\cdots, n - 1$.

这意味着我们将在计算后得到 $(w_n^0, f(w_n^0)), (w_n^1, f(w_n^1)), \cdots(w_n^{n - 1}, f(w_n^{n - 1}))$ 的值。

下面考虑在 $O(\log n)$ 的时间复杂度内计算 $f(x_k)$

按照次数的奇偶分类，一定存在 $f_1(x)$ 和 $f_2(x)$ 满足

$$
f(x) = f_1(x^2) + xf_2(x^2)
$$

其中 $f_1(x)$ 与 $f_2(x)$ 是 $2^{t - 1}$ 次多项式

$\forall x_k, 0\leq k < \frac{n}{2}$ 有:

$$
\begin{aligned}
f(x_k) &= f_1(x_k^2) + x_kf_2(x_k^2)\\
       &= f_1(\omega_n^{2k}) + \omega_n^kf_2(\omega_n^{2k})\\
       &= f_1(\omega_{\frac{n}{2}}^k) + \omega_n^kf_2(\omega_{\frac{n}{2}}^k)
\end{aligned}
$$

$$
\begin{aligned}
f(x_{k + \frac{n}{2}}) &= f_1(x_{k + \frac{n}{2}}^2) + x_{k + \frac{n}{2}}f_2(x_{k + \frac{n}{2}}^2)\\
       &= f_1(\omega_n^{2k+n}) + \omega_n^{k + \frac{n}{2}}f_2(\omega_n^{2k+n})\\
       &= f_1(\omega_n^{2k}) - \omega_n^kf_2(\omega_n^{2k})\\
       &= f_1(\omega_{\frac{n}{2}}^k) - \omega_n^kf_2(\omega_{\frac{n}{2}}^k)
\end{aligned}
$$

由于 $k < \frac{n}{2}, \partial\degree(f_1) = \partial\degree(f_2) = \frac{n}{2}$，$f_2(\omega_{\frac{n}{2}}^k)$ 与 $f_1(\omega_{\frac{n}{2}}^k)$ 的值可以通过递归计算 $f_1$ 与 $f_2$ 来得出。

由于合并的时间复杂度为 $O(1)$，所以计算 $f(x_k)$ 的时间复杂度为递归的深度 $O(\log n)$，类似于快速幂

#### 伪代码

```cpp
void fft(int n, complex* f) {
    if(n == 1) return;
    complex f1[n >> 1], f2[n >> 1];
    for(int i = 0; i < n; i += 2) {
        f1[i >> 1] = f[i];
        f2[i >> 1] = f[i | 1];
    }
    fft(n >> 1, f1), fft(n >> 1, f2);
    complex w(1, 0), wn(cos(2pi/n), sin(2pi/n));
    for(int k = 0; k < (n >> 1); ++k, w = w * wn) {
        f[k] = f1[k] + w * f2[k];
        f[k + (n >> 1)] = f1[j] - w * f2[k];
    }
}
```

### 将点转化为系数

注意到，卷积后的多项式 $fg(x) = a_0 + a_1x + \cdots + a_{N-1}x^{N-1}, N = n + m - 2$ 满足:

为了看起来舒服，下面将 $N$ 用 $n$ 表示

$$
\left[
\begin{matrix}
        \omega_n^0 & \omega_n^0 & \omega_n^0 & \cdots & \omega_n^0 \\
        \omega_n^0 & \omega_n^1 & \omega_n^2 & \cdots & \omega_n^{n-1} \\
        \omega_n^0 & \omega_n^2 & \omega_n^4 & \cdots & \omega_n^{2(n-1)} \\
        \vdots     & \vdots     & \vdots     & \ddots & \vdots\\
        \omega_n^0 & \omega_n^{n-1} & \omega_n^{2(n-1)} & \cdots & \omega_n^{(n-1)(n-1)} \\
\end{matrix}
\right]
\left[
\begin{matrix}
    a_0\\ a_1\\ a_2\\ \vdots\\ a_{n-1}
\end{matrix}
\right] =
\left[
\begin{matrix}
    y_0\\ y_1\\ y_2\\ \vdots\\ y_{n-1}
\end{matrix}
\right]
$$

显然，这是一个范德蒙德矩阵，一定有它的逆矩阵 $M^{-1}$，我们需要做的就是求出 $A = M^{-1}Y$

再次注意到:

$$
M^{-1}=\frac{1}{n}
\left[
\begin{matrix}
        \omega_n^0 & \omega_n^0 & \omega_n^0 & \cdots & \omega_n^0 \\
        \omega_n^0 & \omega_n^{-1} & \omega_n^{-2} & \cdots & \omega_n^{-(n-1)} \\
        \omega_n^0 & \omega_n^{-2} & \omega_n^{-4} & \cdots & \omega_n^{-2(n-1)} \\
        \vdots     & \vdots     & \vdots     & \ddots & \vdots\\
        \omega_n^0 & \omega_n^{-(n-1)} & \omega_n^{-2(n-1)} & \cdots & \omega_n^{-(n-1)(n-1)} \\
\end{matrix}
\right]
$$

验证交给读者，于是有:

$$
\begin{aligned}
    a_k &= \frac{1}{n}\sum_{i = 0}^{n-1} y_i\omega_n^{-ik}\\
        &= \frac{1}{n}\sum_{i = 0}^{n-1} y_i(\omega_n^{-k})^i\\
        &= \frac{1}{n}h(\omega_n^{-k})
\end{aligned}
$$

这条式子非常眼熟，因为上面我们求的就是 $f(w_n^k)$，现在只差一个负号，用相同的方法解决就行了。

## 优化

可以发现，FFT的过程中出现了许多的递归操作，这种做法对算法的复杂度危害极大。可以将其优化成不使用递归的方式

先按照分治的流程操作观察一下递归的迭代过程

![Nanako7_ix](image1.png)

然后计算过程如下，对应上方代码 FFT 递归后的计算（计算过程与中的 $a_i$ 下标无关，仅展示树的结构）

### 计算向上迭代

假设我们的数组已经经过交换变成了 $[a_0, a_4, a_2, a_6, a_1, a_5, a_3, a_7]$ 的顺序，那么我们需要做的就是向上计算出结果。

<img src="image2.png" alt="Nanako7_ix" style="zoom: 33%;" />

先枚举子区间的长度，然后枚举父区间左端点的位置，最后枚举父区间的每一个元素进行更换。观察迭代的方向:

- 旧数组的第 $1$ 和第 $3$ 个元素结合，计算出第 $1$ 和第 $3$ 个元素
- 旧数组的第 $2$ 和第 $4$ 个元素结合，计算出第 $2$ 和第 $4$ 个元素

你还会发现，这甚至可以滚动数组!? 按照这种规律，就可以写出以下代码:

- 枚举子区间的长度
- 枚举父区间的左端点
- 枚举子区间中的每一个位置，用来计算父区间的值

```cpp
for(int i = 1; i < n; i <<= 1) {
    complex<double> wn(cos(pi/i), sgn * sin(pi/i));
    for(int j = 0; j < n; j += (i << 1)) {
        complex<double> w(1, 0);
        for(int k = 0; k < i; ++k, w = w * wn) {
            auto x = f[j + k], y = w * f[j + i + k];
            f[j + k] = x + y, f[j + i + k] = x - y;
        }
    }
}
```

### 初始改变顺序

然后 FFT 递归前的操作是改变 `f[]` 的顺序。从上文的图片可以看出，数组从 $[a_0, a_1, a_2, a_3, a_4, a_5, a_6, a_7]$ 变成了 $[a_0, a_4, a_2, a_6, a_1, a_5, a_3, a_7]$

下面研究如何将初始数组转化成后面这个看起来没有什么规律的数组

可以发现，交换的两个数的二进制是 `reverse` 的关系

$$
0 \leftrightarrow 0 : 000 \leftrightarrow 000\\
1 \leftrightarrow 4 : 001 \leftrightarrow 100\\
2 \leftrightarrow 2 : 010 \leftrightarrow 010\\
3 \leftrightarrow 6 : 011 \leftrightarrow 110\\
4 \leftrightarrow 1 : 100 \leftrightarrow 001\\
5 \leftrightarrow 5 : 101 \leftrightarrow 101\\
6 \leftrightarrow 3 : 110 \leftrightarrow 011\\
7 \leftrightarrow 7 : 111 \leftrightarrow 111
$$

所以开始的时候需要预处理出一个 `rev[]` 用来表示数字与 $i$ 对应的数字，计算方法如下，其中 $t$ 表示 $n = 2^t$:

```cpp
rev[i] = (rev[i >> 1] >> 1) | ((i & 1) << (t - 1))
```

每次在 FFT 之前，先对原来的数组进行交换即可

### 利用虚数域上多项式环 $\mathbb{C}[x]$ 优化

上面的情况需要两个多项式做操作，也需要三次 FFT。如果令 $F(x) = f(x) + ig(x)$，那么 $F^2(x) = f^2(x) - g^2(x) + 2ifg(x)$。容易发现虚部多项式正是卷积结果的两倍。所以我们只需要用 $F(x)$ 和自己做一次卷积，然后取卷积结果的虚部的一半，就是 $f(x)$ 和 $g(x)$ 的卷积。这样只需要做两次 FFT，优化还是很明显的。

## 收尾工作及细节提醒

经过以下操作后，$f$ 和 $g$ 卷积后的多项式 $fg$ 已经得到，它的系数储存在 `vector<complex<double>> a` 中了

如果不出错的话，$a_i$ 的虚部应该是 $0$ 或因为精度非常接近于 $0$，所以只需要对实部**除以 $N$** （因为 $a_k = \frac{1}{n}h(\omega_n^{-k})$）后四舍五入即可四舍五入的操作是 `ans[i] = LL(a[i].real() / N + 0.5);`

如果（大部分情况）是进行高精度乘法，还需要对 `ans[]` 进行进一步处理，也就是进位与处理 $N$ 的值

## 代码

[洛谷P1919](https://www.luogu.com.cn/problem/P1919)

```cpp
inline int Convolution(const string& a, const string& b, vector<LL>& res) {
    constexpr long double pi = acos(-1.0);
    int n = a.size(), m = b.size();
    int N = 1, t = 0;
    while(N < n + m) N <<= 1, ++t;
    vector<complex<double>> f(N);
    for(int i = 0; i < n; ++i) f[i].real(a[n - i - 1] ^ 48);
    for(int i = 0; i < m; ++i) f[i].imag(b[m - i - 1] ^ 48);

    vector<int> rev(N);
    for(int i = 0; i < N; ++i)
        rev[i] = (rev[i >> 1] >> 1) | ((i & 1) << (t - 1));
    auto fft = [&](int sgn) {
        for(int i = 0; i < N; ++i)
            if(i < rev[i]) swap(f[i], f[rev[i]]);
        for(int i = 1; i < N; i <<= 1) {
            complex<double> wn(cos(pi / i), sgn * sin(pi / i));
            for(int j = 0; j < N; j += (i << 1)) {
                complex<double> w(1, 0);
                for(int k = 0; k < i; ++k, w *= wn) {
                    auto x = f[j + k], y = w * f[i + j + k];
                    f[j + k] = x + y, f[i + j + k] = x - y;
                }
            }
        }
    };
    
    fft(1);
    for(int i = 0; i < N; ++i) f[i] *= f[i];
    fft(-1);
    res.resize(N + 3);
    for(int i = 0; i < N; ++i) {
        res[i] = LL(f[i].imag() / (2 * N) + 0.5);
    }
    return N;
}

void solve() {
    string s1, s2;
    cin >> s1 >> s2;
    vector<LL> ans;
    int N = Convolution(s1, s2, ans);
    for(int i = 0; i <= N; ++i) {
        if(ans[i] >= 10) {
            ans[i + 1] += ans[i] / 10;
            ans[i] %= 10;
            N += (i == N);
        }
    }
    while(!ans[N] && N >= 1) N--;
    N++;
    while(--N >= 0) cout << ans[N];
}
```

## Ending!

至此 FFT 就结束了，数学含量挺高的。内容看起来很难但其实很多都是计算过程，这里把不影响阅读的证明都删除了 (~~因为不会~~)。
