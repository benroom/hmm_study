# データ作成
n_1 <- 500
n_2 <- 500
n_3 <- 500
x <- c(rpois(n_1, 3), rpois(n_2, 10), rpois(n_3, 20))
x <- c(1, 1, 1, 5, 5, 5, 10, 10, 10)

# 想定するクラス
K <- 3

# 長さ(書籍に倣いNで書く)
N <- length(x)

# prior/posteriorの型
## 1.目的変数の分布に関わるもの
## λ~Ga(a,b)
## y~po(λ)
## ガンマ分布のパラメータ(ポアソン分布の事前分布)
## 各潜在クラスの次元ごとのポアソン分布の平均
a <- a_update <- rep(2, K)
b <- b_update <- rep(0.5, K)

## 2.初期状態に関わるもの
## ディリクレ分布のパラメータ(カテゴリ分布の事前分布)
## aはn=1時点での状態：K次元(潜在クラスの数)
## π~Dir(α)
## s1~Cat(π)
alpha <- rep(20, K) + c(10, 0, -10)
alpha <- alpha_update <- rep(20, K) + c(10, 0,-10)

## 3.遷移行列に関わるもの
## betaはn>1時点でのAのパラメータ：K×K次元(遷移行列)
## A_{i,j}=p(s_{n-1}=j →s_{n}=i)
## p(A)=Dir(A|β)
beta <- matrix(10, nrow = K, ncol = K) + diag(10, K)
beta <- beta_update <- matrix(10, nrow = K, ncol = K) + diag(10, K)

# このあと計算に使うもの
## :q(s)どのクラスにいるか
s_expectedValue <- matrix(0, nrow = N, ncol = K)
## p(x_n|s_n)の配列、N行K列
p_xn_sn_all <- matrix(0, nrow = N, ncol = K)
## f(s)の配列、N行K列
p_forward <- matrix(0, nrow = N, ncol = K)
## b(s)の配列、N行K列
p_backward <- matrix(0, nrow = N, ncol = K)

# 推定
iter_max <- 200
for (i in 1:iter_max) {
  # 1. 更新式に使用するパラメータの期待値を計算する
  # # 求めたい値と型(勝手に宣言されるので初期化しなくてOK)
  # #〈λ〉
  # lambda_expectedValue <- rep(0, K)
  # #〈lnλ〉
  # ln_lambda_expectedValue <- rep(0, K)
  # #〈lnπ_i〉
  # ln_pi_expectedValue <- rep(0, K)
  # ##〈lnA_{j,i}〉
  # ln_A_expectedValue <- matrix(0, nrow = K , ncol = K)
  
  ## ポアソン分布
  ## (5.105)<λ>
  lambda_expectedValue <- a_update / b_update
  ## (5.106)<lnλ>
  ln_lambda_expectedValue <- digamma(a_update) - log(b_update)
  ## カテゴリ分布
  ## (5.107)<lnπ>
  ln_pi_expectedValue <-
    digamma(alpha_update) - digamma(sum(alpha_update))
  ## (5.108)<lnA>
  ## 列に対して正規化
  ln_A_expectedValue <-
    t(t(digamma(beta_update)) - digamma(colSums(beta_update)))
  
  ## 実際に使用する遷移行列A=p'(s_n|s_{n-1}) （※時点ごとに変化しない)
  ## (5.115)→(5.74)式
  p_sn_sn_1 <- exp(ln_A_expectedValue)
  p_sn_sn_1 <- t(t(p_sn_sn_1) / colSums(p_sn_sn_1))
  
  # 2.forward/backwardで周辺確率を計算
  ## 初期化
  p_forward <- matrix(0, nrow = N, ncol = K)
  ##forwardループ
  for (n in 1:N) {
    # 正規化して更新する
    #　時点nでのp'(x_n|s_n)
    ## 以降、書くのがめんどくさいのでチルダを'で書く
    ##p'(x_n|s_n):(5.116)→(5.95)式
    ln_xn_sn <-
      x[n] * ln_lambda_expectedValue - lambda_expectedValue
    ## 正規化して計算を求める(5.116)式のexpを計算
    p_xn_sn <- exp(ln_xn_sn)
    ## 総和で割って確率に直す
    p_xn_sn <- p_xn_sn / sum(p_xn_sn)
    ## n時点目に代入
    p_xn_sn_all[n, ] <- p_xn_sn
    
    ## (5.122)の計算
    if (n == 1) {
      ## n=1時点での遷移確率
      p_s1 <- exp(ln_pi_expectedValue)
      p_s1 <- p_s1 / sum(p_s1)
      ## p'(x_1|s_1)p'(s_1)
      p_forward[n, ] <- p_xn_sn * p_s1
    } else{
      ## Σ_{s_{n-1}}p'(s_n|s_{n-1})f(s_{n-1})
      p_forward[n, ] <-
        p_xn_sn * rowSums(t(p_forward[n - 1, ] * t(p_sn_sn_1)))
    }
    ##時点ごとに確率となるよう行方向に正規化する
    p_forward[n, ] = p_forward[n, ] / sum(p_forward[n, ])
  }
  
  ##初期化
  p_backward <- matrix(0, nrow = N, ncol = K)
  ##backwardループ
  for (n in N:1) {
    if (n == N) {
      p_backward[n, ] <- rep(1, K)
    } else{
      ## 変数名使いまわしているけど実際はn+1時点のp'(x_{n+1}|s_{n+1})
      ln_xn_sn <-
        x[n + 1] * ln_lambda_expectedValue - lambda_expectedValue
      p_xn_sn <- exp(ln_xn_sn)
      p_xn_sn <- p_xn_sn / sum(p_xn_sn)
      ## Σ_{s_{n+1}}p'(x_{n+1}|s_{n+1})p'(s_n|s_n-1)b(s_{n+1})
      p_backward[n, ] <-
        colSums(p_backward[n + 1, ] * p_sn_sn_1 * p_xn_sn)
    }
    ## 正規化
    p_backward[n, ] <-
      p_backward[n, ] / sum(p_backward[n,])
  }
  
  # q(sn):=f(sn)b(sn)
  s_expectedValue <- p_forward * p_backward
  ##行ごとに和が1になるように正規化
  s_expectedValue <- s_expectedValue / rowSums(s_expectedValue)
  
  # 3.パラメータ更新
  ## 潜在クラスごとのポアソン分布の平均パラメータ更新
  ##(5.86)でaを更新
  a_update <- colSums(x * s_expectedValue) + a
  ##(5.86)でbを更新
  b_update <- colSums(s_expectedValue) + b
  
  ## n=1の初期位置に関するパラメータ更新
  ## (5.90)でαを更新
  alpha_update <- s_expectedValue[1, ] + alpha
  
  ## n>1に関するパラメータ更新
  ##(5.93)でβを更新
  ##前処理：〈s_{n-1,i}s_{n,j}〉∝f_n b_n’ ・ p(s_n|s_{n-1}) p(x_n|s_n)
  outer_temp <- matrix(0, nrow = K , ncol = K)
  for (n in 1:(N - 1)) {
    temp <-
      matrix(p_backward[n, ], nrow = K, ncol = 1) %*%
      matrix(p_forward[n + 1, ], nrow = 1, ncol = K) * p_sn_sn_1 * p_xn_sn_all[n + 1]
    temp <- t(t(temp) / colSums(temp))
    outer_temp <- outer_temp + temp + beta
  }
  # β更新
  beta_update <- outer_temp
}