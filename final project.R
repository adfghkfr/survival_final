rm(list = ls(all = TRUE))
data <- read.table("/Users/guoyixuan/Documents/北大統研/存活分析/final/DATA_4.txt", head = TRUE)
head(data)
library(survival)
# X1-發生肝纖維化(肝硬化)的時間(時間一般很長) X2:發生肝癌的時間(survival time一般只有5.6年)
# 一般先發生肝纖維化再發生肝癌(反過來很少見，而且肝癌survival time較短，常常等不到發生肝硬化)
# Health ------ popuulation_T1 ------ population_T2
# Health ---------------- population_T2
#-------------------------------------------

# 競爭狀態：使用multistat model
p.time <- seq(1, 30, by = 0.5)
T1 <- data$X1
T2 <- data$X2
d2 <- data$D # 現在這個人是不是有得到肝癌 # delta_2
d1 <- (T1<T2) # TRUE/FALSE # delta_1 (X=T1)
#X1=26.827737, X2=26.82774, D=0：還沒得到肝癌就離開實驗了(健康的人)(#9177)
#X1=18.8,X2=19.1,D=1：在X2得到肝癌，也得到肝硬化
#X1=19.1,X2=19.1,D=1：沒有肝硬化直接肝癌
Z <- data$S
X=cbind(data$age, data$sex, data$smoking, data$alcohol, data$alt1)
summary(X)
m = length(T1)
m # 5000 資料中，censor rate 高達 0.9998 (很正常)
mean(d2) # 0.0112 # 肝癌
mean(d1) # 0.0304 # 肝硬化
# 約有1%得到肝癌，有3%得到肝硬化

# 研究目的：肝癌是很重要的國民疾病，今天搜集到一個資料，希望研究是否肝癌的發生率在給定有無肝硬化的情況下會不同，
# 而資料中有肝硬化的時間、C型肝炎的表面抗原、年齡、抽菸(0:無)、飲酒習慣(0:無)、肝指數(0:正常 1:輕度過標)
# 資料名稱：R資料
MM = as.matrix(cbind(Z, X))
Lev_M = apply(MM, 2, median)
BB1 = coxph(Surv(T1, d1)~MM)
names(BB1$coef) = names(data)[4:9]
names(BB1$coef)[1] = "HCV"
BB1
# 根據表格，C型肝原帶原者比起非帶原者得到肝硬化風險的指數比增高了7.3倍
# 另一個顯著的因子是肝發炎指數，可知每增加一個單位的肝指數發炎狀況，會比前一個單位多2.103倍的風險得到肝硬化
# 此外，性別、抽菸、喝酒從表格來看並無統計顯著

st1 = sort(T1, index.return = TRUE)
sd1 = d1[st1$ix]
MMst1 = MM[st1$ix, ] # 排序後的covariate
st1 = st1$x # 排序後的T1

exB1 = exp(apply(t(MMst1) * BB1$coef, 2, sum))
exB1
dLA1 = sd1/(cumsum(exB1[m:1])[m:1]) # lambda_1(t_j)(見照片)
LA1 = cumsum(dLA1) # 累加所有的lambda_d(t_j)

# LA1_i = pLA(T1, st1, LA1)

# plot survival curve
S1 = exp(-LA1)
S1_HBV = exp(-LA1*exp(BB1$coef[1]))
S1_age = exp(-LA1*exp(BB1$coef[2]*48)) # coef * median
S1_alt = exp(-LA1*exp(BB1$coef[6]))
plot(st1, S1, type = "s", ylim=c(0.96, 1)) # health
points(st1, S1_HBV, type = "s", col = "2") 
points(st1, S1_age, type = "s", col = "3")
points(st1, S1_alt, type = "s", col = "4")

# ----------------------------------------------------
# 沒有經過肝硬化直接從健康到得到癌症 ((1-delta_1)*delta_2)
# 只看罹癌者
T02 = ifelse(d1 == 1, T1, T2) # 把經過肝硬化的人宣告成censor
D02 = ifelse(d1 == 1, 0, d2) # 新的censor indicator
BB2 = coxph(Surv(T02, D02) ~ MM) 
names(BB2$coef) = names(data)[4:9]
names(BB2$coef)[1] = "HCV"
BB2
# 顯著的剩下年齡
# 對於是否得癌症，年齡是顯著的有害因子
# C型帶原者相較於沒帶原者罹癌風險增加5.4倍
# 每增加一歲，相較於前一歲罹癌的風險增加1.09倍
# 其他因子對罹癌風險的影響則無統計顯著

st2 = sort(T02, index.return = TRUE)
sd2 = D02[st2$ix]
MMst2 = MM[st2$ix, ]
st2 = st2$x


exB2 = exp(apply(t(MMst2) * BB2$coef, 2, sum))
exB2
dLA2 = sd2/(cumsum(exB2[m:1])[m:1]) # lambda_1(t_j)(見照片)
LA2 = cumsum(dLA2) # 累加所有的lambda_d(t_j)


pLA = function(yy, tt, LL){
  pla = function(yi){
    loc = sum(yi >= tt)
    if(loc==0)ans=0 # ans- hazard:0, survival:1
    if(loc>0)ans=LL[loc]
    ans
  }
  apply(matrix(yy), 1, pla)
}

LA2 = pLA(p.time, st2, LA2)
# plot survival curve
S2 = exp(-LA2)
S2_HBV = exp(-LA2*exp(BB2$coef[1]))
S2_age = exp(-LA2*exp(BB2$coef[2]*48)) # coef * median

plot(p.time, S2, type = "s", ylim=c(0.99, 1)) # health
points(p.time, S2_HBV, type = "s", col = "2") 
points(p.time, S2_age, type = "s", col = "3")

#####
# left truncation的問題
LN1 = T2 > T1 # 一定抓到肝硬化的人
M3 = cbind(Z, X)[LN1, ]
D03 = (d1*d2)[LN1]
T3 = T2[LN1]
R3 = T1[LN1] # 修正得到肝硬化的時間點
BB3 = coxph(Surv(R3, T3, D03) ~ M3)
BB3

st3 = sort(T2, index.return = TRUE)
sd3 = (d1*d2)[st3$ix]
sr3 = T1[st3$ix]
MMst3 = MM[st3$ix, ]
st3 = st3$x
m3 = length(st3)

exB3 = exp(apply(t(MMst3) * BB3$coef, 2, sum))
R3m = matrix(sr3, m3, m3)
T3m = matrix(st3, m3, m3, byrow = TRUE)
Risk3 


# 見照片
#---------------------

dLA3 = sd3/(cumsum(exB3[m:1])[m:1]) # lambda_1(t_j)(見照片)
LA3 = cumsum(dLA3) # 累加所有的lambda_d(t_j)
LA3 = pLA(p.time, st3, LA3)
# plot survival curve
S3 = exp(-LA3)
S3_HBV = exp(-LA3*exp(BB3$coef[1]))
S3_age = exp(-LA3*exp(BB3$coef[2]*48)) # coef * median

plot(p.time, S3, type = "s", ylim=c(0.99, 1)) # health
points(p.time, S3_HBV, type = "s", col = "2") 


