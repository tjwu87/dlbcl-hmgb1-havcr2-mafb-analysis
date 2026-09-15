## ============================================================
## 分析3 路线a：per-SD Cox 敏感性分析（cutoff不依赖性检验）
## 输入：All_cohorts_risk_scores.csv（与脚本同目录）
## 输出：Table S20 数据（analysis3_perSD_Cox.csv）
## ============================================================

library(survival)

df <- read.csv("All_cohorts_risk_scores.csv", stringsAsFactors = FALSE)
df$Cohort <- trimws(df$Cohort)

cohorts <- c("GSE10846", "GSE87371", "GSE11318", "GSE181063", "TCGA-DLBC")

## 论文已发表的 median-split HR（Table S14），用于自检对照
published_hr <- c("GSE10846" = 2.702, "GSE87371" = 2.367,
                  "GSE11318" = 1.872, "GSE181063" = 1.418,
                  "TCGA-DLBC" = 2.166)

## ------------------------------------------------------------
## 自检A：各队列 n / events（先确认数据加载和分组没出错）
## ------------------------------------------------------------
cat("==== 自检A：队列规模 ====\n")
for (co in cohorts) {
  d <- subset(df, Cohort == co)
  cat(sprintf("%-10s  n = %3d   events = %3d\n",
              co, nrow(d), sum(d$OS_event)))
}
## 预期：GSE10846 n≈412/163，GSE87371 n=111，GSE11318 n≈200±，
##       GSE181063 n≈441，TCGA n=45/9 —— 与论文Table S14核对

## ------------------------------------------------------------
## 自检B：median-split Cox（用现成RiskGroup复现论文HR）
## ------------------------------------------------------------
cat("\n==== 自检B：median-split HR vs 论文Table S14 ====\n")
chk <- data.frame()
for (co in cohorts) {
  d <- subset(df, Cohort == co)
  d$RiskGroup <- relevel(factor(d$RiskGroup), ref = "Low")  # HR = High vs Low
  s <- summary(coxph(Surv(OS_years, OS_event) ~ RiskGroup, data = d))
  chk <- rbind(chk, data.frame(
    cohort = co,
    HR_reproduced = s$conf.int[1],
    HR_published  = published_hr[[co]]
  ))
}
print(chk, digits = 3)
## 判定标准：HR_reproduced 与 HR_published 之差 < 0.05 即通过。
## 全部通过 ⇒ 您的RiskScore/临床数据与当年发表管线一致，继续往下跑。
## 若某个队列对不上 ⇒ 停，把输出发回来，别用这套结果。

## ------------------------------------------------------------
## 核心：路线a —— 各队列内标准化（mean 0, SD 1）后做 per-SD Cox
## ------------------------------------------------------------
res <- list()
for (co in cohorts) {
  d <- subset(df, Cohort == co)
  d$risk_sd <- as.numeric(scale(d$RiskScore))   # 队列内z-score，与Methods口径一致
  s <- summary(coxph(Surv(OS_years, OS_event) ~ risk_sd, data = d))
  res[[co]] <- data.frame(
    cohort    = co,
    n         = nrow(d),
    events    = sum(d$OS_event),
    HR_per_SD = s$conf.int[1],
    CI_lower  = s$conf.int[3],
    CI_upper  = s$conf.int[4],
    p_value   = s$coefficients["risk_sd", "Pr(>|z|)"]
  )
}
tab <- do.call(rbind, res)
tab$significant <- tab$CI_upper < 1 | tab$CI_lower > 1   # CI不含1
cat("\n==== 路线a核心结果：per-SD Cox（Table S20）====\n")
print(tab, digits = 3)

## ------------------------------------------------------------
## 加分项：五个队列 per-SD logHR 的逆方差合并（meta估计）
## 不需要额外安装包
## ------------------------------------------------------------
tab$logHR <- log(tab$HR_per_SD)
tab$se    <- (log(tab$CI_upper) - log(tab$CI_lower)) / (2 * 1.96)
w   <- 1 / tab$se^2
pool_fix <- exp(sum(w * tab$logHR) / sum(w))
Q   <- sum(w * (tab$logHR - sum(w * tab$logHR)/sum(w))^2)
pool_rnd_se <- sqrt(1/sum(w) + sum(Q - (length(cohorts)-1)) /
                    (sum(w)^2 * (sum(1/ w) - sum(1/w^2)/sum(w))) )  # DerSimonian-Laird
## 随机效应权重（手写DL法，稳妥起见两种都报）
tau2 <- max(0, (Q - (length(cohorts)-1)) /
            (sum(w) - sum(w^2)/sum(w)))
w_r  <- 1 / (tab$se^2 + tau2)
pool_rnd <- exp(sum(w_r * tab$logHR) / sum(w_r))
se_rnd   <- sqrt(1 / sum(w_r))

cat(sprintf("\n合并估计（固定效应）HR = %.3f\n", pool_fix))
cat(sprintf("合并估计（随机效应）HR = %.3f  (95%% CI %.3f–%.3f)\n",
            pool_rnd, exp(log(pool_rnd) - 1.96*se_rnd), exp(log(pool_rnd) + 1.96*se_rnd)))

## ------------------------------------------------------------
## 导出 Table S20
## ------------------------------------------------------------
write.csv(tab[, c("cohort","n","events","HR_per_SD","CI_lower","CI_upper","p_value","significant")],
          "analysis3_perSD_Cox.csv", row.names = FALSE)
cat("\n已导出：analysis3_perSD_Cox.csv → 整理进补充表格即 Table S20\n")

