# ====== 路径（按需修改） ======
pgs_path  <- "/Users/li/Desktop/UTAustin/Lab_Rotation/Arbel_Lab/Project/PGS_with_pop_trait.tsv"
hdl_csv <- "/Users/li/Desktop/UTAustin/Lab_Rotation/Arbel_Lab/Project/trait/NCD_RisC_Nature_2020_Cholesterol_age_standardised_countries.csv"
map_path  <- "/Users/li/Desktop/UTAustin/Lab_Rotation/Arbel_Lab/Project/trait/pop_to_country.tsv"
out_dir   <- "/Users/li/Desktop/UTAustin/Lab_Rotation/Arbel_Lab/Project/results_hdl_pgs"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(stringr)
  library(ggplot2); library(broom); library(purrr)
})
library(ggrepel)


hdl_raw <- read_csv(hdl_csv, show_col_types = FALSE)

hdl_tidy <- hdl_raw %>%
  rename(Country = `Country/Region/World`) %>%
  filter(Year == 2013) %>%
  select(Country, Sex, HDL = `Mean HDL cholesterol (mmol/L)`) %>%
  mutate(
    Sex = recode(Sex, "Men" = "Boys", "Women" = "Girls")
  )

# ====== 2) 读取群体-国家映射 ======
pop_map <- read_tsv(map_path, show_col_types = FALSE) %>%
  rename(Population = Population, SuperPop = SuperPop, Country = Country)

# ====== 3) 读取 PGS：筛 Height，按群体取均值（PGS 是个体层面的） ======
pgs_raw <- read_tsv(pgs_path, show_col_types = FALSE)
pgs_height <- pgs_raw %>%
  filter(TraitName == "Height") %>%
  # 你的文件中：Population=群体代码（ACB/CHB...），SuperPop=AFR/EAS...，PGS=个体分值
  group_by(Population, SuperPop) %>%
  summarise(Mean_PGS = mean(PGS, na.rm = TRUE), .groups = "drop") %>%
  inner_join(pop_map, by = c("Population", "SuperPop"))  # 加上 Country

# ====== 4) 合并到国家层面的表型；注意：PGS 本身不分性别 → 复制给 Boys/Girls 两个性别 ======
pgs_by_sex <- pgs_height %>%
  tidyr::expand_grid(Sex = c("Boys", "Girls"))

dat <- pgs_by_sex %>%
  inner_join(hdl_tidy, by = c("Country", "Sex")) %>%
  select(Sex, SuperPop, Population, Country, Mean_PGS, HDL)

# ====== 5a) 计算全局 z-score（在各自性别内标准化） ======
dat_global <- dat %>%
  group_by(Sex) %>%  # 按性别分组后再标准化
  mutate(PGS_z = as.numeric(scale(Mean_PGS))) %>%
  ungroup()

# ====== 5b) 计算组内 z-score（在每个 super-pop 内、并在各自性别内标准化） ======
dat_within <- dat %>%
  group_by(Sex, SuperPop) %>%  # 按性别+super-pop 同时分组
  mutate(PGS_z = as.numeric(scale(Mean_PGS))) %>%
  ungroup()

# ====== 6) 回归分析 ======
# 定义一个通用回归函数
fit_one <- function(df) {
  model <- lm(HDL ~ PGS_z, data = df)
  s <- summary(model)
  tibble(
    n = nrow(df),
    beta1 = coef(model)[["PGS_z"]],
    se = s$coefficients["PGS_z", "Std. Error"],
    p = s$coefficients["PGS_z", "Pr(>|t|)"],
    r = cor(df$PGS_z, df$HDL)
  )
}

## ----- (A) 全局分析：PGS_z = global z -----
res_global <- dat_global %>%
  group_by(Sex) %>%
  group_modify(~ fit_one(.x)) %>%
  mutate(Level = "Global", SuperPop = "All")

## ----- (B) 组内分析：PGS_z = within z -----
res_within <- dat_within %>%
  group_by(Sex, SuperPop) %>%
  group_modify(~ fit_one(.x)) %>%
  mutate(Level = "Within-superpop")

# 合并结果
res_all <- bind_rows(res_global, res_within) %>%
  select(Level, Sex, SuperPop, n, beta1, se, r, p)

# 保存
write_csv(res_all, file.path(out_dir, "hdl_pgs_regression_results.csv"))


# ====== 7) 作图（全局 vs 各 super-pop） ======
# --- 全局图 ---
# 1) 先按 Sex 计算统计量
# 1) 统计量 + 标注位置（统一放右上角）
stats_global <- dat_global %>%
  dplyr::group_by(Sex) %>%
  dplyr::do({
    df <- .
    m <- lm(HDL ~ PGS_z, data = df)
    s <- summary(m)
    tibble::tibble(
      beta1 = unname(coef(m)["PGS_z"]),
      p = s$coefficients["PGS_z", "Pr(>|t|)"],
      r = cor(df$PGS_z, df$HDL, use = "complete.obs"),
      x = Inf,   # 关键：锚到右上角（两面板一致）
      y = Inf
    )
  }) %>% 
  dplyr::ungroup() %>%
  dplyr::mutate(label = sprintf("β₁ = %.2f\nr = %.2f\np = %.2g", beta1, r, p))


# 2) 画图（仅一条全局回归线）+ 右上角标注
p_global <- dat_global %>%
  ggplot(aes(x = PGS_z, y = HDL, color = SuperPop, label = Population)) +
  geom_point(size = 1.5) +
  geom_text_repel(aes(label = Population), size = 2, max.overlaps = 50, show.legend = FALSE)+
  geom_smooth(aes(color = NULL), method = "lm", se = TRUE,
              color = "black", linewidth = 1) +
  facet_wrap(~ Sex, scales = "free_y") +
  geom_text(
    data = stats_global,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 1.05, vjust = 1.1,
    size = 3, color = "black"
  )+
  labs(x = "PGS (global z-score)", y = "HDL cholesterol(mmol/L)",
       title = "Global regression: HDL cholesterol vs PGS (z, all populations combined)") +
  theme_bw()

ggsave(file.path(out_dir, "hdl_vs_pgs_global_annotated.png"),
       p_global, width = 8, height = 4.5, dpi = 300)

# --- 每个 super-pop 各画一张 ---
# --- 每个 super-pop 各画一张（带 β₁、r、p 标注）---
for (sp in unique(dat_within$SuperPop)) {
  df_sp <- dat_within %>% dplyr::filter(SuperPop == sp)
  
  # 1) 计算该 super-pop 内、按 Sex 的统计量，并把标注锚到右上角
  stats_within <- df_sp %>%
    dplyr::group_by(Sex) %>%
    dplyr::do({
      d <- .
      # 小样本保护：<3个点就跳过
      if (nrow(d) < 3) {
        tibble::tibble(beta1 = NA_real_, p = NA_real_, r = NA_real_, x = Inf, y = Inf)
      } else {
        m <- lm(HDL ~ PGS_z, data = d)
        s <- summary(m)
        tibble::tibble(
          beta1 = unname(coef(m)["PGS_z"]),
          p = s$coefficients["PGS_z", "Pr(>|t|)"],
          r = cor(d$PGS_z, d$HDL, use = "complete.obs"),
          x = Inf,  # 锚到右上角（不受坐标范围影响）
          y = Inf
        )
      }
    }) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(label = sprintf("β₁ = %.2f\nr = %.2f\np = %.2g", beta1, r, p))
  
  # 2) 作图（该 super-pop 内按 Sex 分面）
  p_sp <- df_sp %>%
    ggplot(aes(x = PGS_z, y = HDL, label = Population)) +
    geom_point(size = 1.5) +
    geom_text_repel(aes(label = Population), size = 2, max.overlaps = 50, show.legend = FALSE)+
    geom_smooth(method = "lm", se = TRUE) +
    facet_wrap(~ Sex, scales = "free_y") +
    geom_text(
      data = stats_within,
      aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      hjust = 1.05, vjust = 1.1,
      size = 3, color = "black"
    )+
    labs(x = paste0("PGS within ", sp, " (z-score)"),
         y = "HDL cholesterol(mmol/L)",
         title = paste0("Within ", sp, ": HDL cholesterol vs PGS (z)")) +
    theme_bw()
  
  ggsave(file.path(out_dir, paste0("hdl_vs_pgs_within_", sp, "_annotated.png")),
         p_sp, width = 8, height = 4.5, dpi = 300)
}



