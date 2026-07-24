# Age prediction and uncertainty quantification using quantile regression (QR)
# Models were trained using UK Biobank Olink proteomics data
# Reference (TBA)
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(splines2)
library(ggplot2)
library(survival)
options(scipen = 999)

dir = paste0(".")

## 1. Phenotype data
## 1.1 Chronological age at sample collection and other covarites such as sex
df.age = fread(file = paste0(dir, "/example/example.age.tsv"),
               header = T, sep = "\t", data.table = F, stringsAsFactors = F)

## 1.2 Protein expression matrix
## IID column: unique sample identifier
## remaining columns: protein expression values, with one column per protein and one row per sample
df.prot = fread(file = paste0(dir, "/example/example.proteomics.tsv"),
                header = T, sep = "\t", data.table = F, stringsAsFactors = F)
mt.prot = df.prot %>% select(-IID) %>% scale() %>% as.data.frame()
mt.prot[is.na(mt.prot)] = 0
mt.prot = mt.prot %>% mutate(Intercept = 1)



## 2. Organ age models
## Select the identifier column from the model weight table to 
## match the protein names in expression matrix.
## For example, 
## if proteins are named using the Olink panel name and gene symbol,
## create an identifier column in the model weight table: 
## df.beta = mutate(Predictor = paste0(Protein.Panel, "_", Gene.Symbol))
## By default, 
## the assay target ID is assumed to be the common identifier
## used to map the proteins between the expression matrix and the model weights
df.beta = fread(file = paste0(dir, "/model/qr.beta.proteomic.age.tsv"),
                header = T, sep = "\t", data.table = F, stringsAsFactors = F)
df.beta = df.beta %>% mutate(Predictor = Assay.Target)



## 3. Organ age predictions
df.pred = data.frame()
ogn.pool = df.beta %>% distinct(Organ) %>% pull(Organ)
for (ogn in ogn.pool) {
  #### 3.1 Applying the QR model to the test set
  beta = df.beta %>% 
    filter(Organ == ogn,
           Predictor %in% colnames(mt.prot)) %>%
    select(Predictor, starts_with("QR"))
  pred = as.data.frame(as.matrix(mt.prot %>% select(all_of(beta$Predictor))) %*% as.matrix(beta %>% select(-Predictor))) %>%
    mutate(df.prot %>% select(IID), .before = 1)
  
  #### 3.2 Post hoc sorting
  df.qntl = data.frame(Name = pred %>% select(starts_with("QR")) %>% colnames(), 
                       stringsAsFactors = F) %>%
    mutate(Idx = 1:nrow(.))
  pred = pred %>% 
    pivot_longer(cols = starts_with("QR"), names_to = "Name", values_to = "Value") %>%
    group_by(IID) %>%
    mutate(Idx = row_number()) %>%
    arrange(Value, Idx, .by_group = T) %>%
    mutate(Idx = row_number()) %>%
    ungroup() %>%
    select(-Name) %>%
    left_join(y = df.qntl, by = join_by(Idx)) %>%
    pivot_wider(id_cols = IID, 
                names_from = Name, values_from = Value)
  pred = pred %>% 
    left_join(y = df.age, by = join_by(IID)) %>%
    relocate(Sex, Age, .after = IID)
  
  #### 3.3 Age-calibrated age gap and standardized age gap z-score
  #### age prediction => Age.QR0.5
  #### age-calibrated age gap => Gap.QR0.5
  #### standardized age-calibrated age gap z-score => Gap.Std.QR0.5
  fit = lm(data = pred, formula = QR0.5 ~ Age)
  adj = unname(fit$fitted.values)
  pred = pred %>% 
    mutate(across(.cols = starts_with("QR"),
                  .fns = ~ .x - adj,
                  .names = "Gap.{.col}"),
           across(.cols = starts_with("QR"),
                  .fns = ~ .x - adj,
                  .names = "Gap.Std.{.col}")) %>%
    rename_with(.cols = starts_with("QR"), 
                .fn = ~ paste0("Age.", .x))
  m.pred = pred %>% pull(Gap.Std.QR0.5) %>% mean()
  sd.pred = pred %>% pull(Gap.Std.QR0.5) %>% sd()
  pred = pred %>% 
    mutate(across(.cols = starts_with("Gap.Std"),
                  .fns = ~ (.x - m.pred) / sd.pred))
  
  #### 3.4 Tail probability P(Gap > 0)
  #### tail probability => Tail.Prob
  pred = pred %>% mutate(Tail.Prob = 0)
  th.prob = 0
  n.sample = 10000
  set.seed(seed = 256)
  prob.int = runif(n.sample, 0, 1)
  prob.cum = c(0.025, (1:9) / 10, 0.975)
  set.seed(seed = 256)
  for (idx in 1:nrow(pred)) {
    tmp = pred %>% 
      slice(idx) %>% 
      select(starts_with("Gap.QR")) %>%
      pivot_longer(cols = starts_with("Gap.QR"), 
                   names_to = "Quantile", 
                   values_to = "QR") %>% 
      mutate(Prob = prob.cum)
    fit.spline = lm(data = tmp, formula = QR ~ bSpline(Prob, Boundary.knots = c(0, 1)))
    pred.spline = predict(object = fit.spline, newdata = data.frame(Prob = prob.int))
    pred$Tail.Prob[idx] = sum(pred.spline >= th.prob) / n.sample
  }
  
  df.pred = rbind(df.pred,
                  pred %>% mutate(Organ = ogn, .before = 1))
  print(paste0("==== ", ogn, " ==== ", Sys.time(), " ===="))
}

#### 3.5 Age gap interval length
#### prediction interval length => Interval.Length
df.pred = df.pred %>% mutate(Interval.Length = Gap.QR0.975 - Gap.QR0.025)

#### 3.6 Extreme ager status (aging outlier) by age gap
th.z = 1.5
df.pred = df.pred %>% 
  mutate(Outlier.Gap = case_when(Gap.Std.QR0.5 > th.z ~ "Extreme age",
                                 Gap.Std.QR0.5 < -th.z ~ "Extreme youth",
                                 .default = "Normal"))



## 4. Use age gap, tail probability, and interval length as aging measures
## age gap => Gap.Std.QR0.5
## tail probability => Tail.Prob
## interval length => Interval.Length

ofile = paste0(dir, "/example/example.prediction.tsv")
fwrite(x = df.pred,
       file = ofile,
       quote = F, sep = "\t", row.names = F, col.names = T)



#### 4.1 Incident disease/mortality risk and aging measures
# df.pred = fread(file =  paste0(dir, "/example/example.prediction.tsv"),
#                 header = T, sep = "\t", data.table = F, stringsAsFactors = F)

df.event = fread(file = paste0(dir, "/example/example.disease.tsv"),
                 header = T, sep = "\t", data.table = F, stringsAsFactors = F)

#### exclude prevalent disease cases and retain incident cases
df.event = df.event %>% filter(Duration.Event > 0, !is.na(Status.Event))

#### scale of predictor (standardized z-score, quintile strata)
scl.pool = c("Standardized", "Quintile")

#### predictor (age gap, tail probability, interval length)
prd.pool = c("Gap", "Pro", "Len")

#### disease or mortality
dis.pool = c("Hypertension")
df.hr.dis = data.frame()
for (scl in scl.pool) {
  for (prd in prd.pool) {
    for (dis in dis.pool) {
      for (ogn in ogn.pool) {
        tmp = df.event %>% 
          filter(Event == dis) %>% 
          inner_join(y = df.pred %>% filter(Organ == ogn), 
                     by = join_by(IID))
        #### predictor => standardized values (continuous variable)
        if (scl == "Standardized") {
          tmp = tmp %>%
            mutate(Gap = as.numeric(scale(Gap.Std.QR0.5)),
                   Pro = as.numeric(scale(Tail.Prob)),
                   Len = as.numeric(scale(Interval.Length)))
        }
        #### predictor => quintile strata (ordinal variable)
        if (scl == "Quintile") {
          tmp = tmp %>%
            mutate(Gap = ntile(x = Gap.Std.QR0.5, n = 5),
                   Pro = ntile(x = Tail.Prob, n = 5),
                   Len = ntile(x = Interval.Length, n = 5))
        }
        tmp = tmp %>% filter(!is.na(!!as.symbol(prd)))
        if (nrow(tmp) == 0 | sum(tmp$Status.Event) == 0) {
          df.hr.dis = rbind(df.hr.dis,
                            data.frame(Phenotype = dis, Organ = ogn, Predictor = prd, Scale = scl,
                                       N.Total = nrow(tmp), N.Event = sum(tmp$Status.Event),
                                       HR = NA, CI.Lower = NA, CI.Upper = NA, P = NA,
                                       stringsAsFactors = F))
        }
        #### cox proportional hazards regression
        if (nrow(tmp) > 0 & sum(tmp$Status.Event) > 0) {
          fml = paste0("Surv(time = Duration.Event, event = Status.Event) ~ ", prd, " + Age + Sex")
          fit.cox = coxph(data = tmp, formula = as.formula(fml))
          # stat.cox = broom::tidy(x = fit.cox, conf.int = T, exponentiate = T) %>% filter(term == prd)
          stat.cox = summary(fit.cox)
          df.hr.dis = rbind(df.hr.dis,
                            data.frame(Phenotype = dis, Organ = ogn, Predictor = prd, Scale = scl,
                                       N.Total = stat.cox$n, N.Event = stat.cox$nevent,
                                       HR = stat.cox$coefficients[prd, "exp(coef)"],
                                       CI.Lower = stat.cox$conf.int[prd, "lower .95"],
                                       CI.Upper = stat.cox$conf.int[prd, "upper .95"],
                                       P = stat.cox$coefficients[prd, "Pr(>|z|)"], 
                                       stringsAsFactors = F))
        }
      }
    }
  }
}

df.hr.dis = df.hr.dis %>% 
  mutate(Predictor = case_when(Predictor == "Gap" ~ "Age gap",
                               Predictor == "Pro" ~ "Tail probability",
                               Predictor == "Len" ~ "Interval length"))

ofile = paste0(dir, "/example/example.hr.disease.tsv")
fwrite(x = df.hr.dis,
       file = ofile,
       quote = F, sep = "\t", row.names = F, col.names = T)



#### 4.2 Age prediction performance
#### pearson correlation between proteomic age vs chronological age 
df.eval = df.pred %>% 
  mutate(Is.Covered.rho95 = 1 * (Age >= Age.QR0.025 & Age <= Age.QR0.975),
         Is.Covered.rho80 = 1 * (Age >= Age.QR0.1 & Age <= Age.QR0.9)) %>%
  group_by(Organ) %>%
  summarise(Coverage.rho95 = 100 * sum(Is.Covered.rho95)/n(),
            Coverage.rho80 = 100 * sum(Is.Covered.rho80)/n(),
            PCC = unname(cor.test(x = Age, y = Age.QR0.5)$estimate),
            MAE = mean(abs(Age - Age.QR0.5))) %>%
  ungroup()

df.eval = df.eval %>% 
  left_join(df.beta %>% filter(Predictor != "Intercept") %>% count(Organ) %>% rename(N.Protein = n), 
            by = join_by(Organ))

level.ogn = df.eval %>% arrange(desc(PCC)) %>% pull(Organ)
level.ogn.str = df.eval %>% arrange(desc(PCC)) %>% mutate(Organ.str = paste0(Organ, "\n", "# Proteins = ", N.Protein)) %>% pull(Organ.str)

df.pred = df.pred %>% 
  mutate(Organ.fct = factor(x = Organ, levels = level.ogn),
         Organ.str.fct = factor(x = Organ, levels = level.ogn, labels = level.ogn.str))

df.eval.lab = df.eval %>% 
  mutate(Label = paste0("r = ", formatC(x = PCC, format = "f", digits = 3), "\n",
                        "MAE = ", formatC(x = MAE, format = "f", digits = 2), "\n")) %>%
  mutate(Organ.fct = factor(x = Organ, levels = level.ogn),
         Organ.str.fct = factor(x = Organ, levels = level.ogn, labels = level.ogn.str)) %>% 
  left_join(y = df.pred %>% 
              group_by(Organ) %>%
              summarise(Max.Pred = max(Age.QR0.5),
                        Min.Pred = min(Age.QR0.5)) %>%
              ungroup(), 
            by = join_by(Organ)) %>%
  mutate(X.Label = 40,
         Y.Label = Max.Pred - 0.18 * (Max.Pred - Min.Pred))

g1 = ggplot(data = df.pred, 
            mapping = aes(x = Age, y = Age.QR0.5)) +
  facet_wrap(facets = vars(Organ.str.fct), nrow = 3, scale = "free") +
  geom_point(shape = 20, size = 1, color = "#fb6a4a", fill = "#fb6a4a", alpha = 0.1) +
  geom_smooth(method = "lm", se = F, linewidth = 0.6, colour = "#000000", alpha = 0.5) + 
  geom_text(data = df.eval.lab,
            mapping = aes(x = X.Label, y = Y.Label, label = Label),
            hjust = 0, size = 3) +
  labs(x = "Chronological age", 
       y = "Proteomic age prediction",
       title = paste0("Performance of proteomic aging models")) +
  theme_minimal() +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "#000000"),
        text = element_text(family = "Helvetica"),
        strip.text = element_text(size = 10, color = "#000000"),
        axis.text = element_text(size = 10, color = "#000000"),
        axis.title = element_text(size = 10, color = "#000000"),
        plot.title = element_text(size = 12, color = "#000000"))

ofile = paste0(dir, "/example/fig.age.corr.png")
ggsave(plot = g1, 
       filename = ofile,
       device = "png", bg = "#ffffff", width = 9, height = 6.9)



#### 4.3 Distribution of prediction interval lengths
lvl.ogn.len = df.pred %>%
  group_by(Organ) %>%
  summarise(Len = median(Interval.Length)) %>%
  ungroup() %>%
  arrange(Len)

df.pred = df.pred %>% 
  mutate(Organ.fct = factor(x = Organ, levels = lvl.ogn.len$Organ))

g2 = ggplot(data = df.pred, 
            mapping = aes(x = Organ.fct, y = Interval.Length)) +
  geom_boxplot(outlier.shape = 16, outlier.size = 0.4, outlier.alpha = 0.3,
               width = 0.6, linewidth = 0.2, alpha = 0.7,
               color = "#2171b5", fill = "#2171b5") +
  labs(y = "Length",
       title = paste0("Length of age gap prediction intervals")) +
  theme_minimal() +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "#000000"),
        text = element_text(family = "Helvetica"),
        axis.text.x = element_text(size = 10, color = "#000000", angle = 45, hjust = 1),
        axis.text.y = element_text(size = 10, color = "#000000"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10, color = "#000000"),
        plot.title = element_text(size = 12, color = "#000000"))

ofile = paste0(dir, "/example/fig.dist.length.png")
ggsave(plot = g2, 
       filename = ofile,
       device = "png", bg = "#ffffff", width = 5, height = 4)



#### 4.4 Distributions of age gap and tail probability for extreme agers
#### outlier by age gap
df.extr = df.pred %>% 
  filter(Outlier.Gap != "Normal") %>%
  pivot_longer(cols = c(Gap.Std.QR0.5, Tail.Prob), names_to = "Measure", values_to = "Score") %>%
  mutate(Measure = case_when(Measure == "Gap.Std.QR0.5" ~ "Age gap",
                             Measure == "Tail.Prob" ~ "Tail probability",
                             Measure == "Interval.Length" ~ "Interval length"))

lvl.ogn.extr = df.extr %>%
  filter(Measure == "Tail probability") %>%
  group_by(Organ) %>%
  summarise(SD = sd(Score)) %>%
  ungroup() %>%
  arrange(desc(SD)) %>%
  pull(Organ)

df.extr = df.extr %>% 
  mutate(Intercept = if_else(Measure == "Age gap", 0, 0.5),
         Organ.fct = factor(x = Organ, levels = lvl.ogn.extr),
         Outlier.Gap.fct = factor(x = Outlier.Gap, levels = c("Extreme age", "Extreme youth")))

g3 = ggplot(data = df.extr, 
              mapping = aes(x = Organ.fct, y = Score, 
                            color = Outlier.Gap.fct, fill = Outlier.Gap.fct)) +
  facet_wrap(facets = vars(Measure), ncol = 1, scale = "free_y", strip.position = "left") +
  geom_hline(mapping = aes(yintercept = Intercept), linetype = "longdash", linewidth = 0.2, color = "#969696") +
  geom_boxplot(position = position_dodge(width = 0), 
               outlier.shape = 16, outlier.size = 0.4, outlier.alpha = 0.3,
               width = 0.8, linewidth = 0.4, alpha = 0.7) +
  scale_color_manual(values = c("#fb6a4a", "#6baed6")) +
  scale_fill_manual(values = c("#fb6a4a", "#6baed6")) +
  labs(x = NULL, y = NULL, 
       fill = NULL, color = NULL,
       title = "Age gap and tail probability for extreme agers") +
  theme_minimal() +
  theme(legend.position = "top",
        legend.justification  = "center",
        legend.key.size = unit(x = 10, units = "pt"),
        legend.box.margin = margin(t = 0, b = 0, unit = "pt" ),
        legend.box.spacing = unit(x = 0, units = "pt"),
        panel.border = element_rect(fill = NA, color = "#000000"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.placement = "outside",
        text = element_text(family = "Helvetica"),
        strip.text = element_text(size = 10, color = "#000000"),
        axis.text.x = element_text(size = 10, color = "#000000", angle = 45, hjust = 1),
        axis.text.y = element_text(size = 10, color = "#000000"),
        legend.text = element_text(size = 10, color = "#000000"),
        axis.title = element_blank(),
        plot.title = element_text(size = 10, color = "#000000"))
g3

ofile = paste0(dir, "/example/fig.dist.gap.prob.extreme.png")
ggsave(plot = g3, 
       filename = ofile,
       device = "png", bg = "#ffffff", width = 4, height = 5.283)





























