# Age prediction and uncertainty quantification using quantile regression (QR)
# Train age models using user's own Olink proteomics data
# Reference (TBA)
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(splines2)
library(fsQRPPA)
library(RcppParallel)
options(scipen = 999)

dir = paste0(".")

## 1. Phenotype data
#### 1.1 Chronological age at sample collection and other covarites such as sex
df.age = fread(file = paste0(dir, "/example/example.age.tsv"),
               header = T, sep = "\t", data.table = F, stringsAsFactors = F)

#### 1.2 Protein expression matrix for training set
#### IID column: unique sample identifier
#### remaining columns: protein expression values, with one column per protein and one row per sample
df.prot = fread(file = paste0(dir, "/example/example.proteomics.tsv"),
                header = T, sep = "\t", data.table = F, stringsAsFactors = F)
df.prot = df.age %>% 
  select(IID, Age) %>%
  left_join(y = df.prot, by = join_by(IID))

#### Select training, validation, and test samples based on the user's data
#### The input dataset can be partitioned into training, validation, and test sets 
#### in a single step using the code below. For example, the dataset is randomly partitioned into 
#### three sets, containing 60%, 20%, and 20% of the samples, respectively.
#### Alternatively, 
#### users may split the available dataset into training and validation sets for model development 
#### and reserve a completely independent dataset for subsequent testing. 
set.seed(seed = 256)
sample.tr = df.prot %>% 
  select(IID) %>% 
  slice_sample(prop = 0.6) %>%
  mutate(Set = "training")
sample.va = df.prot %>% 
  select(IID) %>% 
  filter(!IID %in% sample.tr$IID) %>%
  slice_sample(prop = 0.5) %>%
  mutate(Set = "validation")
sample.te = df.prot %>% 
  select(IID) %>% 
  filter(!IID %in% sample.tr$IID) %>%
  filter(!IID %in% sample.va$IID) %>%
  mutate(Set = "test")
df.prot = rbind(sample.tr, sample.va, sample.te) %>%
  left_join(y = df.prot, by = join_by(IID))

#### protein expression values are standardized to z-scores for each protein across all individuals
m.tr = sapply(X = df.prot %>% filter(Set == "training") %>% select(-IID, -Set, -Age), 
              FUN = mean, na.rm = T)
sd.tr = sapply(X = df.prot %>% filter(Set == "training") %>% select(-IID, -Set, -Age), 
               FUN = sd, na.rm = T)
mt.prot = df.prot %>% 
  select(-IID, -Set, -Age) %>% 
  scale(center = m.tr, scale = sd.tr) %>% 
  as.data.frame()
mt.prot[is.na(mt.prot)] = 0
mt.prot = mt.prot %>% mutate(Intercept = 1)
df.prot = df.prot %>% 
  select(IID, Set, Age) %>% 
  mutate(mt.prot)
rm(mt.prot)

#### 1.3. Organ enriched proteins for proteomics profiling quantified using the Olink Explore 3072 platform
#### https://www.nature.com/articles/s41591-025-03798-1
####
#### Select the identifier column from the model weight table to 
#### match the protein names in expression matrix.
#### For example, 
#### if proteins are named using the Olink panel name and gene symbol,
#### create an identifier column in the model weight table: 
#### df.organ = df.organ %>% mutate(Predictor = paste0(Protein.Panel, "_", Gene.Symbol))
####
#### By default, 
#### the assay target ID is assumed to be the common identifier
#### used to map the proteins between the expression matrix and the model weights
####
#### The file of organ enriched proteins includes one row per organ for the intercept term 
#### in addition to the protein predictors. The intercept should be included when performing 
#### model training and calculating predictions.
df.organ = fread(file = paste0(dir, "/model/qr.beta.proteomic.age.tsv"),
                header = T, sep = "\t", data.table = F, stringsAsFactors = F)
df.organ = df.organ %>% mutate(Predictor = Assay.Target)



## 2. Model training
#### 2.1 Model parameters
#### adjust the maximum number of CPU threads based on the user's computing environment
max.nthread = 16
RcppParallel::setThreadOptions(numThreads = max.nthread)

#### fsQRPPA parameters
mu = 1e-5
c_eta = 1e-3
tol_abs = 1e-4
tol_rel = 1e-3
lambda_min_ratio = 0.01
lambda_max_ratio = 1.1
n_lambda = 20
step_LLA = 1 
max_iter = 1e4
min_iter = 10



#### 2.2 training
#### The example code trains all composite and organ-specific models sequentially in a single loop. 
#### Because each organ model is trained independently, users can reduce the runtime by training 
#### multiple organ models in parallel.
tau.pool = c(0.025, (1:9) / 10, 0.975)
ogn.pool = df.organ %>% distinct(Organ) %>% arrange(Organ) %>% pull(Organ)
# pen.pool = c("LASSO", "SCAD", "MCP", "ELASTIC_NET")
pen = "ELASTIC_NET"
a_param = ifelse(pen == "SCAD", 3.7, 3)

df.beta = data.frame() 
df.pred = data.frame()
for (ogn in ogn.pool) {
  prot.pool = df.prot %>% select(any_of(df.organ %>% filter(Organ == ogn) %>% pull(Predictor))) %>% colnames()
  n.group = min(length(prot.pool), max.nthread)
  for (tau in tau.pool) {
    lambda_seq = gen_lambda_seq(X = df.prot %>% 
                                  filter(Set == "training") %>% 
                                  select(all_of(prot.pool)) %>% 
                                  relocate(Intercept, .before = 1) %>% 
                                  as.matrix(), 
                                tau = tau, n_lambda = n_lambda, 
                                lambda_max_ratio = lambda_max_ratio,
                                lambda_min_ratio = lambda_min_ratio,
                                G = n.group, standardized = T)
    lambda.pool = c(lambda_seq, 0)
    zeta.pool = c(1e-2, 1e-4, 1e-6, 0)
    if (pen != "ELASTIC_NET") {
      zeta.pool = 0
    }
    param.grid = expand.grid(lambda = lambda.pool,
                             zeta = zeta.pool) %>% 
      mutate(warmup = lambda + zeta) %>%
      arrange(desc(warmup), desc(lambda), desc(zeta)) %>%
      filter(warmup != 0)
    
    fit = fsQRPPA(X = df.prot %>% 
                    filter(Set == "training") %>% 
                    select(all_of(prot.pool)) %>% 
                    relocate(Intercept, .before = 1) %>% 
                    as.matrix(), 
                   y = df.prot %>% 
                    filter(Set == "training") %>% 
                    pull(Age), 
                   tau = tau, lambda = param.grid$lambda, zeta = param.grid$zeta,
                   penalty = pen, a = a_param, c_eta = c_eta, mu = mu, 
                   max_iter = max_iter, min_iter = min_iter, tol_abs = tol_abs, tol_rel = tol_rel, 
                   step_LLA = step_LLA, G = n.group, standardized = T, verbose = F)
    
    loss.pinball = pinball_loss_path(X = df.prot %>% 
                                       filter(Set == "validation") %>% 
                                       select(all_of(prot.pool)) %>% 
                                       relocate(Intercept, .before = 1) %>% 
                                       as.matrix(),
                                      y = df.prot %>% 
                                       filter(Set == "validation") %>% 
                                       pull(Age),
                                      tau = tau, G = n.group, Beta = fit$beta, standardized = T)
    opt.idx = which.min(loss.pinball)
    opt.lambda = param.grid$lambda[opt.idx]
    opt.zeta = param.grid$zeta[opt.idx]
    opt.beta = fit$beta[, opt.idx]
    opt.iter = fit$n_iter[opt.idx]
    opt.pred = as.matrix(df.prot %>% select(all_of(prot.pool)) %>% relocate(Intercept, .before = 1)) %*% as.matrix(opt.beta)
    
    df.beta = rbind(df.beta,
                    data.frame(Organ = ogn, Tau = tau, 
                               Predictor = prot.pool, Beta = opt.beta))
    df.pred = rbind(df.pred,
                    as.data.frame(opt.pred) %>%
                      mutate(Organ = ogn, Tau = tau) %>%
                      mutate(df.prot %>% select(IID)))
  }
}



## 3. Model and predictions
#### 3.1 Age prediction model coefficients
df.beta = df.beta %>% 
  pivot_wider(id_cols = c(Organ, Predictor), 
              names_from = Tau, names_prefix = "QR", 
              values_from = Beta)

ofile = paste0(dir, "/example/example.model.beta.tsv")
fwrite(x = df.beta,
       file = ofile,
       quote = F, sep = "\t", row.names = F, col.names = T)



#### 3.2 Age prediction results for all samples
#### post hoc sorting
df.qntl = data.frame(Name = df.pred %>% distinct(Tau) %>% arrange(Tau) %>% pull(Tau), 
                     stringsAsFactors = F) %>%
  mutate(Idx = 1:nrow(.), Name = paste0("QR", Name))
df.pred = df.pred %>% 
  group_by(Organ, IID) %>%
  arrange(V1, Tau, .by_group = T) %>%
  mutate(Idx = row_number()) %>%
  ungroup() %>%
  left_join(y = df.qntl, by = join_by(Idx)) %>%
  pivot_wider(id_cols = c(Organ, IID), 
              names_from = Name, 
              values_from = V1)
df.pred = df.pred %>% 
  left_join(y = df.prot %>% select(IID, Set, Age), by = join_by(IID)) %>%
  relocate(Set, Age, .after = IID)

ofile = paste0(dir, "/example/example.model.prediction.tsv")
fwrite(x = df.pred,
       file = ofile,
       quote = F, sep = "\t", row.names = F, col.names = T)



#### 3.3 Age gap and predictive uncertainty for test samples
df.gap = data.frame()
ogn.pool = df.beta %>% distinct(Organ) %>% pull(Organ)
for (ogn in ogn.pool) {
  #### age prediction => Age.QR0.5
  #### age-calibrated age gap => Gap.QR0.5
  #### standardized age-calibrated age gap z-score => Gap.Std.QR0.5
  pred = df.pred %>% 
    filter(Organ == ogn, Set == "test")
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
  
  df.gap = rbind(df.gap, pred)
  # print(paste0("==== ", ogn, " ==== ", Sys.time(), " ===="))
}

#### prediction interval length => Interval.Length
df.gap = df.gap %>% mutate(Interval.Length = Gap.QR0.975 - Gap.QR0.025)

ofile = paste0(dir, "/example/example.model.uncertainty.tsv")
fwrite(x = df.pred,
       file = ofile,
       quote = F, sep = "\t", row.names = F, col.names = T)


