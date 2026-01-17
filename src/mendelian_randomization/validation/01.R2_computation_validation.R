source("../f_statistics.R")


df <- fread("Bouras2022.ST2.CTACK.csv", data.table=F)
N.CTACK <- 3631 # From ST1

# df.1kg.frq <- fread("/data1/sanghyeon/wonlab_contribute/combined/data_common/reference_panel/1kGp3/EUR/reference.1kG.EUR.maf_0.005.geno_0.02.frq",
#                     data.table=F, select=c("SNP", "A1", "MAF")) %>%
#     rename(MAF.1kg=MAF)

# df1 <- df %>%
#     left_join(df.1kg.frq, by=c("rsid"="SNP"))

# View(df1)

df1 <- df %>%
    rowwise() %>%
    mutate(R2.self = compute_R2(BETA_exposure, SE_exposure, MAF, N.CTACK))
View(df1)

# 정확한 MAF은 모르기 때문에 R2가 완전히 동일하지는 않을 것.

df1 <- df1 %>%
    rowwise() %>%
    mutate(F.self = compute_f_stat(R2.self, N.CTACK, 1),
           F.self2 = (1 + (N.CTACK*R2.self)/(1-R2.self)))
