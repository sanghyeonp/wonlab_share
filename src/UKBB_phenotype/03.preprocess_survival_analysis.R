library(data.table)
library(dplyr)

### Follow-up loss: f.191
f.191 <- "UKBB.191.follow_up_loss.rds"
df.followup_loss <- readRDS(f.191) %>%
    mutate(f.191.0.0 = as.character(f.191.0.0)) %>%
    rename(date.follow_up_loss=f.191.0.0)
# df.followup_loss %>% filter(!is.na(date.follow_up_loss))

###  Inpatient record origin: f.40022
# f.40022 <- "UKBB.40022.inpatient_record_origin.rds"
# df.record_origin <- readRDS(f.40022) %>%
#     mutate(date.end_follow=ifelse(f.40022.0.0 %in% c("HES", "PEDW"), "2018-01-31", 
#                                   ifelse(f.40022.0.0 == "SMR", "2016-11-30", NA))) %>%
#     dplyr::select(f.eid, date.end_follow)
# # df.record_origin %>% filter(!is.na(date.end_follow))

### Health record origin
# 현재 (2024.11 기준) 사용하는 tab 파일은: ukb_uni20220812.tab
# 아래 URL에 따르면, 2022.05.19에 release 된 파일은 PEDW (Wales) = 2018년 2월 28일; HES( England/Scotland) = 2021년 9월 30일
# https://web.archive.org/web/20220519123145/https://biobank.ndph.ox.ac.uk/ukb/exinfo.cgi?src=Data_providers_and_dates
f.40022 <- "UKBB.40022.inpatient_record_origin.rds"
df.40022 <- readRDS(f.40022) %>%
    mutate(date.end_follow = ifelse((nchar(f.40022.0.0) == 0 | f.40022.0.0 == "PEDW") &
                                    (nchar(f.40022.0.1) == 0 | f.40022.0.1 == "PEDW") &
                                    (nchar(f.40022.0.2) == 0 | f.40022.0.2 == "PEDW"),
                                    "2018-02-28", "2021-09-30")) %>%
    dplyr::select(f.eid, date.end_follow)
# df.assessment %>% filter(is.na(date.end_follow))

###  Death date: f.40000
# Death date at instance 0
f.40000 <- "UKBB.40000.date_of_death.rds"
df.death <- readRDS(f.40000) %>%
    mutate(f.40000.0.0 = as.character(f.40000.0.0)) %>%
    rename(date.death=f.40000.0.0)
# df.death %>% filter(!is.na(date.death))

#### Follow-up censor date
df.follow_up_censor <- df.followup_loss %>%
    full_join(df.40022, by="f.eid") %>%
    full_join(df.death, by="f.eid")

df.follow_up_censor$date.follow_up_loss <- as.Date(df.follow_up_censor$date.follow_up_loss, format="%Y-%m-%d")
df.follow_up_censor$date.end_follow <- as.Date(df.follow_up_censor$date.end_follow, format="%Y-%m-%d")
df.follow_up_censor$date.death <- as.Date(df.follow_up_censor$date.death, format="%Y-%m-%d")
df.follow_up_censor$date.follow_up_end <- apply(df.follow_up_censor, 1, function(x){min(x[2], x[3], x[4], na.rm=T)})
head(df.follow_up_censor)

nrow(df.follow_up_censor %>% filter(!is.na(date.follow_up_end)))
nrow(df.follow_up_censor %>% filter(is.na(date.follow_up_end)))

saveRDS(df.follow_up_censor, "UKBB.follow_up_censor_date.f191.f40022.f40000.rds")
