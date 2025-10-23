####################################################################################################

module Vars
using ..Paths

####################################################################################################
# sigma experiment
####################################################################################################

const SIG1R_HT_xlsx = joinpath(Paths.XLSX, "ffl_2023-06-14_11_59_sig1r_transcriptoms_hetero_male_female_macro_OneClickMacroV.2.53.1-slice30mins.mac_1.xlsx")
const SIG1R_WT_xlsx = joinpath(Paths.XLSX, "ffl_2023-06-20_12_14_sig1r_transcriptoms_sigko_wt_male_macro_OneClickMacroV.2.53.1-slice30mins.mac_1.xlsx")
const KO_WT_xlsx = joinpath(Paths.XLSX, "ffl_2023-06-27_11_15_ko_wt_female_macro_OneClickMacroV.2.53.1-slice30mins.mac_1.xlsx")
xvars_xlsx = ["DateTime", "DurationMin", "Animal"]

const SIG1R_HT_csv = joinpath(Paths.CSV, "FFL_2023-06-14_11_59_Sig1R_Transcriptoms_Hetero_Male_Female_rt.csv")
const SIG1R_WT_csv = joinpath(Paths.CSV, "FFL_2023-06-20_12_14_Sig1R_Transcriptoms_SigKo_WT_Male_rt.csv")
const KO_WT_csv = joinpath(Paths.CSV, "FFL_2023-06-27_11_15_KO_WT_FEMALE_rt.csv")
xvars_csv = ["Date_Time"]

####################################################################################################
# HMGCR experiement
####################################################################################################

const COHORT1_0120 = joinpath(Paths.HMGCR, "cohort1_2025_01_20_male.xlsx")
const COHORT1_0131 = joinpath(Paths.HMGCR, "cohort1_2025_01_31_male.xlsx")
const COHORT2_0306 = joinpath(Paths.HMGCR, "cohort2_2025-03-06_male.xlsx")
const COHORT2_0314 = joinpath(Paths.HMGCR, "cohort2_2025_03_14_12_male.xlsx")
const COHORT3_0321 = joinpath(Paths.HMGCR, "cohort3_2025_03_21_female.xlsx")
const COHORT4_0430 = joinpath(Paths.HMGCR, "cohort4_2025_04_30_female.xlsx")
const COHORT5_0522 = joinpath(Paths.HMGCR, "cohort5_2025_05_22_male.xlsx")
const COHORT6_0529 = joinpath(Paths.HMGCR, "cohort5_2025-05-29_male.xlsx")

end

####################################################################################################
