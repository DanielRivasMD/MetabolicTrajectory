####################################################################################################

module Vars
using ..Paths

####################################################################################################
# sigma experiment
####################################################################################################

const SMETA_xlsx    = joinpath(Paths.SIGMA, "sigma_groups.xlsx")

const SIG1R_HT_xlsx = joinpath(Paths.SXLSX, "ffl_2023-06-14_11_59_sig1r_transcriptoms_hetero_male_female_macro_OneClickMacroV.2.53.1-slice30mins.mac_1.xlsx")
const SIG1R_WT_xlsx = joinpath(Paths.SXLSX, "ffl_2023-06-20_12_14_sig1r_transcriptoms_sigko_wt_male_macro_OneClickMacroV.2.53.1-slice30mins.mac_1.xlsx")
const KO_WT_xlsx    = joinpath(Paths.SXLSX, "ffl_2023-06-27_11_15_ko_wt_female_macro_OneClickMacroV.2.53.1-slice30mins.mac_1.xlsx")
xvars_xlsx = ["DateTime", "DurationMin", "Animal"]

const SIG1R_HT_csv  = joinpath(Paths.SCSV,  "FFL_2023-06-14_11_59_Sig1R_Transcriptoms_Hetero_Male_Female_rt.csv")
const SIG1R_WT_csv  = joinpath(Paths.SCSV,  "FFL_2023-06-20_12_14_Sig1R_Transcriptoms_SigKo_WT_Male_rt.csv")
const KO_WT_csv     = joinpath(Paths.SCSV,  "FFL_2023-06-27_11_15_KO_WT_FEMALE_rt.csv")
xvars_csv = ["Date_Time"]

####################################################################################################
# HMGCR experiement
####################################################################################################

const C1_0120_xlsx  = joinpath(Paths.HXLSX, "cohort1_2025_01_20_male.xlsx")
const C1_0131_xlsx  = joinpath(Paths.HXLSX, "cohort1_2025_01_31_male.xlsx")
const C2_0306_xlsx  = joinpath(Paths.HXLSX, "cohort2_2025-03-06_male.xlsx")
const C2_0314_xlsx  = joinpath(Paths.HXLSX, "cohort2_2025_03_14_12_male.xlsx")
const C3_0321_xlsx  = joinpath(Paths.HXLSX, "cohort3_2025_03_21_female.xlsx")
const C4_0430_xlsx  = joinpath(Paths.HXLSX, "cohort4_2025_04_30_female.xlsx")
const C5_0522_xlsx  = joinpath(Paths.HXLSX, "cohort5_2025_05_22_male.xlsx")
const C6_0529_xlsx  = joinpath(Paths.HXLSX, "cohort5_2025-05-29_male.xlsx")

####################################################################################################

end

####################################################################################################
