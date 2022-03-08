# Complex_V = c("ATP5F1A","ATP5F1C","ATP5F1D","ATP5F1E","ATP5PB",
#               "ATP5MC1","ATP5MC2","ATP5MC3","ATP5PD","ATP5ME","ATP5PF",
#               "ATP5MF","ATP5MG","ATP5O")
Complex_V = c("ATP5A1(ATP5F1A)","ATP5C1(ATP5F1C)","ATP5D(ATP5F1D)","ATP5E(ATP5F1E)",
              "ATP5F1(ATP5PB)","ATP5G1(ATP5MC1)","ATP5G2(ATP5MC2)","ATP5G3(ATP5MC3)",
              "ATP5H(ATP5PD)","ATP5I(ATP5ME)","ATP5J(ATP5PF)","ATP5J2(ATP5MF)",
              "ATP5L(ATP5MG)","ATP5O")
Complex_V = sapply(Complex_V, function(x){
   a = str_split_fixed(x, "\\(", 2)[1,1]
   return(a)
})
Complex_V = as.character(Complex_V)
Complex_III = c("COQ9","UQCR10","UQCR11","UQCRB","UQCRC1","UQCRC2","UQCRFS1","UQCRQ")
Complex_IV = c("COX10","COX11","COX15","COX17","COX4I1","COX4I2","COX5A","COX5B","COX6A1",
               "COX6A2","COX6B1","COX6B2","COX6C","COX7A1","COX7A2","COX7A2L",
               "COX7B","COX7C","COX8A")
Complex_I = c("NDUFA1","NDUFA10","NDUFA11","NDUFA12","NDUFA13","NDUFA2","NDUFA3",
              "NDUFA4","NDUFA4L2","NDUFA5","NDUFA6","NDUFA7","NDUFA8","NDUFA9",
              "NDUFAB1","NDUFB1","NDUFB10","NDUFB11","NDUFB2","NDUFB3","NDUFB4",
              "NDUFB5","NDUFB6","NDUFB7","NDUFB8","NDUFB9","NDUFC1","NDUFC2",
              "NDUFS1","NDUFS2","NDUFS3","NDUFS4","NDUFS5","NDUFS6","NDUFS7",
              "NDUFS8","NDUFV1","NDUFV2","NDUFV3")
Complex_II = c("SDHA","SDHB","SDHC","SDHD","SUCLA2","SUCLG1","SUCLG2")
glycolyticAndTCA = c("ACO1","ACO2","ALDOA","ALDOB","ALDOC","BPGM","CS","CYB5R1",
                     "CYC1","DLAT","DLD","DLST","ECI1","ENO1","ENO2","ENO3","ETFB",
                     "FH","GAPDH","GCK","GOT1","GOT2","GPD1","GPD2","GPI","HK1",
                     "HK2","HK3","IDH1","IDH2","IDH3A","IDH3B","IDH3G","LDHA",
                     "LDHB","MDH1","MDH2","ME1","ME2","ME3","MICU1","MPC1","MPC2",
                     "MRPL49","OGDH","OGDHL","PC","PCK1","PCK2","PFKL","PFKM","PFKP",
                     "PGAM1","PGAM2","PGK1","PGK2","PKLR","PPA1","PPA2","SCO1",
                     "SLC25A11","SLC25A13","SLC2A1","SLC2A2","SLC2A3","SLC2A4","TPI1")
MRP_positive = c("MRPL16","MRPL49","MRPS35","MRPL46","MRPL10","MRPL39","MRPL19","MRPL33",
                 "MRPS9","MRPL44","MRPL1","MRPS18C","MRPS30","MRPS36","MRPS18A","MRPS18B",
                 "MRPS28","MRPL50","MRPS2")
MRP_negative = c("MRPL12","MRPL14","MRPL15","MRPL17","MRPL18","MRPL23","MRPL24","MRPL34",
                 "MRPL36","MRPL38","MRPL40","MRPL47","MRPL53","MRPL55","MRPS6","MRPS11",
                 "MRPS15","MRPS16","MRPS25")
mt_gene = c("MT-CO1","MT-CO2","MT-CO3","MT-ATP6","MT-ATP8",
            "MT-ND1","MT-ND3","MT-ND4","MT-ND4L")