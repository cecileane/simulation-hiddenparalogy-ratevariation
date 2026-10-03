# Calibrate the branches of the angiosperm species tree that ASTRAL
# cannot estimate (the 8 external branches and the ingroup stem) in CU
#
# Idea: convert CU to substitutions per site, which we estimated ourselves
# on every branch from the 260 gene trees.
# One ratio k = CU / substitutions, from the stem of the water lilies, gives
# the water-lily tips
# the rest of the tree is made ultrametric and the ingroup stem follows 
# from the substitutions on the path from Pinus to the ingroup.

library(ape)

#------------------------------------------------------------------#
# 1. observations: the 5 internal branches of the ASTRAL tree
#------------------------------------------------------------------#
# ASTRAL coalescent units (data/kew_a353/astral/species.tre, 
# just use 2 digits)
# and substitutions per site on the same branches
# Below coal_obs shows the internal branches estimated from astral 
#            Chl|Lir  Cer|..  Sch|..  Amb|rest NB stem
coal_obs  = c(0.44,   0.95,   0.79,   0.16,    3.77)
subst_u   = c(0.0224, 0.0330, 0.0311, 0.0065,  0.1721) 
branch    = c("Cer|Chl+Lir", "Sch|Cer..", "Nym|Sch..", "Amb|rest", "NB stem")

# substitutions per site on the external branches
# Pinus: the path Pinus tip -> ingroup crown (Pinus branch + ingroup stem):
# only the sum, 0.227 + 0.227 = 0.454.
subst_ext = c(Liriodendron = 0.135, Chloranthus = 0.147, Ceratophyllum = 0.32,
              Schisandra = 0.167, Brasenia = 0.136, Nymphaea = 0.131,
              Amborella = 0.254, Pinus_path = 0.454)

round(coal_obs / subst_u, 1)
# 19.6 28.8 25.4 24.6 21.9 : nearly constant from around 19 to 25

# the least-squares slope through the origin using all five branches. 
k_all= sum(coal_obs * subst_u) / sum(subst_u^2)   # 22.2, all 5

# k_backbone = the same fit using only the four backbone branches, 
# leaving out the water-lily stem. 
# If the water-lily stem behaved differently from the backbone, 
# these two numbers would be far apart.They are not.
k_backbone = sum((coal_obs * subst_u)[1:4]) / sum(subst_u[1:4]^2) # 25.7
fit_log = lm(log(coal_obs) ~ log(subst_u))
coef(fit_log) # slope 0.98: slope 1 on log-log axes = proportional
# All three agree with the simple ratio of the water-lily stem alone, 
# 3.77 / 0.172 = 21.9 (SEE BELOW), so it is safe to use the water-lily stem 
# alone to compute the ratio (see below). 

#------------------------------------------------------------------#
# 2. check that substitutions are also proportional to time on every branch
#------------------------------------------------------------------#
# node ages from Zuntini et al. 2024, old tree and Morris et al.
# 2018 for the root (348 Ma), here Zuntini et al. does not have the root age:
age_root  = 348  # crown seed plants: Pinus | angiosperms
age_crown = 247.0 # angiosperm crown: Amborella | rest
age_Nym  = 245.7 # Nymphaeales split: (Nymphaea, Brasenia) | rest
age_Sch  = 244.5 # Schisandra | rest
age_Cer  = 243.5 # Ceratophyllum | (Chloranthus, Liriodendron)
age_ChlLir = 238.3 # Chloranthus | Liriodendron
age_NB = 126.0 # Nymphaea | Brasenia

# time span of a branch = age of its parent node - age of its child node (mya)
# internal branches, in the order of coal_obs, then the external branches
Ma_branch = c("Cer|Chl+Lir" = age_Cer - age_ChlLir,   # 5.2
              "Sch|Cer.."  = age_Sch - age_Cer,     # 1.0
              "Nym|Sch.."  = age_Nym - age_Sch,     # 1.2
              "Amb|rest"   = age_crown - age_Nym,    # 1.3
              "NB stem"    = age_Nym - age_NB,     # 119.7
              Liriodendron = age_ChlLir,         # 238.3
              Chloranthus  = age_ChlLir,        # 238.3
              Ceratophyllum = age_Cer,      # 243.5
              Schisandra  = age_Sch,      # 244.5
              Brasenia  = age_NB,        # 126.0
              Nymphaea   = age_NB,     # 126.0
              Amborella  = age_crown, # 247.0
              # Pinus path = Pinus tip (348) + ingroup stem (348 - 247)
              # Pinus path means the length from ingroup to the tip of Pinus 
              Pinus_path    = age_root + (age_root - age_crown)) # 449
su_branch = c(subst_u, subst_ext)   # substitutions on the same branches
names(su_branch) = names(Ma_branch)
round(su_branch / Ma_branch, 5)# substitutions per Mya
# Cer|Chl+Lir 0.00431  Sch|Cer.. 0.03300  Nym|Sch.. 0.02592
# Amb|rest 0.00500     NB stem 0.00144
# Liriodendron 0.00057 Chloranthus 0.00062 Ceratophyllum 0.00131
# Schisandra 0.00068   Brasenia 0.00108    Nymphaea 0.00104
# Amborella 0.00103    Pinus path 0.00101

# Some notes from the above substitutions per Mya: 
# On the 9 long branches (NB stem, the tips, the Pinus path) the rate is
# 0.0006-0.0014 substitutions per Mya: constant within a factor of 2.5, which
# is the lineage rate variation (Ceratophyllum fast, Chloranthus and
# Liriodendron slow). So substitutions are a usable clock for the branches
# ASTRAL cannot see.

# On the 4 backbone branches the published ages give 0.004-0.033, i.e. 4 to
# 30 times faster. These branches are short in million years (1-9 Mya in
# total) because treePL from the original paper 
# stacked their nodes against the maximum age imposed
# at the angiosperm crown, so their time spans cannot anchor a rate. 
# This could be why we calibrate against substitutions, which we
# measured ourselves. 

#------------------------------------------------------------------#
# 3. the method: one ratio k, from the NB stem
#------------------------------------------------------------------#
# k from the NB stem: the longest branch measured in both units. The 4
# backbone ratios (20-29) agree with it, which checks the calibration and
# argues against hidden paralogy having inflated the backbone discordance
# (that would push those ratios below 22).
# sensitivity: k = 20 to 29 gives water-lily tips 2.7-3.9, stem 1.3-2.8,
# height 7.8-10.4 CU.
k = 22 # = CU / substitutions of the NB stem = 3.77 / 0.172 = 21.9

coal_ext = k * subst_ext
round(coal_ext, 2)
# 2.97 3.23 7.04 3.67 2.99 2.88 5.59 9.99
# These cannot all be used as they are: an ultrametric tree needs tips at
# the same height, and the substitution lengths vary with the lineage rate
# (Ceratophyllum is fast). So we use the two water-lily tips (the cherry on
# the calibration branch) to set the height, and the other tips follow from
# ultrametricity.
# Using cherry for ultrametrizing the tree is the same procedure as our 
# reptile tree. 

#------------------------------------------------------------------#
# 4. build the tree in coalescent units
#------------------------------------------------------------------#
cu_NBtip   = k * mean(subst_ext[c("Nymphaea", "Brasenia")]) # 2.94, both tips
# height of the ingroup: from the angiosperm crown down to the tips

# height of the ingroup: from the angiosperm crown down to the tips 
H_ingroup  = coal_obs[4] + coal_obs[5] + cu_NBtip # 0.16 + 3.77 + 2.94 = 6.87
# ingroup stem: the path Pinus tip -> ingroup crown is k * 0.454 = 9.99 CU
# and, the tree being ultrametric, Pinus tip = stem + H_ingroup.
# So 2 * stem + H_ingroup = 9.99.

# This calculates the stem height from root to ingroup 
cu_stem    = unname(coal_ext["Pinus_path"] - H_ingroup) / 2 # 1.56 CU
H_total    = H_ingroup + cu_stem  # 8.43 = Pinus tip
# other tips: height minus the depth of their parent node (ASTRAL branches)
cu_Amborella     = H_ingroup   # 6.87
cu_Schisandra    = H_ingroup-0.16-0.79   # 5.92
cu_Ceratophyllum = H_ingroup-0.16-0.79-0.95    # 4.97
cu_ChlLir        = H_ingroup-0.16-0.79-0.95-0.44   # 4.53

tree_cu = sprintf(paste0(
  "(Pinus:%.2f,(((Nymphaea:%.2f,Brasenia:%.2f):3.77,",
  "(Schisandra:%.2f,(Ceratophyllum:%.2f,(Chloranthus:%.2f,Liriodendron:%.2f)",
  ":0.44):0.95):0.79):0.16,Amborella:%.2f):%.2f);"),
  H_total, cu_NBtip, cu_NBtip, cu_Schisandra, cu_Ceratophyllum,
  cu_ChlLir, cu_ChlLir, cu_Amborella, cu_stem)
cat(tree_cu, "\n")
# (Pinus:8.43,(((Nymphaea:2.94,Brasenia:2.94):3.77,(Schisandra:5.92,
#  (Ceratophyllum:4.97,(Chloranthus:4.53,Liriodendron:4.53):0.44):0.95):0.79)
#  :0.16,Amborella:6.87):1.56);
# This was tested by QuartetNetworkGoodnessFit.ultrametrize! 
# see step 5 of scripts/speciestree/speciestree_angiosperm.jl 

writeLines(tree_cu, "data/kew_a353/calibration/speciestree_cu_from_subs.tre") 

#------------------------------------------------------------------#
# 5. generations for SimPhy: same height as the reptile tree (3440)
#------------------------------------------------------------------#
Ne2 = 50 * round(3440 / H_total / 50) # 2Ne = 3440 / 8.43 = 408 -> 400
Ne2
round(c(backbone = sum(coal_obs[1:4]), stem = cu_stem, NB_stem = 3.77,
        NB_tip = cu_NBtip, Pinus = H_total) * Ne2)
# backbone 936, stem 624, NB stem 1508, NB tip 1175, Pinus 3371 generations
# (from unrounded CU; speciestree_angiosperm.jl rounds the tree first:
#  1176, 3372)
# (reptile tree: height 3440, stem 500, internal branches 170-2440)

# substitutions on the stem and the Pinus branch: only their sum (0.454) is
# known. Under this approach the rate is the same on both, so split in
# proportion to CU:
round(0.454 * c(stem = cu_stem, Pinus = H_total) / (cu_stem + H_total), 3)
# stem 0.071, Pinus 0.383 (instead of 0.227 / 0.227 from the least squares)

#------------------------------------------------------------------#
# 6. figure: the fit and the tree
#------------------------------------------------------------------#
png("data/kew_a353/calibration_subs.png", width = 16, height = 7,
    units = "in", res = 130)
layout(matrix(1:2, 1, 2), widths = c(1, 1.3))
par(mar = c(5, 5, 4, 1))
plot(subst_u, coal_obs, pch = 19, cex = 1.6, xlim = c(0, 0.2),
     ylim = c(0, 4.5), xaxs = "i", yaxs = "i", cex.lab = 1.2,
     xlab = "substitutions per site (from the 260 gene trees)",
     ylab = "ASTRAL branch length (coalescent units)")
abline(0, k, lwd = 2, col = "firebrick")
abline(0, k_all, lwd = 1, lty = 2, col = "grey40")
text(subst_u, coal_obs, branch, pos = c(4, 4, 2, 4, 2), cex = 1)
legend("topleft", bty = "n", cex = 1.1,
       legend = c(sprintf("CU = %d x subs (ratio of the NB stem)", k),
                  sprintf("least squares, 5 branches: %.1f", k_all)),
       col = c("firebrick", "grey40"), lty = c(1, 2), lwd = c(2, 1))
title("A. the 5 internal branches: CU vs substitutions", cex.main = 1.3)
par(mar = c(4,1,4,1))
tr = read.tree(text = tree_cu)
plot(tr, edge.width = 2, cex = 1.2, label.offset = 0.1, x.lim = c(0, 10.5))
edgelabels(sprintf("%.2f", tr$edge.length), frame = "none",
           adj = c(0.5, -0.4), cex = 1)
axisPhylo(); mtext("coalescent units", side = 1, line = 2.3)
title(sprintf(paste("B. species tree in CU: water-lily tips = %d x subs,",
                    "stem from the Pinus path, rest ultrametric"), k),
      cex.main = 1.2)
mtext(sprintf("height %.2f CU, stem %.2f CU; 2Ne = %d for 3440 generations",
              H_total, cu_stem, Ne2), side = 3, line = 0.2, cex = 1)
invisible(dev.off())
