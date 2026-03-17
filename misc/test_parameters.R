library(methylKit)

file.list = list(
  "/home/users/allstaff/yan.a/davidson_longread/yan.a/emseq/data/CB_A_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz",
  "/home/users/allstaff/yan.a/davidson_longread/yan.a/emseq/data/CB_V_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz"
)

results = data.frame(mincov=integer(), difference=integer(), sig_regions=integer())

for (minc in c(5, 10)) {
  print(paste("Running with mincov =", minc))
  
  myobj = methRead(
    file.list,
    sample.id = list("CB_A_1", "CB_V_1"),
    assembly = "hg38",
    treatment = c(1, 0),
    context = "CpG",
    pipeline = "bismarkCoverage",
    header = FALSE,
    mincov = minc
  )
  
  meth = unite(myobj, destrand = FALSE, min.per.group = 1L)
  myDiff = calculateDiffMeth(meth, mc.cores = 16)
  
  for (diff_val in c(15, 20, 25)) {
    print(paste("  Filtering with difference =", diff_val))
    myDiff_filtered = getMethylDiff(myDiff, difference=diff_val, qvalue=0.01)
    
    # getData can return NULL if 0 rows, so we handle that safely
    dat = getData(myDiff_filtered)
    n_sig = if(is.null(dat)) 0 else nrow(dat)
    
    results = rbind(results, data.frame(mincov=minc, difference=diff_val, sig_regions=n_sig))
  }
}

print("Final Results:")
print(results)
write.csv(results, "/home/users/allstaff/yan.a/davidson_longread/yan.a/emseq/data/parameter_testing_results.csv", row.names=FALSE)
