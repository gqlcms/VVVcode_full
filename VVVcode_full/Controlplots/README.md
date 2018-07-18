cd ./PDFs
root
root [0] gSystem->AddIncludePath("-I$ROOFITSYS/include")
root [1] .L PdfDiagonalizer.cc++
root [2] .q

root
root [0] gSystem->AddIncludePath("-I$ROOFITSYS/include")
root [1] .L RooRelBWRunningWidth.cxx++
root [2] .q

root
root [0] gSystem->AddIncludePath("-I$ROOFITSYS/include")
root [1] .L HWWLVJRooPdfs.cxx++
root [2] .q

root
root [0] gSystem->AddIncludePath("-I$ROOFITSYS/include")
root [1] .L Util.cxx++
root [2] .q


python g1_exo_doFit.py --control -c el
python g1_exo_doFit.py --control -c mu

