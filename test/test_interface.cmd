# help
spower -h
spower LOGIT -h
spower PAR -h
spower BLNR -h
spower LNR -h
spower ELNR -h
spower show -h
spower execute -h
# 
spower LOGIT data/EA.txt -a 2 --def_valid_locus 2 1000 --sample_size 20000 --alpha 2.5E-6 -v1 -o tmp
spower LOGIT data/EA.txt -a 2 --def_valid_locus 2 1000 --power 0.8 --alpha 2.5E-6 -v1 -o tmp
spower LOGIT data/EA.txt -a 2 --def_valid_locus 2 1000 --sample_size 20000 --alpha 2.5E-6 -P 0.5 -r 100 -v1 -o tmp
spower LOGIT data/EA.txt -a 2 --def_valid_locus 2 1000 --power 0.8 --alpha 2.5E-6 -P 0.5 -r 100 -v1 -o tmp
# 
spower PAR data/EA.txt -a .05 --def_valid_locus 2 1000 --sample_size 20000 --alpha 2.5E-6 -v1 -o tmp
spower PAR data/EA.txt -a .05 --def_valid_locus 2 1000 --power 0.8 --alpha 2.5E-6 -v1 -o tmp
spower PAR data/EA.txt -a .05 --def_valid_locus 2 1000 --sample_size 20000 --alpha 2.5E-6 -P 0.5 -r 100 -v1 -o tmp
spower PAR data/EA.txt -a .05 --def_valid_locus 2 1000 --power 0.8 --alpha 2.5E-6 -P 0.5 -r 100 -v1 -o tmp
# 
spower LOGIT data/EA.txt --resampling -a 2 --def_valid_locus 2 1000 --sample_size 200 --alpha 2.5E-6 -v2 -o tmp -r 50 -m CFisher 
spower LOGIT data/EA.txt --resampling -a 2 --def_valid_locus 2 1000 --sample_size 200 --alpha 2.5E-6 -P 0.5 -r 50 -v2 -o tmp -m CFisher
# 
spower PAR data/EA.txt --resampling -a .05 --def_valid_locus 2 1000 --sample_size 200 --alpha 2.5E-6 -v2 -o tmp -r 50 -m CFisher 
spower PAR data/EA.txt --resampling -a .05 --def_valid_locus 2 1000 --sample_size 200 --alpha 2.5E-6 -P 0.5 -r 50 -v2 -o tmp -m CFisher
# 
spower BLNR data/EA.txt -a .05 --def_valid_locus 2 1000 --sample_size 200 --alpha 2.5E-6 -v2 -o tmp -r 50 -m CFisher 
spower BLNR data/EA.txt -a .05 --def_valid_locus 2 1000 --sample_size 200 --alpha 2.5E-6 -P 0.5 -r 50 -v2 -o tmp -m CFisher
# 
spower LNR data/EA.txt -a 1 --def_valid_locus 2 1000 --sample_size 20000 --alpha 2.5E-6 -v1 -o tmp
spower LNR data/EA.txt -a 1 --def_valid_locus 2 1000 --power 0.8 --alpha 2.5E-6 -v1 -o tmp
spower LNR data/EA.txt -a 1 --def_valid_locus 2 1000 --sample_size 20000 --alpha 2.5E-6 -P 0.5 -r 100 -v1 -o tmp
spower LNR data/EA.txt -a 1 --def_valid_locus 2 1000 --power 0.8 --alpha 2.5E-6 -P 0.5 -r 100 -v1 -o tmp
# 
spower ELNR data/EA.txt -a .05 --def_valid_locus 2 1000 --sample_size 200 --alpha 2.5E-6 -v2 -o tmp -r 50 -m BurdenQt 
spower ELNR data/EA.txt -a .05 --def_valid_locus 2 1000 --sample_size 200 --alpha 2.5E-6 -v2 -o tmp -r 50 -m BurdenQt --p1 0.5 
spower ELNR data/EA.txt -a .05 --def_valid_locus 2 1000 --sample_size 200 --alpha 2.5E-6 -P 0.5 -r 50 -v2 -o tmp -m BurdenQt
# TDT 
# spower M8 data/EA.txt -a 2 --def_valid_locus 2 1000 --sample_size 20000 --alpha 2.5E-6 -v1 -o tmp
# spower M8 data/EA.txt -a 2 --def_valid_locus 2 1000 --power 0.8 --alpha 2.5E-6 -v1 -o tmp
# spower M8 data/EA.txt -a 2 --def_valid_locus 2 1000 --sample_size 20000 --alpha 2.5E-6 -P 0.5 -r 100 -v1 -o tmp
# spower M8 data/EA.txt -a 2 --def_valid_locus 2 1000 --power 0.8 --alpha 2.5E-6 -P 0.5 -r 100 -v1 -o tmp
# show
spower show tests
spower show test CFisher
# RV tests
# spower LOGIT data/EA.txt -a 2 --def_valid_locus 2 1000 --sample_size 20000 --alpha 0.01 -v1 -o tmp -m "SSeq_rare --MAF 0.01"
# GDAT test
spower LOGIT data/boyko1.gdat --resampling -a 2 --def_valid_locus 2 1000 --sample_size 200 --alpha 2.5E-6 -v2 -o tmp -r 50 -m CFisher 
