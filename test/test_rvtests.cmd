spower LNR data/EA.txt -a .5 --def_valid_locus 2 1000 --sample_size 1000 --alpha 0.05 -v1 -o QT -r 100 -m "BurdenQt" "RVTests --name ZegginiTest" "RVTests --name SkatTest" "RVTests --name CMCTest" -j6
spower show QT.csv name method *power*
spower LOGIT data/EA.txt -a2 --def_valid_locus 2 1000 --sample_size 1000 --alpha 0.05 -v1 -o BT -r 100 -m "BurdenBt" "CFisher" "RVTests --name ZegginiTest --is-binary" "RVTests --name SkatTest --is-binary" "RVTests --name CMCTest --is-binary" "RVTests --name CMCFisherExactTest --is-binary"  -j6
spower show BT.csv name method *power*