fsz=8
fszxtick=7
fszytick=7
fszaxlab=9
lw=0.7
mw=0.5
msz=3.2
lgdw=1.0
lgdfsz=7

# python burkelab_SimScripts/simulateIDT_Shao_1x4_ESSCI.py \
# --figwidth 6.5 --figheight 2.5 --fsz $fszaxlab --fszxtick $fszxtick --fszytick $fszytick --fszaxlab $fszaxlab \
# --lw $lw --mw $mw --msz $msz --lgdw $lgdw --lgdfsz 6 --gridsz 9 &

# python burkelab_SimScripts/simulateIDT_Shao_4x1_PCI.py \
# --figwidth 2.5 --figheight 6.66667 --fsz $fszaxlab --fszxtick $fszxtick --fszytick $fszytick --fszaxlab $fszaxlab \
# --lw $lw --mw $mw --msz $msz --lgdw $lgdw --lgdfsz 7 --gridsz 9 &

# # python burkelab_SimScripts/simulateJSR_H2O_1x3_ESSCI.py \
# # --figwidth 6.5 --figheight 2.5 --fsz $fsz --fszxtick $fszxtick --fszytick $fszytick --fszaxlab $fszaxlab \
# # --lw $lw --mw $mw --msz $msz --lgdw $lgdw --lgdfsz 6 --gridsz 50 &

# python burkelab_SimScripts/simulateJSR_H2O_3x1_PCI.py \
# --figwidth 2.5 --figheight 5 --fsz $fszxtick --fszxtick $fszxtick --fszytick $fszytick --fszaxlab $fszaxlab \
# --lw $lw --mw $mw --msz $msz --lgdw $lgdw --lgdfsz 7 --gridsz 50 &

# # python burkelab_SimScripts/simulateJSR_NH3_1x3_ESSCI.py \
# # --figwidth 6.5 --figheight 2.5 --fsz $fsz --fszxtick $fszxtick --fszytick $fszytick --fszaxlab $fszaxlab \
# # --lw $lw --mw $mw --msz $msz --lgdw $lgdw --lgdfsz $lgdfsz --gridsz 50 &

# python burkelab_SimScripts/simulateJSR_NH3_3x1_PCI.py \
# --figwidth 2.5 --figheight 5 --fsz $fsz --fszxtick $fszxtick --fszytick $fszytick --fszaxlab $fszaxlab \
# --lw $lw --mw $mw --msz $msz --lgdw $lgdw --lgdfsz 7 --gridsz 50 &

# python burkelab_SimScripts/simulateshocktubeShao_1x1_ESSCI.py \
# --figwidth 5.5 --figheight 2.5 --fsz $fsz --fszxtick $fszxtick --fszytick $fszytick --fszaxlab $fszaxlab \
# --lw $lw --mw $mw --msz $msz --lgdw $lgdw --lgdfsz $lgdfsz --gridsz 14 &

python burkelab_SimScripts/simulateshocktubeShao_1x1_PCI.py \
--figwidth 2.5 --figheight 2.5 --fsz $fszaxlab --fszxtick $fszxtick --fszytick $fszytick --fszaxlab $fszaxlab \
--lw $lw --mw $mw --msz $msz --lgdw $lgdw --lgdfsz 7 --gridsz 3 &


# # # python burkelab_SimScripts/simulateIDT_Keromnes.py \
# # # --figwidth 3.5 --figheight 6.66667 --fsz $fsz --fszxtick $fszxtick --fszytick $fszytick --fszaxlab $fszaxlab \
# # # --lw $lw --mw $mw --msz $msz --lgdw $lgdw --lgdfsz 7 --gridsz 9 &

# # ## FLAME SPEED CALCULATIONS
# # # python burkelab_SimScripts/simulateflamespeedRonney_NH3_H2.py --gridsz 30 --date 'May28' --slopeVal 0.01 --curveVal 0.01 --transport 'multicomponent'
# # # python burkelab_SimScripts/simulateflamespeedRonney_NH3_H2.py --gridsz 30 --date 'May28' --slopeVal 0.05 --curveVal 0.05 --transport 'multicomponent'
# python burkelab_SimScripts/simulateflamespeedBurke.py --gridsz 15 --date 'Jun02' --slopeVal 0.05 --curveVal 0.05 --transport 'multicomponent'
# # # python burkelab_SimScripts/simulateflamespeedBurke.py --gridsz 40 --date 'May28' --slopeVal 0.05 --curveVal 0.05 --transport 'multicomponent'
# # # python burkelab_SimScripts/simulateflamespeedRonney_NH3_H2.py --gridsz 40 --date 'May31' --slopeVal 0.05 --curveVal 0.05 --transport 'multicomponent'
# # # python burkelab_SimScripts/simulateflamespeedRonney_NH3only.py --gridsz 50 --date 'May30' --slopeVal 0.05 --curveVal 0.05 --transport 'multicomponent'

# # # ## PLOT FLAME SPEED DATA
# python burkelab_SimScripts/simulateflamespeedBurke_FromData_ESSCI.py \
# --figwidth 5.5 --figheight 2.5 --fsz $fszaxlab --fszxtick $fszxtick --fszytick $fszytick --fszaxlab $fszaxlab \
# --lw $lw --mw $mw --msz $msz --lgdw $lgdw --lgdfsz 7 --date 'Jun02' --slopeVal 0.05 --curveVal 0.05 &

# python burkelab_SimScripts/simulateflamespeedBurke_FromData_PCI.py \
# --figwidth 2.5 --figheight 2.5 --fsz $fszaxlab --fszxtick $fszxtick --fszytick $fszytick --fszaxlab $fszaxlab \
# --lw $lw --mw $mw --msz $msz --lgdw $lgdw --lgdfsz 7 --date 'May31' --slopeVal 0.01 --curveVal 0.01 --title 'Burke/Song (slope=0.01 curve=0.01)' &

# python burkelab_SimScripts/simulateflamespeedRonney_0.6NH3_0.4H2_FromData.py \
# --figwidth 6.5 --figheight 2.72 --fsz $fszaxlab --fszxtick $fszxtick --fszytick $fszytick --fszaxlab $fszaxlab \
# --lw $lw --mw $mw --msz $msz --lgdw $lgdw --lgdfsz $lgdfsz --date 'May31' --slopeVal 0.05 --curveVal 0.05 --paper 'ESSCI' &

# python burkelab_SimScripts/simulateflamespeedRonney_0.6NH3_0.4H2_FromData.py \
# --figwidth 2.5 --figheight 4.5 --fsz $fszaxlab --fszxtick $fszxtick --fszytick $fszytick --fszaxlab $fszaxlab \
# --lw $lw --mw $mw --msz $msz --lgdw $lgdw --lgdfsz 7 --date 'May31' --slopeVal 0.05 --curveVal 0.05 --paper 'PCI' &

# python burkelab_SimScripts\\ChemicalMechanismCalculationScripts\\AllArrheniusFits.py \
# --figwidth 6.5 --figheight 6 --fsz 7 --fszxtick $fszxtick --fszytick $fszytick --fszaxlab $fszaxlab \
# --lw $lw --mw $mw --msz $msz --lgdw $lgdw --lgdfsz 7 &

# python burkelab_SimScripts/simulateflamespeedRonney_NH3_H2_FromData.py \
# --figwidth 6.5 --figheight 4 --fsz 5 --fszxtick $fszxtick --fszytick $fszytick --fszaxlab $fszaxlab \
# --lw $lw --mw $mw --msz $msz --lgdw $lgdw --lgdfsz 5 --date 'May28' --slopeVal 0.03 --curveVal 0.06 --title 'Ronney (slope=0.03 curve=0.06)' &

# python burkelab_SimScripts/simulateflamespeedRonney_NH3only_FromData.py \
# --figwidth 3.5 --figheight 2.72222 --fsz 5 --fszxtick $fszxtick --fszytick $fszytick --fszaxlab $fszaxlab \
# --lw $lw --mw $mw --msz $msz --lgdw $lgdw --lgdfsz $lgdfsz --date 'May30' --slopeVal 0.05 --curveVal 0.05 --paper 'ESSCI' &

# python burkelab_SimScripts/simulateflamespeedRonney_NH3only_FromData.py \
# --figwidth 3.5 --figheight 2.72222 --fsz 5 --fszxtick $fszxtick --fszytick $fszytick --fszaxlab $fszaxlab \
# --lw $lw --mw $mw --msz $msz --lgdw $lgdw --lgdfsz $lgdfsz --date 'May30' --slopeVal 0.05 --curveVal 0.05 --paper 'PCI' &


## USEFUL COMMANDS
# # Convert a cti to a yaml
# python interfaces\\cython\\cantera\\cti2yaml.py "G:\\Mon disque\\Columbia\\Burke Lab\\07 Mechanisms\\09 Nitrogen\\Shrestha\\shrestha2018.cti" 

# # To make this file executable:
# chmod +x burkelab_SimScripts/burkelabSims.sh

# # To run this file:
# ./burkelab_SimScripts/burkelabSims.sh
