
root -q -l -b finalPlot2016.C+'(1,1,"Emulated p_{T}^{miss} [GeV]","","inputs_zh_2016/histo_1.root","wz_fakemet_allcuts_postfit", 0,"",1,0,"WZ CR","","",1)';
root -q -l -b finalPlot2016.C+'(1,1,"Emulated p_{T}^{miss} [GeV]","","inputs_zh_2016/histo_0.root","zz_fakemet_allcuts_postfit", 0,"H(125)",1,0,"ZZ CR","","",1)';
root -q -l -b finalPlot2016.C+'(2,1,"p_{T}^{miss} [GeV]","","inputs_zh_2016/histo_2.root","fullsel_met_ll_postfit",1,"ZH(inv.)",1,0,"ee+#mu#mu","inputs_zh_2016/histo_2.root","#splitline{DM, Vector coupling}{m_{DM} = 150 GeV, m_{med} = 500 GeV}",1)';

sed -i 's/double SFBinWidth = 1/double SFBinWidth = 0.05/' StandardPlot2016.C;

root -q -l -b finalPlot2016.C+'(1,1,"BDT classifier","","inputs_zh_2016/histo_1_bdt.root","fullsel_bdt_wz_postfit",0,"ZH(125)",1,0,"WZ CR","","",1,0.2,.899)';
root -q -l -b finalPlot2016.C+'(1,1,"BDT classifier","","inputs_zh_2016/histo_0_bdt.root","fullsel_bdt_zz_postfit",0,"ZH(125)",1,0,"ZZ CR","","",1,0.2,.899)';
root -q -l -b finalPlot2016.C+'(2,1,"BDT classifier","","inputs_zh_2016/histo_2_bdt.root","fullsel_bdt_ll_postfit",0,"ZH(125)",1,0,"ee+#mu#mu","","",1,0.2,.899)';

sed -i 's/double SFBinWidth = 0.05/double SFBinWidth = 1/' StandardPlot2016.C;

root -l -b -q PlotLimit_ZHinv.C+'("inputs_zh_2016/ana_zhinv_bdt_nj.txt" ,"ana_hzinv_met_nj" ,"35.9 fb^{-1} (13 TeV)",110,1000,1,0,"ZH #rightarrow 2l+p_{T}^{miss} + #leq 1 jet",1,5,"pdf")';

root -l -b -q PlotLimit_ZHinv.C+'("inputs_zh_2016/ana_zhinv_bdt_nj_ratios.txt" ,"ana_hzinv_met_nj_ratios" ,"35.9 fb^{-1} (13 TeV)",110,300,0,0,"qq #rightarrow ZH #rightarrow 2l+p_{T}^{miss} + #leq 1 jet",1,4,"pdf")';
