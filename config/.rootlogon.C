{
  // TMVA Plotting
  TMVA::gConfig().fVariablePlotting.fUsePaperStyle=0; //make png files
  
  // Panda Objects
  gSystem->Load("libPandaAnalysisFlat.so");
  gSystem->Load("libPandaCoreDrawers.so");
  gSystem->Load("libPandaCoreLearning.so");
  gSystem->Load("libPandaCoreStatistics.so");
  gSystem->Load("libPandaCoreTools.so");
  gSystem->Load("libPandaTreeFramework.so");
  gSystem->Load("libPandaTreeObjects.so");
  gSystem->Load("libPandaTreeUtils.so");

}
