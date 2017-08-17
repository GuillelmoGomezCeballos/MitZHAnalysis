#include <TString.h>
#include <TFile.h>
class FlatFile {
  public:
    int category;
    TString name;
    int signalIndex;
    TFile *file;
    FlatFile(TString name_, int category_, int signalIndex_=-1);
    FlatFile();
    ~FlatFile();
    //FlatFile(const FlatFile& theFlatFile) { name=theFlatFile.name; category=theFlatFile.category; signalIndex=theFlatFile.signalIndex; }
    TFile* Open() { file=TFile::Open(name, "READ"); return file; }
    void Close() { if(file) file->Close(); }
};
FlatFile::FlatFile(TString name_, int category_, int signalIndex_) {
  name=name_; category=category_; signalIndex=signalIndex_;
}
FlatFile::FlatFile() {
  name=""; category=-1; signalIndex=-999;
}
FlatFile::~FlatFile() {}

