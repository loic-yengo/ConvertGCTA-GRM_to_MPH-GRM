// STD include
#include <string.h>
#include <math.h>
#include <cstdio>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdlib>
#include <vector>

using namespace std;

int main(int argc, char *argv[]){
  int i,j;
  string sw;
  
  if(argc==1){
    cerr<<"\tArguments must be specified. Type --help for more details."<<endl;
    exit(1);
  }
  
  // Read arguments
  sw = argv[1];
  if (sw == "--help"){
    cerr<<"\t--grm     : Prefix for genetic relationship matrices (GRM)."<<endl;
    cerr<<"\t--out     : Specify prefix for output GRM."<<endl;
    exit(1);
  }else{
    if (argc == 1) {
      cerr<<"\tArguments must be specified. Type --help for more details."<<endl;
      exit(1);
    }
  }
  
  string grmBinFile;
  string grmIdFile;
  
  string grmInPrefix    = "";
  string grmOutPrefix = "";

  for(i = 1; i<argc;i++){
    sw = argv[i];
    if (sw == "--grm"){
      grmInPrefix = argv[i + 1];
    }
    if (sw == "--out"){
      grmOutPrefix = argv[i + 1];
    }
  }
  

  string line = "";
  int N       =  0;
  grmIdFile   = grmInPrefix+".grm.id";
  ifstream idStream;
  idStream.open(grmIdFile.c_str());

  string fid, iid;
  string grmOutIdFile = grmOutPrefix+".grm.iid";
  ofstream fileId(grmOutIdFile.c_str());
  while(idStream){
     getline(idStream,line);
     if(line!=""){
        stringstream ss;
        ss << line;
        ss >> fid;
	ss >> iid;
	fileId<<iid<<endl;
        N++;
     }
   }
   idStream.close();
   fileId.close();

   cout<<N<<" samples detected."<<endl;

  // Read GRM now
  int NP = N * (N + 1) / 2;
  float *GRM = new float[NP];

  //MatrixXf GRM = MatrixXf::Zero(N,N);  
  int size = sizeof (float);
  grmBinFile = grmInPrefix+".grm.bin";
  ifstream binData(grmBinFile.c_str(), ios::in | ios::binary);
  if(!binData){
    cerr << "[readGRM] Error reading file "<<grmBinFile<<endl;
    exit(1);
  }
  
  cout << "Reading the GRM from [" + grmBinFile + "]." << endl;

  float Nf     = (float) N;
  float Npairs = Nf * (Nf + 1.) / 2.;
  float g_buf  = 0.0;
  double Mo    = 0.;
  double Mo2   = 0.;
  double Md    = 0.;
  double Md2   = 0.;
  int index;
  for (i = 0; i < N; i++) {
    for (j = 0; j <= i; j++) {
      if (!(binData.read((char*) &g_buf, size))) throw ("\tError: the size of the [" + grmBinFile + "] file is incomplete?");
      //GRM(i,j) = g_buf;
      index = i * (i + 1) / 2 + j;
      GRM[index] = g_buf;
      if(j<i){
        Mo  += (double) g_buf;
        Mo2 += ((double) g_buf) * ((double) g_buf);
      }else{
        Md  += (double) g_buf;
        Md2 += ((double) g_buf) * ((double) g_buf);
      }
    }
  }
  

  Mo  = Mo  / Npairs;
  Mo2 = Mo2 / Npairs;
  Md  = Md  / N;
  Md2 = Md2 / N;
  double varOffDiag = Mo2 - Mo * Mo;
  double varDiag = Md2 - Md * Md;
    
  double Me = 2.0 / varOffDiag;
  cout<<"\nMean of diagonal elements: "<<Md<<".\n";
  cout<<"Mean of off-diagonal elements: "<<Mo<<".\n";
  cout<<"Variance of off-diagonal elements: "<<varOffDiag<<".\n";
  cout<<"Variance of diagonal elements: "<<varDiag<<".\n";
  cout<<"\nEstimated Me = "<<Me<<"."<<endl;
    
  binData.close();

  // Output GRM to MPH format
  float f_buf  = 0.;
  float sum2pq = 1.;
  string grmOutBinFile = grmOutPrefix + ".grm.bin";
  fstream A_Bin(grmOutBinFile.c_str(), ios::out | ios::binary);
  if (!A_Bin) throw ("Error: can not open the file [" + grmOutBinFile + "] to write.");
  A_Bin.write((char*) &N, sizeof(int));
  A_Bin.write((char*) &sum2pq, sizeof(float));

  for (i=0;i<N;i++){
    for (j=i;j<N;j++){
      index = j * (j + 1)/2 + i;	    
      f_buf = GRM[index];
      A_Bin.write((char*) &f_buf, size);
    }
  }
  A_Bin.close();
  cout << "GRM of " << N << " individuals has been saved in the file [" + grmOutBinFile + "] (in binary format)." << endl;

  delete [] GRM;
  return EXIT_SUCCESS;
}

