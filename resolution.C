const int N_CH = 2692;
double fadc[N_CH];
double sum[5000];
double ch2keV = 5500./200000.; // 5.5 MeV alpha gives 200k channels in integral

//TH1F *h_reso;

void clear_fadc(){
  for(int ii=0;ii<N_CH;ii++)
    fadc[ii]=0;
}

double mean(int Nev){
  double m=0;
  for(int i=0;i<Nev;i++)
    m += sum[i];
  return ( m / double(Nev) );
}

double sdv(int Nev, double mean){
  double s = 0;
  for(int i=0;i<Nev;i++)
    s += pow( sum[i]-mean , 2 );
  return ( sqrt( s / double(Nev-1) ) );
}

void resolution(int gate=150 , int Nev=1000){

  int pos = 0;
  //h_reso = new TH1F("h_reso","");

  std::ifstream fDATA("./Dump.txt" , std::ios::in);
//  std::ifstream fDATA("/home/geant4/Data/Noise/MERGED/noise.data", std::ios::in);
  for(int ev=0;ev<Nev;ev++){
    clear_fadc();

    for(int ch=0;ch<N_CH;ch++){
      fDATA >> fadc[ch];
    }

//    if(ev==0){
//      for(int c=0;c<N_CH;c++)
//         std::cout << fadc[c] << " " ;
//      std::cout << "\n";
//    }

    pos = int( gRandom->Rndm()*(N_CH - gate - 10 ) );

    sum[ev]=0;
    for(int p=pos;p<pos+gate;p++)
      sum[ev] += fadc[p];

//    std::cout << "  --> " << sum[ev] << "\n";
  }

  double mn = mean(Nev    );
  double sd = sdv (Nev, mn);

  std::cout << " Nev  = " << Nev         << " events\n";
  std::cout << " GATE = " << gate * 0.04 << " us\n";
  std::cout << " MEAN = " << mn          << " a.u.\n";
  std::cout << " SDV. = " << sd          << " a.u.\n";
  std::cout << " RESO = " << sd*ch2keV   << " keV\n";

  gSystem->Exit(1);
}
