#define TRACELENGTH 2692
#define CHANNELS 70

class tpcEvent
{
  public:
    tpcEvent();
   ~tpcEvent();

    void SetChannel(int anod, double x_val, double y_val );
    void CheckAll();
    void Reset();
    void CorrectBaseLevel(int anod);
    void EvalBaseDev(int anod);
    void FindTime(int anod );
    void FindEnergyDeposit(int anod );

    int  GetEventNumber();
    int  GetWindowWidth();
    double  GetThresholdFactor();
    double  GetBaseLevel(int anod);

    void DumpSpectrum(int anod, int length_ch);

    void SetEventNumber(int N );
    void SetWindowWidth(int w );
    void SetThresholdFactor(double t );

    void Print();

    void PrepareAnimation(int anod );
    void AnimatedDisplay( );

    void DrawHist(int anod );
    void DrawDisplay( );
    void DrawAll( );
    void DrawTwo( int anod1, int anod2 );
    void DrawSector(int sector);

  private:

    int peakPosition[CHANNELS];
    int timeStamp[CHANNELS];
    double energyDeposit[CHANNELS];
    double baseLevel[CHANNELS];
    double baseDev[CHANNELS];
    double FADC[CHANNELS][TRACELENGTH];
    int  color[CHANNELS];
    int ev;
    bool fired[CHANNELS];
    TText  ttext[CHANNELS];
    TText ev_ttext;
    TText an_ttext;
    TText tm_ttext;
    TText pk_ttext;
    TText en_ttext;
    TLine ll;
    TText  sec_ttext[16];
    int window    = 10;
    double threshold = 1.02;
    int n_anods = CHANNELS;
    int n_ch    = TRACELENGTH;
};

tpcEvent::tpcEvent(){
  for(int i = 0; i < CHANNELS; i++){
    fired[i]=false;
    baseLevel[i]=0;
    baseDev[i]=0;
    peakPosition[i]=0;
    timeStamp[i]=0;
    for(int j =0;j<n_ch;j++)
      FADC[i][j] = 0;
  }

  int dd=1;
  for(int sector=0;sector<16;sector++){
    for(int dist=0;dist<4;dist++){
      if(dist==0) { color[dd] = 209; }
      if(dist==1) { color[dd] =   4; }
      if(dist==2) { color[dd] =   6; }
      if(dist==3) { color[dd] =  46; }
      dd++;
    }
  }
  color[65]=1; color[66] = 2;
  window    = 10;
  threshold = 1.02;
  n_anods = CHANNELS;
  n_ch    = TRACELENGTH;
};

tpcEvent::~tpcEvent(){ };

void tpcEvent::Reset(){
  for(int i = 0; i < CHANNELS; i++){
    fired[i]=false;
    baseLevel[i]=0;
    peakPosition[i]=0;
    timeStamp[i]=0;
    for(int j =0;j<n_ch;j++)
      FADC[i][j] = 0;
  }

  int dd=1;
  for(int sector=0;sector<16;sector++){
    for(int dist=0;dist<4;dist++){
      if(dist==0) { color[dd] = 209; }
      if(dist==1) { color[dd] =   4; }
      if(dist==2) { color[dd] =   6; }
      if(dist==3) { color[dd] =  46; }
      dd++;
    }
  }
  color[65]=1; color[66] = 2;
  window    = 10;
  threshold = 1.02;
  n_anods = CHANNELS;
  n_ch    = TRACELENGTH;
};

int    tpcEvent::GetEventNumber(){ return ev; };

int    tpcEvent::GetWindowWidth(){ return window; };

double tpcEvent::GetThresholdFactor(){return threshold;};

double tpcEvent::GetBaseLevel(int anod){return baseLevel[anod];};

void   tpcEvent::SetEventNumber(int N){ ev = N; };

void   tpcEvent::SetWindowWidth(int w){ window = w; };

void   tpcEvent::SetThresholdFactor(double t){ threshold = t; };

void   tpcEvent::DumpSpectrum(int anod, int length_ch){
  if( !fired[anod] ){
    ofstream dump_file;
    dump_file.open("dump.txt", std::ios_base::app);
    for(int ch=0;ch<length_ch;ch++){
      dump_file << FADC[anod][ch] ;
      if( (ch+1) != length_ch) dump_file << ",";
      else dump_file << "\n";
    }
    dump_file.close();
  }
};

void   tpcEvent::Print(){
  cout << "Event     : " << ev        << endl;
  cout << "Window    : " << window    << endl;
  cout << "Threshold : " << threshold << endl;
};
void tpcEvent::SetChannel(int anod, double x_val, double y_val ){
  int ch = (x_val-1.);
  FADC[anod][ch] = y_val;
};

// DYSPLAY
void tpcEvent::DrawDisplay(){
  TCanvas* cnvd = new TCanvas("cnvd","display",900,900);
  TCrown* cr[67];
  for(int ii=0;ii<67;ii++){
    cr[ii] = new TCrown(0.5,0.5,0,0);
    cr[ii]->SetLineColor( color[ii] );
    if( fired[ii] ){
      cr[ii]->SetFillStyle(3004);
      cr[ii]->SetFillColor( color[ii] );
    }
  }
  char *txt = new char[256];
  double D2R = 3.14159265 / 180.;
  double cx = 0.5; double cy = 0.5;
  double r_off = 0.003; double phi_off = 0.5;
  int num=1;  double phi_min, phi_max, r_in, r_out, r_av, phi_av;
  for(int sector=0;sector<16;sector++){
    phi_min = 90.+(360./16.)*sector     + phi_off;
    phi_max = 90.+(360./16.)*(sector+1) - phi_off;
    phi_av = 0.5*(phi_min + phi_max);
    for(int dist=0;dist<4;dist++){
      if(dist<3){ r_in = 0.1 + 0.1*dist + r_off; r_out = 0.1 + 0.10*(dist+1)    - r_off; }
      else      { r_in = 0.1 + 0.1*dist + r_off; r_out = 0.1 + 0.10*dist + 0.05 - r_off; }
      sprintf(   txt, "%d", num);
      r_av   = 0.5*(r_out + r_in) ;
      ttext[num].SetText( cx + r_av*cos(phi_av*D2R) , cy + r_av*sin(phi_av*D2R) , txt );
      ttext[num].SetTextSize(0.014);

 //     cout << " " << num << " " << color[num] << " " << fired[num] << "\n";
  //    if( fired[num] ){
  //      cr[num]->SetFillStyle(3004);
 //       cr[num]->SetFillColor( color[num] );
 //     }

      cr[num]->DrawCrown(cx, cy, r_in, r_out, phi_min, phi_max);
      ttext[num].Draw();
      num++;
    }

    r_av = 0.475;
    sprintf(   txt, "%d", sector+1);
    sec_ttext[sector].SetText( cx + r_av*cos(phi_av*D2R) , cy + r_av*sin(phi_av*D2R) , txt );
    sec_ttext[sector].SetTextSize(0.03);
    sec_ttext[sector].Draw();
  }
  cr[65]->DrawCrown(cx, cy, 0.00           , 0.05 - r_off/2. );
  cr[66]->DrawCrown(cx, cy, 0.05 + r_off/2., 0.10 - r_off/2. );

  ttext[65].SetText( 0.51 , 0.51 , "65" );
  ttext[65].SetTextSize(0.014);
  ttext[66].SetText( 0.43 , 0.52 , "66" );
  ttext[66].SetTextSize(0.014);
  ttext[65].Draw();
  ttext[66].Draw();

  sprintf(   txt, "EVENT %d", ev);
  ev_ttext.SetText(0.01, 0.95, txt);
  ev_ttext.Draw();

  cnvd->Print("display.png");
  cnvd->Close();
};

//ANIMATED
void tpcEvent::AnimatedDisplay(){

  TString fig_name;
  for(int frame=0;frame<1000;frame++){

  if(              frame<  10 ) fig_name.Form("disp_000%d.png",frame);
  if( frame>  9 && frame< 100 ) fig_name.Form("disp_00%d.png" ,frame);
  if( frame> 99 && frame<1000 ) fig_name.Form("disp_0%d.png"  ,frame);
  if( frame>999               ) fig_name.Form("disp_%d.png"   ,frame);

  TCanvas* cnvd = new TCanvas("cnvd","display",900,900);
  cnvd->SetGrayscale();
  TCrown* cr[67];
  for(int ii=0;ii<67;ii++){
    cr[ii] = new TCrown(0.5,0.5,0,0);
    if(ii==66) cout <<  632 +  int( 10.*(FADC[ii][frame] - baseLevel[ii])/baseDev[ii]) << endl;
    cr[ii]->SetLineColor(  632 +  int( 10.*(FADC[ii][frame] - baseLevel[ii])/baseDev[ii]) );
    if( fired[ii] ){
      cr[ii]->SetFillColor( 632 +  int( 10.*(FADC[ii][frame] - baseLevel[ii])/baseDev[ii]) );
    }
  }
  char *txt = new char[256];
  double D2R = 3.14159265 / 180.;
  double cx = 0.5; double cy = 0.5;
  double r_off = 0.003; double phi_off = 0.5;
  int num=1;  double phi_min, phi_max, r_in, r_out, r_av, phi_av;
  for(int sector=0;sector<16;sector++){
    phi_min = 90.+(360./16.)*sector     + phi_off;
    phi_max = 90.+(360./16.)*(sector+1) - phi_off;
    phi_av = 0.5*(phi_min + phi_max);
    for(int dist=0;dist<4;dist++){
      if(dist<3){ r_in = 0.1 + 0.1*dist + r_off; r_out = 0.1 + 0.10*(dist+1)    - r_off; }
      else      { r_in = 0.1 + 0.1*dist + r_off; r_out = 0.1 + 0.10*dist + 0.05 - r_off; }
      sprintf(   txt, "%d", num);
      r_av   = 0.5*(r_out + r_in) ;
      ttext[num].SetText( cx + r_av*cos(phi_av*D2R) , cy + r_av*sin(phi_av*D2R) , txt );
      ttext[num].SetTextSize(0.014);

 //     cout << " " << num << " " << color[num] << " " << fired[num] << "\n";
  //    if( fired[num] ){
  //      cr[num]->SetFillStyle(3004);
 //       cr[num]->SetFillColor( color[num] );
 //     }

      cr[num]->DrawCrown(cx, cy, r_in, r_out, phi_min, phi_max);
      ttext[num].Draw();
      num++;
    }

    r_av = 0.475;
    sprintf(   txt, "%d", sector+1);
    sec_ttext[sector].SetText( cx + r_av*cos(phi_av*D2R) , cy + r_av*sin(phi_av*D2R) , txt );
    sec_ttext[sector].SetTextSize(0.03);
    sec_ttext[sector].Draw();
  }
  cr[65]->DrawCrown(cx, cy, 0.00           , 0.05 - r_off/2. );
  cr[66]->DrawCrown(cx, cy, 0.05 + r_off/2., 0.10 - r_off/2. );

  ttext[65].SetText( 0.51 , 0.51 , "65" );
  ttext[65].SetTextSize(0.014);
  ttext[66].SetText( 0.43 , 0.52 , "66" );
  ttext[66].SetTextSize(0.014);
  ttext[65].Draw();
  ttext[66].Draw();

  sprintf(   txt, "TIME %d", frame);
  ev_ttext.SetText(0.01, 0.95, txt);
  ev_ttext.Draw();

  cnvd->Print(fig_name);
  cnvd->Close();
  }
};


void tpcEvent::DrawHist( int anod ){
  TCanvas* cnv = new TCanvas("cnvh","histo",900,450);
  TH1F* hF = new TH1F("hF",";ch;val",n_ch,0.5,n_ch+0.5);
  TH1F* hB = new TH1F("hB",";ch;val",n_ch,0.5,n_ch+0.5);
  hB->SetLineColor(2);

  double minY = 1000000;
  double maxY = 0;

  for(int ch=0;ch<n_ch;ch++){
    hF->SetBinContent(ch+1,FADC[anod][ch]);
    hF->SetBinError(ch+1,0);
    hB->SetBinContent(ch+1,baseLevel[anod]);
    hB->SetBinError(ch+1,0);
    if( FADC[anod][ch] < minY ) minY = FADC[anod][ch] ;
    if( FADC[anod][ch] > maxY ) maxY = FADC[anod][ch] ;

  }
  double step = (maxY - minY) / 10.;
  hF->GetYaxis()->SetRangeUser(minY-30, maxY+30);
  hF->Draw();
  hB->Draw("same");
  ll.SetX1( timeStamp[anod]                );
  ll.SetY1( baseLevel[anod]                );
  ll.SetX2( timeStamp[anod]                );
  ll.SetY2( FADC[anod][peakPosition[anod]] );
  ll.SetLineColor(2);
  ll.SetLineStyle(7);
  ll.Draw();
  char *txt = new char[256];
  sprintf(   txt, "EVENT  %d", ev);
  ev_ttext.SetText(100, maxY - 1.*step, txt);
  ev_ttext.Draw();
  sprintf(   txt, "ANOD   %d", anod);
  an_ttext.SetText(100, maxY - 2.*step, txt);
  an_ttext.Draw();
  sprintf(   txt, "TIME   %d", timeStamp[anod]);
  tm_ttext.SetText(100, maxY - 3.*step, txt);
  tm_ttext.Draw();
  sprintf(   txt, "PEAK   %d", peakPosition[anod]);
  pk_ttext.SetText(100, maxY - 4.*step, txt);
  pk_ttext.Draw();
  sprintf(   txt, "ENERGY %d", energyDeposit[anod]);
  en_ttext.SetText(100, maxY - 5.*step, txt);
  en_ttext.Draw();
  cnv->Print("anode.png");
  cnv->Close();
};

void tpcEvent::PrepareAnimation( int anod ){
  TString fig_name;
  for(int frame=0;frame<240;frame++){

  if(              frame<  10 ) fig_name.Form("temp_000%d.png",frame);
  if( frame>  9 && frame< 100 ) fig_name.Form("temp_00%d.png" ,frame);
  if( frame> 99 && frame<1000 ) fig_name.Form("temp_0%d.png"  ,frame);
  if( frame>999               ) fig_name.Form("temp_%d.png"   ,frame);


  TCanvas* cnv = new TCanvas("cnvh","histo",900,450);
  TH1F* hF = new TH1F("hF",";ch;val",n_ch,0.5,n_ch+0.5);
  TH1F* hB = new TH1F("hB",";ch;val",n_ch,0.5,n_ch+0.5);
  hB->SetLineColor(2);

  double minY = 1000000;
  double maxY = 0;

  for(int ch=0;ch<n_ch;ch++){
    hF->SetBinContent(ch+1,FADC[anod][ch]);
    hF->SetBinError(ch+1,0);
    hB->SetBinContent(ch+1,baseLevel[anod]);
    hB->SetBinError(ch+1,0);
    if( FADC[anod][ch] < minY ) minY = FADC[anod][ch] ;
    if( FADC[anod][ch] > maxY ) maxY = FADC[anod][ch] ;

  }
  double step = (maxY - minY) / 10.;
  hF->GetYaxis()->SetRangeUser(minY-30, maxY+30);
  hF->Draw();
  hB->Draw("same");
  ll.SetX1( frame                          );
  ll.SetY1( baseLevel[anod]                );
  ll.SetX2( frame                          );
  ll.SetY2( FADC[anod][peakPosition[anod]] );
  ll.SetLineColor(2);
  ll.SetLineStyle(7);
  ll.Draw();
  char *txt = new char[256];
  sprintf(   txt, "EVENT  %d", ev);
  ev_ttext.SetText(100, maxY - 1.*step, txt);
  ev_ttext.Draw();
  sprintf(   txt, "ANOD   %d", anod);
  an_ttext.SetText(100, maxY - 2.*step, txt);
  an_ttext.Draw();
  sprintf(   txt, "TIME   %d", frame );
  tm_ttext.SetText(100, maxY - 3.*step, txt);
  tm_ttext.Draw();
  sprintf(   txt, "PEAK   %d", peakPosition[anod]);
  pk_ttext.SetText(100, maxY - 4.*step, txt);
  pk_ttext.Draw();
  sprintf(   txt, "ENERGY %d", energyDeposit[anod]);
  en_ttext.SetText(100, maxY - 5.*step, txt);
  en_ttext.Draw();
  cnv->Print(fig_name);
  cnv->Close();
  }

};


void tpcEvent::DrawTwo( int anod1, int anod2 ){
  TCanvas* cnv = new TCanvas("cnvh","histo",900,450);
  TH1F* hF1 = new TH1F("hF1",";ch;val",n_ch,0.5,n_ch+0.5);
  TH1F* hB1 = new TH1F("hB1",";ch;val",n_ch,0.5,n_ch+0.5);
  TH1F* hF2 = new TH1F("hF2",";ch;val",n_ch,0.5,n_ch+0.5);
  TH1F* hB2 = new TH1F("hB2",";ch;val",n_ch,0.5,n_ch+0.5);
  hF1->SetLineColor(1);
  hB1->SetLineColor(1);
  hF2->SetLineColor(2);
  hB2->SetLineColor(2);
  hB2->SetLineStyle(7);
  hB1->SetLineStyle(7);
  TLine l1,l2;
  for(int ch=0;ch<n_ch;ch++){
    hF1->SetBinContent(ch+1,FADC[anod1][ch]);
    hF1->SetBinError(ch+1,0);
    hB1->SetBinContent(ch+1,baseLevel[anod1]);
    hB1->SetBinError(ch+1,0);
    hF2->SetBinContent(ch+1,FADC[anod2][ch]);
    hF2->SetBinError(ch+1,0);
    hB2->SetBinContent(ch+1,baseLevel[anod2]);
    hB2->SetBinError(ch+1,0);
  }
  hF1->Draw();
  hB1->Draw("same");
  hF2->Draw("same");
  hB2->Draw("same");
  cnv->Print("two_anodes.png");
  cnv->Close();
};

void tpcEvent::DrawSector( int sector ){
  TCanvas* cnvd = new TCanvas("cnvh","display",900,450);

  TH1F* hF1 = new TH1F("hF1",";ch;val",n_ch,0.5,n_ch+0.5);
  TH1F* hF2 = new TH1F("hF2",";ch;val",n_ch,0.5,n_ch+0.5);
  TH1F* hF3 = new TH1F("hF3",";ch;val",n_ch,0.5,n_ch+0.5);
  TH1F* hF4 = new TH1F("hF4",";ch;val",n_ch,0.5,n_ch+0.5);
  TH1F* hF5 = new TH1F("hF5",";ch;val",n_ch,0.5,n_ch+0.5);
  TH1F* hF6 = new TH1F("hF6",";ch;val",n_ch,0.5,n_ch+0.5);
  hF1->SetLineColor(1);
  hF2->SetLineColor(2);
  hF3->SetLineColor(209);
  hF4->SetLineColor(4);
  hF5->SetLineColor(6);
  hF6->SetLineColor(46);

  double minY = 1000000;
  double maxY = 0;

  for(int ch=0;ch<n_ch;ch++){
    hF1->SetBinContent(ch+1,FADC[65][ch]);
    hF2->SetBinContent(ch+1,FADC[66][ch]);
    hF3->SetBinContent(ch+1,FADC[ (sector-1)*4 + 1 ][ch]);
    hF4->SetBinContent(ch+1,FADC[ (sector-1)*4 + 2 ][ch]);
    hF5->SetBinContent(ch+1,FADC[ (sector-1)*4 + 3 ][ch]);
    hF6->SetBinContent(ch+1,FADC[ (sector-1)*4 + 4 ][ch]);

    if( FADC[ 65               ][ch] < minY ) minY = FADC[ 65               ][ch] ;
    if( FADC[ 66               ][ch] < minY ) minY = FADC[ 66               ][ch] ;
    if( FADC[ (sector-1)*4 + 1 ][ch] < minY ) minY = FADC[ (sector-1)*4 + 1 ][ch] ;
    if( FADC[ (sector-1)*4 + 2 ][ch] < minY ) minY = FADC[ (sector-1)*4 + 2 ][ch] ;
    if( FADC[ (sector-1)*4 + 3 ][ch] < minY ) minY = FADC[ (sector-1)*4 + 3 ][ch] ;
    if( FADC[ (sector-1)*4 + 4 ][ch] < minY ) minY = FADC[ (sector-1)*4 + 4 ][ch] ;

    if( FADC[ 65               ][ch] > maxY ) maxY = FADC[ 65               ][ch] ;
    if( FADC[ 66               ][ch] > maxY ) maxY = FADC[ 66               ][ch] ;
    if( FADC[ (sector-1)*4 + 1 ][ch] > maxY ) maxY = FADC[ (sector-1)*4 + 1 ][ch] ;
    if( FADC[ (sector-1)*4 + 2 ][ch] > maxY ) maxY = FADC[ (sector-1)*4 + 2 ][ch] ;
    if( FADC[ (sector-1)*4 + 3 ][ch] > maxY ) maxY = FADC[ (sector-1)*4 + 3 ][ch] ;
    if( FADC[ (sector-1)*4 + 4 ][ch] > maxY ) maxY = FADC[ (sector-1)*4 + 4 ][ch] ;

    hF1->SetBinError(ch+1,0);
    hF2->SetBinError(ch+1,0);
    hF3->SetBinError(ch+1,0);
    hF4->SetBinError(ch+1,0);
    hF5->SetBinError(ch+1,0);
    hF6->SetBinError(ch+1,0);
  }

  hF1->GetYaxis()->SetRangeUser(minY-30, maxY+30);
  hF2->GetYaxis()->SetRangeUser(minY-30, maxY+30);
  hF3->GetYaxis()->SetRangeUser(minY-30, maxY+30);
  hF4->GetYaxis()->SetRangeUser(minY-30, maxY+30);
  hF5->GetYaxis()->SetRangeUser(minY-30, maxY+30);
  hF6->GetYaxis()->SetRangeUser(minY-30, maxY+30);

  hF1->Draw();
  hF2->Draw("same");
  hF3->Draw("same");
  hF4->Draw("same");
  hF5->Draw("same");
  hF6->Draw("same");
  cnvd->Print("sector.png");
  cnvd->Close();
};

void tpcEvent::DrawAll( ){
  TCanvas* cnv_all = new TCanvas("cnv_all","all",1200,900);
  cnv_all->Divide(8,9);
  TH1F* hF[67];
  TH1F* hB[67];
  for(int anod=1;anod<67;anod++){
    cnv_all->cd(anod);
    char *hff = new char[256];
    char *hbb = new char[256];
    sprintf( hff, "hF_%d", anod);
    sprintf( hbb, "hB_%d", anod);
    hF[anod] = new TH1F( hff,";ch;val",n_ch,0.5,n_ch+0.5);
    hB[anod] = new TH1F( hbb,";ch;val",n_ch,0.5,n_ch+0.5);
    hB[anod]->SetLineColor(2);
    //hF[anod]->GetYaxis()->SetRangeUser(8000,12000);
    for(int ch=0;ch<n_ch;ch++){
      hF[anod]->SetBinContent(ch+1,FADC[anod][ch]);
      hF[anod]->SetBinError(ch+1,0);
      hB[anod]->SetBinContent(ch+1,threshold*baseLevel[anod]);
      hB[anod]->SetBinError(ch+1,0);
    }
    hF[anod]->Draw();
    hB[anod]->Draw("same");
  }
};
void tpcEvent::CheckAll( ){

  //cout << "Window    : " << window    << endl;
  //cout << "Threshold : " << threshold << endl;
  for(int anod=1;anod<67;anod++){
    baseLevel[anod] = 0;
    fired[anod] = false;
  }
  for(int anod=1;anod<67;anod++){
    double sum;
    for(int ch=0;ch<n_ch;ch++){
      baseLevel[anod] = baseLevel[anod] + FADC[anod][ch]/ n_ch;
    }
    sum = 0.;
    for(int ch = 0; ch < window; ch++){
      sum = sum + FADC[anod][ch] / window;
    }
    if( sum > threshold*baseLevel[anod] )
      fired[anod] = true;
    for(int ch=0; ch < n_ch - window; ch++ ){
      sum = sum - FADC[anod][ch] / window;
      sum = sum + FADC[anod][ch+window] / window;
//      if(anod==16 && ch>1300 && ch<1500)
//         cout << " Channel: " << ch << "\t" << baseLevel[anod]
//              << "\t" << sum << "\t" << threshold*baseLevel[anod] << "\n";

      if( sum > threshold*baseLevel[anod] )
        fired[anod] = true;
    }
    if(fired[anod]){
      CorrectBaseLevel(anod);
      EvalBaseDev(anod);
      FindTime(anod);
      FindEnergyDeposit(anod);
      //cout << " Anod : " << anod << "\t   time: " << timeStamp[anod]
      //                           << "\t   peak: " << peakPosition[anod]
      //                           << "\t energy: " << energyDeposit[anod] << "\n";
    }
  }
};

void tpcEvent::CorrectBaseLevel( int anod ){
  if( fired[anod] ){
    int peak_pos = 0; double peak_level = FADC[anod][peak_pos];
    for(int ch=1;ch<TRACELENGTH;ch++){
      if( FADC[anod][ch] > peak_level ){
        peak_pos   = ch;
        peak_level = FADC[anod][ch];
      }
    }
    peakPosition[anod] = peak_pos;

    double sum = 0 ; double nch = 0 ;
    int veto_width = 150;

    // if peak is close left edge
    if(peak_pos<1500){
      for(int ch=0;ch<TRACELENGTH;ch++){
        if( ch < peak_pos - veto_width || ch > peak_pos + (2.*veto_width) ){
          sum = sum + FADC[anod][ch];
          nch = nch + 1.;
        }
      }
    }
    else{ // if it's as usual take only left part
      for(int ch=0;ch<TRACELENGTH;ch++){
        if( ch < peak_pos - veto_width ){
          sum = sum + FADC[anod][ch];
          nch = nch + 1.;
        }
      }
    }

    baseLevel[anod] = sum / nch ;

  }
};

void tpcEvent::EvalBaseDev( int anod ){
  if( fired[anod] ){

    double max_dev = 0 ;
    int veto_width = 150;

    // if peak is close left edge
    if(peakPosition[anod]<1500){
      for(int ch=0;ch<TRACELENGTH;ch++){
        if( ch < peakPosition[anod] - veto_width || ch > peakPosition[anod] + (2.*veto_width) ){
          if( sqrt(pow(baseLevel[anod] - FADC[anod][ch],2))>max_dev ){
            max_dev = sqrt(pow(baseLevel[anod] - FADC[anod][ch],2));
          }
        }
      }
    }
    else{ // if it's as usual take only left part
      for(int ch=0;ch<TRACELENGTH;ch++){
        if( ch < peakPosition[anod] - veto_width ){
          if( sqrt(pow(baseLevel[anod] - FADC[anod][ch],2))>max_dev ){
            max_dev = sqrt(pow(baseLevel[anod] - FADC[anod][ch],2));
          }
        }
      }
    }

    baseDev[anod] = max_dev;

  }
};



void tpcEvent::FindTime( int anod ){
  const int gate  = 80;
  int shift = 40;
  double fact = 0.5;
  if( fired[anod] ){
    if( peakPosition[anod]>gate ){
      double part[gate];
      for(int ii =0; ii<gate; ii++){
        part[ii] = ( FADC[anod][ peakPosition[anod]-gate + ii ] - baseLevel[anod] ) ;
        part[ii] = part[ii] - fact*( FADC[anod][ peakPosition[anod]-gate+ii+shift ] - baseLevel[anod] ) ;
        if( ii>0 ){
          if( part[ii-1]<0  && part[ii]>=0 )
            timeStamp[anod] = peakPosition[anod]-gate + ii ;
        }
      }
    }
  }
};
void tpcEvent::FindEnergyDeposit( int anod ){
  energyDeposit[anod] = 0.;
  if( fired[anod] ){
    bool go_left  = true;
    bool go_right = true;
    int pos = peakPosition[anod];
    while(go_left){
      pos--;
      if( FADC[anod][pos] > baseLevel[anod] )
        energyDeposit[anod] += FADC[anod][pos];
      else
        go_left = false;
    }
    pos = peakPosition[anod]-1;
    while(go_right){
      pos++;
      if( FADC[anod][pos] > baseLevel[anod] )
        energyDeposit[anod] += FADC[anod][pos];
      else
        go_right = false;
    }
  }
};

