#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "TFile.h"
#include <TGraph.h>
#include "TH1.h"

//
// class declaration
//

class DemoAnalyzer6 : public edm::EDAnalyzer {
   public:
      explicit DemoAnalyzer6(const edm::ParameterSet&);
      ~DemoAnalyzer6();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:

	// configuration parameters

	edm::EDGetTokenT<std::vector<reco::Track>> tracksToken_;
	edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_;
         //edm::EDGetTokenT<edm::View<reco::Vertex>> primary_verticesToken_;  
         edm::EDGetTokenT<std::vector<reco::Vertex> > vertexToken_;


      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
  
// histograms for tracks in 0.1,0.2,0.5 windows
  TH1D* histo1;
  TH1D* histo2;
  TH1D* histo5;
 
// histogram for the number of tracks in the events
  TH1D* histoNo;
  
// histograms for the number of tracks in the event after 1st and 2nd selections
  TH1D* histoSl1;
  TH1D* histoSl2;

// histograms for the track impact parameter with respect to (0,0,0)
  TH1D* histod0;
  TH1D* histodz;

// histograms for the track impact parameter with respect to the vertex with maximum tracks
  TH1D* histod0pv; 
  TH1D* histodzpv; 

// histograms for the track impact parameter with respect to beamspot
  TH1D* histod0bs; 
  TH1D* histodzbs; 

// histogram for the number of valid hits of the tracks
TH1D* histohn;


  //TGraph* g;
  //TFile* f;

int i =0; // i for the event number
int j =0; // j for the number of tracks after the 1st selection
int p=0; // p for the number of tracks after the 2nd selection

int u=0; // u for the number of vertices in the event

double dxxy=0.0; double dzz=0.0; // dxxy and dzz are the position of the vertex with the maximum tracks


int k1 =0; // used as a comparsion to get the maxium number of tracks in 0.1 window 
int y2 =0; // used as a comparsion to get the maxium number of tracks in 0.2 window
int m5 =0; // used as a comparsion to get the maxium number of tracks in 0.5 window

int l =0; // the maximum number of tracks in 0.1 window
int z =0; // the maximum number of tracks in 0.2 window
int s =0; // the maximum number of tracks in 0.5 window

//TGraph* graphno;

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
DemoAnalyzer6::DemoAnalyzer6(const edm::ParameterSet& iConfig)
{
	tracksToken_ = mayConsume<std::vector<reco::Track> >(iConfig.getParameter<edm::InputTag> ("trck"));
	vertexToken_ = mayConsume<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("vrtx"));
	beamSpotToken_ = consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("bmspt"));
}


DemoAnalyzer6::~DemoAnalyzer6()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
DemoAnalyzer6::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
i = i+1;

   using namespace edm;
   using namespace std;
//*******************************************************************************************************************************//
cout <<"************************** begining of the event no( " << i << " )*****************************" << endl;

// handles
Handle<vector<reco::Track>> tracks;
Handle<vector<reco::Vertex>> primary_vertices;
Handle<reco::BeamSpot> beamSpot;
//********************************************************************************************************************************//
// vector defined to get the vertex with maxium number of tracks
vector <const reco::Vertex*> vertexvector;

//********************************************************************************************************************************//

// read beamspot 
	iEvent.getByToken(beamSpotToken_, beamSpot);
// read tracks
         iEvent.getByToken(tracksToken_, tracks);
// read primary vertices 
        iEvent.getByToken(vertexToken_, primary_vertices);

//*********************************************************************************************************************************//
// vx is the first vertex in primary_vertices
const reco::Vertex *vx = &(*primary_vertices ->begin());

//*********************************************************************************************************************************//

//for loop to arrange the vertices again and fill vertexvector 
 for (vector<reco::Vertex>::const_iterator vertices1 = primary_vertices->begin() ; vertices1 != primary_vertices->end(); vertices1 ++ )
{
if (vertices1 ->tracksSize()> vx ->tracksSize())
{
vertexvector.push_back(&*vertices1);
}
if (vertices1 ->tracksSize()== vx ->tracksSize())
{
vertexvector.push_back(&*vertices1);
}}

//***********************************************************************************************************************************//
// vxx is the vertex with the maxium track density
const reco::Vertex *vxx = *vertexvector.begin();

//**********************************************************************************************************************************//
// to print the number of vertices in vertexvector
cout <<"the number of vertices after selection " << vertexvector.size()<<endl;

//***********************************************************************************************************************************//

dxxy = vx->position().Rho();
dzz = vx ->position().Z();
cout <<"***********************************************************************"<<endl;
cout << "the 1st primary vertex has "<<vx ->tracksSize() << " number of tracks" <<endl;
cout <<"position of the 1st primary vertex in z "<<vx ->position().Z()<<" && in rho "<<vx ->position().Rho() <<endl; 
cout <<"***********************************************************************"<<endl;

//************************************************************************************************************************************//

  // print the number of primary vertex in this event to screen
cout <<"----------------------------------------------------------------------------------" <<endl;
cout << " ^_^ ^_^ the number of primary vertices in this event is = " << primary_vertices -> size () << endl;
   
//*************************************************************************************************************************************//

   // to print the number of tracks in the primary vertices. also print the vertex position in z and rho
   for (vector<reco::Vertex>::const_iterator vertices = primary_vertices->begin() ; vertices != primary_vertices->end(); vertices ++ )
   {
 if(!vertices->isFake()) {
u = u+1;
cout <<"*************************** "<<"the number of this vertex is "<< u <<" ********************************"<<endl;
   cout << "  => Vertex has nTracks: " << vertices ->nTracks() << " and tracksSize:" << vertices ->tracksSize() << endl;
   cout << "  => position of vertex in z = " << vertices ->position().z() << " && in rho = " << vertices -> position().Rho() << endl;
cout <<"----------------------------------------------------------------------------------" <<endl; 
//cout << "    ===> position of the vertex with maximum no of tracks in z = " << dzz << " && in rho = " << dxxy << endl;

if (vertices ->tracksSize() >= vx ->tracksSize())
{
dxxy = vertices ->position().Rho();
dzz = vertices ->position().Z();
}
}
}

//**********************************************************************************************************************************//

cout <<"$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"<<endl;
 cout << "    ===> position of the vertex with maximum no of tracks in z = " << dzz << " && in rho = " << dxxy << endl;
cout <<"      ===> position of the vertex with maxium no of tracks after selction in Z = " <<vxx->position().Z() <<" && in rho= "<<vxx->position().Rho() <<endl; 
cout <<"$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"<<endl;

//***********************************************************************************************************************************//

for (reco::TrackCollection::const_iterator track2 = tracks->begin() ; track2 != tracks->end(); track2 ++ )
   {
if (track2 -> quality(reco::Track::qualityByName("highPurity"))                  // 1- high quality     
     && fabs(track2 -> eta()) < 2.4                                              // 2- |η| < 2.4
     && track2 ->pt() >0.4  && track2 ->pt() <5.0                                // 3- pt >0.4 GeV/c.
     && track2 ->ptError()/track2 ->pt() < 0.1                                   // 4- σ(pT)/pT < 10 %                                
     && fabs(track2 ->dz()/ track2 ->dzError()) < 3.0                            // 5- dz/ σ(dz) < 3
      && track2 ->numberOfValidHits() >5                                         // 5- number of valid hits >5
     && fabs(track2 ->dxy()/track2 ->dxyError()) < 3.0    )                      // 6- dxy/σ(dxy)
j = j+1;
}

//**********************************************************************************************************************************//

for (reco::TrackCollection::const_iterator track3 = tracks->begin() ; track3 != tracks->end(); track3 ++ )
   {

if (track3 -> quality(reco::Track::qualityByName("highPurity"))                  // 1- high quality     
     && fabs(track3 -> eta()) < 2.4                                              // 2- |η| < 2.4
     &&track3 ->pt() >0.4   &&track3 ->pt() <5.0                                 // 3- pt >0.4 GeV/c. && pt<5.0 GeV/c
     && track3 ->ptError()/track3 ->pt() < 0.1                                   // 4- σ(pT)/pT < 10                                
     && track3 ->numberOfValidHits() >5                                          // 5- number of valid hits >5
         // 5- dz/ σ(dz) < 3
     && fabs(track3 ->dz(vxx ->position())/sqrt(pow(track3->dzError(),2)+pow(vxx->zError(),2))) < 3.0        
        // 6- dxy/ σ(dxy) < 3
     && fabs(track3 ->dxy(vxx ->position())/sqrt(pow(track3->dxyError(),2)+pow(vxx->xError(),2)+pow(vxx->yError(),2)) < 3.0)) 
{
//fill histograms after tracks selection
histod0 ->Fill(track3 ->dxy());
histodz ->Fill(track3 ->dz());

histod0bs ->Fill(track3 ->dxy(beamSpot->position()));
histodzbs ->Fill(track3 ->dz(beamSpot->position()));

histod0pv ->Fill(track3 ->dxy(vxx ->position()));
histodzpv ->Fill(track3 ->dz(vxx ->position()));

histohn ->Fill(track3 -> numberOfValidHits());

p = p+1;
}}


//************************************************************************************************************************************//
   
   for (reco::TrackCollection::const_iterator track = tracks->begin() ; track != tracks->end(); track ++ )
   {
       // tracks selection for track
if (track -> quality(reco::Track::qualityByName("highPurity"))                  // 1- high quality     
     && fabs(track -> eta()) < 2.4                                              // 2- |η| < 2.4
     && track ->pt() >0.4 &&track ->pt() <5.0                                   // 3-pt >0.4 GeV/c && pt < 5.0 GeV/c. 
     && track ->ptError()/track ->pt() < 0.1                                    // 4- σ(pT)/pT < 10                  
      && track ->numberOfValidHits() >5                                          // 5- number of valid hits >5
     // 5- dz/ σ(dz) < 3              
     && fabs(track ->dz(vxx ->position())/sqrt(pow(track->dzError(),2)+pow(vxx->zError(),2))) < 3.0 
     //6-dxy/ σ(dxy)<3
     && fabs(track ->dxy(vxx ->position())/sqrt(pow(track->dxyError(),2)+pow(vxx->xError(),2)+pow(vxx->yError(),2)) < 3.0)) 


{
   for (reco::TrackCollection::const_iterator track1 = tracks->begin() ; track1 != tracks->end(); track1 ++ )
   {
      // tracks selection for track1

if (track1 -> quality(reco::Track::qualityByName("highPurity"))                 // 1- high quality     
     && fabs(track1 -> eta()) < 2.4                                             // 2- |η| < 2.4
     &&track1 ->pt() >0.4 &&track1 ->pt() <5.0                                  // 3- pt  0.1 GeV/c && pt <5.0 GeV/c.
     && track1 ->ptError()/track1 ->pt() < 0.1                                  // 4- σ(pT)/pT < 10                                
 && track1 ->numberOfValidHits() >5                                          // 5- number of valid hits >5     
      // 5- dz/ σ(dz) < 3
     && fabs(track1 ->dz(vxx ->position())/sqrt(pow(track1->dzError(),2)+pow(vxx->zError(),2))) < 3.0         
           //6-dxy/σ(dxy)
     && fabs(track1 ->dxy(vxx ->position())/sqrt(pow(track1->dxyError(),2)+pow(vxx->xError(),2)+pow(vxx->yError(),2)) < 3.0)) 
{                                                                           
       if (abs(track ->eta() - track1 ->eta()) <= 0.1)
     { k1 =k1+1;}
  
       if (abs(track ->eta() - track1 ->eta()) <= 0.2)
     { y2 =y2+1;}   
      
       if (abs(track ->eta() - track1 ->eta()) <= 0.5)
     { m5 =m5+1;}   
 }}  
   }
//*************************************************************************************************************************//
   // k is for 0.1 for loop and l is the maxmium for 0.1
      if(k1 > l)
               { l =k1;}
   // z is for 0.2 for loop and z is the maxmium for 0.2
       if(y2 > z)
               { z =y2;} 
   // m is for 0.1 for loop and s is the maxmium for 0.5
       if(m5 > s)
               { s =m5;}
   // after fishing for loop and give numbers for l, z and s make k, m and y zero again for the second loop;
       k1 =0; y2 =0; m5 =0; }//}
   
//****************************************************************************************************************************// 
//fill the histogram with the event tracks and print the event tracks number
histoNo ->Fill(tracks ->size());
cout <<"----------------------------------------------------------------------------------" <<endl;
cout <<"*** the number of tracks in this event = " << tracks ->size()<<endl;
//**************************************************************************************************************************//
 //  g ->SetPoint(1,i,tracks ->size());
//*************************************************************************************************************************//

//fill the histogram histosl1 the tracks passes the first selections
histoSl1 ->Fill(j);
cout <<"*** the numer of tracks 1st after selection in this event = " <<j <<endl;
j=0;
//fill the histogram histosl2 the tracks passes the second selections
histoSl2 ->Fill(p);
cout <<"*** the numer of tracks  2nd after selection in this event = " <<p <<endl;
p=0;

//***************************************************************************************************************************//
   histo1 ->Fill(l);
cout <<"======================================================================================"<<endl;
   cout <<"         the maximum no of tracks with respect to 0.1 window  = "<< l << endl;
l =0;
   histo2 ->Fill(z);
   cout <<"         the maximum no of tracks with respect to 0.2 window  = "<< z << endl;
z=0;
   histo5 ->Fill(s);
   cout <<"         the maximum no of tracks with respect to 0.5 window  = "<< s << endl;
s=0;  u =0;dxxy =0.0 ; dzz =0.0;
//*******************************************************************************************************************************//

cout <<"======================================================================================"<<endl;
cout <<"********************************** ending of the event *******************************" <<endl;

//********************************************************************************************************************************//
#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif

}

// ------------ method called once each job just before starting event loop  ------------
void 
DemoAnalyzer6::beginJob()
{



  edm::Service<TFileService> fs;

//histograms
  
//******************************************************************************************************************************//

// histograms for number of tracks in 0.1,0.2,0.5 windows
   histo1 = fs->make<TH1D>("eta(0.1)" , "eta(0.1)" , 100 , 0 , 300 );
   histo2 = fs->make<TH1D>("eta(0.2)" , "eta(0.2)" , 100 , 0 , 600 );
   histo5 = fs->make<TH1D>("eta(0.5)" , "eta(0.5)" , 100 , 0 , 1000 );   

//******************************************************************************************************************************//

// histogram for the number of tracks in the event 
   histoNo = fs->make<TH1D>("tracksno." , "tracksno." , 200 , 0 , 16000 );

//******************************************************************************************************************************//

// histogram for the number of tracks in the event after the 1st and 2nd selsection  
   histoSl1 = fs->make<TH1D>("trackssl1." , "tracksSl1." , 100 , 0 , 1000 );
   histoSl2 = fs->make<TH1D>("trackssl2." , "tracksSl2." , 200 , 0 , 4000 );

//******************************************************************************************************************************//

// histograms for impact paramter with respect to (0,0,0)
   histod0 = fs->make<TH1D>("tracksd0" , "tracksd0" , 100 , -3 , 3 );
   histodz = fs->make<TH1D>("tracksdz" , "tracksdz" , 100 , -10 ,10 );

//******************************************************************************************************************************//

// histograms for impact paramter with respect to the vertex with maxium tracks
   histod0pv = fs->make<TH1D>("tracksd0pv" , "tracksd0pv" , 100 , -3 , 3 );
   histodzpv = fs->make<TH1D>("tracksdzpv" , "tracksdzpv" , 100 , -10, 10) ;
//******************************************************************************************************************************//

// histograms for impact paramter with respect to beamspot
   histod0bs = fs->make<TH1D>("tracksd0bs" , "tracksd0bs" , 100 , -3 , 3 );
   histodzbs = fs->make<TH1D>("tracksdzbs" , "tracksdzbs" , 100 , -10, 10) ;
//******************************************************************************************************************************//
  histohn = fs->make<TH1D>("trackshitsno" , "trackshitsno" , 200 , 0, 100) ;

  //graphno = fs->make<TH1D>("gtracksno." , "gtracksno." , 1000 , 0 , 15000 );
//f= new TFile ("test.root","RECREATE");
//g = new TGraph();
//int i =0;

}

// ------------ method called once each job just after ending the event loop  ------------
void 
DemoAnalyzer6::endJob() 
{
//g->Write();
//f->Close();

}

// ------------ method called when starting to processes a run  ------------
/*
void 
DemoAnalyzer6::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
DemoAnalyzer6::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
DemoAnalyzer6::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
DemoAnalyzer6::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
DemoAnalyzer6::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DemoAnalyzer6);
