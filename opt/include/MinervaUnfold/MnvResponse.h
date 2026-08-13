#ifndef MNV_MnvResponse_h
#define MNV_MnvResponse_h 1

//Root Stuff
#include "TObject.h"
#include "TNamed.h"
#include "TFile.h"
#include "TError.h"
#include "PlotUtils/MnvH2D.h"
#include "PlotUtils/MnvH3D.h"
#include "PlotUtils/HistWrapper.h"
#include "PlotUtils/Hist2DWrapper.h"
#include "PlotUtils/ChainWrapper.h"
#include "PlotUtils/AnaBinning.h"
#include "PlotUtils/HistogramUtils.h" //for NotPhysicsShift

#include "RooUnfold/RooUnfoldResponse.h"
#include <map>
#include <iterator>

using namespace PlotUtils;

namespace MinervaUnfold{

  class MnvResponse
  {
    public:
      //! Default constructor
      MnvResponse();
      //***********************************************/
      //***** one-dimensional constructors*************/
      //***********************************************/
      //! MCReco and MCTruth just used for binning, name and title
      MnvResponse(const char* name, const char* title, TH1D* MCReco, TH1D* MCTruth);

      MnvResponse(const char* name, const char* title, PlotUtils::axis_binning reco, PlotUtils::axis_binning truth);

      //Adding for CVUniverse...
      MnvResponse(const char* name,const char* title,PlotUtils::axis_binning reco,PlotUtils::axis_binning truth,std::map<const  std::string,  int> & error_bands);

      MnvResponse(const char* name, const char* title, const Int_t nbins_reco, const Double_t* bins_reco, const Int_t nbins_truth, const Double_t* bins_truth);

      // new HMS to allow one to handle particular error bands.  This already existed for 2D.

       MnvResponse(const char* name,const char* title,const Int_t nbins_reco, const Double_t* bins_reco, const Int_t nbins_truth, const Double_t* bins_truth ,std::map< std::string, int > &error_bands);
      //***********************************************/
      // *******two-dimensional constructors.**********/
      //***********************************************/
      //! Construct with MC Reco and True binning from TH2Ds
      //! MC Reco and MC True only used for shape
      MnvResponse(const char* name, const char* title, TH2D* MCReco, TH2D* MCTruth);

      //! Construct with MC Reco and True binning from axis_binning
      MnvResponse(const char* name, const char* title, PlotUtils::axis_binning x_reco, PlotUtils::axis_binning y_reco, PlotUtils::axis_binning x_truth, PlotUtils::axis_binning y_truth);


      //! From doubles*
      MnvResponse(const char* name, const char* title, const Int_t nbinsx_reco, const Double_t* xbins_reco, const Int_t nbinsy_reco, const Double_t* ybins_reco, const Int_t nbinsx_truth, const Double_t* xbins_truth, const Int_t nbinsy_truth, const Double_t* ybins_truth);

     //!For NSF

     //lets try one more time...
      template <typename T>
      MnvResponse(const char* name, const char* title, PlotUtils::axis_binning x_reco, PlotUtils::axis_binning y_reco, PlotUtils::axis_binning x_truth, PlotUtils::axis_binning y_truth,std::map< std::string, T*> & error_bands );


      MnvResponse(const char* name, const char* title, PlotUtils::axis_binning x_reco, PlotUtils::axis_binning y_reco, PlotUtils::axis_binning x_truth, PlotUtils::axis_binning y_truth,std::map<const std::string,  int > & error_bands );

      MnvResponse(const char* name, const char* title, const Int_t nbinsx_reco, const Double_t* xbins_reco, const Int_t nbinsy_reco, const Double_t* ybins_reco, const Int_t nbinsx_truth, const Double_t* xbins_truth, const Int_t nbinsy_truth, const Double_t* ybins_truth,std::map <const std::string, int >& error_bands);
      MnvResponse(const char* name, const char* title, std::vector<double> xbins_vector, std::vector<double> ybins_vector ,std::map <const std::string, int >& error_bands);


      //! From doubles*

    //  MnvResponse(const char* name, const char* title, const Int_t nbinsx_reco, const Double_t* xbins_reco, const Int_t nbinsy_reco, const Double_t* ybins_reco, const Int_t nbinsx_truth, const Double_t* xbins_truth, const Int_t nbinsy_truth, const Double_t* ybins_truth,std::map<const  std::string, int> & error_bands );

      //***********************************************/
      //******* three-dimensional constructors********/
      //***********************************************/
      //! MC Reco and MC True only used for shape
      MnvResponse(const char* name, const char* title, TH3D* MCReco, TH3D* MCTruth);

      //! Construct with MC Reco and True binning from axis_binning
      MnvResponse(const char* name, const char* title, PlotUtils::axis_binning x_reco, PlotUtils::axis_binning y_reco, PlotUtils::axis_binning z_reco, PlotUtils::axis_binning x_truth, PlotUtils::axis_binning y_truth, PlotUtils::axis_binning z_truth);

      //! From doubles*
      MnvResponse(const char* name, const char* title, const Int_t nbinsx_reco, const Double_t* xbins_reco, const Int_t nbinsy_reco, const Double_t* ybins_reco, const Int_t nbinsz_reco, const Double_t* zbins_reco, const Int_t nbinsx_truth, const Double_t* xbins_truth, const Int_t nbinsy_truth, const Double_t* ybins_truth, const Int_t nbinsz_truth, const Double_t* zbins_truth);

      //! Default destructor
      virtual ~MnvResponse();

      //! 1D name/binning Setup
      void Setup(const char* name, const char* title, TH1D* MCReco, TH1D* MCTruth);
      void Setup(const char* name, const char* title, PlotUtils::axis_binning reco, PlotUtils::axis_binning truth);
      void Setup(const char* name, const char* title, const Int_t nbins_reco, const Double_t* bins_reco, const Int_t nbins_truth, const Double_t* bins_truth);

  // new HMS  need to set it up
      void Setup(const char* name, const char* title, const Int_t nbins_reco, const Double_t* bins_reco, const Int_t nbins_truth, const Double_t* bins_truth,std::map< std::string,  int >  error_bands);


      ///////////////////
      //! 2D name/binning Setup
      ///////////////////

      void Setup(const char* name, const char* title, TH2D* MCReco, TH2D* MCTruth);
      void Setup(const char* name, const char* title, PlotUtils::axis_binning x_reco, PlotUtils::axis_binning y_reco, PlotUtils::axis_binning x_truth, PlotUtils::axis_binning y_truth);
      void Setup(const char* name, const char* title, const Int_t nbinsx_reco, const Double_t* xbins_reco, const Int_t nbinsy_reco, const Double_t* ybins_reco, const Int_t nbinsx_truth, const Double_t* xbins_truth, const Int_t nbinsy_truth, const Double_t* ybins_truth);
      void Setup(const char* name, const char* title, const Int_t nbinsx_reco, const Double_t* xbins_reco, const Int_t nbinsy_reco, const Double_t* ybins_reco, const Int_t nbinsx_truth, const Double_t* xbins_truth, const Int_t nbinsy_truth, const Double_t* ybins_truth , std::map<const std::string,  int > & error_bands);

      void Setup(const char* name, const char* title, PlotUtils::axis_binning reco, PlotUtils::axis_binning truth,std::map<const  std::string,  int > & error_bands);
      void Setup(const char* name, const char* title, PlotUtils::axis_binning x_reco, PlotUtils::axis_binning y_reco, PlotUtils::axis_binning x_truth, PlotUtils::axis_binning y_truth, std::map<const std::string,  int > & error_bands);

      //Implementations for the NSF....
      //! 1D name/binning Setup
    //   template <class T>
    //  void Setup(const char* name, const char* title, TH1D* MCReco, TH1D* MCTruth,std::map<const  std::string, const int > & error_bands);
     //  template <class T>


    //  template <class T>
    //  void Setup(const char* name, const char* title, const Int_t nbins_reco, const Double_t* bins_reco, const Int_t nbins_truth, const Double_t* bins_truth,std::map<const std::string, const int > & error_bands);

      //!2D Binning Setup for NSF
    //  template <class T>
    //  void Setup(const char* name, const char* title, TH2D* MCReco, TH2D* MCTruth,std::map<const std::string, const int > error_bands);
     // template <class T>
     // template <class T>
    //  void Setup(const char* name, const char* title, const Int_t nbinsx_reco, const Double_t* xbins_reco, const Int_t nbinsy_reco, const Double_t* ybins_reco, const Int_t nbinsx_truth, const Double_t* xbins_truth, const Int_t nbinsy_truth,const Double_t* ybins_truth,std::map<const  std::string, const int >  error_bands);

    //  template <typename T>
    //  void Setup(const char* name, const char* title, PlotUtils::axis_binning x_reco, PlotUtils::axis_binning y_reco, PlotUtils::axis_binning x_truth, PlotUtils::axis_binning y_truth,std::map<std::string, std::vector<T*> > & error_bands);



    ///////////////////
    //! 3D name/binning Setup
    ///////////////////

      void Setup(const char* name, const char* title, TH3D* MCReco, TH3D* MCTruth);
      void Setup(const char* name, const char* title, PlotUtils::axis_binning x_reco, PlotUtils::axis_binning y_reco, PlotUtils::axis_binning z_reco, PlotUtils::axis_binning x_truth, PlotUtils::axis_binning y_truth, PlotUtils::axis_binning z_truth);
      void Setup(const char* name, const char* title, const Int_t nbinsx_reco, const Double_t* xbins_reco, const Int_t nbinsy_reco, const Double_t* ybins_reco, const Int_t nbinsz_reco, const Double_t* zbins_reco, const Int_t nbinsx_truth, const Double_t* xbins_truth, const Int_t nbinsy_truth, const Double_t* ybins_truth, const Int_t nbinsz_truth, const Double_t* zbins_truth);

      Int_t Fill( double meas, double truth, double w = 1.0 );

     // HMS once again, allows you do access particular error bands  lik in the 2D case
      Int_t Fill( double meas_x,  double truth_x, const std::string name, const int univ, double w = 1.0 );

      //! Fill 2D histograms. Entries outside the histogram range are brought within range using \ref BringWithinHistLimits()
      Int_t Fill( double meas_x, double meas_y, double truth_x, double truth_y, double w = 1.0 );
     // template <class T>
      Int_t Fill(double meas_x,double meas_y,double truth_x,double truth_y,const std::string error_name,const int univ,double w = 1.0);
      Int_t Fill( double meas_x, double meas_y, double meas_z, double truth_x, double truth_y, double truth_z, double w = 1.0 );

      //! Check for presence of vertical error band
      bool HasVertErrorBand( const std::string& name ) const;

      std::vector<std::string> GetVertErrorBandNames() const;

      //!Handle Vertical Errors
      bool AddVertErrorBand( const std::string& name, const int nhists = -1 );

      // Fill 1D Vertical Errors
      bool FillVertErrorBand( const std::string& name, const double meas, const double truth, const std::vector<double>& weights, const double cvweight  = 1.0, double cvWeightFromMe = 1. );
      bool FillVertErrorBand( const std::string& name, const double meas, const double truth, const double* weights, const double cvweight  = 1.0, double cvWeightFromMe = 1. );
      //! Fill the weights of a MnvVertErrorBand's 2 universes with these 2 weights (must have 2)
      bool FillVertErrorBand( const std::string& name, const double meas, const double truth, const double weightDown, const double weightUp, const double cvweight  = 1.0, double cvWeightFromMe = 1. );
      // Fill 2D Vertical Errors
      //! Fill the weights of a MnvVertErrorBand's universes from a vector
      bool FillVertErrorBand( const std::string& name, const double x_meas, const double y_meas, const double x_truth, const double y_truth, const std::vector<double>& weights, const double cvweight  = 1.0, double cvWeightFromMe = 1. );
      bool FillVertErrorBand( const std::string& name, const double x_meas, const double y_meas, const double x_truth, const double y_truth, const double* weights, const double cvweight  = 1.0, double cvWeightFromMe = 1. );
      //! Fill the weights of an MnvVertErrorBand's 2 universes with these 2 weights (must have 2)
      bool FillVertErrorBand( const std::string& name, const double x_meas, const double y_meas, const double x_truth, const double y_truth, const double weightDown, const double weightUp, const double cvweight, double cvWeightFromMe );
      // Fill 3D Vertical Errors
      //! Fill the weights of a MnvVertErrorBand's universes from a vector
      bool FillVertErrorBand( const std::string& name, const double x_meas, const double y_meas, const double z_meas, const double x_truth, const double y_truth, const double z_truth, const std::vector<double>& weights, const double cvweight  = 1.0, double cvWeightFromMe = 1. );
      bool FillVertErrorBand( const std::string& name, const double x_meas, const double y_meas, const double z_meas, const double x_truth, const double y_truth, const double z_truth, const double* weights, const double cvweight  = 1.0, double cvWeightFromMe = 1. );
      //! Fill the weights of a MnvVertErrorBand's 2 universes with these 2 weights (must have 2)
      bool FillVertErrorBand( const std::string& name, const double x_meas, const double y_meas, const double z_meas, const double x_truth, const double y_truth, const double z_truth, const double weightDown, const double weightUp, const double cvweight, double cvWeightFromMe );

      //! Check for presence of lateral error band
      bool HasLatErrorBand( const std::string& name ) const;

      std::vector<std::string> GetLatErrorBandNames() const;

      //!Handle Lateral Errors
      bool AddLatErrorBand( const std::string& name, const int nhists = -1 );

      // Fill 1D Lateral Errors
      bool FillLatErrorBand( const std::string& name, const double meas, const double truth, const std::vector<double>& meas_shifts, const double cvweight  = 1.0, const bool fillcv = true, const double* weights = 0 );
      //! Fill the weights of an MnvLatErrorBand's universes from array
      bool FillLatErrorBand( const std::string& name, const double meas, const double truth, const double *meas_shifts, const double cvweight  = 1.0, const bool fillcv = true, const double* weights = 0 );
      //! Fill the weights of an MnvLatErrorBand's 2 universes with these 2 weights (must have 2)
      bool FillLatErrorBand( const std::string& name, const double meas, const double truth, const double meas_shiftDown, const double meas_shiftUp, const double cvweight = 1.0, const bool fillcv = true );

      // Fill 2D Lateral Errors
      bool FillLatErrorBand( const std::string& name, const double x_meas, const double y_meas, const double x_truth, const double y_truth, const std::vector<double>& meas_xshifts, const std::vector<double>& meas_yshifts, const double cvweight  = 1.0, const bool fillcv = true, const double* weights = 0 );
      //! Fill the weights of an MnvLatErrorBand's universes from array
      bool FillLatErrorBand( const std::string& name, const double x_meas, const double y_meas, const double x_truth, const double y_truth, const double *meas_xshifts, const double *meas_yshifts, const double cvweight  = 1.0, const bool fillcv = true, const double* weights = 0 );
      //! Fill the weights of an MnvLatErrorBand's 2 universes with these 2 weights (must have 2)
      bool FillLatErrorBand( const std::string& name, const double x_meas, const double y_meas, const double x_truth, const double y_truth, const double meas_xshiftDown, const double meas_xshiftUp, const double meas_yshiftDown, const double meas_yshiftUp, const double cvweight = 1.0, const bool fillcv = true );

      // Fill 3D Lateral Errors
      bool FillLatErrorBand( const std::string& name, const double x_meas, const double y_meas, const double z_meas, const double x_truth, const double y_truth, const double z_truth, const std::vector<double>& meas_xshifts, const std::vector<double>& meas_yshifts, const std::vector<double>& meas_zshifts, const double cvweight  = 1.0, const bool fillcv = true, const double* weights = 0 );
      //! Fill the weights of an MnvLatErrorBand's universes from array
      bool FillLatErrorBand( const std::string& name, const double x_meas, const double y_meas, const double z_meas, const double x_truth, const double y_truth, const double z_truth, const double *meas_xshifts, const double *meas_yshifts, const double *meas_zshifts, const double cvweight  = 1.0, const bool fillcv = true, const double* weights = 0 );
      //! Fill the weights of an MnvLatErrorBand's 2 universes with these 2 weights (must have 2)
      bool FillLatErrorBand( const std::string& name, const double x_meas, const double y_meas, const double z_meas, const double x_truth, const double y_truth, const double z_truth, const double meas_xshiftDown, const double meas_xshiftUp, const double meas_yshiftDown, const double meas_yshiftUp, const double meas_zshiftDown, const double meas_zshiftUp, const double cvweight = 1.0, const bool fillcv = true );

      bool GetMigrationObjects(MnvH2D* &h_migration, MnvH3D* &h_reco, MnvH3D* &h_truth);
      bool GetMigrationObjects(MnvH2D* &h_migration, MnvH2D* &h_reco, MnvH2D* &h_truth);
      bool GetMigrationObjects(MnvH2D* &h_migration, MnvH1D* &h_reco, MnvH1D* &h_truth);



      Int_t GetUnfoldingDimension();

      MnvH2D* GetMigrationMatrix() {return fMigration;}
      MnvH2D* GetTruth2D() {return fTruth2D;}
      MnvH2D* GetReco2D() {return fReco2D;}
      // MnvH2D* GetReco1D(){return fReco1D;}

    private:
      //3D case
      std::map<const std::string, int>::iterator itr;
      MnvH3D* fReco3D;
      MnvH3D* fTruth3D;
      //2D case
      MnvH2D* fReco2D;
      MnvH2D* fTruth2D;
     // template <typename T>
     // Hist2DWrapper<T> *WReco2D;
     // template <typename T>
     // Hist2DWrapper<T> *WTruth2D;

      //template <typename T>
      //Hist2DWrapper<T>*WMigration;

      //1D case
      MnvH1D* fReco1D;
      MnvH1D* fTruth1D;
      //HistWrapper Object

      //Migration Histogram is always 2D
      MnvH2D* fMigration;

      bool    fIs1D;
      bool    fIs2D;
      bool    fIs3D;

      /***********************/
      /***Helper Functions****/
      /***********************/
      bool checkIfMnvH1DNull(const MnvH1D* hist1D, const char* name = "MnvH1D Object");
      bool checkIfMnvH2DNull(const MnvH2D* hist2D, const char* name = "MnvH2D Object");
      bool checkIfMnvH3DNull(const MnvH3D* hist2D, const char* name = "MnvH3D Object");
      bool SetupResponse3D( const char* name, const char* title );
      bool SetupResponse2D( const char* name, const char* title );
      bool SetupResponse1D( const char* name, const char* title );
      void is1D();
      void is2D();
      void is3D();
      void Init();

  };
}

#endif
