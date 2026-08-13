#ifndef MNV_MnvUnfold_h 
#define MNV_MnvUnfold_h 1

//RooUnfold dependencies
#include "RooUnfold/RooUnfold.h" //include to get method enum
class RooUnfoldResponse;

//forwards
class TFile;
class TH1D;
class TH2D;
class TH3D;


namespace PlotUtils {
  class MnvH1D;
  class MnvH2D;
  class MnvH3D;
}

namespace MinervaUnfold{

  //! Methods to interface MnvH*Ds and the RooUnfold unfolding framework
  class MnvUnfold
  {
    private: 
    
      //! turn on better statistical error calculation
      bool m_useBetterStatErrorCalc;
      
  
    public:
    
      //! Default constructor
      MnvUnfold();
    
      //! Set better statistical error calculation
      inline void setUseBetterStatErrorCalc(bool useBetterStatErrorCalc) {
        m_useBetterStatErrorCalc = useBetterStatErrorCalc;
      }

      //! singleton gettor
      static MnvUnfold& Get();

      //===================================
      //===================================
      /** @name 1D Unfolding
        Functions to unfold in 1 dimension.
        @{*/

      /*! Unfold one histogram.
        @param[out] h_unfold The resulting unfolded histo.  Creates a *new* histogram if input is null.
        @param[out] cov covariance of unfolding error
        @param[in] h_migration Migration matrix of reco(X) vs true(Y) kinematics
        @param[in] h_mc_reco histogram binned in reco kinematics
        @param[in] h_mc_true histogram binned in true kinematics
        @param[in] h_data the folded histogram you will unfold
        @param[in] method See RooUnfold::Algorithm for options
        @param[in] regparam Regulates the unfolding algorithm
        @return true if all is OK
       */
      bool UnfoldHisto( 
          TH1D*& h_unfold, 
          TMatrixD &cov, 
          const TH2D* h_migration, 
          const TH1D* h_mc_reco, 
          const TH1D* h_mc_true, 
          const TH1D *h_data,
          RooUnfold::Algorithm method = RooUnfold::kBayes, 
          const Double_t regparam = -1e30 
          ) const;

      /*! Unfold one histogram and throw away the covariance matrix
        @param[out] h_unfold The resulting unfolded histo.   Creates a *new* histogram if input is null.
        @param[in] h_migration Migration matrix of reco(X) vs true(Y) kinematics
        @param[in] h_mc_reco histogram binned in reco kinematics
        @param[in] h_mc_true histogram binned in true kinematics
        @param[in] h_data the folded histogram you will unfold
        @param[in] method See RooUnfold::Algorithm for options
        @param[in] regparam Regulates the unfolding algorithm
        @return true if all is OK
       */
      bool UnfoldHisto( 
          TH1D*& h_unfold, 
          const TH2D* h_migration, 
          const TH1D* h_mc_reco, 
          const TH1D* h_mc_true, 
          const TH1D *h_data, 
          RooUnfold::Algorithm method = RooUnfold::kBayes, 
          const Double_t regparam = -1e30 
          ) const;

      /*! Unfold an MnvH1D and possibly its error bands
        @param[out] h_unfold The resulting unfolded histo.  Creates a *new* histogram if input is null.
        @param[in] h_migration Migration matrix of reco(X) vs true(Y) kinematics
        @param[in] h_mc_reco histogram binned in reco kinematics
        @param[in] h_mc_true histogram binned in true kinematics
        @param[in] h_data the folded histogram you will unfold
        @param[in] method See RooUnfold::Algorithm for options
        @param[in] regparam Regulates the unfolding algorithm
        @return true if all is OK
       */
      bool UnfoldHisto( 
          PlotUtils::MnvH1D*& h_unfold, 
          const TH2D* h_migration, 
          const TH1D* h_mc_reco, 
          const TH1D* h_mc_true,
          const  PlotUtils::MnvH1D *h_data, 
          RooUnfold::Algorithm method = RooUnfold::kBayes, 
          const Double_t regparam = -1e30, 
          bool addSystematics = true 
          ) const;

      /*! Unfold aTH1D
        @param[out] h_unfold The resulting unfolded histo.  Creates a *new* histogram if input is null.
        @param[in] h_migration Migration matrix of reco(X) vs true(Y) kinematics
        @param[in] h_data the folded histogram you will unfold
        @param[in] method See RooUnfold::Algorithm for options
        @param[in] regparam Regulates the unfolding algorithm
        @return true if all is OK
       */
      bool UnfoldHisto( 
          TH1D*& h_unfold,
          const TH2D* h_migration,
          const TH1D *h_data,
          RooUnfold::Algorithm method = RooUnfold::kBayes,
          const Double_t regparam = -1e30
          ) const;


      /*! Unfold an MnvH1D and possibly its error bands
        @param[out] h_unfold The resulting unfolded histo.  Creates a *new* histogram if input is null.
        @param[in] h_migration Migration matrix of reco(X) vs true(Y) kinematics
        @param[in] h_data the folded histogram you will unfold
        @param[in] method See RooUnfold::Algorithm for options
        @param[in] regparam Regulates the unfolding algorithm
        @param[in] addSystematics Unfold the error bands to create systematics
        @return true if all is OK
       */
      bool UnfoldHisto( 
          PlotUtils::MnvH1D*& h_unfold, 
          const TH2D* h_migration, 
          const PlotUtils::MnvH1D *h_data, 
          RooUnfold::Algorithm method = RooUnfold::kBayes, 
          const Double_t regparam = -1e30, 
          bool addSystematics = true 
          ) const;

      /*! Unfold an MnvH1D and possibly its error bands using a specially formatted file (formatted how?)
        @param[out] h_unfold The resulting unfolded histo.  Creates a *new* histogram if input is null.
        @param[in] f_MigrationFile A file that holds migration matrices of reco(X) vs true(Y) kinematics
        @param[in] h_nobck the folded histogram you will unfold
        @param[in] method See RooUnfold::Algorithm for options
        @param[in] regparam Regulates the unfolding algorithm
        @param[in] name_migration Name of the migration matrix
        @param[in] addSystematics Unfold the error bands to create systematics
        @param[in] useSysVariatedMigrations Use a special migration matrix for each universe for each error band
        @return true if all is OK
       */
      bool UnfoldHisto( 
          PlotUtils::MnvH1D*& h_unfold, 
          TFile *f_MigrationFile, 
          const PlotUtils::MnvH1D* h_nobck,
          RooUnfold::Algorithm method = RooUnfold::kBayes, 
          const Double_t regparam = -1e30, 
          const std::string name_migration = "h_enu_migration_qe", 
          bool addSystematics = true, 
          bool useSysVariatedMigrations = false 
          ) const;

      /*! Unfold one histogram using MnvH2D as migration.  Creates a *new* histogram if input is null.
        @param[out] h_unfold The resulting unfolded histo
        @param[out] cov covariance of unfolding error
        @param[in] h_migration Migration matrix of reco(X) vs true(Y) kinematics
        @param[in] h_data the folded histogram you will unfold
        @param[in] method See RooUnfold::Algorithm for options
        @param[in] regparam Regulates the unfolding algorithm
        @param[in] addSystematics Unfold the error bands to create systematics
        @param[in] useSysVariatedMigrations Use a special migration matrix for each universe for each error ba
        @return true if all is OK
       */
      bool UnfoldHisto( 
          PlotUtils::MnvH1D*& h_unfold, 
          TMatrixD& covmx,
          const PlotUtils::MnvH2D *h_migration, 
          const PlotUtils::MnvH1D *h_data, 
          RooUnfold::Algorithm method = RooUnfold::kBayes, 
          const Double_t regparam = -1e30, 
          bool addSystematics = true, 
          bool useSysVariatedMigrations = true
          );
      //@}
      //=========== END OF 1D UNFOLDING =============//
      //=============================================//


      //===================================
      //
      /** @name 1D Folding
        Functions to fold/smear a generated distribution.
        @*/
      /*! Unfold one 1D histogram
        @param[out] h_folded The result of folding the true distribution with detector smearing
        @param[in] h_migration Migration matrix of reco(X) vs true(Y) kinematics
        @param[in] h_generated The generated distribution binned in true kinematics
        */
      bool FoldHisto(
          TH1D* &h_folded,
          const TH2D* h_migration,
          const TH1D* h_generated
          ) const;

      /*! Unfold one MnvH1D histogram and its error bands
        @param[out] h_folded The result of folding the true distribution with detector smearing
        @param[in] h_migration Migration matrix of reco(X) vs true(Y) kinematics
        @param[in] h_generated The generated distribution binned in true kinematics
        @param[in] addSystematics Do you want to fold the systematics too?
        */
      bool FoldHisto(
          PlotUtils::MnvH1D* &h_folded,
          const PlotUtils::MnvH2D* h_migration,
          const PlotUtils::MnvH1D* h_generated,
          bool addSystematics = true
          ) const;
      //@}
      //=========== END OF 1D FOLDING =============//
      //=============================================//



      //===================================
      //===================================
      /** @name 2D Unfolding
        Functions to unfold in 2 dimensions.
        DocDB 8687 and 8773 have more detail on 2D unfolding. 

        2D unfolding is handled via the \ref MnvResponse class. Create
        a MnvResponse, fill it, and get the necessary histograms out
        with \ref MnvResponse::GetMigrationObjects().

        @{*/

      /*! Unfold one 2D histogram.
        @param[out] h_unfold The resulting unfolded histo.  Creates a *new* histogram if input is null.
        @param[out] cov covariance of unfolding error
        @param[in] h_migration Migration matrix of reco(X) vs true(Y) kinematics
        @param[in] h_mc_reco histogram binned in reco kinematics
        @param[in] h_mc_true histogram binned in true kinematics
        @param[in] h_data the folded histogram you will unfold
        @param[in] regparam Regulates the unfolding algorithm (which is iterative Bayesian)
        @return true if all is OK
       */
      bool UnfoldHisto2D(
          TH2D* &h_unfold, 
          TMatrixD &cov, 
          const TH2D* h_migration, 
          const TH2D* h_mc_reco, 
          const TH2D* h_mc_true, 
          const TH2D *h_data, 
          const Double_t regparam = 4. 
          ) ;

      /*! Unfold one 2D histogram and ignore unfolding covariance.
        @param[out] h_unfold The resulting unfolded histo.  Creates a *new* histogram if input is null.
        @param[in] h_migration Migration matrix of reco(X) vs true(Y) kinematics
        @param[in] h_mc_reco histogram binned in reco kinematics
        @param[in] h_mc_true histogram binned in true kinematics
        @param[in] h_data the folded histogram you will unfold
        @param[in] regparam Regulates the unfolding algorithm (which is iterative Bayesian)
        @return true if all is OK
       */
      bool UnfoldHisto2D(
          TH2D* &h_unfold, 
          const TH2D* h_migration, 
          const TH2D* h_mc_reco, 
          const TH2D* h_mc_true, 
          const TH2D *h_data, 
          const Double_t regparam = 4. 
          ) ;

      /*! Unfold a TH2D using a specially formatted file (formatted how?)
        @todo this belongs in the ccqe package
        @param[out] h_unfold The resulting unfolded histo.  Creates a *new* histogram if input is null.
        @param[in] f_EventSelection A file that holds information needed for unfolding ccqe distrubutions
        @param[in] h_nobck the folded histogram you will unfold
        @param[in] regparam Regulates the unfolding algorithm
        @param[in] name_migration_enu_q2 Name of the migration matrix
        @param[in] name_migration_q2_mcreco name of hist with qe events in reco Q2 bins
        @param[in] name_migration_q2_mctrue name of hist with qe events in true Q2 bins
        @return true if all is OK
       */
      bool UnfoldHisto2D(
          TH2D* &h_unfold, 
          TFile *f_EventSelection, 
          const TH2D* h_nobck, 
          const Double_t regparam = 4., 
          const std::string& name_migration_enu_q2 = "h_enu_q2_migration_qe", 
          const std::string& name_enu_q2_mcreco = "h_enu_q2_qe", 
          const std::string& name_enu_q2_mctrue = "h_enu_q2_generated_qe" 
          ) ;

      /*! Unfold one 2D histogram and possibly add systematic uncertainties
        @param[out] h_unfold The resulting unfolded histo.  Creates a *new* histogram if input is null.
        @param[in] h_migration Migration matrix of reco(X) vs true(Y) kinematics
        @param[in] h_mc_reco histogram binned in reco kinematics
        @param[in] h_mc_true histogram binned in true kinematics
        @param[in] h_data the folded histogram you will unfold
        @param[in] regparam Regulates the unfolding algorithm (which is iterative Bayesian)
        @param[in] addSystematics Do you want to add many universe systematics?
        @return true if all is OK
       */      
      bool UnfoldHisto2D(
          PlotUtils::MnvH2D* &h_unfold, 
          const TH2D *h_migration, 
          const TH2D *h_mc_reco, 
          const TH2D *h_mc_true, 
          const PlotUtils::MnvH2D *h_data, 
          const Double_t regparam = 4., 
          bool addSystematics = true 
          ) ;

      /*! Unfold one 2D histogram and possibly add systematic uncertainties with many migration matrices
        @param[out] h_unfold The resulting unfolded histo.  Creates a *new* histogram if input is null.
        @param[in] h_migration Migration matrix of many universes of reco(X) vs true(Y) kinematics
        @param[in] h_mc_reco histogram binned in reco kinematics
        @param[in] h_mc_true histogram binned in true kinematics
        @param[in] h_data the folded histogram you will unfold
        @param[in] regparam Regulates the unfolding algorithm (which is iterative Bayesian)
        @param[in] addSystematics Do you want to add many universe systematics?
        @param[in] useSysVariatedMigrations Do you want to use a different migration matrix for each universe?
        @return true if all is OK
       */
      bool UnfoldHisto2D(
          PlotUtils::MnvH2D* &h_unfold, 
          const PlotUtils::MnvH2D *h_migration, 
          const PlotUtils::MnvH2D *h_mc_reco, 
          const PlotUtils::MnvH2D *h_mc_true, 
          const PlotUtils::MnvH2D *h_data, 
          const Double_t regparam = 4., 
          bool addSystematics = true, 
          bool useSysVariatedMigrations = true 
          ) ;

      /*! Unfold one 2D histogram using RooUnfoldResponse
        @param[out] h_unfold The resulting unfolded histo.  Creates a *new* histogram if input is null.
        @param[out] cov covariance of unfolding error
        @param[in] response Detector response to migration
        @param[in] h_data the folded histogram you will unfold
        @param[in] regparam Regulates the unfolding algorithm (which is iterative Bayesian)
        @return true if all is OK
       */
      bool UnfoldHisto2D(
          TH2D* &h_unfold, 
          TMatrixD &cov, 
          RooUnfoldResponse* response, 
          const TH2D *h_data, 
          const Double_t regparam = 4. 
          ) ;

      /*! Unfold one 2D histogram using RooUnfoldResponse and ignore unfolding covariance.
        @param[out] h_unfold The resulting unfolded histo.  Creates a *new* histogram if input is null.
        @param[in] response Detector response to migration
        @param[in] h_data the folded histogram you will unfold
        @param[in] regparam Regulates the unfolding algorithm (which is iterative Bayesian)
        @return true if all is OK
       */
      bool UnfoldHisto2D(
          TH2D* &h_unfold, 
          RooUnfoldResponse* response, 
          const TH2D *h_data, 
          const Double_t regparam = 4. 
          ) ;


      //@}
      //=========== END OF 2D UNFOLDING =============//
      //=============================================//

      //===================================
      //===================================
      /** @name 3D Unfolding
        Functions to unfold in 3 dimensions. 
        @{*/

      //3D unfolding
      bool UnfoldHisto3D(
          TH3D* &h_unfold, 
          TMatrixD &cov, 
          const TH2D* h_migration, 
          const TH3D* h_mc_reco, 
          const TH3D* h_mc_true, 
          const TH3D *h_data, 
          const Double_t regparam = 4. 
          ) const;

      bool UnfoldHisto3D(
          TH3D* &h_unfold, 
          const TH2D* h_migration, 
          const TH3D* h_mc_reco, 
          const TH3D* h_mc_true, 
          const TH3D *h_data, 
          const Double_t regparam = 4. 
          ) const;

      bool UnfoldHisto3D(
          PlotUtils::MnvH3D* &h_unfold, 
          const TH2D *h_migration, 
          const TH3D *h_mc_reco, 
          const TH3D *h_mc_true, 
          const PlotUtils::MnvH3D *h_data, 
          const Double_t regparam = 4., 
          bool addSystematics = true 
          ) const;

      //! Unfold3D using MnvResponse  migration, reco and true outputs
      bool UnfoldHisto3D(
          PlotUtils::MnvH3D* &h_unfold, 
          const PlotUtils::MnvH2D *h_migration, 
          const PlotUtils::MnvH3D *h_mc_reco, 
          const PlotUtils::MnvH3D *h_mc_true, 
          const PlotUtils::MnvH3D *h_data, 
          const Double_t regparam = 4., 
          bool addSystematics = true, 
          bool useSysVariatedMigrations = true 
	  ) const;

      //=========== END OF 3D UNFOLDING =============//


    template <class T>
    bool UnfoldHistoWithFakes(T* h_unfold, TMatrixD &covmx, const TH2* const h_migration, const T* const h_data, const T* const  h_model_reco=nullptr, const T* const h_model_truth=nullptr,  const T* const h_model_background = nullptr, double regparam=4) const;

    template <class T>
    bool UnfoldVertErrorBands(T* h_unfold, const PlotUtils::MnvH2D* const h_migration, const T* const h_data, const T* const  h_model_reco, const T* const h_model_truth,  const T* const h_model_background,std::string errorband_name, double regparam,bool useSysVariatedMigrations) const;
    template <class T>
    bool UnfoldLatErrorBands(T* h_unfold, const PlotUtils::MnvH2D* const h_migration, const T* const h_data, const T* const  h_model_reco, const T* const h_model_truth,  const T* const h_model_background, std::string errorband_name,double regparam,bool useSysVariatedMigrations) const;
    template <class T>
    bool UnfoldUncorrError(T* h_unfold, const PlotUtils::MnvH2D* const h_migration, const T* const h_data, const T* const  h_model_reco, const T* const h_model_truth,  const T* const h_model_background, std::string errorband_name,double regparam,bool useSysVariatedMigrations) const;


    bool UnfoldHistoWithFakes(PlotUtils::MnvH1D* &h_unfold, TMatrixD &covmx, const PlotUtils::MnvH2D* const h_migration, const PlotUtils::MnvH1D* const h_data, const PlotUtils::MnvH1D* const  h_model_reco=nullptr, const PlotUtils::MnvH1D* const h_model_truth=nullptr, const PlotUtils::MnvH1D* const h_model_background = nullptr, double regparam=4,bool addSystematics = false,bool useSysVariatedMigrations = false) const;

    bool UnfoldHistoWithFakes(PlotUtils::MnvH2D* &h_unfold, TMatrixD &covmx, const PlotUtils::MnvH2D* const h_migration, const PlotUtils::MnvH2D* const h_data, const PlotUtils::MnvH2D* const  h_model_reco=nullptr, const PlotUtils::MnvH2D* const h_model_truth=nullptr, const PlotUtils::MnvH2D* const h_model_background = nullptr, double regparam=4,bool addSystematics = false,bool useSysVariatedMigrations = false) const;

    bool UnfoldHistoWithFakes(PlotUtils::MnvH3D* &h_unfold, TMatrixD &covmx, const PlotUtils::MnvH2D* const h_migration, const PlotUtils::MnvH3D* const h_data, const PlotUtils::MnvH3D* const  h_model_reco=nullptr, const PlotUtils::MnvH3D* const h_model_truth=nullptr, const PlotUtils::MnvH3D* const h_model_background = nullptr, double regparam=4,bool addSystematics = false,bool useSysVariatedMigrations = false) const;


    
  private:
    //! Does \a h have content in the under/overflow bins? If so, the 2D unfolding won't work
    bool HasOverflowContent(const TH2D* h) const;
  };

}

#endif
